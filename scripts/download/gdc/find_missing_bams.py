#!/usr/bin/env python3
"""
find_missing_bams.py  —  Audit all MS-searched TMT plexes for missing GDC BAMs.

For every real patient channel across all plexes in plex_list.txt, checks whether
that (case_submitter_id, sample_type) has a UUID in gdc_meta_matched.tsv. If not,
queries the GDC API to locate the BAM and outputs a download manifest + metadata
rows ready to append to gdc_meta_matched.tsv.

Outputs (all in OUT_DIR = DATA/MS_SEARCH/missing_bams/):
  missing_manifest.tsv  — GDC-format manifest ready for download.py / gdc-client
  missing_metadata.tsv  — rows to append to gdc_meta_matched.tsv
  missing_report.tsv    — full audit (found, not_found, skipped)

Usage (on cluster from repo root):
  python3 scripts/download/gdc/find_missing_bams.py

Environment overrides:
  REPO=/path/to/repo  DATA=/path/to/scratch
"""

import os, re, json, time, logging
import urllib.parse, urllib.request
import pandas as pd
from collections import defaultdict

# ── CONFIG ─────────────────────────────────────────────────────────────────────
REPO      = os.environ.get("REPO", "/home/leduc.an/AAS_Evo_project/AAS_Evo")
DATA      = os.environ.get("DATA", "/scratch/leduc.an/AAS_Evo")
TMT_MAP   = f"{REPO}/metadata/PDC_meta/pdc_file_tmt_map.tsv"
GDC_META  = f"{REPO}/metadata/GDC_meta/gdc_meta_matched.tsv"
PLEX_LIST = f"{DATA}/MS_SEARCH/plex_list.txt"
OUT_DIR   = f"{DATA}/MS_SEARCH/missing_bams"
GDC_API   = "https://api.gdc.cancer.gov/files"
API_DELAY = 0.5   # seconds between GDC API calls to avoid rate limiting

# CPTAC patient ID pattern — C3L-XXXXX or C3N-XXXXX
CPTAC_CASE_RE = re.compile(r'^C3[LN]-\d{5}$')

# case_submitter_id values that are never real patients
SKIP_CASE_LOWER = {"ref", "reference", "pooled", "pooled sample", "pool",
                   "nan", "not reported", ""}

# sample_types we intentionally exclude from BAM processing
SKIP_SAMPLE_TYPES = {"Blood Derived Normal", "Not Reported", "Cell Lines"}

TMT_CHANNEL_MAP = {
    "tmt_126": "126",  "tmt_127n": "127N", "tmt_127c": "127C",
    "tmt_128n": "128N","tmt_128c": "128C", "tmt_129n": "129N",
    "tmt_129c": "129C","tmt_130n": "130N", "tmt_130c": "130C",
    "tmt_131":  "131N","tmt_131c": "131C",
}
# ───────────────────────────────────────────────────────────────────────────────

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


def gdc_query(case_id: str, sample_type: str) -> list:
    """Query GDC API for WXS BAM files matching case + sample_type.
    Returns list of file-metadata dicts from the GDC hits array."""
    filt = json.dumps({
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.submitter_id",
                                    "value": case_id}},
            {"op": "=", "content": {"field": "experimental_strategy",
                                    "value": "WXS"}},
            {"op": "=", "content": {"field": "data_format",
                                    "value": "BAM"}},
            {"op": "=", "content": {"field": "cases.samples.sample_type",
                                    "value": sample_type}},
        ]
    })
    fields = ",".join([
        "file_id", "file_name", "md5sum", "file_size", "state",
        "cases.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.tissue_type",
        "cases.project.primary_site",
        "cases.project.disease_type",
        "cases.project.project_id",
    ])
    url = GDC_API + "?" + urllib.parse.urlencode({
        "filters": filt,
        "fields":  fields,
        "size":    "10",
        "format":  "json",
    })
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            data = json.loads(r.read())
        return data["data"]["hits"]
    except Exception as e:
        log.warning(f"    API error for {case_id} ({sample_type}): {e}")
        return []


def extract_meta(hit: dict, case_id: str, sample_type: str) -> dict:
    """Flatten a GDC API hit into a gdc_meta_matched.tsv-compatible row."""
    cases   = hit.get("cases", [{}])
    c       = cases[0] if cases else {}
    samples = c.get("samples", [])
    sample  = next(
        (s for s in samples if s.get("sample_type") == sample_type),
        samples[0] if samples else {}
    )
    proj = c.get("project", {})
    return {
        "file_id":           hit["file_id"],
        "file_name":         hit.get("file_name", ""),
        "case_submitter_id": case_id,
        "sample_type":       sample_type,
        "tissue_type":       sample.get("tissue_type", ""),
        "primary_site":      proj.get("primary_site", ""),
        "disease_type":      proj.get("disease_type", ""),
        "project_id":        proj.get("project_id", ""),
    }


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load inputs ───────────────────────────────────────────────────────────
    log.info("Loading TMT map and GDC metadata...")
    tmt = pd.read_csv(TMT_MAP, sep="\t")
    gdc = pd.read_csv(GDC_META, sep="\t")
    if "file_id" in gdc.columns and "gdc_file_id" not in gdc.columns:
        gdc = gdc.rename(columns={"file_id": "gdc_file_id"})

    # Build set of (case_submitter_id, sample_type) already in our metadata
    have_meta = set(zip(gdc["case_submitter_id"], gdc["sample_type"]))

    with open(PLEX_LIST) as f:
        plex_ids = [l.strip() for l in f if l.strip()]
    log.info(f"Scanning {len(plex_ids):,} plexes...")

    # ── Scan all plexes ───────────────────────────────────────────────────────
    # Collect unique (case, sample_type) pairs missing from gdc_meta,
    # tracking which plexes they appear in (for reporting).
    missing_pairs: dict[tuple, set] = defaultdict(set)   # (case, stype) → plexes
    skipped: list[dict] = []

    for plex_id in plex_ids:
        plex_tmt = (
            tmt[tmt["run_metadata_id"] == plex_id]
            [["tmt_channel", "case_submitter_id", "sample_type"]]
            .drop_duplicates()
        )
        plex_tmt["channel"] = plex_tmt["tmt_channel"].map(TMT_CHANNEL_MAP)

        for _, row in plex_tmt.iterrows():
            case  = str(row.get("case_submitter_id", "")).strip()
            stype = str(row.get("sample_type", "")).strip()

            # Non-patient channels (pooled, reference, etc.)
            if case.lower() in SKIP_CASE_LOWER:
                skipped.append({"case": case, "sample_type": stype,
                                 "plex_id": plex_id, "reason": "non_patient"})
                continue

            # Non-CPTAC IDs: GTEx samples, NX### external references, QC lanes
            if not CPTAC_CASE_RE.match(case):
                skipped.append({"case": case, "sample_type": stype,
                                 "plex_id": plex_id, "reason": "non_cptac_id"})
                continue

            # Sample types we intentionally skip
            if stype in SKIP_SAMPLE_TYPES:
                skipped.append({"case": case, "sample_type": stype,
                                 "plex_id": plex_id, "reason": "skip_sample_type"})
                continue

            # Check if this (case, sample_type) is already in our metadata
            if (case, stype) not in have_meta:
                missing_pairs[(case, stype)].add(plex_id)

    log.info(f"Found {len(missing_pairs):,} unique (case, sample_type) pairs "
             f"missing from gdc_meta_matched.tsv")

    # ── Query GDC API ─────────────────────────────────────────────────────────
    report_rows:   list[dict] = []
    manifest_rows: list[dict] = []
    metadata_rows: list[dict] = []

    seen_file_ids: set[str] = set()

    for i, ((case, stype), plexes) in enumerate(sorted(missing_pairs.items()), 1):
        log.info(f"  [{i}/{len(missing_pairs)}] {case}  {stype}")
        hits = gdc_query(case, stype)
        time.sleep(API_DELAY)

        if not hits:
            report_rows.append({
                "case":        case,
                "sample_type": stype,
                "plexes":      ";".join(sorted(plexes)),
                "status":      "not_found_in_gdc",
                "file_id":     "",
                "file_name":   "",
            })
            continue

        for hit in hits:
            fid = hit["file_id"]
            report_rows.append({
                "case":        case,
                "sample_type": stype,
                "plexes":      ";".join(sorted(plexes)),
                "status":      "found_in_gdc",
                "file_id":     fid,
                "file_name":   hit.get("file_name", ""),
            })
            if fid not in seen_file_ids:
                seen_file_ids.add(fid)
                manifest_rows.append({
                    "id":       fid,
                    "filename": hit.get("file_name", ""),
                    "md5":      hit.get("md5sum", ""),
                    "size":     hit.get("file_size", ""),
                    "state":    hit.get("state", "released"),
                })
                metadata_rows.append(extract_meta(hit, case, stype))

    # ── Write outputs ─────────────────────────────────────────────────────────
    report_df = pd.DataFrame(report_rows)
    report_df.to_csv(f"{OUT_DIR}/missing_report.tsv", sep="\t", index=False)
    log.info(f"Report → {OUT_DIR}/missing_report.tsv")

    if manifest_rows:
        manifest_df = pd.DataFrame(manifest_rows).drop_duplicates(subset=["id"])
        manifest_df.to_csv(f"{OUT_DIR}/missing_manifest.tsv", sep="\t", index=False)
        log.info(f"Manifest ({len(manifest_df)} BAMs) → {OUT_DIR}/missing_manifest.tsv")

    if metadata_rows:
        meta_df = pd.DataFrame(metadata_rows).drop_duplicates(subset=["file_id"])
        meta_df.to_csv(f"{OUT_DIR}/missing_metadata.tsv", sep="\t", index=False)
        log.info(f"Metadata ({len(meta_df)} rows) → {OUT_DIR}/missing_metadata.tsv")

    # ── Summary ───────────────────────────────────────────────────────────────
    n_found   = (report_df["status"] == "found_in_gdc").sum()    if len(report_df) else 0
    n_missing = (report_df["status"] == "not_found_in_gdc").sum() if len(report_df) else 0
    by_reason = pd.Series([s["reason"] for s in skipped]).value_counts()

    print(f"""
╔══════════════════════════════════════════════════════╗
║              Missing BAM Audit — Summary              ║
╠══════════════════════════════════════════════════════╣
  Plexes scanned:               {len(plex_ids):>6,}
  Missing (case, sample_type):  {len(missing_pairs):>6,}
    ✓ Found in GDC:             {n_found:>6,}
    ✗ Not found in GDC:         {n_missing:>6,}
  Skipped channels:             {len(skipped):>6,}
    · non_patient (pooled/ref): {by_reason.get("non_patient", 0):>6,}
    · non_cptac_id (GTEx/NX):   {by_reason.get("non_cptac_id", 0):>6,}
    · skip_sample_type (blood): {by_reason.get("skip_sample_type", 0):>6,}
╠══════════════════════════════════════════════════════╣
  Next steps:
  1. Download missing BAMs:
     python3 scripts/download/gdc/download.py \\
         --manifest {OUT_DIR}/missing_manifest.tsv \\
         --token metadata/GDC_meta/.gdc-user-token.txt

  2. Append metadata to gdc_meta_matched.tsv:
     tail -n +2 {OUT_DIR}/missing_metadata.tsv \\
         >> metadata/GDC_meta/gdc_meta_matched.tsv

  3. Run variant calling on downloaded BAMs:
     find {DATA}/BAMS -name "*.bam" > {DATA}/bam_list_missing.txt
     CHUNK_NAME=missing BAMS_LIST={DATA}/bam_list_missing.txt \\
         bash scripts/proc_bams/run_pipeline.sh variant-call missing

  4. Run VEP, then regenerate mutant FASTAs for new UUIDs.
╚══════════════════════════════════════════════════════╝
""")


if __name__ == "__main__":
    main()
