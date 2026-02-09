#!/usr/bin/env python3
"""
Fetch PDC metadata via the GraphQL API.

Queries the PDC API for each study in studies.tsv and writes per-tissue CSVs
in the same format as the PDC portal exports:
  - PDC_file_manifest.csv    (file download info + signed URLs)
  - PDC_study_biospecimen.csv (sample/patient metadata)
  - PDC_study_experimental.csv (TMT plex layout)

These CSVs are consumed by consolidate_metadata.py downstream.

Usage:
    python3 scripts/setup/fetch_pdc_metadata.py

    # Specific tissues only:
    python3 scripts/setup/fetch_pdc_metadata.py --tissues CCRCC LUAD

    # Skip signed URL fetching (faster, URLs added later via refresh_urls.py):
    python3 scripts/setup/fetch_pdc_metadata.py --skip-urls
"""

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import ssl
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import META_DIR

PDC_GRAPHQL_URL = "https://pdc.cancer.gov/graphql"
STUDIES_FILE = META_DIR / "studies.tsv"
PDC_META_DIR = META_DIR / "PDC_meta"

REQUEST_TIMEOUT = 120
RETRY_ATTEMPTS = 3
REQUEST_DELAY = 1.0  # seconds between API calls

TMT_CHANNELS = [
    "tmt_126", "tmt_127n", "tmt_127c", "tmt_128n", "tmt_128c",
    "tmt_129n", "tmt_129c", "tmt_130n", "tmt_130c", "tmt_131", "tmt_131c",
]


def _make_ssl_context():
    """Create SSL context, falling back to unverified if certs unavailable."""
    try:
        ctx = ssl.create_default_context()
        # Test that it works
        return ctx
    except Exception:
        ctx = ssl._create_unverified_context()
        return ctx


def query_pdc(graphql_query, variables=None):
    """Send a GraphQL query to PDC API with retries."""
    payload = json.dumps({"query": graphql_query, "variables": variables or {}})
    ctx = _make_ssl_context()

    for attempt in range(RETRY_ATTEMPTS):
        try:
            req = urllib.request.Request(
                PDC_GRAPHQL_URL,
                data=payload.encode("utf-8"),
                headers={"Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT, context=ctx) as resp:
                result = json.loads(resp.read())

            if "errors" in result:
                print(f"  GraphQL errors: {result['errors'][0]['message']}")
                return None
            return result

        except Exception as e:
            if attempt < RETRY_ATTEMPTS - 1:
                wait = 2 ** attempt
                print(f"  Retry {attempt + 1}/{RETRY_ATTEMPTS} after error: {e}")
                time.sleep(wait)
            else:
                print(f"  ERROR: API request failed after {RETRY_ATTEMPTS} attempts: {e}")
                return None

    return None


# ---------------------------------------------------------------------------
# GraphQL queries
# ---------------------------------------------------------------------------

def fetch_files(pdc_study_id, include_urls=True):
    """Fetch all RAW files for a study via filesPerStudy."""
    url_field = "signedUrl { url }" if include_urls else ""
    query = f"""
    query {{
        filesPerStudy(pdc_study_id: "{pdc_study_id}", data_category: "Raw Mass Spectra") {{
            file_id
            file_name
            file_type
            data_category
            md5sum
            file_size
            file_location
            study_submitter_id
            pdc_study_id
            study_id
            {url_field}
        }}
    }}
    """
    result = query_pdc(query)
    if result and "data" in result:
        return result["data"].get("filesPerStudy", [])
    return []


def fetch_biospecimen(pdc_study_id):
    """Fetch biospecimen data for a study via biospecimenPerStudy."""
    query = f"""
    query {{
        biospecimenPerStudy(pdc_study_id: "{pdc_study_id}") {{
            aliquot_id
            aliquot_submitter_id
            sample_id
            sample_submitter_id
            case_id
            case_submitter_id
            project_name
            sample_type
            primary_site
            disease_type
            aliquot_is_ref
            aliquot_status
        }}
    }}
    """
    result = query_pdc(query)
    if result and "data" in result:
        return result["data"].get("biospecimenPerStudy", [])
    return []


def fetch_experimental_design(pdc_study_id):
    """Fetch experimental design (TMT plex layout) for a study."""
    tmt_fields = "\n".join(
        f"{ch} {{ aliquot_id aliquot_submitter_id aliquot_run_metadata_id label }}"
        for ch in TMT_CHANNELS
    )
    query = f"""
    query {{
        studyExperimentalDesign(pdc_study_id: "{pdc_study_id}") {{
            study_run_metadata_id
            study_run_metadata_submitter_id
            experiment_type
            experiment_number
            plex_dataset_name
            analyte
            acquisition_type
            polarity
            {tmt_fields}
        }}
    }}
    """
    result = query_pdc(query)
    if result and "data" in result:
        return result["data"].get("studyExperimentalDesign", [])
    return []


# ---------------------------------------------------------------------------
# File-to-run matching
# ---------------------------------------------------------------------------

def derive_run_metadata_id(file_name, known_run_ids):
    """Derive Run Metadata ID from RAW file name.

    RAW files follow: {plex}CPTAC_{tissue}_W_{inst}_{date}_{instr}_f{frac}.raw
    Run IDs follow:   {plex}CPTAC_{tissue}_Proteome_{inst}_{date}

    Strategy: replace _W_ with _Proteome_, truncate after 8-digit date,
    then verify against known run IDs.
    """
    # Try pattern: prefix_W_inst_YYYYMMDD_...
    m = re.match(r"(.+?)_W_(.+?_)(\d{8})", file_name)
    if m:
        candidate = m.group(1) + "_Proteome_" + m.group(2) + m.group(3)
        if candidate in known_run_ids:
            return candidate

    # Fallback: match by plex number prefix and date
    plex_match = re.match(r"(\d+)", file_name)
    date_match = re.search(r"(\d{8})", file_name)
    if plex_match and date_match:
        plex_num = plex_match.group(1)
        date = date_match.group(1)
        for run_id in known_run_ids:
            if run_id.startswith(plex_num) and date in run_id:
                return run_id

    return ""


# ---------------------------------------------------------------------------
# CSV writers (portal-compatible format)
# ---------------------------------------------------------------------------

def write_file_manifest(files, runs, tissue_dir, study_name, pdc_study_id,
                        include_urls=True):
    """Write portal-format file manifest CSV."""
    known_run_ids = {
        r["study_run_metadata_submitter_id"]
        for r in runs
        if r.get("study_run_metadata_submitter_id")
    }

    output = tissue_dir / "PDC_file_manifest.csv"
    headers = [
        "File ID", "File Name", "Run Metadata ID", "Protocol", "Study Name",
        "PDC Study ID", "PDC Study Version", "Study ID", "Project Name",
        "Data Category", "File Type", "Access", "File Size (in bytes)",
        "Md5sum", "Downloadable", "File Download Link",
    ]

    matched = 0
    with open(output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for file in files:
            file_name = file.get("file_name", "")
            run_id = derive_run_metadata_id(file_name, known_run_ids)
            if run_id:
                matched += 1

            signed_url = ""
            if include_urls:
                su = file.get("signedUrl")
                if su:
                    signed_url = su.get("url", "")

            writer.writerow([
                file.get("file_id", ""),
                file_name,
                run_id,
                "",  # Protocol (not available via API)
                study_name,
                pdc_study_id,
                "",  # PDC Study Version
                file.get("study_id", ""),
                "",  # Project Name
                file.get("data_category", ""),
                file.get("file_type", ""),
                "Open",
                file.get("file_size", ""),
                file.get("md5sum", ""),
                "Yes",
                signed_url,
            ])

    unmatched = len(files) - matched
    if unmatched > 0:
        print(f"    WARNING: {unmatched}/{len(files)} files could not be matched to a run")

    return output


def write_biospecimen(biospecimens, tissue_dir):
    """Write portal-format biospecimen CSV."""
    output = tissue_dir / "PDC_study_biospecimen.csv"
    headers = [
        "Aliquot ID", "Aliquot Submitter ID", "Sample ID",
        "Sample Submitter ID", "Case ID", "Case Submitter ID",
        "Project Name", "Sample Type", "Primary Site", "Disease Type",
        "Aliquot Is Ref", "Aliquot Status",
    ]

    with open(output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for b in biospecimens:
            writer.writerow([
                b.get("aliquot_id", ""),
                b.get("aliquot_submitter_id", ""),
                b.get("sample_id", ""),
                b.get("sample_submitter_id", ""),
                b.get("case_id", ""),
                b.get("case_submitter_id", ""),
                b.get("project_name", ""),
                b.get("sample_type", ""),
                b.get("primary_site", ""),
                b.get("disease_type", ""),
                b.get("aliquot_is_ref", ""),
                b.get("aliquot_status", ""),
            ])

    return output


def write_experimental_design(runs, aliquot_to_case, tissue_dir):
    """Write portal-format experimental design CSV.

    The portal format stores TMT channels as 'CaseID\\nSampleType'.
    The API returns aliquot_submitter_id per channel, so we join
    via the biospecimen data to get case_submitter_id + sample_type.
    """
    output = tissue_dir / "PDC_study_experimental.csv"
    headers = ["Run Metadata ID", "Number of Fractions", "Protocol"] + TMT_CHANNELS

    with open(output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for run in runs:
            run_id = run.get("study_run_metadata_submitter_id", "")

            row = [run_id, "", ""]  # Number of Fractions, Protocol not in API

            for ch in TMT_CHANNELS:
                channel_data = run.get(ch) or []
                if channel_data:
                    # Take first aliquot in channel (usually just one)
                    aliquot = channel_data[0]
                    aliquot_sub_id = aliquot.get("aliquot_submitter_id", "")

                    # Look up case info
                    case_info = aliquot_to_case.get(aliquot_sub_id, {})
                    case_id = case_info.get("case_submitter_id", aliquot_sub_id)
                    sample_type = case_info.get("sample_type", "")

                    row.append(f"{case_id}\n{sample_type}")
                else:
                    row.append("")

            writer.writerow(row)

    return output


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_studies(tissues_filter=None):
    """Load study registry from studies.tsv."""
    studies = []
    with open(STUDIES_FILE, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if tissues_filter and row["tissue_folder"] not in tissues_filter:
                continue
            studies.append(row)
    return studies


def process_study(tissue_folder, pdc_study_id, study_name, include_urls=True):
    """Fetch and write all metadata for one study."""
    tissue_dir = PDC_META_DIR / tissue_folder
    os.makedirs(tissue_dir, exist_ok=True)

    # 1. Experimental design (needed for file-to-run matching)
    print(f"  Fetching experimental design...")
    runs = fetch_experimental_design(pdc_study_id)
    print(f"    {len(runs)} TMT plexes")
    time.sleep(REQUEST_DELAY)

    # 2. Biospecimen (needed for aliquot→case mapping in experimental design)
    print(f"  Fetching biospecimen data...")
    biospecimens = fetch_biospecimen(pdc_study_id)
    print(f"    {len(biospecimens)} biospecimens")
    time.sleep(REQUEST_DELAY)

    # Build aliquot_submitter_id → case info map
    aliquot_to_case = {}
    for b in biospecimens:
        aliquot_sub_id = b.get("aliquot_submitter_id", "")
        if aliquot_sub_id:
            aliquot_to_case[aliquot_sub_id] = {
                "case_submitter_id": b.get("case_submitter_id", ""),
                "sample_type": b.get("sample_type", ""),
            }

    # 3. Files (RAW only, with optional signed URLs)
    print(f"  Fetching RAW file manifest...")
    files = fetch_files(pdc_study_id, include_urls=include_urls)
    print(f"    {len(files)} RAW files")
    time.sleep(REQUEST_DELAY)

    # 4. Write CSVs
    write_file_manifest(files, runs, tissue_dir, study_name, pdc_study_id,
                        include_urls)
    write_biospecimen(biospecimens, tissue_dir)
    write_experimental_design(runs, aliquot_to_case, tissue_dir)

    return len(files), len(biospecimens), len(runs)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch PDC metadata via GraphQL API"
    )
    parser.add_argument(
        "--tissues", nargs="+",
        help="Specific tissue folders to fetch (default: all in studies.tsv)"
    )
    parser.add_argument(
        "--skip-urls", action="store_true",
        help="Skip signed URL fetching (faster; add URLs later via refresh_urls.py)"
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Fetch PDC Metadata")
    print("=" * 60)

    studies = load_studies(set(args.tissues) if args.tissues else None)
    if not studies:
        print("No studies found in studies.tsv")
        sys.exit(1)

    print(f"\nStudies to fetch: {len(studies)}")
    if args.skip_urls:
        print("  (skipping signed URLs)")
    print()

    total_files = 0
    total_bio = 0
    total_plexes = 0

    for i, study in enumerate(studies):
        tissue = study["tissue_folder"]
        pdc_id = study["pdc_study_id"]
        name = study["study_name"]

        print(f"[{i+1}/{len(studies)}] {tissue} ({pdc_id})")
        n_files, n_bio, n_plexes = process_study(
            tissue, pdc_id, name, include_urls=not args.skip_urls
        )
        total_files += n_files
        total_bio += n_bio
        total_plexes += n_plexes
        print()

    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"  Studies:      {len(studies)}")
    print(f"  RAW files:    {total_files}")
    print(f"  Biospecimens: {total_bio}")
    print(f"  TMT plexes:   {total_plexes}")
    print(f"\nOutput: {PDC_META_DIR}")
    print("\nNext step:")
    print("  python3 scripts/download/pdc/consolidate_metadata.py")


if __name__ == "__main__":
    main()
