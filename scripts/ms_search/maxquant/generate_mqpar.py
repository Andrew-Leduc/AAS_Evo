# -*- coding: utf-8 -*-
"""
generate_mqpar.py

Generate MaxQuant mqpar.xml parameter files for per-plex TMT-11 MS searches.

MaxQuant handles decoys internally (reversed sequences) — no separate decoy
step needed. With minPeptidesProtein=1, single-peptide protein identifications
are reported, which is critical for our mutant entries (each has exactly one
tryptic peptide). This avoids the ProteinProphet single-peptide penalty.

Results land in {output_dir}/{plex_id}/combined/

Usage:
    python3 generate_mqpar.py \
        --tmt-map /path/to/pdc_file_tmt_map.tsv \
        --raw-dir /path/to/RAW/ \
        --fasta-dir /path/to/FASTA/per_plex/ \
        -o /path/to/MQ_SEARCH/ \
        [--threads 16]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from xml.etree import ElementTree as ET
from xml.dom import minidom


TMT_CHANNEL_MAP = {
    "tmt_126":  "126",
    "tmt_127n": "127N",
    "tmt_127c": "127C",
    "tmt_128n": "128N",
    "tmt_128c": "128C",
    "tmt_129n": "129N",
    "tmt_129c": "129C",
    "tmt_130n": "130N",
    "tmt_130c": "130C",
    "tmt_131":  "131N",
    "tmt_131c": "131C",
}

# TMT11plex reporter ion masses (monoisotopic, Da)
TMT11_CHANNELS = [
    ("TMT11-126",  126.127726),
    ("TMT11-127N", 127.124761),
    ("TMT11-127C", 127.131081),
    ("TMT11-128N", 128.128116),
    ("TMT11-128C", 128.134436),
    ("TMT11-129N", 129.131471),
    ("TMT11-129C", 129.137791),
    ("TMT11-130N", 130.134826),
    ("TMT11-130C", 130.141146),
    ("TMT11-131N", 131.138181),
    ("TMT11-131C", 131.144501),
]


def load_tmt_map(tmt_map_path):
    plex_to_files    = defaultdict(set)
    plex_to_channels = defaultdict(list)
    with open(tmt_map_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            run_id      = row["run_metadata_id"].strip()
            file_name   = row["file_name"].strip()
            channel     = row.get("tmt_channel", "").strip().lower()
            case_id     = row.get("case_submitter_id", "").strip()
            sample_type = row.get("sample_type", "").strip()
            plex_to_files[run_id].add(file_name)
            plex_to_channels[run_id].append({
                "file_name":   file_name,
                "channel":     channel,
                "case_id":     case_id,
                "sample_type": sample_type,
            })
    return plex_to_files, plex_to_channels


def channel_label(case_id, sample_type):
    if case_id.lower() in ("ref", "reference", "pooled", "pool", ""):
        return "Reference"
    return f"{case_id}_{sample_type}" if sample_type else case_id


def sub(parent, tag, text=None):
    el = ET.SubElement(parent, tag)
    if text is not None:
        el.text = str(text)
    return el


def build_mqpar(raw_files, fasta_path, out_dir, plex_id,
                channel_entries, n_threads):
    """
    Build and return a MaxQuant mqpar.xml ElementTree for one plex.

    Key choices:
      - minPeptidesProtein = 1  → single-peptide proteins reported (our mutants)
      - decoyMode = revert      → MaxQuant generates reversed decoys internally
      - TMT11plex fixed mods    → correct TMT quantification
      - reporterPIF = 0.75      → precursor ion fraction filter (recommended)
    """
    root = ET.Element("MaxQuantParams")

    # ── FASTA ─────────────────────────────────────────────────────────────────
    fastas = sub(root, "fastaFiles")
    fi = sub(fastas, "FastaFileInfo")
    sub(fi, "fastaFilePath",         fasta_path)
    # Standard UniProt parse rule — works for sp|ACC|ENTRY and sp|ACC-swap-hash|ENTRY-mut
    sub(fi, "identifierParseRule",   r">.*\|(.+)\|.*")
    sub(fi, "descriptionParseRule",  r">(.*)")
    sub(fi, "taxonomyParseRule",     "")
    sub(fi, "variationParseRule",    "")
    sub(fi, "modificationParseRule", "")
    sub(fi, "taxonomyId",            "")

    sub(root, "fastaFilesFirstSearch")

    # ── RAW FILES ─────────────────────────────────────────────────────────────
    fp_el = sub(root, "filePaths")
    exp_el = sub(root, "experiments")
    frac_el = sub(root, "fractions")
    ptm_el = sub(root, "ptms")
    pg_el = sub(root, "paramGroupIndices")

    for raw in sorted(raw_files):
        sub(fp_el, "string", raw)
        sub(exp_el, "string", plex_id)
        sub(frac_el, "short", 32767)    # 32767 = no fractionation
        sub(ptm_el, "boolean", "False")
        sub(pg_el, "int", 0)

    sub(root, "referenceChannel", "")

    # ── GLOBAL SETTINGS ───────────────────────────────────────────────────────
    sub(root, "numThreads",            n_threads)
    sub(root, "emailAddress",          "")
    sub(root, "smtpHost",              "")
    sub(root, "emailFromAddress",      "")
    sub(root, "fixedCombinedFolder",   out_dir)
    sub(root, "fullMinMz",             -1.79769313486232e+308)
    sub(root, "fullMaxMz",              1.79769313486232e+308)
    sub(root, "sendEmail",             "False")
    sub(root, "ionCountIntensities",   "False")
    sub(root, "verboseColumnOutput",   "True")
    sub(root, "writeMsScansTable",     "False")
    sub(root, "writeMsmsScansTable",   "True")
    sub(root, "writeMs3ScansTable",    "False")
    sub(root, "writeMsmsTable",        "True")
    sub(root, "writeModificationPositionTable", "False")
    sub(root, "writeUnresolvedMsmsTable",       "False")
    sub(root, "ibaq",                  "True")
    sub(root, "ibaqLogFit",            "False")
    sub(root, "separateLfq",          "False")
    sub(root, "lfqStabilizeLargeRatios", "True")
    sub(root, "lfqRequireMsms",        "True")
    sub(root, "decoyMode",             "revert")
    sub(root, "includeContaminants",   "True")
    sub(root, "maxPeptideMass",        4600.0)
    sub(root, "epsilonMutationScore",  0.01)
    sub(root, "mutatedPeptidesSeparately", "True")
    sub(root, "useMultModifications",  "True")
    sub(root, "maxMissedCleavages",    2)

    # FDR thresholds
    sub(root, "psmFdr",                0.01)
    sub(root, "proteinFdr",            0.01)
    sub(root, "siteFdr",               0.01)
    sub(root, "minPeptidesProtein",    1)      # ← key: allow single-peptide hits
    sub(root, "minScoreForMBR",        0)
    sub(root, "useNormRatiosForPeptideAndProteinResults", "True")
    sub(root, "minUniquePeptides",     0)      # 0 = no unique-peptide requirement
    sub(root, "minRazorPeptides",      1)

    # Quantification
    sub(root, "quantMode",             1)      # 1 = reporter ion (TMT/iTRAQ)
    sub(root, "reporterMassTolerance", 0.003)  # 3 mDa (high-res MS2)
    sub(root, "reporterPIF",           0.75)
    sub(root, "filterPIF",             "True")
    sub(root, "lcmsRunType",           "Reporter ion MS2")

    # ── TMT CHANNEL LABELS ────────────────────────────────────────────────────
    # Build channel_key → label from channel_entries (dedup per channel)
    seen = {}
    for entry in channel_entries:
        ch = entry["channel"]
        if ch not in seen:
            seen[ch] = channel_label(entry["case_id"], entry["sample_type"])

    iso_el = sub(root, "isobaricLabels")
    for internal_label, mass in TMT11_CHANNELS:
        # Map TMT11-127N → 127N → look up in seen via TMT_CHANNEL_MAP values
        ch_key = internal_label.replace("TMT11-", "").replace("TMT11-", "")
        # find which tmt_channel maps to this channel string
        sample_name = next(
            (seen[k] for k, v in TMT_CHANNEL_MAP.items()
             if v == ch_key and k in seen),
            ch_key)   # fall back to channel name if not in annotation

        ili = sub(iso_el, "IsobaricLabelInfo")
        sub(ili, "internalLabel",        internal_label)
        sub(ili, "terminalLabel",        internal_label)
        sub(ili, "correctionFactorM2",   0)
        sub(ili, "correctionFactorM1",   0)
        sub(ili, "correctionFactorP1",   0)
        sub(ili, "correctionFactorP2",   0)
        sub(ili, "tmtLike",              "True")
        sub(ili, "reporterMass",         mass)
        sub(ili, "name",                 sample_name)

    # ── PARAMETER GROUP (enzyme, mods, mass tolerances) ───────────────────────
    pgs = sub(root, "parameterGroups")
    pg  = sub(pgs, "parameterGroup")

    # MS1 mass tolerance
    sub(pg, "maxCharge",        7)
    sub(pg, "msInstrument",     0)   # 0 = default (Orbitrap)
    sub(pg, "massTolerancePpm", 20)  # 20 ppm MS1

    # Fixed modifications: carbamidomethylation + TMT11 on K and peptide N-term
    fixed = sub(pg, "fixedModifications")
    sub(fixed, "string", "Carbamidomethyl (C)")
    sub(fixed, "string", "TMT11plex (K)")
    sub(fixed, "string", "TMT11plex (Peptide N-term)")

    # Variable modifications
    var = sub(pg, "variableModifications")
    sub(var, "string", "Oxidation (M)")
    sub(var, "string", "Acetyl (Protein N-term)")

    sub(pg, "useEntorhinal",                 "False")
    sub(pg, "horizon",                        100)
    sub(pg, "centroidMatchTol",               8)
    sub(pg, "centroidHalfWidth",              35)
    sub(pg, "valleyFactor",                   1.4)
    sub(pg, "isotopeValleyFactor",            1.2)
    sub(pg, "advancedPeakSplitting",          "False")
    sub(pg, "customProteinQuantification",    "False")
    sub(pg, "customProteinQuantificationFile","")
    sub(pg, "minScanNumber",                  0)
    sub(pg, "maxScanNumber",                  2147483647)

    # Enzyme
    enzymes = sub(pg, "enzymes")
    sub(enzymes, "string", "Trypsin/P")
    sub(pg, "enzymesFirstSearch")
    sub(pg, "allowedMissedCleavages", 2)

    # MS2 mass tolerance (Orbitrap)
    sub(pg, "msmsTolerancePpm", 20)

    # Peptide length and mass range
    sub(pg, "minPeptideLen",     7)
    sub(pg, "maxPeptideLen",     60)

    # Match between runs
    sub(pg, "matchBetweenRuns", "False")  # off for TMT (fractions handle it)

    sub(pg, "precursorMassTolerancePpm", 20)
    sub(pg, "msmsTolerancePpmForCalib",  20)

    return root


def prettify(root):
    raw = ET.tostring(root, encoding="unicode")
    parsed = minidom.parseString(raw)
    return parsed.toprettyxml(indent="  ", encoding=None)


def main():
    parser = argparse.ArgumentParser(
        description="Generate MaxQuant mqpar.xml files for per-plex TMT searches."
    )
    parser.add_argument("--tmt-map", required=True,
                        help="TMT mapping file (pdc_file_tmt_map.tsv)")
    parser.add_argument("--raw-dir", required=True,
                        help="Directory containing RAW files")
    parser.add_argument("--fasta-dir", required=True,
                        help="Directory containing per-plex FASTAs")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for MQ_SEARCH setup")
    parser.add_argument("--threads", type=int, default=16,
                        help="CPU threads per search (default: 16)")
    args = parser.parse_args()

    for path, label in [(args.tmt_map,  "TMT map"),
                        (args.raw_dir,  "RAW directory"),
                        (args.fasta_dir,"FASTA directory")]:
        if not os.path.exists(path):
            sys.exit(f"Error: {label} not found: {path}")

    mqpar_dir  = os.path.join(args.output, "mqpar")
    plex_list_path = os.path.join(args.output, "plex_list.txt")
    os.makedirs(mqpar_dir, exist_ok=True)

    print(f"Loading TMT map: {args.tmt_map}")
    plex_to_files, plex_to_channels = load_tmt_map(args.tmt_map)
    print(f"  {len(plex_to_files)} plexes")

    available_raws = {f for f in os.listdir(args.raw_dir)
                      if f.lower().endswith(".raw")}
    available_fastas = {f.replace(".fasta", "")
                        for f in os.listdir(args.fasta_dir)
                        if f.endswith(".fasta")}
    print(f"  {len(available_raws)} RAW files")
    print(f"  {len(available_fastas)} per-plex FASTAs")

    plex_list = []
    n_skip_fasta = n_skip_raw = 0

    for plex_id in sorted(plex_to_files.keys()):
        if plex_id not in available_fastas:
            n_skip_fasta += 1
            continue

        found_raws = [
            os.path.join(args.raw_dir, f)
            for f in sorted(plex_to_files[plex_id])
            if f in available_raws
        ]
        if not found_raws:
            n_skip_raw += 1
            continue

        fasta_path = os.path.abspath(
            os.path.join(args.fasta_dir, f"{plex_id}.fasta"))
        out_dir = os.path.join(args.output, "results", plex_id)
        os.makedirs(out_dir, exist_ok=True)

        root = build_mqpar(
            raw_files=found_raws,
            fasta_path=fasta_path,
            out_dir=os.path.abspath(out_dir),
            plex_id=plex_id,
            channel_entries=plex_to_channels[plex_id],
            n_threads=args.threads,
        )

        mqpar_path = os.path.join(mqpar_dir, f"{plex_id}.xml")
        with open(mqpar_path, "w") as f:
            f.write(prettify(root))

        plex_list.append(plex_id)
        print(f"  {plex_id}: {len(found_raws)} RAW files → {mqpar_path}")

    with open(plex_list_path, "w") as f:
        for p in plex_list:
            f.write(p + "\n")

    print(f"\nPlexes ready:             {len(plex_list)}")
    print(f"Skipped (no FASTA):       {n_skip_fasta}")
    print(f"Skipped (no RAW files):   {n_skip_raw}")
    print(f"mqpar files:              {mqpar_dir}")
    print(f"Plex list:                {plex_list_path}")
    print(f"\nNext: bash scripts/ms_search/maxquant/run_maxquant.sh")


if __name__ == "__main__":
    main()
