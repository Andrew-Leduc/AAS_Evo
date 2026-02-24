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
    "tmt_126":  "126C",
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

# MaxQuant internal modification names for TMT11plex.
# Names must exactly match titles in modifications.xml.
# Key naming rules (verified from modifications.xml):
#   - Channel 126:  no N/C suffix  → "TMT10plex-Lys126"  / "TMT10plex-Nter126"
#   - Channels 127–130: N/C suffix → "TMT10plex-Lys127N" / "TMT10plex-Nter127N" etc.
#   - Channel 131N: no N/C suffix  → "TMT10plex-Lys131"  / "TMT10plex-Nter131"
#   - Channel 131C: TMT11plex      → "TMT11plex-Lys131C" / "TMT11plex-Nter131C"
# Tuple: (internalLabel for Lys, terminalLabel for N-term, channel_key)
TMT11_CHANNELS = [
    ("TMT10plex-Lys126",   "TMT10plex-Nter126",   "126C"),
    ("TMT10plex-Lys127N",  "TMT10plex-Nter127N",  "127N"),
    ("TMT10plex-Lys127C",  "TMT10plex-Nter127C",  "127C"),
    ("TMT10plex-Lys128N",  "TMT10plex-Nter128N",  "128N"),
    ("TMT10plex-Lys128C",  "TMT10plex-Nter128C",  "128C"),
    ("TMT10plex-Lys129N",  "TMT10plex-Nter129N",  "129N"),
    ("TMT10plex-Lys129C",  "TMT10plex-Nter129C",  "129C"),
    ("TMT10plex-Lys130N",  "TMT10plex-Nter130N",  "130N"),
    ("TMT10plex-Lys130C",  "TMT10plex-Nter130C",  "130C"),
    ("TMT10plex-Lys131",   "TMT10plex-Nter131",   "131N"),
    ("TMT11plex-Lys131C",  "TMT11plex-Nter131C",  "131C"),
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

    Structure validated against a working TMT11 GUI-generated mqpar (v2.4.7).
    Key choices:
      - minPeptides = 1         → single-peptide proteins reported (our mutants)
      - decoyMode = revert      → MaxQuant generates reversed decoys internally
      - isobaricLabels          → TMT11 channels using MaxQuant internal mod names
      - writeSdrf = False       → disable SDRF writer (crashes when channels empty)
    """
    root = ET.Element("MaxQuantParams",
                      attrib={
                          "xmlns:xsd": "http://www.w3.org/2001/XMLSchema",
                          "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
                      })

    # ── FASTA ─────────────────────────────────────────────────────────────────
    fastas = sub(root, "fastaFiles")
    fi = sub(fastas, "FastaFileInfo")
    sub(fi, "fastaFilePath",         fasta_path)
    sub(fi, "identifierParseRule",   r">([^\s]*)")
    sub(fi, "descriptionParseRule",  r">(.*)")
    sub(fi, "taxonomyParseRule",     "")
    sub(fi, "variationParseRule",    "")
    sub(fi, "modificationParseRule", "")
    sub(fi, "taxonomyId",            "")

    sub(root, "fixedSearchFolder",   "")

    # ── GLOBAL SETTINGS (matching reference mqpar field order) ────────────────
    sub(root, "andromedaCacheSize",  350000)
    sub(root, "advancedRatios",      "True")
    sub(root, "pvalThres",           0.005)
    sub(root, "rtShift",             "False")
    sub(root, "separateLfq",         "False")
    sub(root, "lfqStabilizeLargeRatios", "True")
    sub(root, "lfqRequireMsms",      "True")
    sub(root, "lfqBayesQuant",       "False")
    sub(root, "decoyMode",           "revert")
    sub(root, "includeContaminants", "True")
    sub(root, "maxPeptideMass",      4600)
    sub(root, "minDeltaScoreUnmodifiedPeptides", 0)
    sub(root, "minDeltaScoreModifiedPeptides",   6)
    sub(root, "minScoreUnmodifiedPeptides",      0)
    sub(root, "minScoreModifiedPeptides",        40)
    sub(root, "secondPeptide",       "True")
    sub(root, "matchBetweenRuns",    "False")
    sub(root, "matchUnidentifiedFeatures", "False")
    sub(root, "matchBetweenRunsFdr", "False")
    sub(root, "dependentPeptides",   "False")
    sub(root, "dependentPeptideFdr", 0)
    sub(root, "dependentPeptideMassBin", 0)
    sub(root, "dependentPeptidesBetweenRuns", "False")
    sub(root, "dependentPeptidesWithinExperiment", "False")
    sub(root, "dependentPeptidesWithinParameterGroup", "False")
    sub(root, "dependentPeptidesRestrictFractions", "False")
    sub(root, "dependentPeptidesFractionDifference", 0)
    sub(root, "ibaq",               "False")
    sub(root, "top3",               "False")
    sub(root, "independentEnzymes", "False")
    sub(root, "useDeltaScore",      "False")
    sub(root, "splitProteinGroupsByTaxonomy", "False")
    sub(root, "taxonomyLevel",      "Species")
    sub(root, "avalon",             "False")
    sub(root, "nModColumns",        3)
    sub(root, "ibaqLogFit",         "False")
    sub(root, "razorProteinFdr",    "True")
    sub(root, "deNovoSequencing",   "False")
    sub(root, "massDifferenceSearch", "False")
    sub(root, "minPepLen",          7)
    sub(root, "psmFdrCrosslink",    0.01)
    sub(root, "peptideFdr",         1)    # 100% at peptide level; protein FDR controls
    sub(root, "proteinFdr",         1)    # use minPeptides=1 filter instead
    sub(root, "siteFdr",            0.01)
    sub(root, "minPeptides",        1)    # ← key: allow single-peptide protein IDs
    sub(root, "minRazorPeptides",   1)
    sub(root, "minUniquePeptides",  0)
    sub(root, "useCounterparts",    "False")
    sub(root, "quantMode",          1)    # 1 = reporter ion (TMT)
    sub(root, "mainSearchMaxCombinations", 200)
    sub(root, "writeMsScansTable",  "True")
    sub(root, "writeMsmsScansTable","True")
    sub(root, "writeMs3ScansTable", "False")
    sub(root, "writeAllPeptidesTable", "True")
    sub(root, "writeMzTab",         "False")
    sub(root, "writeSdrf",          "False")   # ← CRITICAL: prevents SDRF crash
    sub(root, "useSeriesReporters", "False")
    sub(root, "numThreads",         n_threads)
    sub(root, "emailAddress",       "")
    sub(root, "smtpHost",           "")
    sub(root, "emailFromAddress",   "")
    sub(root, "fixedCombinedFolder",out_dir)
    sub(root, "fullMinMz",          -1.79589544172745e+308)
    sub(root, "fullMaxMz",           1.79589544172745e+308)
    sub(root, "sendEmail",          "False")
    sub(root, "ionCountIntensities","False")

    # ── RAW FILES ─────────────────────────────────────────────────────────────
    fp_el   = sub(root, "filePaths")
    exp_el  = sub(root, "experiments")
    frac_el = sub(root, "fractions")
    ptm_el  = sub(root, "ptms")
    pg_el   = sub(root, "paramGroupIndices")

    for raw in sorted(raw_files):
        sub(fp_el,   "string",  raw)
        sub(exp_el,  "string",  plex_id)
        sub(frac_el, "short",   32767)    # 32767 = no fractionation
        sub(ptm_el,  "boolean", "False")
        sub(pg_el,   "int",     0)

    ref_ch = sub(root, "referenceChannel")
    for _ in sorted(raw_files):
        sub(ref_ch, "string", "")

    # ── PARAMETER GROUP ───────────────────────────────────────────────────────
    pgs = sub(root, "parameterGroups")
    pg  = sub(pgs, "parameterGroup")

    sub(pg, "msInstrument",  0)    # 0 = Orbitrap
    sub(pg, "maxCharge",     7)
    sub(pg, "lcmsRunType",   "Reporter MS2")   # inside parameterGroup
    sub(pg, "maxMissedCleavages", 2)
    sub(pg, "multiplicity",  1)
    sub(pg, "complementaryReporterType", 0)
    sub(pg, "reporterNormalization",     0)

    # Fixed mods: only Carbamidomethyl — TMT is handled via isobaricLabels
    fixed = sub(pg, "fixedModifications")
    sub(fixed, "string", "Carbamidomethyl (C)")

    # Enzyme
    enzymes = sub(pg, "enzymes")
    sub(enzymes, "string", "Trypsin/P")
    sub(pg, "enzymesFirstSearch")
    sub(pg, "useVariableModificationsFirstSearch", "False")

    # Variable modifications
    var = sub(pg, "variableModifications")
    sub(var, "string", "Oxidation (M)")
    sub(var, "string", "Acetyl (Protein N-term)")

    # TMT11 isobaric labels — use MaxQuant's internal modification names.
    # These names must exactly match entries in MaxQuant's modifications.xml.
    # Channels 126–131N: TMT10plex prefix; channel 131C: TMT11plex prefix.
    iso_el = sub(pg, "isobaricLabels")
    for internal_lys, terminal_nter, ch_key in TMT11_CHANNELS:
        ili = sub(iso_el, "IsobaricLabelInfo")
        sub(ili, "internalLabel",      internal_lys)
        sub(ili, "terminalLabel",      terminal_nter)
        sub(ili, "correctionFactorM2", 0)
        sub(ili, "correctionFactorM1", 0)
        sub(ili, "correctionFactorP1", 0)
        sub(ili, "correctionFactorP2", 0)
        sub(ili, "tmtLike",            "True")

    sub(pg, "neucodeLabels")
    sub(pg, "variableModificationsFirstSearch")

    # Mass tolerances (inside parameterGroup per reference)
    sub(pg, "firstSearchTol",        20)   # MS1 first pass (ppm)
    sub(pg, "mainSearchTol",         4.5)  # MS1 main search (ppm)
    sub(pg, "searchTolInPpm",        "True")
    sub(pg, "isotopeMatchTol",       2)
    sub(pg, "isotopeMatchTolInPpm",  "True")
    sub(pg, "reporterMassTolerance", 0.003)  # 3 mDa reporter ion (high-res MS2)
    sub(pg, "reporterPif",           0.75)   # precursor ion fraction filter
    sub(pg, "filterPif",             "True")
    sub(pg, "isobaricSumOverWindow", "True")
    sub(pg, "isobaricWeightExponent","0.75")

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

        # Create symlinks for RAW files inside out_dir.
        # MaxQuant determines combined/ location from where the RAW files live
        # (it creates combined/ as a subdirectory of the RAW parent dir).
        # Symlinking RAW files into out_dir makes MaxQuant place
        # combined/ at MQ_SEARCH/results/{plex_id}/combined/ as desired.
        linked_raws = []
        for raw in found_raws:
            raw_name = os.path.basename(raw)
            link = os.path.join(out_dir, raw_name)
            if not os.path.islink(link) and not os.path.exists(link):
                os.symlink(raw, link)
            linked_raws.append(link)

        root = build_mqpar(
            raw_files=linked_raws,
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
