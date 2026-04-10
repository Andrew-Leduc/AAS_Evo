# -*- coding: utf-8 -*-
"""
generate_mqpar.py

Generate MaxQuant mqpar.xml parameter files for per-plex TMT-11 MS searches
by patching a working template mqpar (sample_mqpar.xml) with per-plex values.

All fixed parameters (TMT11 isobaric labels, search tolerances, instrument
settings, etc.) live in the template. Only the per-plex fields are modified:
  - fastaFilePath
  - filePaths  (RAW file paths, hardlinked into results/{plex_id}/)
  - experiments (file basenames without .raw)
  - fractions   (all 32767 — no fractionation)
  - ptms        (all False)
  - paramGroupIndices (all 0)
  - referenceChannel  (all empty — use raw TMT intensities)
  - numThreads / name

MaxQuant places combined/ in the directory containing the RAW files.
RAW files are hardlinked into results/{plex_id}/ so combined/ lands there.

Usage:
    python3 generate_mqpar.py \\
        --tmt-map /path/to/pdc_file_tmt_map.tsv \\
        --raw-dir /path/to/RAW/ \\
        --fasta-dir /path/to/FASTA/per_plex/ \\
        -o /path/to/MQ_SEARCH/ \\
        [--template scripts/ms_search/maxquant/sample_mqpar.xml] \\
        [--threads 16]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from xml.etree import ElementTree as ET

# Preserve namespace declarations when round-tripping through ElementTree
ET.register_namespace("xsd", "http://www.w3.org/2001/XMLSchema")
ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")


def load_tmt_map(path):
    """Return {run_metadata_id: set(file_name)} from the TMT map."""
    plex_to_files = defaultdict(set)
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            plex_to_files[row["run_metadata_id"].strip()].add(
                row["file_name"].strip()
            )
    return plex_to_files


def patch_template(template_path, raw_files, fasta_path, plex_id, n_threads):
    """
    Parse template_path, patch per-plex fields, return XML string.
    """
    tree = ET.parse(template_path)
    root = tree.getroot()

    # --- scalar fields --------------------------------------------------
    root.find(".//fastaFilePath").text = fasta_path
    root.find("name").text            = plex_id
    root.find("numThreads").text      = str(n_threads)

    # --- per-file array fields ------------------------------------------
    sorted_raws = sorted(raw_files)
    n           = len(sorted_raws)
    basenames   = [os.path.splitext(os.path.basename(r))[0] for r in sorted_raws]

    arrays = [
        ("filePaths",         "string",  sorted_raws),
        ("experiments",       "string",  basenames),
        ("fractions",         "short",   ["32767"] * n),
        ("ptms",              "boolean", ["False"]  * n),
        ("paramGroupIndices", "int",     ["0"]      * n),
        ("referenceChannel",  "string",  [""]       * n),
    ]
    for tag, child_tag, values in arrays:
        el = root.find(tag)
        if el is None:
            sys.exit(f"ERROR: <{tag}> not found in template — wrong template?")
        for child in list(el):
            el.remove(child)
        for val in values:
            ET.SubElement(el, child_tag).text = val

    return ET.tostring(root, encoding="unicode")


def main():
    script_dir   = os.path.dirname(os.path.abspath(__file__))
    default_tmpl = os.path.join(script_dir, "sample_mqpar.xml")

    parser = argparse.ArgumentParser(
        description="Generate per-plex MaxQuant mqpar.xml files from a working template."
    )
    parser.add_argument("--tmt-map",  required=True,
                        help="pdc_file_tmt_map.tsv")
    parser.add_argument("--raw-dir",  required=True,
                        help="Directory containing RAW files")
    parser.add_argument("--fasta-dir", required=True,
                        help="Directory with per-plex FASTA files")
    parser.add_argument("-o", "--output", required=True,
                        help="MQ_SEARCH output directory")
    parser.add_argument("--template", default=default_tmpl,
                        help=f"Template mqpar.xml (default: {default_tmpl})")
    parser.add_argument("--threads", type=int, default=16,
                        help="CPU threads per search (default: 16)")
    args = parser.parse_args()

    for path, label in [
        (args.tmt_map,   "TMT map"),
        (args.raw_dir,   "RAW directory"),
        (args.fasta_dir, "FASTA directory"),
        (args.template,  "template mqpar.xml"),
    ]:
        if not os.path.exists(path):
            sys.exit(f"Error: {label} not found: {path}")

    mqpar_dir      = os.path.join(args.output, "mqpar")
    plex_list_path = os.path.join(args.output, "plex_list.txt")
    os.makedirs(mqpar_dir, exist_ok=True)

    print(f"Template:        {args.template}")
    print(f"TMT map:         {args.tmt_map}")
    plex_to_files = load_tmt_map(args.tmt_map)
    print(f"  {len(plex_to_files)} plexes found")

    available_raws   = {f for f in os.listdir(args.raw_dir)
                        if f.lower().endswith(".raw")}
    available_fastas = {os.path.splitext(f)[0]
                        for f in os.listdir(args.fasta_dir)
                        if f.endswith(".fasta")}
    print(f"  {len(available_raws)} RAW files in {args.raw_dir}")
    print(f"  {len(available_fastas)} per-plex FASTAs in {args.fasta_dir}")

    plex_list = []
    n_skip_fasta = n_skip_raw = 0

    for plex_id in sorted(plex_to_files):
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

        # Hardlink RAW files into results/{plex_id}/.
        # MaxQuant creates combined/ next to the RAW files, so this puts
        # combined/ at MQ_SEARCH/results/{plex_id}/combined/ automatically.
        # Hardlinks (not symlinks) are required — MaxQuant's Thermo reader
        # silently reads 0 scans from symlink paths on Linux.
        linked_raws = []
        for raw in found_raws:
            link = os.path.join(out_dir, os.path.basename(raw))
            if not os.path.exists(link):
                try:
                    os.link(raw, link)
                except OSError:
                    os.symlink(raw, link)   # fallback for cross-device mounts
            linked_raws.append(link)

        xml_str    = patch_template(args.template, linked_raws,
                                    fasta_path, plex_id, args.threads)
        mqpar_path = os.path.join(mqpar_dir, f"{plex_id}.xml")
        with open(mqpar_path, "w") as f:
            f.write('<?xml version="1.0" encoding="utf-8"?>\n')
            f.write(xml_str)
            f.write("\n")

        plex_list.append(plex_id)
        print(f"  {plex_id}: {len(linked_raws)} RAW files → {mqpar_path}")

    with open(plex_list_path, "w") as f:
        f.write("\n".join(plex_list) + "\n")

    print(f"\nPlexes ready:           {len(plex_list)}")
    print(f"Skipped (no FASTA):     {n_skip_fasta}")
    print(f"Skipped (no RAW files): {n_skip_raw}")
    print(f"mqpar files:            {mqpar_dir}")
    print(f"Plex list:              {plex_list_path}")


if __name__ == "__main__":
    main()
