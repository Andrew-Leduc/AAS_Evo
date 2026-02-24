# -*- coding: utf-8 -*-
"""
generate_manifests.py

Generate FragPipe manifests and TMT annotations for per-plex MS searches.

Each TMT plex is searched independently against its custom FASTA database
(reference proteome + plex-specific mutant + compensatory entries).

Reads the TMT mapping file to group RAW files by plex, then generates:
  - Per-plex FragPipe manifest (.fp-manifest)
  - Per-plex TMT channel annotation
  - Per-plex workflow (from template, with FASTA path patched)
  - Plex list for SLURM array submission

Usage:
    python3 generate_manifests.py \
        --tmt-map /path/to/pdc_file_tmt_map.tsv \
        --raw-dir /path/to/RAW/ \
        --fasta-dir /path/to/FASTA/per_plex/ \
        --workflow-template /path/to/template.workflow \
        -o /path/to/MS_SEARCH/
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


# TMT channel name mapping: pdc_file_tmt_map format -> FragPipe format
TMT_CHANNEL_MAP = {
    "tmt_126": "126",
    "tmt_127n": "127N",
    "tmt_127c": "127C",
    "tmt_128n": "128N",
    "tmt_128c": "128C",
    "tmt_129n": "129N",
    "tmt_129c": "129C",
    "tmt_130n": "130N",
    "tmt_130c": "130C",
    "tmt_131": "131N",
    "tmt_131c": "131C",
}


def load_tmt_map(tmt_map_path):
    """
    Load TMT mapping file.

    Returns:
        plex_to_files: dict of run_metadata_id -> set of file_names
        plex_to_channels: dict of run_metadata_id -> list of
            {"file_name", "channel", "case_id", "sample_type"}
    """
    plex_to_files = defaultdict(set)
    plex_to_channels = defaultdict(list)

    with open(tmt_map_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            run_id = row["run_metadata_id"].strip()
            file_name = row["file_name"].strip()
            channel = row.get("tmt_channel", "").strip().lower()
            case_id = row.get("case_submitter_id", "").strip()
            sample_type = row.get("sample_type", "").strip()

            plex_to_files[run_id].add(file_name)
            plex_to_channels[run_id].append({
                "file_name": file_name,
                "channel": channel,
                "case_id": case_id,
                "sample_type": sample_type,
            })

    return plex_to_files, plex_to_channels


def get_unique_channel_annotation(channel_entries):
    """
    Collapse channel entries to unique (channel -> sample) mapping.

    Multiple file_names (fractions) map to the same TMT channels, so
    deduplicate to one row per channel for the annotation file.
    """
    channel_to_sample = {}
    for entry in channel_entries:
        ch = entry["channel"]
        if ch not in channel_to_sample:
            case_id = entry["case_id"]
            sample_type = entry["sample_type"]
            if case_id.lower() in ("ref", "reference", "pooled", "pool", ""):
                label = "Reference"
                condition = "reference"
            else:
                label = f"{case_id}_{sample_type}" if sample_type else case_id
                condition = sample_type if sample_type else "unknown"
            channel_to_sample[ch] = {
                "sample": label,
                "condition": condition,
            }
    return channel_to_sample


def write_manifest(raw_files, plex_id, out_path):
    """Write FragPipe manifest file."""
    with open(out_path, "w") as f:
        for raw_path in sorted(raw_files):
            # experiment = plex_id, bioreplicate = 1, data_type = DDA
            f.write(f"{raw_path}\t{plex_id}\t1\tDDA\n")


def write_annotation(channel_map, plex_id, out_path):
    """
    Write TMT channel annotation file.

    Format: sample  channel  condition
    This can be loaded into FragPipe's TMT-Integrator configuration.
    """
    with open(out_path, "w") as f:
        f.write("sample\tplex\tchannel\tcondition\n")
        for channel_key in sorted(channel_map.keys(),
                                   key=lambda c: list(TMT_CHANNEL_MAP.keys()
                                                      ).index(c)
                                   if c in TMT_CHANNEL_MAP else 999):
            info = channel_map[channel_key]
            fp_channel = TMT_CHANNEL_MAP.get(channel_key, channel_key)
            f.write(f"{info['sample']}\t{plex_id}\t{fp_channel}\t"
                    f"{info['condition']}\n")


def patch_workflow(template_path, fasta_path, out_path):
    """
    Copy workflow template with the database path replaced.

    FragPipe workflow files use 'database.db-path=...' for the FASTA path.
    """
    with open(template_path) as f:
        content = f.read()

    # Replace or add database path
    new_line = f"database.db-path={fasta_path}"
    new_content = re.sub(
        r"^database\.db-path=.*$",
        new_line,
        content,
        flags=re.MULTILINE,
    )
    if new_content == content:
        # Key didn't exist in template — append it
        new_content = content.rstrip("\n") + f"\n{new_line}\n"
    content = new_content

    with open(out_path, "w") as f:
        f.write(content)


def main():
    parser = argparse.ArgumentParser(
        description="Generate FragPipe manifests for per-plex MS searches."
    )
    parser.add_argument("--tmt-map", required=True,
                        help="TMT mapping file (pdc_file_tmt_map.tsv)")
    parser.add_argument("--raw-dir", required=True,
                        help="Directory containing flattened RAW files")
    parser.add_argument("--fasta-dir", required=True,
                        help="Directory containing per-plex FASTAs")
    parser.add_argument("--workflow-template", default=None,
                        help="FragPipe workflow template to patch per plex")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for manifests and configs")

    args = parser.parse_args()

    for path, label in [(args.tmt_map, "TMT map"),
                        (args.raw_dir, "RAW directory"),
                        (args.fasta_dir, "FASTA directory")]:
        if not os.path.exists(path):
            sys.exit(f"Error: {label} not found: {path}")

    if args.workflow_template and not os.path.isfile(args.workflow_template):
        sys.exit(f"Error: Workflow template not found: "
                 f"{args.workflow_template}")

    # Create output directories
    manifest_dir = os.path.join(args.output, "manifests")
    annotation_dir = os.path.join(args.output, "annotations")
    workflow_dir = os.path.join(args.output, "workflows")
    results_dir = os.path.join(args.output, "results")
    for d in [manifest_dir, annotation_dir, results_dir]:
        os.makedirs(d, exist_ok=True)
    if args.workflow_template:
        os.makedirs(workflow_dir, exist_ok=True)

    # Load TMT mapping
    print(f"Loading TMT map: {args.tmt_map}")
    plex_to_files, plex_to_channels = load_tmt_map(args.tmt_map)
    print(f"  {len(plex_to_files)} plexes")

    # Scan for available RAW files
    available_raws = set()
    for fname in os.listdir(args.raw_dir):
        if fname.lower().endswith(".raw"):
            available_raws.add(fname)
    print(f"  {len(available_raws)} RAW files in {args.raw_dir}")

    # Scan for available per-plex FASTAs
    available_fastas = set()
    for fname in os.listdir(args.fasta_dir):
        if fname.endswith(".fasta"):
            plex_id = fname.replace(".fasta", "")
            available_fastas.add(plex_id)
    print(f"  {len(available_fastas)} per-plex FASTAs available")

    # Process each plex
    print(f"\nGenerating manifests...")

    plex_list = []
    summary_rows = []
    total_skipped_no_fasta = 0
    total_skipped_no_raws = 0

    for plex_id in sorted(plex_to_files.keys()):
        expected_files = plex_to_files[plex_id]

        # Check FASTA exists
        fasta_path = os.path.join(args.fasta_dir, f"{plex_id}.fasta")
        if plex_id not in available_fastas:
            total_skipped_no_fasta += 1
            continue

        # Find available RAW files for this plex
        found_raws = []
        missing_raws = []
        for fname in sorted(expected_files):
            raw_path = os.path.join(args.raw_dir, fname)
            if fname in available_raws:
                found_raws.append(raw_path)
            else:
                missing_raws.append(fname)

        if not found_raws:
            total_skipped_no_raws += 1
            continue

        # Write manifest
        manifest_path = os.path.join(manifest_dir,
                                     f"{plex_id}.fp-manifest")
        write_manifest(found_raws, plex_id, manifest_path)

        # Write annotation
        channel_map = get_unique_channel_annotation(
            plex_to_channels[plex_id])
        annotation_path = os.path.join(annotation_dir,
                                       f"{plex_id}_annotation.tsv")
        write_annotation(channel_map, plex_id, annotation_path)

        # Patch workflow if template provided
        if args.workflow_template:
            workflow_path = os.path.join(workflow_dir,
                                         f"{plex_id}.workflow")
            patch_workflow(args.workflow_template,
                          os.path.abspath(fasta_path),
                          workflow_path)

        # Create per-plex results directory
        os.makedirs(os.path.join(results_dir, plex_id), exist_ok=True)

        plex_list.append(plex_id)
        summary_rows.append({
            "plex_id": plex_id,
            "expected_raws": len(expected_files),
            "found_raws": len(found_raws),
            "missing_raws": len(missing_raws),
            "n_channels": len(channel_map),
            "fasta": os.path.basename(fasta_path),
        })

    # Write plex list for SLURM array
    plex_list_path = os.path.join(args.output, "plex_list.txt")
    with open(plex_list_path, "w") as f:
        for plex_id in plex_list:
            f.write(plex_id + "\n")

    # Write summary
    summary_path = os.path.join(args.output, "search_setup_summary.tsv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t",
                                fieldnames=["plex_id", "expected_raws",
                                            "found_raws", "missing_raws",
                                            "n_channels", "fasta"])
        writer.writeheader()
        writer.writerows(summary_rows)

    # Print summary
    total_raws = sum(r["found_raws"] for r in summary_rows)
    print(f"\n{'=' * 50}")
    print("MS Search Setup Summary")
    print(f"{'=' * 50}")
    print(f"Plexes ready:           {len(plex_list)}")
    print(f"Plexes skipped (no FASTA): {total_skipped_no_fasta}")
    print(f"Plexes skipped (no RAWs):  {total_skipped_no_raws}")
    print(f"Total RAW files:        {total_raws}")
    print(f"Plex list:              {plex_list_path}")
    print(f"Manifests:              {manifest_dir}")
    print(f"Annotations:            {annotation_dir}")
    if args.workflow_template:
        print(f"Workflows:              {workflow_dir}")
    print(f"Results dir:            {results_dir}")
    print(f"Summary:                {summary_path}")

    if any(r["missing_raws"] > 0 for r in summary_rows):
        print(f"\nWARNING: Some plexes have missing RAW files.")
        print("  These plexes will be searched with available fractions only.")
        for r in summary_rows:
            if r["missing_raws"] > 0:
                print(f"  {r['plex_id']}: {r['missing_raws']} missing "
                      f"of {r['expected_raws']} expected")


if __name__ == "__main__":
    main()
