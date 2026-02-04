# -*- coding: utf-8 -*-
"""
combine_plex_fastas.py

Combine per-sample mutant FASTAs by TMT plex to create custom MS search
databases. Each plex FASTA contains the full reference proteome plus
deduplicated mutant entries from all samples in that plex.

Linking chain:
    TMT plex (run_metadata_id in pdc_file_tmt_map.tsv)
      -> case_submitter_id (skip "Ref" channels)
      -> GDC UUID (via gdc_meta_matched.tsv)
      -> per-sample mutant FASTA

Usage:
    python3 combine_plex_fastas.py \
        --ref-fasta /path/to/uniprot_human_canonical.fasta \
        --sample-dir /path/to/FASTA/per_sample/ \
        --tmt-map /path/to/pdc_file_tmt_map.tsv \
        --gdc-meta /path/to/gdc_meta_matched.tsv \
        -o /path/to/FASTA/per_plex/
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def load_tmt_map(tmt_map_path):
    """
    Load TMT mapping file.

    Returns dict: run_metadata_id -> list of case_submitter_ids
    (excludes reference/pooled channels)
    """
    plex_to_cases = defaultdict(set)

    with open(tmt_map_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            run_id = row["run_metadata_id"].strip()
            case_id = row["case_submitter_id"].strip()
            sample_type = row.get("sample_type", "").strip().lower()

            # Skip reference/pooled channels
            if case_id.lower() in ("ref", "reference", "pooled", "pool", ""):
                continue
            # Some TMT designs use "Internal Reference" or similar
            if "ref" in sample_type and "reference" not in case_id.lower():
                pass  # keep actual patient samples even if type mentions ref

            plex_to_cases[run_id].add(case_id)

    return plex_to_cases


def load_gdc_meta(gdc_meta_path):
    """
    Load GDC metadata to map case_submitter_id -> file_id (UUID).

    Returns dict: case_submitter_id -> list of file_ids
    (one case may have multiple BAMs: tumor + normal)
    """
    case_to_uuids = defaultdict(list)

    with open(gdc_meta_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            case_id = row["case_submitter_id"].strip()
            file_id = row["file_id"].strip()
            case_to_uuids[case_id].append(file_id)

    return case_to_uuids


def load_mutant_fasta(fasta_path):
    """
    Load a mutant FASTA file. Returns list of (header, sequence) tuples.
    """
    entries = []
    current_header = None
    lines = []

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    entries.append((current_header, "".join(lines)))
                current_header = line
                lines = []
            else:
                lines.append(line.strip())

    if current_header is not None:
        entries.append((current_header, "".join(lines)))

    return entries


def write_fasta(entries, out_path):
    """Write list of (header, sequence) tuples to FASTA file."""
    with open(out_path, "w") as f:
        for header, seq in entries:
            f.write(header + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Combine per-sample mutant FASTAs by TMT plex."
    )
    parser.add_argument("--ref-fasta", required=True,
                        help="UniProt reference proteome FASTA")
    parser.add_argument("--sample-dir", required=True,
                        help="Directory containing per-sample *_mutant.fasta files")
    parser.add_argument("--tmt-map", required=True,
                        help="TMT mapping file (pdc_file_tmt_map.tsv)")
    parser.add_argument("--gdc-meta", required=True,
                        help="GDC matched metadata (gdc_meta_matched.tsv)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for per-plex FASTAs")
    parser.add_argument("--compensatory-fasta", default=None,
                        help="Optional compensatory mutant FASTA to append to all plexes")

    args = parser.parse_args()

    for path, label in [(args.ref_fasta, "Reference FASTA"),
                        (args.sample_dir, "Sample directory"),
                        (args.tmt_map, "TMT map"),
                        (args.gdc_meta, "GDC metadata")]:
        if not os.path.exists(path):
            sys.exit(f"Error: {label} not found: {path}")

    os.makedirs(args.output, exist_ok=True)

    # Load reference proteome as raw text entries
    print(f"Loading reference proteome: {args.ref_fasta}")
    ref_entries = load_mutant_fasta(args.ref_fasta)
    print(f"  {len(ref_entries)} reference proteins")

    # Load compensatory entries if provided
    comp_entries = []
    if args.compensatory_fasta and os.path.isfile(args.compensatory_fasta):
        print(f"Loading compensatory FASTAs: {args.compensatory_fasta}")
        comp_entries = load_mutant_fasta(args.compensatory_fasta)
        print(f"  {len(comp_entries)} compensatory entries")

    # Load TMT mapping
    print(f"Loading TMT map: {args.tmt_map}")
    plex_to_cases = load_tmt_map(args.tmt_map)
    print(f"  {len(plex_to_cases)} plexes")

    # Load GDC metadata
    print(f"Loading GDC metadata: {args.gdc_meta}")
    case_to_uuids = load_gdc_meta(args.gdc_meta)
    print(f"  {len(case_to_uuids)} cases with GDC data")

    # Index available mutant FASTAs
    available_fastas = set()
    for fname in os.listdir(args.sample_dir):
        if fname.endswith("_mutant.fasta"):
            uuid = fname.replace("_mutant.fasta", "")
            available_fastas.add(uuid)
    print(f"  {len(available_fastas)} per-sample mutant FASTAs available")

    # Process each plex
    print(f"\nGenerating per-plex FASTAs...")
    summary_path = os.path.join(args.output, "plex_summary.tsv")

    with open(summary_path, "w") as summary_f:
        summary_f.write("plex_id\tcases_in_plex\tcases_with_gdc\t"
                        "samples_with_fasta\tmutant_entries\tunique_mutations\n")

        total_plexes = 0
        for plex_id in sorted(plex_to_cases.keys()):
            cases = plex_to_cases[plex_id]

            # Collect all mutant entries for this plex
            plex_mutants = []
            seen_mutations = set()  # header-based dedup (mutation identity)
            cases_with_gdc = 0
            samples_with_fasta = 0

            for case_id in cases:
                uuids = case_to_uuids.get(case_id, [])
                if uuids:
                    cases_with_gdc += 1

                for uuid in uuids:
                    if uuid not in available_fastas:
                        continue

                    samples_with_fasta += 1
                    fasta_path = os.path.join(args.sample_dir,
                                              f"{uuid}_mutant.fasta")
                    entries = load_mutant_fasta(fasta_path)

                    for header, seq in entries:
                        # Extract mutation identity from header for dedup
                        # Header format: >mut|P04637|TP53_R273H OS=...
                        parts = header.split("|")
                        if len(parts) >= 3:
                            mut_id = parts[1] + "|" + parts[2].split()[0]
                        else:
                            mut_id = header

                        if mut_id not in seen_mutations:
                            seen_mutations.add(mut_id)
                            plex_mutants.append((header, seq))

            # Write plex FASTA: reference + mutants + compensatory
            out_path = os.path.join(args.output, f"{plex_id}.fasta")
            all_entries = ref_entries + plex_mutants + comp_entries
            write_fasta(all_entries, out_path)

            summary_f.write(f"{plex_id}\t{len(cases)}\t{cases_with_gdc}\t"
                            f"{samples_with_fasta}\t{len(plex_mutants)}\t"
                            f"{len(seen_mutations)}\n")

            total_plexes += 1

    print(f"\n{'='*50}")
    print(f"Plex FASTA Generation Summary")
    print(f"{'='*50}")
    print(f"Plexes generated: {total_plexes}")
    print(f"Reference proteins: {len(ref_entries)}")
    print(f"Summary: {summary_path}")
    print(f"Output:  {args.output}")


if __name__ == "__main__":
    main()
