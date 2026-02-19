# -*- coding: utf-8 -*-
"""
combine_plex_fastas.py

Combine per-sample mutant tryptic peptide FASTAs by TMT plex to create custom
MS search databases. Each plex FASTA contains:
  - Full reference proteome
  - Deduplicated mutant peptides from all samples in that plex (with patient info)
  - Compensatory peptides (for mutations observed in that plex, with patient info)

Header format for mutant peptides:
    >mut|{accession}|{gene}|{swap}|genetic|{patient}|{sample_type}

Example:
    >mut|P04637|TP53|R273H|genetic|C3L-00001|tumor

Header format for compensatory peptides:
    >comp|{accession}|{gene}|{orig_swap}_{comp_swap}|predicted|{patient}|{sample_type}

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

    Returns dict: run_metadata_id -> list of (case_submitter_id, sample_type)
    (excludes reference/pooled channels)
    """
    plex_to_cases = defaultdict(list)

    with open(tmt_map_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            run_id = row["run_metadata_id"].strip()
            case_id = row["case_submitter_id"].strip()
            sample_type = row.get("sample_type", "").strip().lower()

            # Skip reference/pooled channels
            if case_id.lower() in ("ref", "reference", "pooled", "pool", ""):
                continue

            plex_to_cases[run_id].append((case_id, sample_type))

    # Deduplicate within each plex
    for run_id in plex_to_cases:
        plex_to_cases[run_id] = list(set(plex_to_cases[run_id]))

    return plex_to_cases


def load_gdc_meta(gdc_meta_path):
    """
    Load GDC metadata.

    Returns:
        case_to_uuids: dict of case_submitter_id -> list of (file_id, sample_type)
        uuid_to_info: dict of file_id -> {"case_id": str, "sample_type": str}
    """
    case_to_uuids = defaultdict(list)
    uuid_to_info = {}

    with open(gdc_meta_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            case_id = row["case_submitter_id"].strip()
            # Column may be "file_id" or "gdc_file_id" depending on metadata version
            file_id = row.get("gdc_file_id", row.get("file_id", "")).strip()
            sample_type = row.get("sample_type", "").strip()
            case_to_uuids[case_id].append((file_id, sample_type))
            uuid_to_info[file_id] = {
                "case_id": case_id,
                "sample_type": sample_type,
            }

    return case_to_uuids, uuid_to_info


def extract_mutation_identity(header):
    """
    Extract mutation identity from a mutant FASTA header.

    Header format: >mut|P04637|TP53|R273H|genetic
    Returns: (accession, gene, swap) tuple or None
    """
    parts = header.split("|")
    if len(parts) < 5:
        return None

    # parts: [">mut", accession, gene, swap, source]
    accession = parts[1]
    gene = parts[2]
    swap = parts[3]

    return (accession, gene, swap)


def extract_comp_original_mutation(header):
    """
    Extract original mutation from a compensatory FASTA header.

    Header format: >comp|P04637|TP53|R273H_G245S|predicted
    Returns: (accession, gene, orig_swap) tuple or None
    """
    parts = header.split("|")
    if len(parts) < 5:
        return None

    accession = parts[1]
    gene = parts[2]
    swap_combo = parts[3]  # e.g., "R273H_G245S"

    # Extract original swap (before _)
    if "_" in swap_combo:
        orig_swap = swap_combo.split("_")[0]
    else:
        orig_swap = swap_combo

    return (accession, gene, orig_swap)


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
            # For peptides, write on single line; for proteins, wrap at 60 chars
            if len(seq) <= 60:
                f.write(seq + "\n")
            else:
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Combine per-sample mutant tryptic peptide FASTAs by TMT plex."
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
                        help="Optional compensatory mutant FASTA to include")

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
    case_to_uuids, uuid_to_info = load_gdc_meta(args.gdc_meta)
    print(f"  {len(case_to_uuids)} cases with GDC data")

    # Index available mutant FASTAs
    available_fastas = set()
    for fname in os.listdir(args.sample_dir):
        if fname.endswith("_mutant.fasta"):
            uuid = fname.replace("_mutant.fasta", "")
            available_fastas.add(uuid)
    print(f"  {len(available_fastas)} per-sample mutant FASTAs available")

    # Index compensatory entries by their original mutation identity
    comp_by_orig_mut = defaultdict(list)
    if comp_entries:
        for header, seq in comp_entries:
            orig_mut = extract_comp_original_mutation(header)
            if orig_mut:
                comp_by_orig_mut[orig_mut].append((header, seq))
        print(f"  {len(comp_by_orig_mut)} unique original mutations "
              f"with compensatory entries")

    # Process each plex
    print(f"\nGenerating per-plex FASTAs...")
    summary_path = os.path.join(args.output, "plex_summary.tsv")

    with open(summary_path, "w") as summary_f:
        summary_f.write("plex_id\tcases_in_plex\tcases_with_gdc\t"
                        "samples_with_fasta\tmutant_peptides\t"
                        "unique_mutations\tcomp_peptides\n")

        total_plexes = 0
        total_comp_entries = 0

        for plex_id in sorted(plex_to_cases.keys()):
            case_sample_pairs = plex_to_cases[plex_id]

            # Collect all mutant entries for this plex
            plex_mutants = []
            seen_peptides = set()  # (accession, swap, peptide_seq) for dedup
            # Track which mutations come from which samples
            # (accession, gene, swap) -> list of (case_id, sample_type)
            mutation_to_samples = defaultdict(list)
            cases_with_gdc = 0
            samples_with_fasta = 0
            unique_cases = set(cp[0] for cp in case_sample_pairs)

            for case_id, pdc_sample_type in case_sample_pairs:
                uuid_pairs = case_to_uuids.get(case_id, [])
                if uuid_pairs:
                    cases_with_gdc += 1

                for uuid, gdc_sample_type in uuid_pairs:
                    if uuid not in available_fastas:
                        continue

                    samples_with_fasta += 1

                    fasta_path = os.path.join(args.sample_dir,
                                              f"{uuid}_mutant.fasta")
                    entries = load_mutant_fasta(fasta_path)

                    for header, seq in entries:
                        # Extract mutation identity
                        mut_id = extract_mutation_identity(header)
                        if not mut_id:
                            continue

                        accession, gene, swap = mut_id

                        # Track sample provenance
                        mutation_to_samples[mut_id].append(
                            (case_id, gdc_sample_type))

                        # Deduplicate by peptide sequence + mutation
                        pep_key = (accession, swap, seq)
                        if pep_key in seen_peptides:
                            continue
                        seen_peptides.add(pep_key)

                        # Build header with patient info
                        # Format: >mut|accession|gene|swap|genetic|patient|sample_type
                        # Normalize sample_type: FASTA headers treat spaces as
                        # description delimiters, which breaks Philosopher lookup.
                        st = gdc_sample_type.replace(" ", "_")
                        new_header = f">mut|{accession}|{gene}|{swap}|genetic|{case_id}|{st}"
                        plex_mutants.append((new_header, seq))

            # Filter compensatory entries to those whose original mutation
            # is present in this plex, and add patient info
            plex_comp_entries = []
            seen_comp_peptides = set()

            if comp_entries:
                for header, seq in comp_entries:
                    orig_mut = extract_comp_original_mutation(header)
                    if not orig_mut:
                        continue

                    # Check if original mutation is observed in this plex
                    if orig_mut not in mutation_to_samples:
                        continue

                    # Get patient info from samples with this mutation
                    samples = mutation_to_samples[orig_mut]
                    if not samples:
                        continue

                    # Use first sample's info (could also list all)
                    case_id, sample_type = samples[0]

                    # Parse original header to rebuild with patient info
                    parts = header.split("|")
                    if len(parts) < 5:
                        continue

                    accession = parts[1]
                    gene = parts[2]
                    swap_combo = parts[3]

                    # Deduplicate by peptide sequence + mutation combo
                    comp_pep_key = (accession, swap_combo, seq)
                    if comp_pep_key in seen_comp_peptides:
                        continue
                    seen_comp_peptides.add(comp_pep_key)

                    # Build header with patient info
                    # Format: >comp|accession|gene|swaps|predicted|patient|sample_type
                    st = sample_type.replace(" ", "_")
                    new_header = f">comp|{accession}|{gene}|{swap_combo}|predicted|{case_id}|{st}"
                    plex_comp_entries.append((new_header, seq))

            # Write plex FASTA: reference + mutants + compensatory
            out_path = os.path.join(args.output, f"{plex_id}.fasta")
            all_entries = ref_entries + plex_mutants + plex_comp_entries
            write_fasta(all_entries, out_path)

            summary_f.write(f"{plex_id}\t{len(unique_cases)}\t{cases_with_gdc}\t"
                            f"{samples_with_fasta}\t{len(plex_mutants)}\t"
                            f"{len(mutation_to_samples)}\t"
                            f"{len(plex_comp_entries)}\n")

            total_plexes += 1
            total_comp_entries += len(plex_comp_entries)

    print(f"\n{'='*50}")
    print(f"Plex FASTA Generation Summary")
    print(f"{'='*50}")
    print(f"Plexes generated:       {total_plexes}")
    print(f"Reference proteins:     {len(ref_entries)}")
    print(f"Compensatory peptides:  {total_comp_entries} (across all plexes)")
    print(f"Summary: {summary_path}")
    print(f"Output:  {args.output}")


if __name__ == "__main__":
    main()
