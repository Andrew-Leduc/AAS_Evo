# -*- coding: utf-8 -*-
"""
generate_compensatory_fastas.py

Generate tryptic peptide FASTA entries for predicted compensatory mutations.

For each compensatory prediction from coevolution analysis, applies both the
original destabilizing mutation AND the predicted compensatory substitution,
then extracts the tryptic peptide(s) containing the mutations.

Header format:
    >comp|{accession}|{gene}|{orig_swap}_{comp_swap}|predicted

Example:
    >comp|P04637|TP53|R273H_G245S|predicted
    SVTCTYSPALNKMFCQLAK

Usage:
    python3 generate_compensatory_fastas.py \
        --predictions /path/to/compensatory_predictions.tsv \
        --ref-fasta /path/to/uniprot_human_canonical.fasta \
        -o /path/to/FASTA/compensatory/

Optional:
    --ranked-mutations /path/to/top_5000_mutations.tsv
    --min-coupling-score 0.0
    --min-preference-shift 0.0
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def parse_reference_fasta(fasta_path):
    """Parse UniProt FASTA. Returns gene_to_protein and accession/entry maps."""
    gene_to_protein = {}
    accession_to_gene = {}
    entry_to_gene = {}

    current_acc = None
    current_entry = None
    current_gene = None
    lines = []

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_acc is not None:
                    seq = "".join(lines)
                    if current_gene and current_gene not in gene_to_protein:
                        gene_to_protein[current_gene] = (
                            current_acc, current_entry, seq
                        )
                    if current_gene:
                        accession_to_gene[current_acc] = current_gene
                        entry_to_gene[current_entry] = current_gene

                parts = line.split("|")
                if len(parts) >= 3:
                    current_acc = parts[1]
                    current_entry = parts[2].split()[0]
                else:
                    current_acc = line[1:].split()[0]
                    current_entry = current_acc

                gn_match = re.search(r"GN=(\S+)", line)
                current_gene = gn_match.group(1) if gn_match else None
                lines = []
            else:
                lines.append(line.strip())

    if current_acc is not None:
        seq = "".join(lines)
        if current_gene and current_gene not in gene_to_protein:
            gene_to_protein[current_gene] = (current_acc, current_entry, seq)
        if current_gene:
            accession_to_gene[current_acc] = current_gene
            entry_to_gene[current_entry] = current_gene

    return gene_to_protein, accession_to_gene, entry_to_gene


def safe_float(val, default=None):
    """Parse float from string, returning default for missing/invalid."""
    if val in ("", "NA", "-", ".", None):
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def find_tryptic_cleavage_sites(sequence):
    """
    Find tryptic cleavage sites in a protein sequence.

    Trypsin cleaves after K or R, unless followed by P.

    Returns list of cleavage positions (0-indexed, position AFTER which cleavage occurs).
    Includes 0 (start) and len(sequence) (end) as boundaries.
    """
    sites = [0]  # Start of protein
    for i, aa in enumerate(sequence):
        if aa in ('K', 'R'):
            if i + 1 < len(sequence) and sequence[i + 1] == 'P':
                continue  # No cleavage before proline
            sites.append(i + 1)
    sites.append(len(sequence))
    return sites


def extract_tryptic_peptides(sequence, mutation_positions, max_missed_cleavages=0):
    """
    Extract tryptic peptide(s) containing any of the mutation positions.

    Args:
        sequence: Full protein sequence (with mutations already applied)
        mutation_positions: List of 1-based positions
        max_missed_cleavages: Include peptides with up to this many missed cleavages

    Returns:
        List of unique (peptide_seq, positions_in_peptide) tuples
    """
    sites = find_tryptic_cleavage_sites(sequence)
    mut_indices = [p - 1 for p in mutation_positions]  # Convert to 0-based

    peptides = {}  # peptide_seq -> set of mutation positions it contains

    for i in range(len(sites) - 1):
        for missed in range(max_missed_cleavages + 1):
            if i + missed + 1 >= len(sites):
                break

            start = sites[i]
            end = sites[i + missed + 1]

            # Check which mutations this peptide contains
            contained_muts = []
            for mut_idx, mut_pos in zip(mut_indices, mutation_positions):
                if start <= mut_idx < end:
                    contained_muts.append(mut_pos)

            if contained_muts:
                peptide = sequence[start:end]
                # Filter out very short or very long peptides
                if 6 <= len(peptide) <= 50:
                    key = peptide
                    if key not in peptides:
                        peptides[key] = set()
                    peptides[key].update(contained_muts)

    return [(pep, sorted(pos)) for pep, pos in peptides.items()]


# ---------------------------------------------------------------------------
# Prediction loading
# ---------------------------------------------------------------------------

MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def load_predictions(predictions_path):
    """
    Load compensatory predictions TSV from coevolution_analysis.py.

    Returns list of dicts with parsed fields.
    """
    predictions = []

    with open(predictions_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Parse original mutation (e.g., "R273H")
            mut_str = row.get("mutation", "").strip()
            m = MUTATION_PATTERN.match(mut_str)
            if not m:
                continue

            orig_ref = m.group(1)
            orig_pos = int(m.group(2))
            orig_alt = m.group(3)

            comp_pos = row.get("covarying_pos", "").strip()
            try:
                comp_pos = int(comp_pos)
            except ValueError:
                continue

            predictions.append({
                "gene": row.get("gene", "").strip(),
                "uniprot_accession": row.get("uniprot_accession", "").strip(),
                "orig_ref": orig_ref,
                "orig_pos": orig_pos,
                "orig_alt": orig_alt,
                "orig_mutation": mut_str,
                "comp_pos": comp_pos,
                "comp_wt": row.get("wildtype_aa", "").strip(),
                "comp_alt": row.get("predicted_compensatory_aa", "").strip(),
                "coupling_score": safe_float(row.get("coupling_score", "")),
                "conditional_score": safe_float(row.get("conditional_score", "")),
                "preference_shift": safe_float(row.get("preference_shift", "")),
                "n_samples": row.get("n_samples", "0").strip(),
                "neff": row.get("neff", "").strip(),
                "hgvsp": row.get("mutation_hgvsp", "").strip(),
            })

    return predictions


def load_ranked_mutations(ranked_path):
    """
    Load ranked mutations to filter predictions.

    Returns set of (gene, mutation_short) tuples.
    """
    ranked = set()

    with open(ranked_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            symbol = row.get("SYMBOL", "").strip()
            mutation = row.get("mutation", "").strip()
            if symbol and mutation:
                ranked.add((symbol, mutation))

    return ranked


# ---------------------------------------------------------------------------
# FASTA generation
# ---------------------------------------------------------------------------

def apply_mutations(seq, mutations):
    """
    Apply a list of (position_1based, alt_aa) mutations to a sequence.

    Returns mutated sequence or None if any position is out of range.
    """
    seq_list = list(seq)
    for pos, alt_aa in mutations:
        idx = pos - 1  # Convert to 0-based
        if idx < 0 or idx >= len(seq_list):
            return None
        seq_list[idx] = alt_aa
    return "".join(seq_list)


def main():
    parser = argparse.ArgumentParser(
        description="Generate tryptic peptide FASTAs with compensatory mutations"
    )
    parser.add_argument("--predictions", required=True,
                        help="Compensatory predictions TSV from coevolution analysis")
    parser.add_argument("--ref-fasta", required=True,
                        help="UniProt reference proteome FASTA")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory")
    parser.add_argument("--ranked-mutations", default=None,
                        help="Optional: limit to mutations in this ranked list")
    parser.add_argument("--min-coupling-score", type=float, default=0.0,
                        help="Minimum coupling score to include (default: 0.0)")
    parser.add_argument("--min-preference-shift", type=float, default=0.0,
                        help="Minimum preference shift to include (default: 0.0)")

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Validate inputs
    if not os.path.isfile(args.predictions):
        print(f"ERROR: Predictions file not found: {args.predictions}",
              file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.ref_fasta):
        print(f"ERROR: Reference FASTA not found: {args.ref_fasta}",
              file=sys.stderr)
        sys.exit(1)

    # Load reference proteome
    print("Loading reference proteome...")
    gene_to_protein, _, _ = parse_reference_fasta(args.ref_fasta)
    print(f"  {len(gene_to_protein)} genes in reference")

    # Load predictions
    print("Loading compensatory predictions...")
    predictions = load_predictions(args.predictions)
    print(f"  {len(predictions)} predictions loaded")

    # Filter by coupling score and preference shift
    if args.min_coupling_score > 0 or args.min_preference_shift > 0:
        before = len(predictions)
        predictions = [
            p for p in predictions
            if (p["coupling_score"] is not None
                and p["coupling_score"] >= args.min_coupling_score)
            and (p["preference_shift"] is not None
                 and p["preference_shift"] >= args.min_preference_shift)
        ]
        print(f"  {before - len(predictions)} filtered by score thresholds")

    # Filter by ranked mutations if provided
    if args.ranked_mutations:
        print(f"Loading ranked mutations: {args.ranked_mutations}")
        ranked = load_ranked_mutations(args.ranked_mutations)
        print(f"  {len(ranked)} ranked mutations loaded")
        before = len(predictions)
        predictions = [
            p for p in predictions
            if (p["gene"], p["orig_mutation"]) in ranked
        ]
        print(f"  {before - len(predictions)} filtered (not in ranked list)")

    print(f"  {len(predictions)} predictions to process")

    # Group predictions by gene
    by_gene = defaultdict(list)
    for pred in predictions:
        by_gene[pred["gene"]].append(pred)

    # Generate FASTA entries
    stats = {
        "genes_processed": 0,
        "genes_skipped_no_ref": 0,
        "peptides_written": 0,
        "entries_skipped_pos_mismatch": 0,
        "entries_deduplicated": 0,
    }
    summary_rows = []

    consolidated_path = os.path.join(args.output, "all_compensatory.fasta")
    consolidated_f = open(consolidated_path, "w")

    for gene in sorted(by_gene.keys()):
        preds = by_gene[gene]

        if gene not in gene_to_protein:
            stats["genes_skipped_no_ref"] += 1
            continue

        accession, entry, ref_seq = gene_to_protein[gene]

        # Deduplicate: (orig_mutation, comp_pos, comp_alt) -> keep highest coupling score
        seen = {}
        for pred in preds:
            key = (pred["orig_mutation"], pred["comp_pos"], pred["comp_alt"])
            if key not in seen:
                seen[key] = pred
            else:
                existing = seen[key]
                if (pred["coupling_score"] or 0) > (existing["coupling_score"] or 0):
                    seen[key] = pred
                else:
                    stats["entries_deduplicated"] += 1

        gene_peptides = 0
        seen_peptides = set()  # Deduplicate peptide sequences

        gene_fasta_path = os.path.join(args.output, f"{gene}_compensatory.fasta")

        with open(gene_fasta_path, "w") as gene_f:
            for key, pred in sorted(seen.items()):
                orig_pos = pred["orig_pos"]
                comp_pos = pred["comp_pos"]
                comp_wt = pred["comp_wt"]
                comp_alt = pred["comp_alt"]
                orig_ref = pred["orig_ref"]
                orig_alt = pred["orig_alt"]

                # Validate: original ref AA matches reference sequence
                if orig_pos < 1 or orig_pos > len(ref_seq):
                    stats["entries_skipped_pos_mismatch"] += 1
                    continue
                if ref_seq[orig_pos - 1] != orig_ref:
                    stats["entries_skipped_pos_mismatch"] += 1
                    continue

                # Validate: compensatory wildtype AA matches reference
                if comp_pos < 1 or comp_pos > len(ref_seq):
                    stats["entries_skipped_pos_mismatch"] += 1
                    continue
                if ref_seq[comp_pos - 1] != comp_wt:
                    stats["entries_skipped_pos_mismatch"] += 1
                    continue

                # Skip if compensatory substitution is same as wildtype
                if comp_wt == comp_alt:
                    continue

                # Apply both mutations to reference sequence
                mutant_seq = apply_mutations(ref_seq, [
                    (orig_pos, orig_alt),
                    (comp_pos, comp_alt),
                ])
                if mutant_seq is None:
                    stats["entries_skipped_pos_mismatch"] += 1
                    continue

                # Extract tryptic peptide(s) containing either mutation
                peptides = extract_tryptic_peptides(
                    mutant_seq,
                    [orig_pos, comp_pos],
                    max_missed_cleavages=0
                )

                if not peptides:
                    continue

                # Header format: >comp|P04637|TP53|R273H_G245S|predicted
                orig_swap = f"{orig_ref}{orig_pos}{orig_alt}"
                comp_swap = f"{comp_wt}{comp_pos}{comp_alt}"

                for peptide, positions in peptides:
                    # Deduplicate peptide sequences per gene
                    pep_key = (accession, orig_swap, comp_swap, peptide)
                    if pep_key in seen_peptides:
                        continue
                    seen_peptides.add(pep_key)

                    header = f">comp|{accession}|{gene}|{orig_swap}_{comp_swap}|predicted"

                    gene_f.write(header + "\n")
                    gene_f.write(peptide + "\n")
                    consolidated_f.write(header + "\n")
                    consolidated_f.write(peptide + "\n")
                    gene_peptides += 1

        if gene_peptides == 0:
            # Remove empty per-gene FASTA
            os.remove(gene_fasta_path)
        else:
            stats["genes_processed"] += 1
            stats["peptides_written"] += gene_peptides
            summary_rows.append({
                "gene": gene,
                "accession": accession,
                "n_peptides": gene_peptides,
                "protein_length": len(ref_seq),
            })

    consolidated_f.close()

    # Write summary TSV
    summary_path = os.path.join(args.output, "compensatory_summary.tsv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene", "accession", "n_peptides", "protein_length"])
        for row in summary_rows:
            writer.writerow([
                row["gene"], row["accession"],
                row["n_peptides"], row["protein_length"],
            ])

    # Print summary
    print(f"\n{'=' * 50}")
    print("Compensatory Tryptic Peptide FASTA Generation Summary")
    print(f"{'=' * 50}")
    print(f"Genes processed:         {stats['genes_processed']}")
    print(f"Genes skipped (no ref):  {stats['genes_skipped_no_ref']}")
    print(f"Peptides written:        {stats['peptides_written']}")
    print(f"Entries deduplicated:    {stats['entries_deduplicated']}")
    print(f"Entries skipped (pos):   {stats['entries_skipped_pos_mismatch']}")
    print(f"Output:                  {consolidated_path}")
    print(f"Summary:                 {summary_path}")


if __name__ == "__main__":
    main()
