# -*- coding: utf-8 -*-
"""
filter_and_rank.py

Filter and rank missense mutations by pathogenicity for downstream analysis.

Reads the consolidated missense TSV (with AlphaMissense + gnomAD scores),
computes per-protein mutation burden per patient, ranks mutations by a
composite harmfulness score, and outputs the top N for MSA generation and
coevolution analysis.

Usage:
    python3 filter_and_rank.py \
        --vep-tsv /path/to/all_missense_mutations.tsv \
        --ref-fasta /path/to/uniprot_human_canonical.fasta \
        -o /path/to/ANALYSIS/

Optional:
    --top-n 5000                   Number of top mutations to keep
    --min-vaf-for-counting 0.3     VAF threshold for burden counting
    --max-gnomad-af 0.01           Exclude common polymorphisms
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}

HGVSP_PATTERN = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


# ---------------------------------------------------------------------------
# Helper functions (duplicated from coevolution_analysis.py per project convention)
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


def parse_mutation(hgvsp, amino_acids, protein_position):
    """Parse missense mutation from VEP output fields."""
    if hgvsp and hgvsp not in ("", "-", "NA"):
        m = HGVSP_PATTERN.search(hgvsp)
        if m:
            ref_1 = AA3_TO_1.get(m.group(1))
            alt_1 = AA3_TO_1.get(m.group(3))
            if ref_1 and alt_1:
                return ref_1, int(m.group(2)), alt_1

    if amino_acids and "/" in amino_acids and protein_position:
        parts = amino_acids.split("/")
        if len(parts) == 2 and len(parts[0]) == 1 and len(parts[1]) == 1:
            try:
                return parts[0], int(protein_position), parts[1]
            except ValueError:
                pass

    return None


def safe_float(val, default=None):
    """Parse float from string, returning default for missing/invalid."""
    if val in ("", "NA", "-", ".", None):
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def load_mutations(vep_tsv_path):
    """
    Load all missense mutation observations from consolidated VEP TSV.

    Returns list of dicts with parsed fields.
    """
    observations = []

    with open(vep_tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            symbol = row.get("SYMBOL", "").strip()
            if not symbol:
                continue

            hgvsp = row.get("HGVSp", "").strip()
            amino_acids = row.get("Amino_acids", "").strip()
            protein_position = row.get("Protein_position", "").strip()

            mutation = parse_mutation(hgvsp, amino_acids, protein_position)
            if mutation is None:
                continue

            ref_aa, pos, mut_aa = mutation
            vaf = safe_float(row.get("VAF", "NA"))
            gnomad_af = safe_float(row.get("gnomADe_AF", "NA"))
            am_path = safe_float(row.get("am_pathogenicity", "NA"))
            am_class = row.get("am_class", "").strip()
            sample_id = row.get("sample_id", "").strip()

            observations.append({
                "sample_id": sample_id,
                "symbol": symbol,
                "ref_aa": ref_aa,
                "pos": pos,
                "mut_aa": mut_aa,
                "hgvsp": hgvsp,
                "amino_acids": amino_acids,
                "chrom": row.get("CHROM", ""),
                "genomic_pos": row.get("POS", ""),
                "ref": row.get("REF", ""),
                "alt": row.get("ALT", ""),
                "gene_id": row.get("Gene", ""),
                "vaf": vaf,
                "gnomad_af": gnomad_af,
                "am_pathogenicity": am_path,
                "am_class": am_class,
                "dp": safe_float(row.get("DP", "NA")),
            })

    return observations


def compute_mutation_burden(observations, min_vaf):
    """
    Count mutations per protein per patient where VAF > min_vaf.

    Returns dict: (sample_id, symbol) -> count
    """
    burden = defaultdict(int)

    for obs in observations:
        if obs["vaf"] is not None and obs["vaf"] > min_vaf:
            burden[(obs["sample_id"], obs["symbol"])] += 1

    return burden


def aggregate_mutations(observations, burden, max_gnomad_af, min_vaf_for_ranking):
    """
    Aggregate observations by unique mutation (gene + position + AA change).

    Returns list of aggregated mutation dicts.
    """
    # Group by unique mutation
    groups = defaultdict(list)
    for obs in observations:
        key = (obs["symbol"], obs["pos"], obs["ref_aa"], obs["mut_aa"])
        groups[key].append(obs)

    aggregated = []
    for key, obs_list in groups.items():
        symbol, pos, ref_aa, mut_aa = key

        # Collect per-observation data
        vafs = [o["vaf"] for o in obs_list if o["vaf"] is not None]
        gnomad_values = [o["gnomad_af"] for o in obs_list if o["gnomad_af"] is not None]
        am_values = [o["am_pathogenicity"] for o in obs_list if o["am_pathogenicity"] is not None]

        gnomad_af = gnomad_values[0] if gnomad_values else None
        am_path = am_values[0] if am_values else None
        am_class = obs_list[0]["am_class"]

        # Pre-filter: exclude common polymorphisms
        if gnomad_af is not None and gnomad_af > max_gnomad_af:
            continue

        # Pre-filter: require at least one observation with VAF above threshold
        if min_vaf_for_ranking is not None:
            if not any(v >= min_vaf_for_ranking for v in vafs):
                continue

        n_samples = len(obs_list)
        mean_vaf = sum(vafs) / len(vafs) if vafs else 0.0
        max_vaf = max(vafs) if vafs else 0.0

        # Mean mutation burden for patients carrying this mutation
        burden_values = []
        for obs in obs_list:
            b = burden.get((obs["sample_id"], obs["symbol"]), 0)
            burden_values.append(b)
        mean_burden = sum(burden_values) / len(burden_values) if burden_values else 0.0

        # Pick representative HGVSp (prefer non-empty)
        hgvsp = ""
        for obs in obs_list:
            if obs["hgvsp"] and obs["hgvsp"] not in ("", "-", "NA"):
                hgvsp = obs["hgvsp"]
                break

        # Pick representative genomic coordinates
        representative = obs_list[0]

        aggregated.append({
            "symbol": symbol,
            "ref_aa": ref_aa,
            "pos": pos,
            "mut_aa": mut_aa,
            "mutation_short": f"{ref_aa}{pos}{mut_aa}",
            "hgvsp": hgvsp,
            "chrom": representative["chrom"],
            "genomic_pos": representative["genomic_pos"],
            "ref": representative["ref"],
            "alt": representative["alt"],
            "gene_id": representative["gene_id"],
            "n_samples": n_samples,
            "mean_vaf": mean_vaf,
            "max_vaf": max_vaf,
            "gnomad_af": gnomad_af,
            "am_pathogenicity": am_path,
            "am_class": am_class,
            "mean_burden": mean_burden,
        })

    return aggregated


def score_mutations(aggregated, weight_am, weight_gnomad, weight_recurrence,
                    weight_burden):
    """
    Compute composite harmfulness score for each aggregated mutation.

    Components:
        - AlphaMissense pathogenicity (0-1, missing = 0.5)
        - gnomAD rarity (1 - AF, missing = 1.0)
        - Sample recurrence (min(1, n_samples/20))
        - Protein mutation burden context (min(1, mean_burden/5))
    """
    for mut in aggregated:
        am = mut["am_pathogenicity"] if mut["am_pathogenicity"] is not None else 0.5
        gnomad = 1.0 - (mut["gnomad_af"] if mut["gnomad_af"] is not None else 0.0)
        recurrence = min(1.0, mut["n_samples"] / 20.0)
        burden_ctx = min(1.0, mut["mean_burden"] / 5.0)

        mut["score"] = (
            weight_am * am +
            weight_gnomad * gnomad +
            weight_recurrence * recurrence +
            weight_burden * burden_ctx
        )


def write_mutation_burden(burden, observations, output_path):
    """Write per-protein-per-patient mutation burden TSV."""
    # Collect mutation lists per (sample, gene)
    mut_lists = defaultdict(list)
    for obs in observations:
        if obs["vaf"] is not None and obs["vaf"] > 0.3:
            key = (obs["sample_id"], obs["symbol"])
            mut_lists[key].append(f"{obs['ref_aa']}{obs['pos']}{obs['mut_aa']}")

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample_id", "SYMBOL", "n_mutations_vaf_gt_0.3", "mutations_list"
        ])

        for (sample_id, symbol), count in sorted(burden.items()):
            muts = mut_lists.get((sample_id, symbol), [])
            writer.writerow([
                sample_id, symbol, count, ";".join(muts)
            ])


def write_ranked_mutations(aggregated, output_path):
    """Write ranked mutations TSV."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "rank", "SYMBOL", "mutation", "HGVSp",
            "CHROM", "POS", "REF", "ALT", "Gene",
            "n_samples", "mean_vaf", "max_vaf",
            "gnomADe_AF", "am_pathogenicity", "am_class",
            "mean_protein_burden", "composite_score",
        ])

        for i, mut in enumerate(aggregated, 1):
            writer.writerow([
                i,
                mut["symbol"],
                mut["mutation_short"],
                mut["hgvsp"],
                mut["chrom"],
                mut["genomic_pos"],
                mut["ref"],
                mut["alt"],
                mut["gene_id"],
                mut["n_samples"],
                f"{mut['mean_vaf']:.4f}" if mut["mean_vaf"] else "NA",
                f"{mut['max_vaf']:.4f}" if mut["max_vaf"] else "NA",
                f"{mut['gnomad_af']:.6f}" if mut["gnomad_af"] is not None else "NA",
                f"{mut['am_pathogenicity']:.4f}" if mut["am_pathogenicity"] is not None else "NA",
                mut["am_class"] or "NA",
                f"{mut['mean_burden']:.2f}",
                f"{mut['score']:.4f}",
            ])


def write_summary(aggregated, top_n, total_obs, n_unique_before_filter,
                   n_genes, output_path):
    """Write human-readable summary."""
    with open(output_path, "w") as f:
        f.write("=== Filter and Rank Summary ===\n\n")
        f.write(f"Total mutation observations:  {total_obs}\n")
        f.write(f"Unique mutations (pre-filter): {n_unique_before_filter}\n")
        f.write(f"Unique mutations (post-filter): {len(aggregated)}\n")
        f.write(f"Top N selected:              {min(top_n, len(aggregated))}\n")
        f.write(f"Unique genes in top N:       {n_genes}\n\n")

        if aggregated:
            scores = [m["score"] for m in aggregated]
            f.write(f"Score range: {min(scores):.4f} - {max(scores):.4f}\n")
            f.write(f"Score median: {sorted(scores)[len(scores)//2]:.4f}\n\n")

            am_count = sum(1 for m in aggregated if m["am_pathogenicity"] is not None)
            f.write(f"AlphaMissense coverage: {am_count}/{len(aggregated)} "
                    f"({100*am_count/len(aggregated):.1f}%)\n")

            gnomad_count = sum(1 for m in aggregated if m["gnomad_af"] is not None)
            f.write(f"gnomAD coverage:        {gnomad_count}/{len(aggregated)} "
                    f"({100*gnomad_count/len(aggregated):.1f}%)\n\n")

        # Score component weights
        f.write("Scoring weights:\n")
        f.write("  AlphaMissense pathogenicity: 0.50\n")
        f.write("  gnomAD rarity:               0.20\n")
        f.write("  Sample recurrence:           0.15\n")
        f.write("  Protein mutation burden:     0.15\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter and rank missense mutations by pathogenicity"
    )
    parser.add_argument("--vep-tsv", required=True,
                        help="Consolidated missense mutations TSV")
    parser.add_argument("--ref-fasta", required=True,
                        help="UniProt reference proteome FASTA")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory")
    parser.add_argument("--top-n", type=int, default=5000,
                        help="Number of top mutations to select (default: 5000)")
    parser.add_argument("--min-vaf-for-counting", type=float, default=0.3,
                        help="VAF threshold for mutation burden counting (default: 0.3)")
    parser.add_argument("--min-vaf-for-ranking", type=float, default=None,
                        help="Minimum VAF for including in ranking (default: no filter)")
    parser.add_argument("--max-gnomad-af", type=float, default=0.01,
                        help="Maximum gnomAD AF to include (default: 0.01)")
    parser.add_argument("--weight-am", type=float, default=0.5,
                        help="Weight for AlphaMissense score (default: 0.5)")
    parser.add_argument("--weight-gnomad", type=float, default=0.2,
                        help="Weight for gnomAD rarity (default: 0.2)")
    parser.add_argument("--weight-recurrence", type=float, default=0.15,
                        help="Weight for sample recurrence (default: 0.15)")
    parser.add_argument("--weight-burden", type=float, default=0.15,
                        help="Weight for protein mutation burden (default: 0.15)")

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Validate input
    if not os.path.isfile(args.vep_tsv):
        print(f"ERROR: VEP TSV not found: {args.vep_tsv}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.ref_fasta):
        print(f"ERROR: Reference FASTA not found: {args.ref_fasta}", file=sys.stderr)
        sys.exit(1)

    # Step 1: Load reference proteome (for gene validation)
    print("Loading reference proteome...")
    gene_to_protein, _, _ = parse_reference_fasta(args.ref_fasta)
    print(f"  {len(gene_to_protein)} genes in reference")

    # Step 2: Load all mutation observations
    print("Loading mutations...")
    observations = load_mutations(args.vep_tsv)
    total_obs = len(observations)
    print(f"  {total_obs} mutation observations loaded")

    # Step 3: Compute mutation burden per protein per patient
    print("Computing mutation burden per protein per patient...")
    burden = compute_mutation_burden(observations, args.min_vaf_for_counting)
    print(f"  {len(burden)} (sample, gene) pairs with mutations above "
          f"VAF {args.min_vaf_for_counting}")

    # Write burden table
    burden_path = os.path.join(args.output, "mutation_burden.tsv")
    write_mutation_burden(burden, observations, burden_path)
    print(f"  Written: {burden_path}")

    # Step 4: Aggregate and filter
    print("Aggregating mutations...")

    # Count unique mutations before filtering (for summary)
    unique_keys = set()
    for obs in observations:
        unique_keys.add((obs["symbol"], obs["pos"], obs["ref_aa"], obs["mut_aa"]))
    n_unique_before = len(unique_keys)
    print(f"  {n_unique_before} unique mutations (pre-filter)")

    aggregated = aggregate_mutations(
        observations, burden,
        max_gnomad_af=args.max_gnomad_af,
        min_vaf_for_ranking=args.min_vaf_for_ranking,
    )
    print(f"  {len(aggregated)} unique mutations (post-filter, "
          f"gnomAD AF <= {args.max_gnomad_af})")

    # Step 5: Score and rank
    print("Scoring mutations...")
    score_mutations(
        aggregated,
        weight_am=args.weight_am,
        weight_gnomad=args.weight_gnomad,
        weight_recurrence=args.weight_recurrence,
        weight_burden=args.weight_burden,
    )

    # Sort by score descending
    aggregated.sort(key=lambda m: m["score"], reverse=True)

    # Write all ranked mutations
    ranked_path = os.path.join(args.output, "ranked_mutations.tsv")
    write_ranked_mutations(aggregated, ranked_path)
    print(f"  Written: {ranked_path}")

    # Step 6: Select top N
    top_n = min(args.top_n, len(aggregated))
    top_mutations = aggregated[:top_n]
    print(f"Selected top {top_n} mutations")

    # Write top N
    top_path = os.path.join(args.output, f"top_{args.top_n}_mutations.tsv")
    write_ranked_mutations(top_mutations, top_path)
    print(f"  Written: {top_path}")

    # Step 7: Extract gene list for MSA generation
    top_genes = sorted(set(m["symbol"] for m in top_mutations))
    # Filter to genes present in reference proteome
    valid_genes = [g for g in top_genes if g in gene_to_protein]
    missing_genes = [g for g in top_genes if g not in gene_to_protein]

    gene_list_path = os.path.join(args.output, "gene_list_for_msa.txt")
    with open(gene_list_path, "w") as f:
        for gene in valid_genes:
            f.write(gene + "\n")
    print(f"  {len(valid_genes)} genes for MSA generation")
    if missing_genes:
        print(f"  {len(missing_genes)} genes not in reference "
              f"(skipped): {', '.join(missing_genes[:10])}"
              + ("..." if len(missing_genes) > 10 else ""))
    print(f"  Written: {gene_list_path}")

    # Step 8: Write summary
    summary_path = os.path.join(args.output, "ranking_summary.txt")
    write_summary(aggregated, top_n, total_obs, n_unique_before,
                  len(valid_genes), summary_path)
    print(f"  Written: {summary_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
