#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
msa_quality.py

Compute quality statistics for an MSA (A3M format).

Metrics computed:
  - n_sequences: Total number of sequences
  - neff: Number of effective sequences (accounts for redundancy)
  - query_length: Length of the query sequence
  - avg_coverage: Average fraction of query covered by each sequence
  - avg_identity: Average pairwise identity to query
  - column_occupancy: Per-position fraction of non-gap residues
  - min/max/median column occupancy

Usage:
    python3 msa_quality.py input.a3m [--output stats.json]

If --output is not specified, prints JSON to stdout.
"""

import argparse
import json
import sys
from collections import Counter


def parse_a3m(filepath):
    """
    Parse A3M format MSA.

    A3M is like FASTA but lowercase letters indicate insertions
    relative to the query (first sequence). We ignore insertions
    for MSA statistics.

    Returns:
        query_seq: The first (query) sequence (uppercase only)
        aligned_seqs: List of aligned sequences (uppercase only, with gaps)
    """
    sequences = []
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                current_seq = []
            else:
                # Remove lowercase (insertions) and keep uppercase + gaps
                filtered = ''.join(c for c in line if c.isupper() or c == '-')
                current_seq.append(filtered)

        if current_seq:
            sequences.append(''.join(current_seq))

    if not sequences:
        return None, []

    return sequences[0], sequences[1:] if len(sequences) > 1 else []


def compute_sequence_weights(sequences, threshold=0.8):
    """
    Compute sequence weights for Neff calculation.

    Each sequence is weighted by 1 / (number of sequences with >= threshold identity).
    This down-weights redundant sequences.

    Returns list of weights (one per sequence).
    """
    n = len(sequences)
    if n == 0:
        return []

    seq_len = len(sequences[0])
    weights = [1.0] * n

    # For efficiency, compute pairwise identities
    # For large MSAs, we use a faster approximation
    if n > 1000:
        # Sample-based approximation for large MSAs
        import random
        sample_size = min(500, n)
        sample_indices = random.sample(range(n), sample_size)

        for i in range(n):
            count = 1  # self
            seq_i = sequences[i]
            for j in sample_indices:
                if i == j:
                    continue
                # Compute identity
                matches = sum(1 for a, b in zip(seq_i, sequences[j])
                             if a == b and a != '-')
                aligned = sum(1 for a, b in zip(seq_i, sequences[j])
                             if a != '-' or b != '-')
                if aligned > 0:
                    identity = matches / aligned
                    if identity >= threshold:
                        count += 1
            # Scale by sampling ratio
            weights[i] = 1.0 / (count * n / sample_size)
    else:
        # Exact computation for smaller MSAs
        for i in range(n):
            count = 0
            seq_i = sequences[i]
            for j in range(n):
                # Compute identity
                matches = sum(1 for a, b in zip(seq_i, sequences[j])
                             if a == b and a != '-')
                aligned = sum(1 for a, b in zip(seq_i, sequences[j])
                             if a != '-' or b != '-')
                if aligned > 0:
                    identity = matches / aligned
                    if identity >= threshold:
                        count += 1
            weights[i] = 1.0 / count if count > 0 else 0

    return weights


def compute_column_occupancy(query_seq, aligned_seqs):
    """
    Compute per-column occupancy (fraction of non-gap residues).

    Returns list of occupancy values (0-1) for each query position.
    """
    if not aligned_seqs:
        return [1.0] * len(query_seq)

    n_seqs = len(aligned_seqs) + 1  # Include query
    occupancy = []

    for i in range(len(query_seq)):
        # Query always contributes 1
        count = 1
        for seq in aligned_seqs:
            if i < len(seq) and seq[i] != '-':
                count += 1
        occupancy.append(count / n_seqs)

    return occupancy


def compute_identity_to_query(query_seq, aligned_seqs):
    """
    Compute pairwise identity of each sequence to the query.

    Returns list of identity values (0-1).
    """
    if not aligned_seqs:
        return []

    identities = []
    query_len = len(query_seq)

    for seq in aligned_seqs:
        matches = 0
        aligned = 0
        for i in range(min(query_len, len(seq))):
            if query_seq[i] != '-':
                aligned += 1
                if seq[i] == query_seq[i]:
                    matches += 1

        identity = matches / aligned if aligned > 0 else 0
        identities.append(identity)

    return identities


def compute_coverage(query_seq, aligned_seqs):
    """
    Compute coverage of query by each aligned sequence.

    Coverage = fraction of query positions covered by non-gap residues.
    """
    if not aligned_seqs:
        return []

    query_len = sum(1 for c in query_seq if c != '-')
    coverages = []

    for seq in aligned_seqs:
        covered = 0
        for i in range(min(len(query_seq), len(seq))):
            if query_seq[i] != '-' and seq[i] != '-':
                covered += 1

        coverage = covered / query_len if query_len > 0 else 0
        coverages.append(coverage)

    return coverages


def compute_msa_stats(filepath):
    """
    Compute comprehensive MSA statistics.

    Returns dict with all metrics.
    """
    query_seq, aligned_seqs = parse_a3m(filepath)

    if query_seq is None:
        return {"error": "Failed to parse MSA"}

    n_sequences = len(aligned_seqs) + 1  # Include query
    query_length = sum(1 for c in query_seq if c != '-')

    # Compute Neff
    all_seqs = [query_seq] + aligned_seqs
    weights = compute_sequence_weights(all_seqs)
    neff = sum(weights) if weights else 1.0

    # Column occupancy
    occupancy = compute_column_occupancy(query_seq, aligned_seqs)

    # Identity to query
    identities = compute_identity_to_query(query_seq, aligned_seqs)

    # Coverage
    coverages = compute_coverage(query_seq, aligned_seqs)

    # Aggregate stats
    stats = {
        "n_sequences": n_sequences,
        "neff": round(neff, 1),
        "query_length": query_length,
        "avg_identity": round(sum(identities) / len(identities), 3) if identities else None,
        "avg_coverage": round(sum(coverages) / len(coverages), 3) if coverages else None,
        "min_coverage": round(min(coverages), 3) if coverages else None,
        "max_coverage": round(max(coverages), 3) if coverages else None,
        "column_occupancy": {
            "min": round(min(occupancy), 3) if occupancy else None,
            "max": round(max(occupancy), 3) if occupancy else None,
            "mean": round(sum(occupancy) / len(occupancy), 3) if occupancy else None,
            "median": round(sorted(occupancy)[len(occupancy) // 2], 3) if occupancy else None,
        },
        "quality_flags": []
    }

    # Quality flags
    if neff < 30:
        stats["quality_flags"].append("LOW_NEFF")
    if neff < 100:
        stats["quality_flags"].append("MARGINAL_NEFF")

    if stats["column_occupancy"]["min"] is not None:
        if stats["column_occupancy"]["min"] < 0.3:
            stats["quality_flags"].append("LOW_OCCUPANCY_COLUMNS")

    if stats["avg_coverage"] is not None and stats["avg_coverage"] < 0.5:
        stats["quality_flags"].append("LOW_COVERAGE")

    if stats["avg_identity"] is not None and stats["avg_identity"] > 0.9:
        stats["quality_flags"].append("LOW_DIVERSITY")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Compute MSA quality statistics"
    )
    parser.add_argument("input", help="Input A3M file")
    parser.add_argument("--output", "-o", help="Output JSON file (default: stdout)")
    parser.add_argument("--tsv", action="store_true",
                        help="Output single-line TSV instead of JSON")

    args = parser.parse_args()

    try:
        stats = compute_msa_stats(args.input)
    except Exception as e:
        stats = {"error": str(e)}

    if args.tsv:
        # Single-line TSV output for easy aggregation
        cols = [
            stats.get("n_sequences", "NA"),
            stats.get("neff", "NA"),
            stats.get("query_length", "NA"),
            stats.get("avg_identity", "NA"),
            stats.get("avg_coverage", "NA"),
            stats["column_occupancy"].get("mean", "NA") if "column_occupancy" in stats else "NA",
            stats["column_occupancy"].get("min", "NA") if "column_occupancy" in stats else "NA",
            ",".join(stats.get("quality_flags", [])) or "OK"
        ]
        output = "\t".join(str(c) for c in cols)
    else:
        output = json.dumps(stats, indent=2)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output + "\n")
    else:
        print(output)


if __name__ == "__main__":
    main()
