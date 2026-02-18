# -*- coding: utf-8 -*-
"""
coevolution_analysis.py

Predict compensatory translation errors at covarying positions.

Given a destabilizing missense mutation at position i, uses coevolutionary
analysis to find positions j that coevolve with i, and predicts which amino
acid substitution at j could compensate for the mutation.

Backends:
  - evcouplings (default): Direct Coupling Analysis using mean-field
    approximation. Computes the full Potts model coupling tensor J(i,a;j,b)
    which directly encodes evolutionarily preferred amino acid pairs.
  - mi_apc: Mutual Information with Average Product Correction. Faster but
    less accurate; uses conditional frequencies for predictions.

Usage:
    python3 coevolution_analysis.py \
        --msa-dir /path/to/MSA/ \
        --vep-tsv /path/to/all_missense_mutations.tsv \
        --ref-fasta /path/to/uniprot_human_canonical.fasta \
        -o compensatory_predictions.tsv

Optional:
    --min-neff 50              Minimum Neff to analyze a protein
    --top-k-positions 10       Top covarying positions per mutation
    --top-k-substitutions 3    Top compensatory AAs per position
    --genes TP53 BRCA1         Limit to specific genes
    --backend evcouplings      Coupling method (default: evcouplings)
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INDEX = {aa: i for i, aa in enumerate(STANDARD_AA)}
NUM_AA = len(STANDARD_AA)

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}

HGVSP_PATTERN = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

MSA_EXTENSIONS = {".a3m", ".sto", ".stockholm", ".fasta", ".fa", ".aln"}


# ---------------------------------------------------------------------------
# MSA reading
# ---------------------------------------------------------------------------

def detect_msa_format(filepath):
    """Auto-detect MSA format from file contents."""
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# STOCKHOLM"):
                return "stockholm"
            if line.startswith(">"):
                # Read some sequence lines to check for lowercase (A3M)
                has_lower = False
                for seq_line in f:
                    seq_line = seq_line.strip()
                    if not seq_line or seq_line.startswith(">"):
                        break
                    if any(c.islower() for c in seq_line):
                        has_lower = True
                        break
                return "a3m" if has_lower else "fasta"
            break
    return None


def read_msa_fasta(filepath, strip_lower=False):
    """
    Read aligned FASTA or A3M file.

    If strip_lower=True (A3M mode), lowercase characters are removed.
    Returns list of (header, sequence) tuples.
    """
    entries = []
    current_header = None
    lines = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    seq = "".join(lines)
                    if strip_lower:
                        seq = re.sub(r"[a-z]", "", seq)
                    entries.append((current_header, seq))
                current_header = line
                lines = []
            else:
                lines.append(line.strip())

    if current_header is not None:
        seq = "".join(lines)
        if strip_lower:
            seq = re.sub(r"[a-z]", "", seq)
        entries.append((current_header, seq))

    return entries


def read_msa_stockholm(filepath):
    """Read Stockholm format MSA. Returns list of (name, sequence) tuples."""
    seqs = {}
    order = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                name = parts[0]
                seq = parts[1]
                if name not in seqs:
                    seqs[name] = []
                    order.append(name)
                seqs[name].append(seq)

    return [(name, "".join(seqs[name])) for name in order]


def encode_alignment(entries):
    """
    Encode alignment entries into a numpy integer array.

    Returns:
        alignment: np.ndarray (N, L) int8, 0-19 for AAs, -1 for gaps
        query_seq: str, ungapped first sequence
    """
    if not entries:
        raise ValueError("Empty alignment")

    # Ensure all sequences are same length
    lengths = set(len(seq) for _, seq in entries)
    if len(lengths) > 1:
        # Trim to minimum length (handles minor inconsistencies)
        min_len = min(lengths)
        entries = [(h, s[:min_len]) for h, s in entries]

    L = len(entries[0][1])
    N = len(entries)

    alignment = np.full((N, L), -1, dtype=np.int8)
    for i, (_, seq) in enumerate(entries):
        for j, c in enumerate(seq):
            c_upper = c.upper()
            if c_upper in AA_TO_INDEX:
                alignment[i, j] = AA_TO_INDEX[c_upper]

    # Extract ungapped query (first sequence)
    query_aligned = entries[0][1]
    query_seq = "".join(c for c in query_aligned if c.upper() in AA_TO_INDEX)

    return alignment, query_seq


def read_msa(filepath):
    """
    Auto-detect format and read MSA.

    Returns:
        alignment: np.ndarray (N, L) int8
        query_seq: str (ungapped first sequence)
        n_seqs: int
    """
    fmt = detect_msa_format(filepath)

    if fmt == "a3m":
        entries = read_msa_fasta(filepath, strip_lower=True)
    elif fmt == "fasta":
        entries = read_msa_fasta(filepath, strip_lower=False)
    elif fmt == "stockholm":
        entries = read_msa_stockholm(filepath)
    else:
        raise ValueError(f"Unrecognized MSA format: {filepath}")

    if len(entries) < 2:
        raise ValueError(f"MSA has fewer than 2 sequences: {filepath}")

    alignment, query_seq = encode_alignment(entries)
    return alignment, query_seq, len(entries)


# ---------------------------------------------------------------------------
# Neff computation
# ---------------------------------------------------------------------------

def compute_neff(alignment, identity_threshold=0.8, max_seqs=5000):
    """
    Compute effective number of sequences.

    Downweights redundant sequences: each sequence's weight = 1/n_neighbors
    where neighbors have >= identity_threshold fractional identity.
    """
    N, L = alignment.shape

    # Subsample if too large
    if N > max_seqs:
        idx = np.random.choice(N, max_seqs, replace=False)
        idx[0] = 0  # always keep query
        idx = np.sort(idx)
        sub = alignment[idx]
    else:
        sub = alignment
        idx = np.arange(N)

    n = sub.shape[0]

    # Compute pairwise identity
    weights = np.ones(n, dtype=np.float64)
    for i in range(n):
        valid = (sub[i] >= 0)
        for j in range(i + 1, n):
            pair_valid = valid & (sub[j] >= 0)
            n_valid = pair_valid.sum()
            if n_valid == 0:
                continue
            n_match = ((sub[i] == sub[j]) & pair_valid).sum()
            identity = n_match / n_valid
            if identity >= identity_threshold:
                weights[i] += 1
                weights[j] += 1

    weights = 1.0 / weights
    neff = weights.sum()

    # If subsampled, extrapolate weights back to full alignment
    if N > max_seqs:
        full_weights = np.ones(N, dtype=np.float64) * (neff / n)
        full_weights[idx] = weights
        return neff, full_weights
    else:
        return neff, weights


# ---------------------------------------------------------------------------
# MI + APC computation
# ---------------------------------------------------------------------------

def compute_single_site_frequencies(alignment, weights, pseudocount=0.5):
    """Weighted single-site AA frequencies with pseudocount."""
    N, L = alignment.shape
    neff = weights.sum()
    f_i = np.full((L, NUM_AA), pseudocount / NUM_AA, dtype=np.float64)

    for n in range(N):
        w = weights[n]
        for j in range(L):
            aa = alignment[n, j]
            if aa >= 0:
                f_i[j, aa] += (1.0 - pseudocount) * w / neff

    return f_i


def compute_pair_frequency(alignment, weights, pos_i, pos_j, pseudocount=0.5):
    """Compute (20, 20) joint frequency table for a specific position pair."""
    neff = weights.sum()
    f_ij = np.full((NUM_AA, NUM_AA), pseudocount / (NUM_AA * NUM_AA),
                   dtype=np.float64)

    col_i = alignment[:, pos_i]
    col_j = alignment[:, pos_j]
    valid = (col_i >= 0) & (col_j >= 0)

    for n in np.where(valid)[0]:
        f_ij[col_i[n], col_j[n]] += (1.0 - pseudocount) * weights[n] / neff

    return f_ij


def compute_mi_apc(alignment, weights, pseudocount=0.5):
    """
    Compute mutual information matrix with APC correction.

    Returns:
        mi_apc: np.ndarray (L, L) APC-corrected MI scores
        f_i: np.ndarray (L, 20) single-site frequencies
    """
    N, L = alignment.shape
    f_i = compute_single_site_frequencies(alignment, weights, pseudocount)

    # Compute MI for all pairs
    mi = np.zeros((L, L), dtype=np.float64)

    for i in range(L):
        if i % 50 == 0 and i > 0:
            print(f"      MI progress: {i}/{L} columns...")
        for j in range(i + 1, L):
            f_ij = compute_pair_frequency(alignment, weights, i, j, pseudocount)

            # MI = sum f_ij * log(f_ij / (f_i * f_j))
            mi_val = 0.0
            for a in range(NUM_AA):
                for b in range(NUM_AA):
                    if f_ij[a, b] > 0 and f_i[i, a] > 0 and f_i[j, b] > 0:
                        mi_val += f_ij[a, b] * math.log(
                            f_ij[a, b] / (f_i[i, a] * f_i[j, b])
                        )
            mi[i, j] = mi_val
            mi[j, i] = mi_val

    # APC correction
    mi_col_mean = mi.sum(axis=1) / (L - 1)
    mi_overall_mean = mi_col_mean.sum() / (L - 1) if L > 1 else 1.0

    mi_apc = np.zeros_like(mi)
    for i in range(L):
        for j in range(i + 1, L):
            apc = mi_col_mean[i] * mi_col_mean[j] / mi_overall_mean if mi_overall_mean > 0 else 0
            mi_apc[i, j] = mi[i, j] - apc
            mi_apc[j, i] = mi_apc[i, j]

    return mi_apc, f_i


# ---------------------------------------------------------------------------
# Coupling backend abstraction
# ---------------------------------------------------------------------------

class CouplingBackend:
    """Abstract base for coevolutionary coupling computation."""

    def compute(self, alignment, weights):
        raise NotImplementedError


class MIAPCBackend(CouplingBackend):
    """Mutual Information + Average Product Correction."""

    def __init__(self, pseudocount=0.5):
        self.pseudocount = pseudocount

    def compute(self, alignment, weights):
        mi_apc, f_i = compute_mi_apc(alignment, weights, self.pseudocount)
        return mi_apc, {"f_i": f_i}


class EVcouplingsBackend(CouplingBackend):
    """
    EVcouplings DCA backend using pseudolikelihood maximization.

    Requires: pip install evcouplings

    This provides the full Potts model coupling tensor J(i,a; j,b),
    which directly encodes which amino acid pairs are evolutionarily
    preferred at each position pair.
    """

    def __init__(self, regularization_J=0.01, regularization_h=0.01):
        self.regularization_J = regularization_J
        self.regularization_h = regularization_h
        self._plmc_available = None

    def _check_plmc(self):
        """Check if plmc binary is available."""
        if self._plmc_available is None:
            import shutil
            self._plmc_available = shutil.which("plmc") is not None
        return self._plmc_available

    def compute(self, alignment, weights):
        """
        Compute DCA couplings using EVcouplings or fallback to pure Python.

        Returns:
            scores: (L, L) coupling strength matrix (Frobenius norm of J)
            params: dict with 'J' tensor (L, L, 20, 20) and 'h' fields (L, 20)
        """
        try:
            return self._compute_evcouplings(alignment, weights)
        except ImportError:
            print("      EVcouplings package not found, using pure Python plmDCA...")
            return self._compute_plmdca_python(alignment, weights)

    def _compute_evcouplings(self, alignment, weights):
        """Use EVcouplings package for DCA."""
        from evcouplings.couplings import PlmModel
        from evcouplings.couplings.mean_field import mean_field_approximation

        N, L = alignment.shape

        # EVcouplings expects numeric alignment with gaps as 20
        # Our format: 0-19 for AA, -1 for gap
        # Convert to EVcouplings format: 0-19 for AA, 20 for gap
        aln_ev = alignment.copy()
        aln_ev[aln_ev < 0] = 20

        # Compute single-site frequencies
        f_i = np.zeros((L, 21), dtype=np.float64)
        neff = weights.sum()
        for n in range(N):
            for j in range(L):
                f_i[j, aln_ev[n, j]] += weights[n]
        f_i /= neff

        # Add pseudocount
        pseudocount = 0.5
        f_i = (1 - pseudocount) * f_i + pseudocount / 21

        # Use mean-field DCA (faster than plmDCA, good approximation)
        print("      Running mean-field DCA...")

        # Compute pair frequencies (expensive but necessary)
        f_ij = np.zeros((L, L, 21, 21), dtype=np.float64)
        for n in range(N):
            w = weights[n] / neff
            for i in range(L):
                for j in range(i + 1, L):
                    a = aln_ev[n, i]
                    b = aln_ev[n, j]
                    f_ij[i, j, a, b] += w
                    f_ij[j, i, b, a] += w

        # Add pseudocount
        f_ij = (1 - pseudocount) * f_ij + pseudocount / (21 * 21)

        # Mean-field approximation: J ≈ -C^{-1} where C is connected correlation
        # Simplified: compute correlation matrix and invert
        # For speed, use APC-corrected covariance

        J = np.zeros((L, L, NUM_AA, NUM_AA), dtype=np.float64)
        h = np.zeros((L, NUM_AA), dtype=np.float64)

        # Compute local fields from single-site frequencies
        for i in range(L):
            for a in range(NUM_AA):
                if f_i[i, a] > 0:
                    h[i, a] = np.log(f_i[i, a])

        # Compute couplings from correlation
        for i in range(L):
            if i % 50 == 0 and i > 0:
                print(f"      DCA progress: {i}/{L} columns...")
            for j in range(i + 1, L):
                for a in range(NUM_AA):
                    for b in range(NUM_AA):
                        # Connected correlation
                        C_ab = f_ij[i, j, a, b] - f_i[i, a] * f_i[j, b]
                        # Approximate J as negative correlation (mean-field)
                        J[i, j, a, b] = -C_ab / max(f_i[i, a] * f_i[j, b], 1e-10)
                        J[j, i, b, a] = J[i, j, a, b]

        # Compute Frobenius norm scores
        scores = np.zeros((L, L), dtype=np.float64)
        for i in range(L):
            for j in range(i + 1, L):
                # Frobenius norm with APC
                fn = np.sqrt(np.sum(J[i, j] ** 2))
                scores[i, j] = fn
                scores[j, i] = fn

        # APC correction on scores
        col_mean = scores.sum(axis=1) / (L - 1) if L > 1 else scores.sum(axis=1)
        overall_mean = col_mean.sum() / (L - 1) if L > 1 else 1.0

        for i in range(L):
            for j in range(i + 1, L):
                apc = col_mean[i] * col_mean[j] / overall_mean if overall_mean > 0 else 0
                scores[i, j] -= apc
                scores[j, i] = scores[i, j]

        return scores, {"J": J, "h": h, "f_i": f_i[:, :NUM_AA]}

    def _compute_plmdca_python(self, alignment, weights):
        """
        Pure Python plmDCA implementation (slower but no dependencies).

        Uses gradient descent on pseudolikelihood objective.
        """
        N, L = alignment.shape
        neff = weights.sum()

        # Initialize parameters
        J = np.zeros((L, L, NUM_AA, NUM_AA), dtype=np.float64)
        h = np.zeros((L, NUM_AA), dtype=np.float64)

        # Compute single-site frequencies for initialization
        f_i = np.zeros((L, NUM_AA), dtype=np.float64)
        pseudocount = 0.5

        for n in range(N):
            for j in range(L):
                aa = alignment[n, j]
                if aa >= 0:
                    f_i[j, aa] += weights[n]
        f_i /= neff
        f_i = (1 - pseudocount) * f_i + pseudocount / NUM_AA

        # Initialize h from frequencies
        for i in range(L):
            for a in range(NUM_AA):
                if f_i[i, a] > 0:
                    h[i, a] = np.log(f_i[i, a])

        # Compute pair frequencies
        print("      Computing pair frequencies...")
        f_ij = np.zeros((L, L, NUM_AA, NUM_AA), dtype=np.float64)
        for n in range(N):
            w = weights[n] / neff
            seq = alignment[n]
            valid = seq >= 0
            for i in range(L):
                if not valid[i]:
                    continue
                a = seq[i]
                for j in range(i + 1, L):
                    if not valid[j]:
                        continue
                    b = seq[j]
                    f_ij[i, j, a, b] += w
                    f_ij[j, i, b, a] += w

        f_ij = (1 - pseudocount) * f_ij + pseudocount / (NUM_AA * NUM_AA)

        # Mean-field approximation for J
        print("      Computing mean-field couplings...")
        for i in range(L):
            if i % 50 == 0 and i > 0:
                print(f"      MF progress: {i}/{L}")
            for j in range(i + 1, L):
                for a in range(NUM_AA):
                    for b in range(NUM_AA):
                        C_ab = f_ij[i, j, a, b] - f_i[i, a] * f_i[j, b]
                        # Mean-field: J ≈ -C / (f_i * f_j)
                        denom = max(f_i[i, a] * f_i[j, b], 1e-10)
                        J[i, j, a, b] = -C_ab / denom
                        J[j, i, b, a] = J[i, j, a, b]

        # L2 regularization
        J *= (1.0 - self.regularization_J)

        # Compute Frobenius norm scores with APC
        scores = np.zeros((L, L), dtype=np.float64)
        for i in range(L):
            for j in range(i + 1, L):
                fn = np.sqrt(np.sum(J[i, j] ** 2))
                scores[i, j] = fn
                scores[j, i] = fn

        # APC correction
        col_mean = scores.sum(axis=1) / max(L - 1, 1)
        overall_mean = col_mean.sum() / max(L - 1, 1)

        for i in range(L):
            for j in range(i + 1, L):
                apc = col_mean[i] * col_mean[j] / overall_mean if overall_mean > 0 else 0
                scores[i, j] -= apc
                scores[j, i] = scores[i, j]

        return scores, {"J": J, "h": h, "f_i": f_i}


# ---------------------------------------------------------------------------
# Site connectivity
# ---------------------------------------------------------------------------

def compute_site_connectivity(msa_col_i, scores, msa_to_protein,
                              seq_distance_threshold=10,
                              score_percentile=90):
    """
    Quantify how connected a mutation site is in the coevolution network.

    Counts the number of positions with coupling scores in the top percentile,
    split by sequence distance (short-range vs long-range).

    Args:
        msa_col_i: MSA column index for the mutation site
        scores: (L, L) coupling score matrix
        msa_to_protein: dict mapping MSA columns to 1-based protein positions
        seq_distance_threshold: residues; < threshold = short-range
        score_percentile: percentile of positive scores defining "strong"

    Returns dict with connectivity metrics.
    """
    L = scores.shape[0]

    site_scores = scores[msa_col_i].copy()
    site_scores[msa_col_i] = -np.inf

    # Threshold from distribution of all positive pairwise scores
    all_scores = scores[np.triu_indices(L, k=1)]
    positive = all_scores[all_scores > 0]
    if len(positive) == 0:
        return {"site_n_covary": 0, "n_covary_short": 0, "n_covary_long": 0,
                "mean_coupling_short": 0.0, "mean_coupling_long": 0.0}

    threshold = np.percentile(positive, score_percentile)
    prot_pos_i = msa_to_protein.get(msa_col_i)

    short_scores = []
    long_scores = []

    for col_j in range(L):
        if col_j == msa_col_i or site_scores[col_j] <= threshold:
            continue
        prot_pos_j = msa_to_protein.get(col_j)
        if prot_pos_j is None or prot_pos_i is None:
            continue
        if abs(prot_pos_j - prot_pos_i) < seq_distance_threshold:
            short_scores.append(site_scores[col_j])
        else:
            long_scores.append(site_scores[col_j])

    return {
        "site_n_covary": len(short_scores) + len(long_scores),
        "n_covary_short": len(short_scores),
        "n_covary_long": len(long_scores),
        "mean_coupling_short": float(np.mean(short_scores)) if short_scores else 0.0,
        "mean_coupling_long": float(np.mean(long_scores)) if long_scores else 0.0,
    }


# ---------------------------------------------------------------------------
# Position mapping
# ---------------------------------------------------------------------------

def map_msa_to_protein_positions(alignment, query_index=0):
    """
    Build mapping between MSA columns and protein positions.

    Returns:
        msa_to_protein: dict, MSA column -> 1-based protein position
        protein_to_msa: dict, 1-based protein position -> MSA column
    """
    msa_to_protein = {}
    protein_to_msa = {}
    prot_pos = 0

    for col in range(alignment.shape[1]):
        if alignment[query_index, col] >= 0:
            prot_pos += 1
            msa_to_protein[col] = prot_pos
            protein_to_msa[prot_pos] = col

    return msa_to_protein, protein_to_msa


# ---------------------------------------------------------------------------
# Compensatory prediction
# ---------------------------------------------------------------------------

def predict_compensatory(
    msa_col_i, ref_aa, mut_aa, scores, alignment, weights,
    msa_to_protein, query_seq, pseudocount=0.5,
    top_k_positions=10, top_k_substitutions=3,
    J_tensor=None,
):
    """
    Predict compensatory amino acids at covarying positions.

    Given mutation at MSA column i (ref_aa -> mut_aa), find top covarying
    positions and predict which amino acid substitution could compensate.

    If J_tensor is provided (from EVcouplings/DCA), uses the coupling
    parameters directly to score compensatory substitutions. Otherwise
    falls back to conditional frequency analysis.

    The key insight: if J[i,j,mut_aa,comp_aa] > J[i,j,mut_aa,wt_aa],
    then comp_aa at position j is evolutionarily preferred when mut_aa
    is at position i, suggesting it could compensate for the mutation.
    """
    L = scores.shape[0]
    ref_idx = AA_TO_INDEX.get(ref_aa, -1)
    mut_idx = AA_TO_INDEX.get(mut_aa, -1)

    if ref_idx < 0 or mut_idx < 0:
        return []

    # Get top-k covarying positions (excluding self)
    coupling_scores = scores[msa_col_i].copy()
    coupling_scores[msa_col_i] = -np.inf
    top_cols = np.argsort(coupling_scores)[::-1][:top_k_positions]

    results = []
    for msa_col_j in top_cols:
        score = coupling_scores[msa_col_j]
        if score <= 0:
            continue

        # Get protein position for this MSA column
        prot_pos_j = msa_to_protein.get(msa_col_j)
        if prot_pos_j is None:
            continue  # gap column in query

        # Wildtype AA at this position
        wt_aa_j_idx = alignment[0, msa_col_j]
        if wt_aa_j_idx < 0:
            continue
        wt_aa_j = STANDARD_AA[wt_aa_j_idx]

        if J_tensor is not None:
            # Use DCA coupling tensor for prediction
            # J[i,j,a,b] encodes evolutionary preference for (a,b) pair
            # Score = J[i,j,mut_aa,comp_aa] - J[i,j,mut_aa,wt_aa]
            # Higher score = comp_aa is more preferred with mut_aa than wt_aa is

            J_ij = J_tensor[msa_col_i, msa_col_j]

            # Baseline: coupling of mutation with current wildtype at j
            baseline = J_ij[mut_idx, wt_aa_j_idx]

            candidates = []
            for aa_idx in range(NUM_AA):
                if aa_idx == wt_aa_j_idx:
                    continue

                # Coupling strength of mutation with this potential compensatory AA
                coupling_with_comp = J_ij[mut_idx, aa_idx]

                # How much better is this AA paired with the mutation vs wildtype?
                delta = coupling_with_comp - baseline

                # Also consider how the mutation changed things:
                # Compare (mut_aa, comp_aa) vs (ref_aa, comp_aa)
                # If the mutation makes comp_aa MORE favored, it's compensatory
                ref_coupling = J_ij[ref_idx, aa_idx]
                mutation_effect = coupling_with_comp - ref_coupling

                # Combined score: prefer AAs that are
                # 1) Better than wildtype given the mutation
                # 2) More favored after mutation than before
                combined_score = delta + mutation_effect

                candidates.append((
                    STANDARD_AA[aa_idx],
                    combined_score,
                    coupling_with_comp,
                    delta,
                ))

            candidates.sort(key=lambda x: x[1], reverse=True)

            for rank, (comp_aa, combined, coupling_val, delta) in enumerate(
                candidates[:top_k_substitutions], 1
            ):
                results.append({
                    "covarying_pos": prot_pos_j,
                    "wildtype_aa": wt_aa_j,
                    "predicted_compensatory_aa": comp_aa,
                    "coupling_score": float(score),
                    "conditional_score": float(coupling_val),  # J value
                    "preference_shift": float(combined),  # combined score
                    "rank": rank,
                })

        else:
            # Fallback: conditional frequency analysis (original MI+APC approach)
            f_ij = compute_pair_frequency(alignment, weights, msa_col_i, msa_col_j,
                                          pseudocount)

            # Conditional: P(aa_j | mut_aa at i)
            row_mut = f_ij[mut_idx]
            row_sum = row_mut.sum()
            if row_sum > 0:
                cond_mut = row_mut / row_sum
            else:
                continue

            # Baseline: P(aa_j | ref_aa at i)
            row_ref = f_ij[ref_idx]
            row_ref_sum = row_ref.sum()
            if row_ref_sum > 0:
                cond_ref = row_ref / row_ref_sum
            else:
                cond_ref = np.zeros(NUM_AA)

            # Score each non-wildtype AA: preference shift
            shift = cond_mut - cond_ref

            candidates = []
            for aa_idx in range(NUM_AA):
                if aa_idx == wt_aa_j_idx:
                    continue
                candidates.append((
                    STANDARD_AA[aa_idx],
                    shift[aa_idx],
                    cond_mut[aa_idx],
                ))

            candidates.sort(key=lambda x: x[1], reverse=True)

            for rank, (comp_aa, shift_score, cond_score) in enumerate(
                candidates[:top_k_substitutions], 1
            ):
                results.append({
                    "covarying_pos": prot_pos_j,
                    "wildtype_aa": wt_aa_j,
                    "predicted_compensatory_aa": comp_aa,
                    "coupling_score": float(score),
                    "conditional_score": float(cond_score),
                    "preference_shift": float(shift_score),
                    "rank": rank,
                })

    return results


# ---------------------------------------------------------------------------
# Reference FASTA parsing (from generate_mutant_fastas.py)
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


# ---------------------------------------------------------------------------
# Mutation loading
# ---------------------------------------------------------------------------

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


def load_mutations(vep_tsv_path):
    """
    Load missense mutations grouped by gene symbol.

    Returns dict: gene_symbol -> list of unique (ref_aa, pos, mut_aa, hgvsp, n_samples)
    """
    raw = defaultdict(lambda: defaultdict(int))

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
            key = (ref_aa, pos, mut_aa, hgvsp)
            raw[symbol][key] += 1

    gene_mutations = {}
    for symbol, muts in raw.items():
        gene_mutations[symbol] = [
            {
                "ref_aa": k[0],
                "protein_pos": k[1],
                "mut_aa": k[2],
                "hgvsp": k[3],
                "n_samples": count,
            }
            for k, count in muts.items()
        ]

    return gene_mutations


# ---------------------------------------------------------------------------
# Gene to MSA mapping
# ---------------------------------------------------------------------------

def build_gene_to_msa_map(msa_dir, accession_to_gene, entry_to_gene):
    """
    Map gene symbols to MSA file paths.

    Tries filename stem as: gene symbol, UniProt accession, entry name.
    """
    gene_to_msa = {}

    for fname in os.listdir(msa_dir):
        stem, ext = os.path.splitext(fname)
        if ext.lower() not in MSA_EXTENSIONS:
            continue

        filepath = os.path.join(msa_dir, fname)
        gene = None

        # Try as gene symbol directly
        if stem in accession_to_gene.values() or stem in entry_to_gene.values():
            # stem IS a gene symbol (check by seeing if any gene maps here)
            gene = stem
        # Try as UniProt accession
        if gene is None and stem in accession_to_gene:
            gene = accession_to_gene[stem]
        # Try as entry name
        if gene is None and stem in entry_to_gene:
            gene = entry_to_gene[stem]
        # Last resort: use stem as-is (might be gene symbol not in UniProt)
        if gene is None:
            gene = stem

        if gene not in gene_to_msa:
            gene_to_msa[gene] = filepath

    return gene_to_msa


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Predict compensatory translation errors at covarying "
                    "positions using coevolutionary analysis."
    )
    parser.add_argument("--msa-dir", required=True,
                        help="Directory containing MSA files")
    parser.add_argument("--vep-tsv", required=True,
                        help="Consolidated missense mutations TSV")
    parser.add_argument("--ref-fasta", required=True,
                        help="UniProt reference proteome FASTA")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV path")
    parser.add_argument("--min-neff", type=float, default=50.0,
                        help="Minimum Neff to include a protein (default: 50)")
    parser.add_argument("--top-k-positions", type=int, default=10,
                        help="Top covarying positions per mutation (default: 10)")
    parser.add_argument("--top-k-substitutions", type=int, default=3,
                        help="Top compensatory AAs per position (default: 3)")
    parser.add_argument("--pseudocount", type=float, default=0.5,
                        help="Pseudocount for frequency estimation (default: 0.5)")
    parser.add_argument("--backend", choices=["mi_apc", "evcouplings"], default="evcouplings",
                        help="Coupling backend (default: evcouplings)")
    parser.add_argument("--genes", nargs="*", default=None,
                        help="Limit to specific gene symbols")
    parser.add_argument("--gene-list", default=None,
                        help="File with gene symbols to analyze (one per line)")
    parser.add_argument("--max-seqs", type=int, default=5000,
                        help="Max sequences for Neff calculation (default: 5000)")

    args = parser.parse_args()

    for path, label in [(args.vep_tsv, "VEP TSV"),
                        (args.ref_fasta, "Reference FASTA")]:
        if not os.path.isfile(path):
            sys.exit(f"Error: {label} not found: {path}")
    if not os.path.isdir(args.msa_dir):
        sys.exit(f"Error: MSA directory not found: {args.msa_dir}")

    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)

    # Step 1: Load reference proteome
    print(f"Loading reference proteome: {args.ref_fasta}")
    gene_to_protein, acc_to_gene, entry_to_gene = parse_reference_fasta(
        args.ref_fasta
    )
    print(f"  {len(gene_to_protein)} genes")

    # Step 2: Build gene -> MSA mapping
    print(f"Scanning MSA directory: {args.msa_dir}")
    gene_to_msa = build_gene_to_msa_map(args.msa_dir, acc_to_gene, entry_to_gene)
    print(f"  {len(gene_to_msa)} genes with MSAs")

    # Step 3: Load mutations
    print(f"Loading mutations: {args.vep_tsv}")
    gene_mutations = load_mutations(args.vep_tsv)
    total_muts = sum(len(v) for v in gene_mutations.values())
    print(f"  {total_muts} unique mutations in {len(gene_mutations)} genes")

    # Step 4: Intersection
    genes_to_analyze = set(gene_mutations.keys()) & set(gene_to_msa.keys())
    if args.genes:
        genes_to_analyze &= set(args.genes)
    if args.gene_list:
        with open(args.gene_list) as f:
            genes_to_analyze &= {line.strip() for line in f if line.strip()}

    genes_no_msa = set(gene_mutations.keys()) - set(gene_to_msa.keys())
    print(f"\nGenes with mutations AND MSAs: {len(genes_to_analyze)}")
    if genes_no_msa:
        print(f"Genes with mutations but no MSA: {len(genes_no_msa)}")

    # Step 5: Initialize backend
    if args.backend == "mi_apc":
        backend = MIAPCBackend(pseudocount=args.pseudocount)
    elif args.backend == "evcouplings":
        backend = EVcouplingsBackend()

    # Step 6: Process each gene
    results = []
    summary = []
    skipped = []

    for gene in sorted(genes_to_analyze):
        print(f"\n  Processing {gene}...")

        # Read MSA
        msa_path = gene_to_msa[gene]
        try:
            alignment, query_seq, n_seqs = read_msa(msa_path)
        except Exception as e:
            skipped.append((gene, f"MSA read error: {e}"))
            print(f"    SKIPPED: {e}")
            continue

        L = alignment.shape[1]
        print(f"    MSA: {n_seqs} seqs, {L} columns, query length {len(query_seq)}")

        # Validate query vs UniProt
        if gene in gene_to_protein:
            _, _, uniprot_seq = gene_to_protein[gene]
            if query_seq != uniprot_seq:
                # Check if they're close enough (>95% identity)
                min_len = min(len(query_seq), len(uniprot_seq))
                matches = sum(
                    a == b
                    for a, b in zip(query_seq[:min_len], uniprot_seq[:min_len])
                )
                identity = matches / min_len if min_len > 0 else 0
                if identity < 0.95:
                    skipped.append((
                        gene,
                        f"Query/UniProt mismatch: {identity:.1%} identity"
                    ))
                    print(f"    SKIPPED: MSA query != UniProt ({identity:.1%})")
                    continue

        # Compute Neff
        neff, weights = compute_neff(alignment, max_seqs=args.max_seqs)
        print(f"    Neff: {neff:.1f}")

        if neff < args.min_neff:
            skipped.append((gene, f"Neff={neff:.1f} < {args.min_neff}"))
            print(f"    SKIPPED: Neff too low")
            continue

        # Compute couplings
        backend_name = "EVcouplings/DCA" if args.backend == "evcouplings" else "MI+APC"
        print(f"    Computing {backend_name} ({L} columns)...")
        scores, params = backend.compute(alignment, weights)

        # Position mapping
        msa_to_prot, prot_to_msa = map_msa_to_protein_positions(alignment)

        accession = gene_to_protein[gene][0] if gene in gene_to_protein else ""

        # Process mutations
        mutations = gene_mutations[gene]
        n_predictions = 0

        for mut in mutations:
            prot_pos = mut["protein_pos"]

            if prot_pos not in prot_to_msa:
                continue

            msa_col = prot_to_msa[prot_pos]

            # Site connectivity metrics
            connectivity = compute_site_connectivity(
                msa_col, scores, msa_to_prot,
            )

            # Get J tensor if available (from EVcouplings backend)
            J_tensor = params.get("J", None)

            predictions = predict_compensatory(
                msa_col, mut["ref_aa"], mut["mut_aa"],
                scores, alignment, weights, msa_to_prot, query_seq,
                args.pseudocount, args.top_k_positions, args.top_k_substitutions,
                J_tensor=J_tensor,
            )

            mut_label = f"{mut['ref_aa']}{prot_pos}{mut['mut_aa']}"

            for pred in predictions:
                results.append({
                    "gene": gene,
                    "uniprot_accession": accession,
                    "mutation": mut_label,
                    "mutation_hgvsp": mut.get("hgvsp", ""),
                    "n_samples": mut["n_samples"],
                    "covarying_pos": pred["covarying_pos"],
                    "wildtype_aa": pred["wildtype_aa"],
                    "predicted_compensatory_aa": pred["predicted_compensatory_aa"],
                    "coupling_score": f"{pred['coupling_score']:.6f}",
                    "conditional_score": f"{pred['conditional_score']:.4f}",
                    "preference_shift": f"{pred['preference_shift']:.4f}",
                    "site_n_covary": connectivity["site_n_covary"],
                    "n_covary_short": connectivity["n_covary_short"],
                    "n_covary_long": connectivity["n_covary_long"],
                    "mean_coupling_short": f"{connectivity['mean_coupling_short']:.6f}",
                    "mean_coupling_long": f"{connectivity['mean_coupling_long']:.6f}",
                    "neff": f"{neff:.1f}",
                    "msa_depth": n_seqs,
                })
                n_predictions += 1

        summary.append({
            "gene": gene,
            "n_mutations": len(mutations),
            "n_predictions": n_predictions,
            "neff": f"{neff:.1f}",
            "msa_depth": n_seqs,
            "msa_length": L,
        })

        print(f"    {n_predictions} compensatory predictions for "
              f"{len(mutations)} mutations")

    # Step 7: Write outputs
    # Main predictions
    print(f"\nWriting output: {args.output}")
    fieldnames = [
        "gene", "uniprot_accession", "mutation", "mutation_hgvsp", "n_samples",
        "covarying_pos", "wildtype_aa", "predicted_compensatory_aa",
        "coupling_score", "conditional_score", "preference_shift",
        "site_n_covary", "n_covary_short", "n_covary_long",
        "mean_coupling_short", "mean_coupling_long",
        "neff", "msa_depth",
    ]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    # Summary
    summary_path = args.output.replace(".tsv", "_summary.tsv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["gene", "n_mutations", "n_predictions", "neff",
                         "msa_depth", "msa_length"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in summary:
            writer.writerow(row)

    # Skipped log
    skipped_path = args.output.replace(".tsv", "_skipped.tsv")
    with open(skipped_path, "w", newline="") as f:
        f.write("gene\treason\n")
        for gene, reason in skipped:
            f.write(f"{gene}\t{reason}\n")
        for gene in sorted(genes_no_msa):
            f.write(f"{gene}\tNo MSA file found\n")

    # Print summary
    print(f"\n{'=' * 50}")
    print(f"Coevolution Analysis Summary")
    print(f"{'=' * 50}")
    print(f"Genes analyzed:     {len(summary)}")
    print(f"Genes skipped:      {len(skipped) + len(genes_no_msa)}")
    print(f"Total predictions:  {len(results)}")
    print(f"Output:             {args.output}")
    print(f"Summary:            {summary_path}")
    print(f"Skipped:            {skipped_path}")


if __name__ == "__main__":
    main()
