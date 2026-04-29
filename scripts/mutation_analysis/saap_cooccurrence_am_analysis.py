#!/usr/bin/env python3
"""
AM score distribution: co-occurring missense vs background
==========================================================
Compares AlphaMissense pathogenicity scores of missense mutations that
co-occur with a SAAP (in the same gene, same TMT set) against the
background distribution of all missense in those same genes.

Outputs
-------
am_distribution_plot.pdf   — histogram + CDF comparison
am_stats.tsv               — summary statistics
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parents[2]

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cooccurrence",
                    default="/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence/saap_missense_cooccurrence.tsv")
    ap.add_argument("--missense",
                    required=True,
                    help="all_missense_with_spurs.tsv")
    ap.add_argument("-o", "--out-dir",
                    default="/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence")
    return ap.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load co-occurrence results ────────────────────────────────────────────
    print("Loading co-occurrence table...", flush=True)
    cooc = pd.read_csv(args.cooccurrence, sep="\t")
    cooc["missense_am"] = pd.to_numeric(cooc["missense_am"], errors="coerce")
    cooc["missense_spurs_ddg"] = pd.to_numeric(cooc["missense_spurs_ddg"], errors="coerce")
    print(f"  Co-occurrence records: {len(cooc):,}", flush=True)

    # Genes that appear in any co-occurrence (the gene universe for background)
    cooc_genes = set(cooc["gene"].dropna().unique())
    print(f"  Unique genes with co-occurrences: {len(cooc_genes)}", flush=True)

    # ── Load background: all missense in those genes ──────────────────────────
    print("Loading missense table...", flush=True)
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["SYMBOL","am_pathogenicity","spurs_ddg"])
    miss["am_pathogenicity"] = pd.to_numeric(miss["am_pathogenicity"], errors="coerce")
    miss["spurs_ddg"]        = pd.to_numeric(miss["spurs_ddg"], errors="coerce")
    print(f"  Total missense rows: {len(miss):,}", flush=True)

    background = miss[miss["SYMBOL"].isin(cooc_genes)]["am_pathogenicity"].dropna()
    print(f"  Background AM scores (same genes): {len(background):,}", flush=True)

    # Co-occurring AM scores
    cooc_am = cooc["missense_am"].dropna()
    print(f"  Co-occurring AM scores: {len(cooc_am):,}", flush=True)

    # ── Statistics ────────────────────────────────────────────────────────────
    u, p_mw = stats.mannwhitneyu(cooc_am, background, alternative="two-sided")
    ks, p_ks = stats.ks_2samp(cooc_am, background)

    print(f"\n── AM SCORE STATISTICS ──────────────────────────────────────────────")
    print(f"Co-occurring:  n={len(cooc_am):,}  mean={cooc_am.mean():.3f}  "
          f"median={cooc_am.median():.3f}")
    print(f"Background:    n={len(background):,}  mean={background.mean():.3f}  "
          f"median={background.median():.3f}")
    print(f"Mann-Whitney U (two-sided): p={p_mw:.2e}")
    print(f"KS test:                   D={ks:.3f}  p={p_ks:.2e}")

    # AM class breakdown
    bins_labels = [(0.0, 0.34, "likely_benign"),
                   (0.34, 0.564, "ambiguous"),
                   (0.564, 1.01, "likely_pathogenic")]
    print(f"\n{'Category':<22} {'Co-occurring':>14} {'Background':>12}")
    for lo, hi, label in bins_labels:
        n_c = ((cooc_am >= lo) & (cooc_am < hi)).sum()
        n_b = ((background >= lo) & (background < hi)).sum()
        pct_c = 100 * n_c / len(cooc_am)
        pct_b = 100 * n_b / len(background)
        print(f"{label:<22} {n_c:>6} ({pct_c:4.1f}%)   {n_b:>7} ({pct_b:4.1f}%)")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel 1: overlapping histograms (density)
    ax = axes[0]
    bins = np.linspace(0, 1, 41)
    ax.hist(background, bins=bins, density=True, alpha=0.5, color="#888888",
            label=f"Background\n(all missense, same genes)\nn={len(background):,}  "
                  f"median={background.median():.3f}")
    ax.hist(cooc_am, bins=bins, density=True, alpha=0.7, color="#d62728",
            label=f"Co-occurring\n(with SAAP in same gene/set)\nn={len(cooc_am):,}  "
                  f"median={cooc_am.median():.3f}")
    for threshold, ls, label in [(0.34, "--", "benign/ambiguous"), (0.564, "-.", "ambiguous/pathogenic")]:
        ax.axvline(threshold, color="k", lw=0.8, ls=ls, alpha=0.5)
    ax.set_xlabel("AlphaMissense pathogenicity score", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(f"AM score distribution\nMann-Whitney p={p_mw:.2e} | KS D={ks:.3f} p={p_ks:.2e}",
                 fontsize=9)
    ax.legend(fontsize=8)

    # Panel 2: empirical CDF
    ax2 = axes[1]
    for data, color, label in [
        (background, "#888888", f"Background (n={len(background):,})"),
        (cooc_am,    "#d62728", f"Co-occurring (n={len(cooc_am):,})"),
    ]:
        xs = np.sort(data)
        ys = np.arange(1, len(xs)+1) / len(xs)
        ax2.plot(xs, ys, color=color, lw=1.5, label=label)
    for threshold, ls in [(0.34, "--"), (0.564, "-.")]:
        ax2.axvline(threshold, color="k", lw=0.8, ls=ls, alpha=0.5)
    ax2.set_xlabel("AlphaMissense pathogenicity score", fontsize=10)
    ax2.set_ylabel("Cumulative fraction", fontsize=10)
    ax2.set_title("Empirical CDF", fontsize=9)
    ax2.legend(fontsize=8)

    plt.suptitle("AlphaMissense scores: co-occurring missense vs background\n"
                 "(background = all missense in genes that have any SAAP co-occurrence)",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig(out_dir / "am_distribution_plot.pdf", dpi=150)
    plt.savefig(out_dir / "am_distribution_plot.png", dpi=150)
    print(f"\nPlot saved to {out_dir}/am_distribution_plot.pdf", flush=True)

    # ── Save stats ────────────────────────────────────────────────────────────
    stat_rows = []
    for lo, hi, label in bins_labels:
        n_c = ((cooc_am >= lo) & (cooc_am < hi)).sum()
        n_b = ((background >= lo) & (background < hi)).sum()
        stat_rows.append({
            "am_class": label,
            "cooc_n": n_c, "cooc_pct": round(100*n_c/len(cooc_am), 2),
            "bg_n": n_b,   "bg_pct": round(100*n_b/len(background), 2),
        })
    stat_df = pd.DataFrame(stat_rows)
    stat_df.to_csv(out_dir / "am_stats.tsv", sep="\t", index=False)
    print(f"Stats saved to {out_dir}/am_stats.tsv", flush=True)


if __name__ == "__main__":
    main()
