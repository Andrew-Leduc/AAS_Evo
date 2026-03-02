#!/usr/bin/env python3
"""
plot_bam_support.py

Visualize not-have BAM alt-read support for low-ratio mutant peptides.

For each low-ratio PSM, shows a bar chart: not-have patient channels on the
x-axis, alt_vaf in their BAM on the y-axis.  Bars are coloured red if
alt_vaf > 0.05 (potentially explains denominator inflation).

Usage:
    python3 plot_bam_support.py [options]

Options:
    --bam-support   Path to bam_support.tsv (default: OUT_DIR/bam_support.tsv)
    --targets       Path to targets.tsv       (default: OUT_DIR/targets.tsv)
    --out-dir       Output directory           (default: MS_SEARCH/bam_interrogate)
    --vaf-thresh    Alt-VAF threshold for "positive" call (default: 0.05)
    --max-plots     Max per-peptide bar charts to save (default: 40)
    --data-dir      Root data dir              (default: /scratch/leduc.an/AAS_Evo)
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam-support", default=None)
    p.add_argument("--targets",     default=None)
    p.add_argument("--out-dir",     default=None)
    p.add_argument("--vaf-thresh",  type=float, default=0.05)
    p.add_argument("--max-plots",   type=int,   default=40)
    p.add_argument("--data-dir",    default="/scratch/leduc.an/AAS_Evo")
    return p.parse_args()


def main():
    args = parse_args()

    DATA_DIR = args.data_dir
    OUT_DIR  = args.out_dir or os.path.join(DATA_DIR, "MS_SEARCH", "bam_interrogate")
    os.makedirs(OUT_DIR, exist_ok=True)

    bam_support_path = args.bam_support or os.path.join(OUT_DIR, "bam_support.tsv")
    targets_path     = args.targets     or os.path.join(OUT_DIR, "targets.tsv")

    if not os.path.exists(bam_support_path):
        sys.exit(f"ERROR: {bam_support_path} not found.\n"
                 f"Run interrogate_bams.sh first, then merge with --merge.")

    print(f"Loading {bam_support_path} ...")
    sup = pd.read_csv(bam_support_path, sep="\t", low_memory=False)
    sup["alt_vaf"]     = pd.to_numeric(sup["alt_vaf"],     errors="coerce").fillna(0.0)
    sup["total_depth"] = pd.to_numeric(sup["total_depth"], errors="coerce").fillna(0)
    sup["alt_depth"]   = pd.to_numeric(sup["alt_depth"],   errors="coerce").fillna(0)
    sup["ratio"]       = pd.to_numeric(sup["ratio"],       errors="coerce")
    sup["have_vaf"]    = pd.to_numeric(sup["have_vaf"],    errors="coerce")

    if os.path.exists(targets_path):
        targets = pd.read_csv(targets_path, sep="\t", low_memory=False)
    else:
        targets = None
        print(f"WARNING: {targets_path} not found — summary statistics may be incomplete.")

    VAF_THRESH = args.vaf_thresh
    print(f"  Rows: {len(sup):,}  |  Unique not-have UUIDs: {sup['not_have_uuid'].nunique()}")
    print(f"  Positive threshold: alt_vaf > {VAF_THRESH}")

    # ── 1. Per-peptide bar charts ─────────────────────────────────────────────
    # Group by (gene, swap, peptide, ratio) — one panel per unique PSM/mutation
    group_key = ["gene", "swap", "peptide", "ratio"]
    groups    = list(sup.groupby(group_key, sort=False))
    n_plots   = min(len(groups), args.max_plots)

    # Sort groups by ratio ascending so we see the worst ones first
    groups_sorted = sorted(groups, key=lambda x: x[0][3] if not pd.isna(x[0][3]) else 99)

    print(f"\nGenerating per-peptide bar charts ({n_plots}/{len(groups)} groups)...")
    bar_pdf = os.path.join(OUT_DIR, "bam_support_per_peptide.pdf")

    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(bar_pdf) as pdf:
        for i, ((gene, swap, peptide, ratio), grp) in enumerate(groups_sorted[:n_plots]):
            grp = grp.copy().sort_values("not_have_channel")
            x   = range(len(grp))
            colors = ["#d62728" if v > VAF_THRESH else "#aec7e8"
                      for v in grp["alt_vaf"]]

            fig, ax = plt.subplots(figsize=(max(6, len(grp) * 0.7 + 2), 3.5))
            bars = ax.bar(x, grp["alt_vaf"], color=colors, edgecolor="white", linewidth=0.5)

            # Annotate depth on bars
            for b, (_, row) in zip(bars, grp.iterrows()):
                depth = int(row["total_depth"])
                if depth > 0:
                    ax.text(b.get_x() + b.get_width() / 2,
                            b.get_height() + 0.005,
                            f"d={depth}", ha="center", va="bottom",
                            fontsize=6, color="grey")

            ax.set_xticks(list(x))
            ax.set_xticklabels(
                [f"{r['not_have_channel']}\n{r['not_have_case']}"
                 for _, r in grp.iterrows()],
                fontsize=7, rotation=45, ha="right")
            ax.set_ylabel("Alt VAF in not-have BAM", fontsize=9)
            ax.axhline(VAF_THRESH, color="#d62728", linestyle="--", linewidth=0.8,
                       label=f"threshold ({VAF_THRESH})")
            ax.set_ylim(0, max(grp["alt_vaf"].max() * 1.3, VAF_THRESH * 3))

            have_vaf_str = (f"{grp['have_vaf'].iloc[0]:.3f}"
                            if pd.notna(grp["have_vaf"].iloc[0]) else "n/a")
            ax.set_title(
                f"{gene}  {swap}  |  peptide: {peptide}\n"
                f"MS ratio (have/not-have) = {ratio:.3f}   "
                f"have VAF = {have_vaf_str}",
                fontsize=9)

            red_patch  = mpatches.Patch(color="#d62728", label=f"alt_vaf > {VAF_THRESH}")
            blue_patch = mpatches.Patch(color="#aec7e8", label=f"alt_vaf ≤ {VAF_THRESH}")
            ax.legend(handles=[red_patch, blue_patch], fontsize=7, loc="upper right")

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    print(f"  → Saved: {bar_pdf}")

    # ── 2. Summary: fraction of low-ratio PSMs with ≥1 not-have hit ──────────
    print("\nSummary statistics:")

    per_psm = sup.groupby(["gene", "swap", "peptide", "ratio"]).agg(
        n_not_have        = ("not_have_uuid", "nunique"),
        n_positive        = ("alt_vaf", lambda v: (v > VAF_THRESH).sum()),
        max_alt_vaf       = ("alt_vaf", "max"),
        mean_alt_vaf      = ("alt_vaf", "mean"),
        mean_depth        = ("total_depth", "mean"),
        have_vaf          = ("have_vaf", "first"),
    ).reset_index()

    per_psm["any_positive"] = per_psm["n_positive"] > 0

    n_psms   = len(per_psm)
    n_pos    = per_psm["any_positive"].sum()
    frac_pos = n_pos / n_psms if n_psms else 0

    print(f"  Low-ratio PSMs interrogated: {n_psms}")
    print(f"  PSMs with ≥1 not-have alt_vaf > {VAF_THRESH}: {n_pos} ({frac_pos:.1%})")
    print(f"  Max alt_vaf in any not-have patient: {sup['alt_vaf'].max():.4f}")

    per_psm_path = os.path.join(OUT_DIR, "bam_support_per_psm_summary.tsv")
    per_psm.to_csv(per_psm_path, sep="\t", index=False)
    print(f"  → Per-PSM summary: {per_psm_path}")

    # ── 3. Distribution of alt_vaf ────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    ax.hist(sup["alt_vaf"], bins=50, color="#4878cf", edgecolor="white", linewidth=0.3)
    ax.axvline(VAF_THRESH, color="#d62728", linestyle="--", linewidth=1,
               label=f"threshold ({VAF_THRESH})")
    ax.set_xlabel("Alt VAF in not-have BAM", fontsize=10)
    ax.set_ylabel("Count (positions)", fontsize=10)
    ax.set_title("Distribution of alt_vaf\nacross all not-have positions", fontsize=10)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.hist(per_psm["max_alt_vaf"], bins=40, color="#6acc65", edgecolor="white", linewidth=0.3)
    ax.axvline(VAF_THRESH, color="#d62728", linestyle="--", linewidth=1)
    ax.set_xlabel("Max alt_vaf per PSM (across all not-have patients)", fontsize=10)
    ax.set_ylabel("Count (PSMs)", fontsize=10)
    ax.set_title("Max alt_vaf per low-ratio PSM", fontsize=10)

    fig.suptitle(
        f"Not-have BAM alt-read support | "
        f"{n_pos}/{n_psms} PSMs ({frac_pos:.1%}) have ≥1 not-have with alt_vaf > {VAF_THRESH}",
        fontsize=11)
    fig.tight_layout()
    dist_pdf = os.path.join(OUT_DIR, "bam_support_distribution.pdf")
    fig.savefig(dist_pdf)
    plt.close(fig)
    print(f"  → Distribution plot: {dist_pdf}")

    # ── 4. Scatter: MS ratio vs max_alt_vaf in not-have ──────────────────────
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(
        per_psm["max_alt_vaf"], per_psm["ratio"],
        c=per_psm["have_vaf"], cmap="RdYlGn_r",
        vmin=0, vmax=1,
        s=20, alpha=0.7, edgecolors="none")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Have-patient VAF", fontsize=9)
    ax.axvline(VAF_THRESH, color="#d62728", linestyle="--", linewidth=0.8,
               label=f"not-have VAF threshold ({VAF_THRESH})")
    ax.set_xlabel("Max alt_vaf in not-have BAMs", fontsize=10)
    ax.set_ylabel("MS enrichment ratio (have / not-have)", fontsize=10)
    ax.set_title("Do not-have patients carry the variant at DNA level?", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    scatter_pdf = os.path.join(OUT_DIR, "bam_support_scatter.pdf")
    fig.savefig(scatter_pdf)
    plt.close(fig)
    print(f"  → Scatter plot:      {scatter_pdf}")

    # ── 5. Depth coverage check ───────────────────────────────────────────────
    low_depth = (sup["total_depth"] < 10).sum()
    pct_low   = low_depth / len(sup) * 100 if len(sup) else 0
    print(f"\nDepth check: {low_depth}/{len(sup)} positions ({pct_low:.1f}%) "
          f"have total_depth < 10 (inconclusive)")

    print("\nDone.")
    print(f"  Outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
