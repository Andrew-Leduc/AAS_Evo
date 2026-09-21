#!/usr/bin/env python3
"""
#2  Observed vs expected RT shift for SAAP/BP pairs in the SAME MS run.

Input: contact_saap/precursor_pairs.tsv (from precursor_saap_bp_pairs.py) — each
row is a SAAP and its base peptide (BP) identified in the same run, so offline
fractionation is controlled and the observed RT difference is a within-run,
online-gradient shift.

  observed  = saap_rt - bp_rt            (the rt_shift column, seconds)
  expected  = rp_rank[alt] - rp_rank[wt] (ordinal RP-LC retention change of the
                                          single substituted residue)

The RP retention scale is ORDINAL (elution order W>F>L>I>...>D>K>R), so the test
is a Spearman rank correlation of observed vs expected — exact coefficients don't
matter, only the direction/ordering. A real substitution should move RT in the
predicted direction; artifacts / mislocalized IDs should not.

(Upgrade path: swap `expected` for MSBooster/DIA-NN predicted-RT delta from the
_edited.pin files for a calibrated prediction — see --pred-rt, not yet wired.)
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats

B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"

# ordinal RP-LC retention rank (most retained -> least), low-pH TFA order
RP_RANK = {aa: i for i, aa in enumerate(reversed(
    list("RKDEHNQSGTAPCYVMILFW")))}   # W=19 (most retained) ... R=0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=f"{B}/precursor_pairs.tsv")
    ap.add_argument("--out", default=B)
    a = ap.parse_args()

    df = pd.read_csv(a.pairs, sep="\t")
    df = df.dropna(subset=["rt_shift", "wt", "alt"]).copy()
    df["expected"] = df["alt"].map(RP_RANK) - df["wt"].map(RP_RANK)
    df = df.dropna(subset=["expected"])
    print(f"same-run pairs with RT: {len(df):,}  ({df.groupby(['acc','pos','wt','alt']).ngroups:,} unique SAAPs)")

    # per-pair correlation
    rho, p = stats.spearmanr(df["expected"], df["rt_shift"])
    r, pr = stats.pearsonr(df["expected"], df["rt_shift"])
    print(f"per-pair : Spearman rho={rho:+.3f} (p={p:.1e}) | Pearson r={r:+.3f} (p={pr:.1e})")

    # per-swap-type means (reduces per-PSM noise)
    g = (df.groupby(["wt", "alt"])
           .agg(expected=("expected", "first"),
                mean_obs=("rt_shift", "mean"),
                median_obs=("rt_shift", "median"),
                n=("rt_shift", "size"))
           .reset_index())
    g = g[g["n"] >= 5]
    rho2, p2 = stats.spearmanr(g["expected"], g["mean_obs"])
    print(f"per-swap : Spearman rho={rho2:+.3f} (p={p2:.1e})  (swap types with n>=5: {len(g)})")
    g.sort_values("expected").to_csv(f"{a.out}/rt_shift_by_swaptype.tsv", sep="\t", index=False)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
        ax[0].hexbin(df["expected"], df["rt_shift"], gridsize=30, cmap="viridis",
                     mincnt=1)
        ax[0].axhline(0, color="w", lw=.7); ax[0].axvline(0, color="w", lw=.7)
        ax[0].set_xlabel("expected RT shift (RP rank Δ, alt−wt)")
        ax[0].set_ylabel("observed RT shift (s), saap−bp")
        ax[0].set_title(f"per-pair (n={len(df):,})  ρ={rho:+.3f}")
        sc = ax[1].scatter(g["expected"], g["mean_obs"], s=g["n"] / g["n"].max() * 200 + 8,
                           c=g["median_obs"], cmap="coolwarm")
        for _, row in g.iterrows():
            ax[1].annotate(f"{row['wt']}{row['alt']}", (row["expected"], row["mean_obs"]),
                           fontsize=6, alpha=.6)
        ax[1].axhline(0, color="k", lw=.5); ax[1].axvline(0, color="k", lw=.5)
        ax[1].set_xlabel("expected RT shift (RP rank Δ)")
        ax[1].set_ylabel("mean observed RT shift (s)")
        ax[1].set_title(f"per-swap-type (n≥5)  ρ={rho2:+.3f}")
        plt.tight_layout()
        plt.savefig(f"{a.out}/rt_shift_obs_vs_expected.png", dpi=200, bbox_inches="tight")
        print(f"wrote {a.out}/rt_shift_obs_vs_expected.png + rt_shift_by_swaptype.tsv")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
