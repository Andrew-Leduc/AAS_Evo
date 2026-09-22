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
    df = df[np.isfinite(df["rt_shift"])]
    df["expected"] = df["alt"].map(RP_RANK) - df["wt"].map(RP_RANK)
    df = df.dropna(subset=["expected"])
    print(f"same-run pairs with RT: {len(df):,}  ({df.groupby(['acc','pos','wt','alt']).ngroups:,} unique SAAPs)")

    # per-pair correlation
    rho, p = stats.spearmanr(df["expected"], df["rt_shift"])
    r, pr = stats.pearsonr(df["expected"], df["rt_shift"])
    print(f"per-pair : Spearman rho={rho:+.3f} (p={p:.1e}) | Pearson r={r:+.3f} (p={pr:.1e})")

    # ── facet by SAAP identification confidence ──
    if "saap_prob" in df.columns and df["saap_prob"].notna().any():
        print("\nRT correlation faceted by SAAP confidence (PeptideProphet prob / 1-PEP):")
        conf_bins = [(0.0, 0.5, "low <0.5"), (0.5, 0.9, "mid 0.5-0.9"),
                     (0.9, 0.99, "high 0.9-0.99"), (0.99, 1.01, "top >=0.99")]
        for lo, hi, lab in conf_bins:
            sub = df[(df["saap_prob"] >= lo) & (df["saap_prob"] < hi)]
            if len(sub) >= 20:
                rr, pp = stats.spearmanr(sub["expected"], sub["rt_shift"])
                print(f"  {lab:14s} n={len(sub):6,d}  Spearman rho={rr:+.3f} (p={pp:.1e})")
            else:
                print(f"  {lab:14s} n={len(sub):6,d}  (too few)")
    else:
        conf_bins = None
        print("(no saap_prob column — re-run precursor_saap_bp_pairs.py to add confidence)")

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
        ncol = 3 if conf_bins else 2
        fig, ax = plt.subplots(1, ncol, figsize=(6.5 * ncol, 5.5))
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

        if conf_bins:
            labs, rhos, ns = [], [], []
            for lo, hi, lab in conf_bins:
                sub = df[(df["saap_prob"] >= lo) & (df["saap_prob"] < hi)]
                if len(sub) >= 20:
                    rr, _ = stats.spearmanr(sub["expected"], sub["rt_shift"])
                    labs.append(lab.split()[0]); rhos.append(rr); ns.append(len(sub))
            ax[2].bar(range(len(labs)), rhos, color="#3a8489")
            ax[2].axhline(0, color="k", lw=.6)
            for i, (rr, nn) in enumerate(zip(rhos, ns)):
                ax[2].text(i, rr, f"n={nn}", ha="center",
                           va="bottom" if rr >= 0 else "top", fontsize=8)
            ax[2].set_xticks(range(len(labs))); ax[2].set_xticklabels(labs)
            ax[2].set_xlabel("SAAP confidence bin")
            ax[2].set_ylabel("Spearman ρ (obs vs expected RT)")
            ax[2].set_title("per-pair ρ by ID confidence")
        plt.tight_layout()
        plt.savefig(f"{a.out}/rt_shift_obs_vs_expected.png", dpi=200, bbox_inches="tight")
        print(f"wrote {a.out}/rt_shift_obs_vs_expected.png + rt_shift_by_swaptype.tsv")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
