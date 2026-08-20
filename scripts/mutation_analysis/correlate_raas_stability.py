#!/usr/bin/env python3
"""
Two analyses:
  (1) compare overall RAAS between carrier and non-carrier (per-swap medians)
  (2) correlate the RAAS level against the ESM single-mutant stability score of the
      AAS itself (does a more stable/tolerated substitution accumulate more?)

ESM sign convention (esm-msr repo): score > 0 = stabilizing, < 0 = destabilizing.
Hypothesis: more stable AAS (higher aas ddG) -> higher RAAS (positive correlation).

Inputs:
  per_swap_raas_summary.tsv  (from the notebook RAW RAAS SANITY cell): gene, swap,
      median, carrier_median, noncarrier_median
  eve_double_mutant_request_all.tsv : gene<->uniprot<->swap map
  esm_msr_singles_scores.csv : ESM single-mutant scores (wt_lora_pred)
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, wilcoxon

B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raas", default=f"{B}/per_swap_raas_summary.tsv")
    ap.add_argument("--request", default=f"{B}/eve_double_mutant_request_all.tsv")
    ap.add_argument("--singles", default=f"{B}/esm_msr_singles_scores.csv")
    ap.add_argument("--score-col", default="wt_lora_pred")
    ap.add_argument("--out", default=f"{B}/raas_vs_aas_stability.pdf")
    a = ap.parse_args()

    raas = pd.read_csv(a.raas, sep="\t")
    req = pd.read_csv(a.request, sep="\t")
    g2u = req.drop_duplicates(["gene", "swap"])[["gene", "swap", "uniprot"]]
    raas = raas.merge(g2u, on=["gene", "swap"], how="left")

    sco = pd.read_csv(a.singles).rename(columns={"code": "uniprot"})
    mt = "mut_type_pdb" if "mut_type_pdb" in sco.columns else "mut_type_renumbered"
    sco = sco.rename(columns={mt: "swap"})
    sco = sco[~sco["swap"].astype(str).str.contains(":")]          # singles only
    sco = sco[["uniprot", "swap", a.score_col]].rename(columns={a.score_col: "aas_ddg"})

    m = raas.merge(sco, on=["uniprot", "swap"], how="inner").dropna(subset=["aas_ddg"])
    print(f"{len(m):,} swaps with RAAS + ESM AAS-single score\n")

    # ── (1) overall carrier vs non-carrier RAAS ──────────────────────────────
    both = m.dropna(subset=["carrier_median", "noncarrier_median"])
    print("=== (1) overall RAAS: carrier vs non-carrier (per-swap medians) ===")
    print(f"  non-carrier median log2 RAAS : {both['noncarrier_median'].median():.3f}")
    print(f"  carrier     median log2 RAAS : {both['carrier_median'].median():.3f}")
    if len(both) > 10:
        w, pw = wilcoxon(both["carrier_median"], both["noncarrier_median"])
        print(f"  paired Wilcoxon (carrier vs non): p={pw:.2e} "
              f"(median delta = {(both['carrier_median'] - both['noncarrier_median']).median():+.3f})")

    # ── (2) RAAS level vs ESM AAS-single stability ───────────────────────────
    print("\n=== (2) RAAS level vs ESM AAS-single ddG (>0 = stabilizing) ===")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, col, lab in zip(axes,
                            ["noncarrier_median", "carrier_median", "median"],
                            ["non-carrier RAAS (AAS baseline)", "carrier RAAS", "overall RAAS"]):
        d = m.dropna(subset=[col])
        x, y = d["aas_ddg"], d[col]
        rho, pr = spearmanr(x, y); rr, _ = pearsonr(x, y)
        print(f"  {lab:32s} n={len(d):5d}  Spearman={rho:+.3f} (p={pr:.2e})  Pearson={rr:+.3f}")
        ax.scatter(x, y, s=8, alpha=0.3, color="#3a7")
        ax.axvline(0, color="#bbb", lw=0.5)
        ax.set_xlabel("ESM AAS single ddG  (>0 = stabilizing)")
        ax.set_ylabel(f"{lab}  [log2 RAAS]")
        ax.set_title(f"{lab}\nSpearman={rho:+.3f} (p={pr:.1e})")
    plt.tight_layout()
    plt.savefig(a.out)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
