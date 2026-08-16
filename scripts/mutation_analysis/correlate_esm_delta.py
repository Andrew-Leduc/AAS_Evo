#!/usr/bin/env python3
"""
Correlate ESM-MSR predictions against the observed carrier-vs-non-carrier RAAS
delta. Tests several ESM "predictors", including single-subtracted variants that
use the individual single-mutant scores (singles+doubles mode also scores the
singles -- rows whose mut_type has no ':').

Predictors (per (missense i, AAS j) pair):
  combined_dddg_pred : epistasis  = ddG(double) - ddG(miss) - ddG(aas)   [both singles]
  double_minus_aas   : ddG(double) - ddG(aas single)   = "carrier(both) - AAS-alone"
                       (the user's framing: AAS-bearing protein, missense present vs not)
  double_minus_miss  : ddG(double) - ddG(missense single) = AAS effect in missense bg
  combined_pred      : ddG(double)  (total)

Usage:
  python3 correlate_esm_delta.py [--delta-col delta_within] [--score-col combined_pred]
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

BASE = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--request", default=f"{BASE}/eve_double_mutant_request_all.tsv")
    ap.add_argument("--scores",  default=f"{BASE}/esm_msr_scores_all.csv")
    ap.add_argument("--out",     default=f"{BASE}/esm_vs_delta.pdf")
    ap.add_argument("--delta-col", default="delta")
    ap.add_argument("--score-col", default="combined_pred",
                    help="ddG column used for the double/single subtractions")
    ap.add_argument("--min-carrier", type=int, default=3)
    ap.add_argument("--min-noncarrier", type=int, default=2)
    a = ap.parse_args()

    req = pd.read_csv(a.request, sep="\t")
    sco = pd.read_csv(a.scores).rename(columns={"code": "uniprot"})
    mtcol = "mut_type_pdb" if "mut_type_pdb" in sco.columns else "mut_type_renumbered"
    sco = sco.rename(columns={mtcol: "mut_type"})

    is_double = sco["mut_type"].astype(str).str.contains(":")
    doubles = sco[is_double].copy()
    singles = sco[~is_double].copy()
    sing = singles.set_index(["uniprot", "mut_type"])[a.score_col].to_dict()
    print(f"scores: {len(doubles):,} double rows, {len(singles):,} single rows")

    def parts(mt):
        p = str(mt).split(":")
        return (p[0], p[1]) if len(p) == 2 else (None, None)
    doubles["miss_mt"], doubles["aas_mt"] = zip(*doubles["mut_type"].map(parts))

    def sub(row, which):
        s = sing.get((row["uniprot"], row[which]))
        return (row[a.score_col] - s) if s is not None and pd.notna(s) else np.nan
    doubles["double_minus_aas"]  = doubles.apply(lambda r: sub(r, "aas_mt"),  axis=1)
    doubles["double_minus_miss"] = doubles.apply(lambda r: sub(r, "miss_mt"), axis=1)
    print(f"  aas-single matched: {doubles['double_minus_aas'].notna().sum():,} | "
          f"miss-single matched: {doubles['double_minus_miss'].notna().sum():,}")

    req["mut_type"] = (req["missense_wt"] + req["missense_pos"].astype(int).astype(str)
                       + req["missense_alt"] + ":" + req["aas_wt"]
                       + req["aas_pos"].astype(int).astype(str) + req["aas_alt"])
    m = (req.merge(doubles, on=["uniprot", "mut_type"], how="inner")
            .drop_duplicates(["uniprot", "mut_type"]))
    print(f"merged {len(m):,} of {len(req):,} request pairs\n")

    predictors = [
        ("combined_dddg_pred", "ESM epistasis ddddG (double - both singles)"),
        ("double_minus_aas",   "ESM double - AAS single  (carrier both - AAS alone)"),
        ("double_minus_miss",  "ESM double - missense single"),
        ("combined_pred",      "ESM total double ddG"),
    ]
    reliable = (m["n_carrier"] >= a.min_carrier) & (m["n_noncarrier"] >= a.min_noncarrier)
    d = pd.to_numeric(m[a.delta_col], errors="coerce")

    fig, axes = plt.subplots(len(predictors), 2, figsize=(12, 3.6 * len(predictors)))
    print(f"=== observed {a.delta_col} vs ESM predictors ===")
    for r, (pcol, plabel) in enumerate(predictors):
        for c, (slabel, mask) in enumerate(
                [("all pairs", pd.Series(True, index=m.index)),
                 (f"n_carrier>={a.min_carrier}", reliable)]):
            ax = axes[r, c]
            x = pd.to_numeric(m.loc[mask, pcol], errors="coerce"); y = d[mask]
            ok = x.notna() & y.notna(); x, y = x[ok], y[ok]
            ax.scatter(x, y, s=8, alpha=0.3)
            ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
            ax.set_xlabel(plabel); ax.set_ylabel(f"obs {a.delta_col}")
            if len(x) > 2:
                rho, pr = spearmanr(x, y); rr, _ = pearsonr(x, y)
                ax.set_title(f"{slabel}  n={len(x)}  Spearman={rho:+.3f} (p={pr:.1e})")
                print(f"  {pcol:20s} {slabel:16s} n={len(x):5d}  "
                      f"Spearman={rho:+.3f} (p={pr:.2e})  Pearson={rr:+.3f}")
    plt.tight_layout()
    plt.savefig(a.out)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
