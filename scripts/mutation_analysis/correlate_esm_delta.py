#!/usr/bin/env python3
"""
Correlate ESM-MSR double-mutant predictions against the observed carrier-vs-
non-carrier RAAS delta.

Joins eve_double_mutant_request_all.tsv (observed delta per (missense, AAS) pair)
to esm_msr_scores_all.csv on (uniprot, mut_type), where mut_type is rebuilt exactly
as build_esm_msr_input.py wrote it: "{miss}:{aas}".

ESM columns of interest:
  combined_dddg_pred : ensemble EPISTASIS (ddddG; departure from additivity) -- the
                       analog of "does the missense make the AAS more/less favorable"
  combined_pred      : total double-mutant ddG (additive + epistasis)

Usage:
  python3 correlate_esm_delta.py           # uses default paths + delta column
  python3 correlate_esm_delta.py --delta-col delta_within
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
    ap.add_argument("--delta-col", default="delta",
                    help="delta, delta_within, or delta_pooled")
    ap.add_argument("--min-carrier", type=int, default=3)
    ap.add_argument("--min-noncarrier", type=int, default=2)
    a = ap.parse_args()

    req = pd.read_csv(a.request, sep="\t")
    sco = pd.read_csv(a.scores).rename(columns={"code": "uniprot"})

    req["mut_type"] = (req["missense_wt"] + req["missense_pos"].astype(int).astype(str)
                       + req["missense_alt"] + ":" + req["aas_wt"]
                       + req["aas_pos"].astype(int).astype(str) + req["aas_alt"])

    # ESM may renumber; join on whichever mut_type column matches more rows
    merged, used = None, None
    for mtcol in ("mut_type_pdb", "mut_type_renumbered"):
        if mtcol not in sco.columns:
            continue
        s = sco.rename(columns={mtcol: "mut_type"})
        mm = req.merge(s, on=["uniprot", "mut_type"], how="inner", suffixes=("", "_sco"))
        if merged is None or len(mm) > len(merged):
            merged, used = mm, mtcol
    if merged is None or len(merged) == 0:
        raise SystemExit("no rows joined -- check mut_type / uniprot keys")
    m = merged.drop_duplicates(["uniprot", "mut_type"])
    print(f"merged {len(m):,} of {len(req):,} request pairs (join via {used})")

    d_all = pd.to_numeric(m[a.delta_col], errors="coerce")
    reliable = (m["n_carrier"] >= a.min_carrier) & (m["n_noncarrier"] >= a.min_noncarrier)
    panels = [("combined_dddg_pred", "ESM combined ddddG (epistasis)"),
              ("combined_pred",      "ESM combined ddG (double mutant)")]
    subsets = [("all pairs", pd.Series(True, index=m.index)),
               (f"n_carrier>={a.min_carrier} & n_noncarrier>={a.min_noncarrier}", reliable)]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    print("\n=== correlations (observed %s vs ESM prediction) ===" % a.delta_col)
    for r, (pcol, plabel) in enumerate(panels):
        for c, (slabel, mask) in enumerate(subsets):
            ax = axes[c, r]
            x = pd.to_numeric(m.loc[mask, pcol], errors="coerce")
            y = d_all[mask]
            ok = x.notna() & y.notna()
            x, y = x[ok], y[ok]
            ax.scatter(x, y, s=8, alpha=0.3)
            ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
            ax.set_xlabel(plabel); ax.set_ylabel(f"observed {a.delta_col}")
            if len(x) > 2:
                rho, pr = spearmanr(x, y); rr, pp = pearsonr(x, y)
                ax.set_title(f"{slabel}\nn={len(x)}  Spearman={rho:.3f} (p={pr:.1e})  "
                             f"Pearson={rr:.3f}")
                print(f"  {pcol:20s} | {slabel:45s} n={len(x):5d}  "
                      f"Spearman={rho:+.3f} (p={pr:.2e})  Pearson={rr:+.3f}")
            else:
                ax.set_title(f"{slabel}\nn={len(x)} (too few)")
    plt.tight_layout()
    plt.savefig(a.out)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
