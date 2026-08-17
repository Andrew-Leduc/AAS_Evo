#!/usr/bin/env python3
"""
Pin ESM-MSR's ddG sign convention by correlating its single-mutant predictions
against SPURS ddG (SPURS convention: positive = destabilizing).

  Spearman(ESM, SPURS) < 0  -> ESM higher = more STABLE (opposite to SPURS)
                               => negative double_minus_aas = destabilizing missense
                               => "DESTABILIZING missense -> more AAS in carriers"
  Spearman(ESM, SPURS) > 0  -> ESM higher = more DESTABILIZING (same as SPURS)
                               => "STABILIZING missense -> more AAS in carriers"

Usage:
  python3 check_esm_spurs_sign.py
"""
import argparse
import glob
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--singles-scores", default=f"{BASE}/esm_msr_singles_scores.csv")
    ap.add_argument("--ddg-dir", default="/scratch/leduc.an/AAS_Evo/SPURS/ddg_matrices")
    ap.add_argument("--esm-col", default="wt_lora_pred")
    a = ap.parse_args()

    ss = pd.read_csv(a.singles_scores).rename(columns={"code": "uniprot"})
    mt = "mut_type_pdb" if "mut_type_pdb" in ss.columns else "mut_type_renumbered"
    ss = ss.rename(columns={mt: "mut_type"})
    ss = ss[~ss["mut_type"].astype(str).str.contains(":")].copy()   # singles only

    p = ss["mut_type"].astype(str).str.extract(r"^([A-Z])(\d+)([A-Z])$")
    ss["pos"] = pd.to_numeric(p[1], errors="coerce")
    ss["alt"] = p[2]
    ss = ss.dropna(subset=["pos", "alt"]); ss["pos"] = ss["pos"].astype(int)

    cache = {}
    def get_ddg(acc, pos, alt):
        if acc not in cache:
            f = glob.glob(os.path.join(a.ddg_dir, f"{acc}.*.ddg_matrix.tsv"))
            cache[acc] = (pd.read_csv(f[0], sep="\t").set_index("pos_1based")
                          if f else None)
        d = cache[acc]
        if d is None or pos not in d.index:
            return np.nan
        col = f"to_{alt}"
        return d.loc[pos, col] if col in d.columns else np.nan

    ss["spurs_ddg"] = [get_ddg(u, po, al)
                       for u, po, al in zip(ss["uniprot"], ss["pos"], ss["alt"])]
    ss["esm"] = pd.to_numeric(ss[a.esm_col], errors="coerce")
    ok = ss["spurs_ddg"].notna() & ss["esm"].notna()
    n = int(ok.sum())
    print(f"matched {n:,} singles with both ESM ({a.esm_col}) and SPURS ddG")
    if n < 10:
        raise SystemExit("too few matched -- are the SPURS ddg_matrices present? "
                         f"(ls {a.ddg_dir}/*.ddg_matrix.tsv)")
    rho, pv = spearmanr(ss.loc[ok, "esm"], ss.loc[ok, "spurs_ddg"])
    print(f"Spearman(ESM {a.esm_col}, SPURS ddG) = {rho:+.3f} (p={pv:.2e})\n")
    if rho < -0.05:
        print("ESM higher = MORE STABLE (opposite sign to SPURS).")
        print("  -> negative double_minus_aas = destabilizing missense")
        print("  -> HEADLINE: DESTABILIZING missense -> more AAS in carriers")
    elif rho > 0.05:
        print("ESM higher = MORE DESTABILIZING (same sign as SPURS).")
        print("  -> negative double_minus_aas = stabilizing missense")
        print("  -> HEADLINE: STABILIZING missense -> more AAS in carriers")
    else:
        print("no clear ESM/SPURS relationship -- convention not pinned this way.")


if __name__ == "__main__":
    main()
