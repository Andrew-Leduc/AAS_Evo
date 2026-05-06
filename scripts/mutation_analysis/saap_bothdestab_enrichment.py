#!/usr/bin/env python3
"""
SAAP enrichment in patients with both-destabilizing missense mutations
======================================================================
For each SAAP that co-occurs (same gene, same TMT set) with a missense mutation
predicted destabilizing by BOTH AlphaMissense (AM>=0.564) AND SPURS (ddG>=0.5),
test whether that SAAP is enriched in the missense-carrying patients vs other
patients who also have the SAAP but no such doubly-destabilizing missense.

Enrichment tested two ways:
  1. Detection frequency: fraction of patients in the TMT set with the SAAP
  2. RAAS: relative AAS abundance (mut/canonical peptide ratio)

Outputs
-------
saap_bothdestab_enrichment.tsv   — per (SAAP, Dataset, TMT set) statistics
saap_bothdestab_summary.txt      — console summary
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

REPO_DIR  = Path(__file__).resolve().parents[2]
TSOUR_DIR = REPO_DIR / "metadata" / "Tsour_et_al"

AM_THRESHOLD    = 0.564
SPURS_THRESHOLD = 0.5

COOC_DEFAULT = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence/saap_missense_cooccurrence.tsv"
GDC_META     = str(REPO_DIR / "metadata" / "GDC_meta" / "gdc_meta_matched.tsv")
OUT_DEFAULT  = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cooccurrence", default=COOC_DEFAULT)
    ap.add_argument("--gdc-meta",     default=GDC_META)
    ap.add_argument("-o", "--out-dir", default=OUT_DEFAULT)
    return ap.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load Tsour tables ─────────────────────────────────────────────────────
    print("Loading Tsour et al. tables...", flush=True)
    p2prot = pd.read_csv(TSOUR_DIR / "pep_to_protein.csv",
                         usecols=["SAAP", "name", "protein.position"])
    p2prot = p2prot.rename(columns={"name": "gene", "protein.position": "saap_pos"})

    p2pat = pd.read_csv(TSOUR_DIR / "peptide_to_patient.csv",
                        usecols=["SAAP", "Sample name", "Dataset", "TMT set",
                                 "RAAS", "SAAP abundance", "BP abundance"])
    p2pat["RAAS"] = pd.to_numeric(p2pat["RAAS"], errors="coerce")
    p2pat = p2pat[np.isfinite(p2pat["RAAS"])]

    # Sample name is like "C3L-00001_tumor" — extract case_submitter_id
    p2pat["case_id"] = p2pat["Sample name"].str.rsplit("_", n=1).str[0]

    print(f"  p2prot: {len(p2prot):,} rows | p2pat: {len(p2pat):,} rows", flush=True)

    # ── GDC UUID -> case_submitter_id map ─────────────────────────────────────
    print("Loading GDC metadata...", flush=True)
    gdc = pd.read_csv(args.gdc_meta, sep="\t",
                      usecols=["gdc_file_id", "case_submitter_id"])
    uuid_to_case = gdc.set_index("gdc_file_id")["case_submitter_id"].to_dict()

    # ── Load co-occurrence, filter to both-destab ─────────────────────────────
    print("Loading co-occurrence table...", flush=True)
    cooc = pd.read_csv(args.cooccurrence, sep="\t")
    cooc["missense_am"]        = pd.to_numeric(cooc["missense_am"],        errors="coerce")
    cooc["missense_spurs_ddg"] = pd.to_numeric(cooc["missense_spurs_ddg"], errors="coerce")

    both = cooc[(cooc["missense_am"]        >= AM_THRESHOLD) &
                (cooc["missense_spurs_ddg"] >= SPURS_THRESHOLD)].copy()
    print(f"  Total co-occurrence records:    {len(cooc):,}", flush=True)
    print(f"  Both AM+SPURS destabilizing:    {len(both):,}", flush=True)
    print(f"  Unique (SAAP, Dataset, TMT set): "
          f"{both.groupby(['SAAP','Dataset','TMT_set']).ngroups:,}", flush=True)

    # Map missense_sample_id (GDC UUID) -> case_submitter_id
    both["missense_case_id"] = both["missense_sample_id"].map(uuid_to_case)

    # ── Build enrichment per (SAAP, Dataset, TMT set) ─────────────────────────
    print("\nRunning enrichment analysis...", flush=True)

    # Index p2pat for fast lookup
    p2pat_idx = p2pat.set_index(["SAAP", "Dataset", "TMT set"])

    results = []
    groups = both.groupby(["SAAP", "Dataset", "TMT_set"])
    n_tested = 0

    for (saap, dataset, tmt_set), grp in groups:
        n_tested += 1
        if n_tested % 500 == 0:
            print(f"  {n_tested}/{len(groups)} groups tested...", flush=True)

        # Case IDs that carry the both-destab missense in this SAAP's gene
        missense_cases = set(grp["missense_case_id"].dropna().unique())
        if not missense_cases:
            continue

        # All patients who have this SAAP detected (in this dataset+set)
        try:
            saap_patients = p2pat_idx.loc[(saap, dataset, tmt_set)].reset_index(drop=True)
        except KeyError:
            # SAAP not in peptide_to_patient for this dataset/set
            continue

        if len(saap_patients) == 0:
            continue

        saap_patients = saap_patients.copy()
        saap_patients["has_missense"] = saap_patients["case_id"].isin(missense_cases)

        with_miss    = saap_patients[saap_patients["has_missense"]]["RAAS"].dropna()
        without_miss = saap_patients[~saap_patients["has_missense"]]["RAAS"].dropna()

        # Total patients in this TMT set (denominator for detection frequency)
        all_cases_in_set = set(
            p2pat[(p2pat["Dataset"] == dataset) &
                  (p2pat["TMT set"] == tmt_set)]["case_id"].unique()
        )
        n_set_total = len(all_cases_in_set)

        # Detection frequency: fraction of TMT-set patients who have this SAAP
        n_saap_with    = saap_patients["has_missense"].sum()
        n_saap_without = (~saap_patients["has_missense"]).sum()
        freq_with    = n_saap_with    / max(len(missense_cases),  1)
        freq_without = n_saap_without / max(n_set_total - len(missense_cases), 1)

        # RAAS comparison
        if len(with_miss) >= 2 and len(without_miss) >= 2:
            _, p_raas = stats.mannwhitneyu(with_miss, without_miss, alternative="greater")
        else:
            p_raas = np.nan

        gene = grp["gene"].iloc[0] if "gene" in grp.columns else ""

        results.append({
            "SAAP":                    saap,
            "gene":                    gene,
            "Dataset":                 dataset,
            "TMT_set":                 tmt_set,
            "n_missense_cases":        len(missense_cases),
            "n_saap_with_missense":    int(n_saap_with),
            "n_saap_without_missense": int(n_saap_without),
            "n_set_total":             n_set_total,
            "freq_with_missense":      round(freq_with,    4),
            "freq_without_missense":   round(freq_without, 4),
            "freq_ratio":              round(freq_with / freq_without, 4) if freq_without > 0 else np.nan,
            "mean_raas_with":          round(with_miss.mean(),    4) if len(with_miss) else np.nan,
            "mean_raas_without":       round(without_miss.mean(), 4) if len(without_miss) else np.nan,
            "raas_diff":               round(with_miss.mean() - without_miss.mean(), 4)
                                       if (len(with_miss) and len(without_miss)) else np.nan,
            "p_raas_mw":               p_raas,
            "missense_am":             grp["missense_am"].mean(),
            "missense_spurs_ddg":      grp["missense_spurs_ddg"].mean(),
        })

    out = pd.DataFrame(results)
    print(f"\nTotal (SAAP, Dataset, TMT set) groups tested: {len(out):,}", flush=True)

    if len(out) == 0:
        print("No results — check input files.", flush=True)
        return

    # ── Summary ───────────────────────────────────────────────────────────────
    testable_raas   = out["p_raas_mw"].notna()
    sig_raas        = testable_raas & (out["p_raas_mw"] < 0.05)
    enriched_freq   = out["freq_ratio"].notna() & (out["freq_ratio"] > 1)

    lines = [
        "── SAAP ENRICHMENT IN BOTH-DESTAB MISSENSE PATIENTS ──────────────────",
        f"(SAAP, Dataset, TMT set) groups tested:          {len(out):,}",
        f"With RAAS data for comparison:                   {testable_raas.sum():,}",
        f"Higher RAAS in missense patients (p<0.05 MW):   {sig_raas.sum():,} "
        f"({100*sig_raas.sum()/max(testable_raas.sum(),1):.1f}%)",
        f"Higher detection frequency in missense patients: {enriched_freq.sum():,} "
        f"({100*enriched_freq.sum()/max(len(out),1):.1f}%)",
        "",
        f"Mean RAAS — with missense:    {out['mean_raas_with'].mean():.4f}",
        f"Mean RAAS — without missense: {out['mean_raas_without'].mean():.4f}",
        f"Mean freq ratio (with/without): {out['freq_ratio'].mean():.3f}",
        "",
        "Top genes by mean RAAS enrichment (with vs without missense):",
    ]
    top = (out.dropna(subset=["raas_diff"])
              .groupby("gene")["raas_diff"]
              .mean()
              .sort_values(ascending=False)
              .head(15))
    for gene, diff in top.items():
        lines.append(f"  {gene:<20} mean RAAS diff = {diff:.4f}")

    summary = "\n".join(lines)
    print(f"\n{summary}", flush=True)
    (out_dir / "saap_bothdestab_summary.txt").write_text(summary)

    out.to_csv(out_dir / "saap_bothdestab_enrichment.tsv", sep="\t", index=False)
    print(f"\nOutputs saved to {out_dir}/", flush=True)


if __name__ == "__main__":
    main()
