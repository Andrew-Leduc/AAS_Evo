#!/usr/bin/env python3
"""
SAAP × missense co-occurrence analysis (TMT-set level)
=======================================================
For each (SAAP, TMT set), asks: does any patient in that TMT set carry a
genomic missense mutation in the same gene as the SAAP?

The SAAP swap itself is excluded — if the patient's missense IS the same
substitution at the same position, it's just the DNA mutation expressing
itself, not a co-occurring independent event.

Outputs
-------
saap_missense_cooccurrence.tsv  — one row per (SAAP, TMT set, cooccurring missense)
saap_set_summary.tsv            — one row per (SAAP, TMT set): n_patients_in_set,
                                  n_with_any_missense, n_with_cooccurring_missense
"""

import argparse
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

REPO_DIR  = Path(__file__).resolve().parents[2]
TSOUR_DIR = REPO_DIR / "metadata" / "Tsour_et_al"

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pep-to-protein", default=str(TSOUR_DIR / "pep_to_protein.csv"))
    ap.add_argument("--pep-to-patient", default=str(TSOUR_DIR / "peptide_to_patient.csv"))
    ap.add_argument("--gdc-meta",       default=str(REPO_DIR / "metadata/GDC_meta/gdc_meta_matched.tsv"))
    ap.add_argument("--missense",       required=True)
    ap.add_argument("-o", "--out-dir",  required=True)
    return ap.parse_args()


def make_swap(amino_acids, protein_position):
    """Build swap string e.g. 'R273H' from Amino_acids='R/H' and Protein_position='273'.
    Falls back to None if malformed. Much faster than regex on HGVSp."""
    try:
        ref, alt = str(amino_acids).split("/")
        pos = str(protein_position).split("-")[0]  # handle ranges
        if len(ref) == 1 and len(alt) == 1 and pos.isdigit():
            return f"{ref}{pos}{alt}"
    except Exception:
        pass
    return None


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load Tsour tables ─────────────────────────────────────────────────────
    print("Loading Tsour tables...", flush=True)
    p2prot = pd.read_csv(args.pep_to_protein, index_col=0)
    p2pat  = pd.read_csv(args.pep_to_patient, index_col=0)

    p2prot = p2prot[["SAAP","name","protein.position","fromto"]].drop_duplicates("SAAP")
    p2prot = p2prot.rename(columns={"name":"gene_symbol","protein.position":"saap_pos"})
    p2prot["saap_ref"]  = p2prot["fromto"].str.split(":").str[0]
    p2prot["saap_alt"]  = p2prot["fromto"].str.split(":").str[1]
    p2prot["saap_swap"] = p2prot["saap_ref"] + p2prot["saap_pos"].astype(str) + p2prot["saap_alt"]

    merged = p2pat.merge(p2prot[["SAAP","gene_symbol","saap_pos","saap_swap"]],
                         on="SAAP", how="left")
    merged["case_submitter_id"] = merged["Sample name"].str.rsplit("_", n=1).str[0]
    print(f"  Tsour entries: {len(merged):,} | unique SAAPs: {merged['SAAP'].nunique()}", flush=True)
    print(f"  Unique (Dataset, TMT set) groups: {merged.groupby(['Dataset','TMT set']).ngroups}", flush=True)

    # ── GDC metadata ─────────────────────────────────────────────────────────
    print("Loading GDC metadata...", flush=True)
    gdc = pd.read_csv(args.gdc_meta, sep="\t",
                      usecols=["case_submitter_id","gdc_file_id","sample_type"])
    gdc_tumor = (gdc[gdc["sample_type"] == "Primary Tumor"]
                 .drop_duplicates("case_submitter_id")
                 .rename(columns={"gdc_file_id":"sample_id"}))
    gdc_any = (gdc.drop_duplicates("case_submitter_id")
               .rename(columns={"gdc_file_id":"sample_id"}))
    gdc_map = (pd.concat([gdc_tumor[["case_submitter_id","sample_id"]],
                          gdc_any[["case_submitter_id","sample_id"]]])
               .drop_duplicates("case_submitter_id")
               .set_index("case_submitter_id")["sample_id"].to_dict())

    merged["sample_id"] = merged["case_submitter_id"].map(gdc_map)
    n_matched = merged["sample_id"].notna().sum()
    print(f"  Tsour patients matched to GDC UUID: {n_matched:,} / {len(merged):,}", flush=True)

    # ── Load missense table ───────────────────────────────────────────────────
    print("Loading missense table (this takes ~2-3 min for 20M rows)...", flush=True)
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["sample_id","SYMBOL","HGVSp","Protein_position","Amino_acids","VAF",
                                "am_pathogenicity","spurs_ddg"])
    print(f"  Read {len(miss):,} rows", flush=True)

    print("  Parsing swap strings from Amino_acids + Protein_position...", flush=True)
    miss["swap"] = np.vectorize(make_swap)(miss["Amino_acids"], miss["Protein_position"])
    miss["VAF"]  = pd.to_numeric(miss["VAF"], errors="coerce")
    miss["am_pathogenicity"] = pd.to_numeric(miss["am_pathogenicity"], errors="coerce")
    miss["spurs_ddg"]        = pd.to_numeric(miss["spurs_ddg"], errors="coerce")
    print(f"  Swap parsed: {miss['swap'].notna().sum():,} / {len(miss):,}", flush=True)

    # ── Build index: gene -> list of (sample_id, swap, VAF, am, spurs, HGVSp) ─
    print("Building gene index...", flush=True)
    gene_index = defaultdict(list)
    for row in miss[["sample_id","SYMBOL","swap","VAF","am_pathogenicity","spurs_ddg","HGVSp"]].itertuples(index=False):
        gene_index[row.SYMBOL].append(row)
    print(f"  Indexed {len(gene_index):,} unique genes", flush=True)

    # ── Build set -> frozenset of sample_ids ──────────────────────────────────
    set_patients = (merged[merged["sample_id"].notna()]
                    .groupby(["Dataset","TMT set"])["sample_id"]
                    .apply(frozenset).to_dict())

    # ── CO-OCCURRENCE SEARCH ──────────────────────────────────────────────────
    saap_set_keys = (merged[merged["sample_id"].notna() & merged["gene_symbol"].notna()]
                     [["SAAP","Dataset","TMT set","gene_symbol","saap_swap","saap_pos"]]
                     .drop_duplicates(["SAAP","Dataset","TMT set"]))

    print(f"Searching {len(saap_set_keys):,} (SAAP, TMT set) combinations...", flush=True)

    records = []
    summary = []
    for i, (_, row) in enumerate(saap_set_keys.iterrows()):
        if i % 1000 == 0:
            print(f"  {i:,} / {len(saap_set_keys):,} done, {len(records):,} co-occurrences found",
                  flush=True)

        saap      = row["SAAP"]
        dataset   = row["Dataset"]
        tmt_set   = row["TMT set"]
        gene      = row["gene_symbol"]
        saap_swap = row["saap_swap"]
        saap_pos  = row["saap_pos"]

        set_uuids  = set_patients.get((dataset, tmt_set), frozenset())
        n_patients = len(set_uuids)

        cooc_samples = set()
        for mrow in gene_index.get(gene, []):
            if mrow.sample_id not in set_uuids:
                continue
            if mrow.swap == saap_swap:   # same substitution — exclude
                continue
            cooc_samples.add(mrow.sample_id)
            records.append({
                "SAAP":               saap,
                "Dataset":            dataset,
                "TMT set":            tmt_set,
                "gene":               gene,
                "saap_swap":          saap_swap,
                "saap_pos":           saap_pos,
                "missense_sample_id": mrow.sample_id,
                "missense_swap":      mrow.swap,
                "missense_HGVSp":     mrow.HGVSp,
                "missense_VAF":       mrow.VAF,
                "missense_am":        mrow.am_pathogenicity,
                "missense_spurs_ddg": mrow.spurs_ddg,
            })

        summary.append({
            "SAAP":               saap,
            "Dataset":            dataset,
            "TMT set":            tmt_set,
            "gene":               gene,
            "saap_swap":          saap_swap,
            "n_patients_in_set":  n_patients,
            "n_patients_with_cooccurring_missense": len(cooc_samples),
            "has_cooccurrence":   len(cooc_samples) > 0,
        })

    cooc = pd.DataFrame(records)
    summ = pd.DataFrame(summary)

    n_saap_sets = len(summ)
    n_with_cooc = summ["has_cooccurrence"].sum()
    print(f"\n── RESULTS ──────────────────────────────────────────────────────────", flush=True)
    print(f"Unique (SAAP, TMT set) combinations tested: {n_saap_sets:,}", flush=True)
    print(f"With ≥1 co-occurring missense in same gene:  {n_with_cooc:,} "
          f"({100*n_with_cooc/n_saap_sets:.1f}%)", flush=True)
    print(f"Co-occurrence records (SAAP × missense):     {len(cooc):,}", flush=True)
    if len(cooc):
        print(f"\nTop genes by co-occurrence count:", flush=True)
        print(cooc.groupby("gene")["SAAP"].count().sort_values(ascending=False).head(10).to_string(),
              flush=True)

    summ.to_csv(out_dir / "saap_set_summary.tsv", sep="\t", index=False)
    cooc.to_csv(out_dir / "saap_missense_cooccurrence.tsv", sep="\t", index=False)
    print(f"\nSaved to {out_dir}/", flush=True)


if __name__ == "__main__":
    main()
