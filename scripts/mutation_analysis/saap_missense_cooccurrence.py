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
import re
import pandas as pd
import numpy as np
from pathlib import Path

# ── DEFAULTS ──────────────────────────────────────────────────────────────────
REPO_DIR  = Path(__file__).resolve().parents[2]
TSOUR_DIR = REPO_DIR / "metadata" / "Tsour_et_al"

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pep-to-protein",  default=str(TSOUR_DIR / "pep_to_protein.csv"))
    ap.add_argument("--pep-to-patient",  default=str(TSOUR_DIR / "peptide_to_patient.csv"))
    ap.add_argument("--gdc-meta",        default=str(REPO_DIR / "metadata/GDC_meta/gdc_meta_matched.tsv"))
    ap.add_argument("--missense",        required=True,
                    help="Path to all_missense_with_spurs.tsv (or all_missense_mutations.tsv)")
    ap.add_argument("-o", "--out-dir",   required=True)
    return ap.parse_args()


def parse_hgvsp_swap(hgvsp):
    """'ENSP00000.4:p.Arg273His' -> 'R273H'. Returns None if unparseable."""
    AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
           "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
           "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
    m = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', str(hgvsp))
    if m:
        ref = AA3.get(m.group(1))
        alt = AA3.get(m.group(3))
        if ref and alt:
            return f"{ref}{m.group(2)}{alt}"
    return None


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load Tsour tables ─────────────────────────────────────────────────────
    p2prot = pd.read_csv(args.pep_to_protein, index_col=0)
    p2pat  = pd.read_csv(args.pep_to_patient, index_col=0)

    # Build SAAP -> gene_symbol, protein_position, saap_swap (e.g. "P45N")
    # fromto is "P:N", protein.position is 45
    p2prot = p2prot[["SAAP","name","protein.position","fromto"]].drop_duplicates("SAAP")
    p2prot = p2prot.rename(columns={"name":"gene_symbol", "protein.position":"saap_pos"})
    p2prot["saap_ref"] = p2prot["fromto"].str.split(":").str[0]
    p2prot["saap_alt"] = p2prot["fromto"].str.split(":").str[1]
    p2prot["saap_swap"] = p2prot["saap_ref"] + p2prot["saap_pos"].astype(str) + p2prot["saap_alt"]

    # Join patient table to protein info
    merged = p2pat.merge(p2prot[["SAAP","gene_symbol","saap_pos","saap_swap"]],
                         on="SAAP", how="left")

    # Parse case_submitter_id from Sample name (e.g. "C3L-01603_Tumor" -> "C3L-01603")
    merged["case_submitter_id"] = merged["Sample name"].str.rsplit("_", n=1).str[0]
    merged["tsour_sample_type"] = merged["Sample name"].str.rsplit("_", n=1).str[1]

    print(f"Tsour entries: {len(merged):,} | unique SAAPs: {merged['SAAP'].nunique()}")
    print(f"Unique (Dataset, TMT set) groups: {merged.groupby(['Dataset','TMT set']).ngroups}")

    # ── Load GDC metadata: case_submitter_id -> gdc_file_id ──────────────────
    gdc = pd.read_csv(args.gdc_meta, sep="\t",
                      usecols=["case_submitter_id","gdc_file_id","sample_type"])
    # Keep one UUID per patient (prefer Primary Tumor; any will do for gene lookup)
    gdc_tumor = (gdc[gdc["sample_type"] == "Primary Tumor"]
                 .drop_duplicates("case_submitter_id")
                 .rename(columns={"gdc_file_id":"sample_id"}))
    gdc_any   = (gdc.drop_duplicates("case_submitter_id")
                 .rename(columns={"gdc_file_id":"sample_id"}))
    # Use tumor UUID first, fall back to any
    gdc_map = (pd.concat([gdc_tumor[["case_submitter_id","sample_id"]],
                           gdc_any[["case_submitter_id","sample_id"]]])
               .drop_duplicates("case_submitter_id")
               .set_index("case_submitter_id")["sample_id"]
               .to_dict())

    merged["sample_id"] = merged["case_submitter_id"].map(gdc_map)
    n_matched = merged["sample_id"].notna().sum()
    print(f"Tsour patients matched to GDC UUID: {n_matched:,} / {len(merged):,}")

    # ── Load missense table ───────────────────────────────────────────────────
    print("Loading missense table...")
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["sample_id","SYMBOL","HGVSp","Protein_position","Amino_acids","VAF"])
    miss["VAF"] = pd.to_numeric(miss["VAF"], errors="coerce")
    miss["swap"] = miss["HGVSp"].apply(parse_hgvsp_swap)
    print(f"Missense rows: {len(miss):,}")

    # Build set-level patient -> UUID map
    # A TMT set is identified by (Dataset, TMT set); get all sample_ids in it
    set_patients = (merged[merged["sample_id"].notna()]
                    .groupby(["Dataset","TMT set"])["sample_id"]
                    .apply(set).to_dict())

    # ── CO-OCCURRENCE SEARCH ──────────────────────────────────────────────────
    # For each unique (SAAP, Dataset, TMT set): get gene_symbol + saap_swap,
    # find all missense mutations in that gene for any patient in the set,
    # exclude those that match the SAAP swap exactly.
    saap_set_keys = (merged[merged["sample_id"].notna() & merged["gene_symbol"].notna()]
                     [["SAAP","Dataset","TMT set","gene_symbol","saap_swap","saap_pos"]]
                     .drop_duplicates(["SAAP","Dataset","TMT set"]))

    records = []
    summary = []

    for _, row in saap_set_keys.iterrows():
        saap       = row["SAAP"]
        dataset    = row["Dataset"]
        tmt_set    = row["TMT set"]
        gene       = row["gene_symbol"]
        saap_swap  = row["saap_swap"]
        saap_pos   = row["saap_pos"]

        set_uuids = set_patients.get((dataset, tmt_set), set())
        n_patients = len(set_uuids)

        # All missense in this gene for any patient in the set
        gene_miss = miss[(miss["SYMBOL"] == gene) & (miss["sample_id"].isin(set_uuids))]

        # Exclude exact same swap (same position + same substitution = same as SAAP source)
        gene_miss = gene_miss[gene_miss["swap"] != saap_swap]

        n_with_miss = gene_miss["sample_id"].nunique()

        summary.append({
            "SAAP":          saap,
            "Dataset":       dataset,
            "TMT set":       tmt_set,
            "gene":          gene,
            "saap_swap":     saap_swap,
            "n_patients_in_set":          n_patients,
            "n_patients_with_cooccurring_missense": n_with_miss,
            "has_cooccurrence": n_with_miss > 0,
        })

        for _, mrow in gene_miss.iterrows():
            records.append({
                "SAAP":        saap,
                "Dataset":     dataset,
                "TMT set":     tmt_set,
                "gene":        gene,
                "saap_swap":   saap_swap,
                "saap_pos":    saap_pos,
                "missense_sample_id": mrow["sample_id"],
                "missense_swap":      mrow["swap"],
                "missense_HGVSp":     mrow["HGVSp"],
                "missense_VAF":       mrow["VAF"],
            })

    cooc = pd.DataFrame(records)
    summ = pd.DataFrame(summary)

    n_saap_sets = len(summ)
    n_with_cooc = summ["has_cooccurrence"].sum()
    print(f"\n── RESULTS ──────────────────────────────────────────────────────────")
    print(f"Unique (SAAP, TMT set) combinations tested: {n_saap_sets:,}")
    print(f"With ≥1 co-occurring missense in same gene:  {n_with_cooc:,} "
          f"({100*n_with_cooc/n_saap_sets:.1f}%)")
    print(f"Co-occurrence records (SAAP × missense):     {len(cooc):,}")
    if len(cooc):
        print(f"\nTop genes by co-occurrence count:")
        print(cooc.groupby("gene")["SAAP"].count().sort_values(ascending=False).head(10).to_string())

    summ.to_csv(out_dir / "saap_set_summary.tsv", sep="\t", index=False)
    cooc.to_csv(out_dir / "saap_missense_cooccurrence.tsv", sep="\t", index=False)
    print(f"\nSaved to {out_dir}/")


if __name__ == "__main__":
    main()
