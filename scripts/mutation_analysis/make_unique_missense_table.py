#!/usr/bin/env python3
"""
Collapse all_missense_mutations.tsv to one row per unique mutation (gene + swap),
adding spurs_ddg, plddt, and per-mutation summary stats (n_samples, mean_vaf, gnomAD AF).

Output is much lighter than the full per-sample table and suitable for
downstream analysis and ranking.

Output columns:
  SYMBOL, uniprot_acc, swap, Protein_position, ref_aa, alt_aa,
  HGVSp, gnomADe_AF, am_pathogenicity, am_class,
  spurs_ddg, plddt, n_samples, mean_vaf, median_vaf
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def build_gene_to_ddgfile(ddg_dir: Path):
    out = {}
    for p in ddg_dir.glob("*.ddg_matrix.tsv"):
        parts = p.name.split(".")
        if len(parts) >= 3:
            gene, acc = parts[1], parts[0]
            out[gene] = (p, acc)
    return out


def load_ddg_and_plddt(ddg_path: Path):
    """Returns (ddg_dict, plddt_dict) keyed by (pos, mut_aa) and pos respectively."""
    ddg = {}
    plddt = {}
    with open(ddg_path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        has_plddt = "plddt" in header
        # to_X columns start after pos_1based, wt_aa, [plddt]
        data_start = 3 if has_plddt else 2
        to_aas = [h.replace("to_", "") for h in header[data_start:]]
        for row in r:
            pos = int(row[0])
            if has_plddt:
                try:
                    plddt[pos] = float(row[2])
                except ValueError:
                    pass
            for aa, v in zip(to_aas, row[data_start:]):
                try:
                    ddg[(pos, aa)] = float(v)
                except ValueError:
                    pass
    return ddg, plddt


def make_swap(amino_acids, protein_position):
    try:
        ref, alt = str(amino_acids).split("/")
        pos = str(protein_position).split("-")[0]
        if len(ref) == 1 and len(alt) == 1 and pos.isdigit():
            return ref, int(pos), alt, f"{ref}{pos}{alt}"
    except Exception:
        pass
    return None, None, None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--missense", required=True,
                    help="all_missense_mutations.tsv (with or without spurs_ddg column)")
    ap.add_argument("--spurs-dir", required=True,
                    help="SPURS base dir containing ddg_matrices/")
    ap.add_argument("-o", "--out", required=True,
                    help="Output unique_missense_mutations.tsv")
    args = ap.parse_args()

    ddg_dir = Path(args.spurs_dir) / "ddg_matrices"
    gene2ddg = build_gene_to_ddgfile(ddg_dir)
    print(f"ddg_matrix files available: {len(gene2ddg)}", flush=True)

    # First pass: collect per-mutation stats
    # key: (SYMBOL, swap)  value: dict of aggregated fields
    mutations = {}
    ddg_cache = {}

    print("Reading missense table...", flush=True)
    n_rows = 0
    with open(args.missense, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n_rows += 1
            if n_rows % 1_000_000 == 0:
                print(f"  {n_rows:,} rows processed, {len(mutations):,} unique mutations", flush=True)

            gene = (row.get("SYMBOL") or "").strip()
            ref_aa, pos, alt_aa, swap = make_swap(
                row.get("Amino_acids", ""), row.get("Protein_position", ""))
            if not swap or not gene:
                continue

            key = (gene, swap)
            if key not in mutations:
                # Look up SPURS ddg + pLDDT on first encounter
                spurs_ddg = plddt_val = None
                uniprot_acc = ""
                if gene in gene2ddg:
                    ddg_path, uniprot_acc = gene2ddg[gene]
                    if gene not in ddg_cache:
                        ddg_cache[gene] = load_ddg_and_plddt(ddg_path)
                    ddg_map, plddt_map = ddg_cache[gene]
                    spurs_ddg  = ddg_map.get((pos, alt_aa))
                    plddt_val  = plddt_map.get(pos)

                mutations[key] = {
                    "SYMBOL":          gene,
                    "uniprot_acc":     uniprot_acc,
                    "swap":            swap,
                    "Protein_position": pos,
                    "ref_aa":          ref_aa,
                    "alt_aa":          alt_aa,
                    "HGVSp":           row.get("HGVSp", ""),
                    "gnomADe_AF":      row.get("gnomADe_AF", ""),
                    "am_pathogenicity": row.get("am_pathogenicity", ""),
                    "am_class":        row.get("am_class", ""),
                    "spurs_ddg":       f"{spurs_ddg:.6f}" if spurs_ddg is not None else "NA",
                    "plddt":           f"{plddt_val:.2f}" if plddt_val is not None else "NA",
                    "_n_samples":      0,
                    "_vafs":           [],
                }

            mutations[key]["_n_samples"] += 1
            try:
                mutations[key]["_vafs"].append(float(row.get("VAF", "")))
            except (ValueError, TypeError):
                pass

    print(f"Total rows: {n_rows:,}  Unique mutations: {len(mutations):,}", flush=True)

    # Write output
    fieldnames = [
        "SYMBOL", "uniprot_acc", "swap", "Protein_position", "ref_aa", "alt_aa",
        "HGVSp", "gnomADe_AF", "am_pathogenicity", "am_class",
        "spurs_ddg", "plddt", "n_samples", "mean_vaf", "median_vaf",
    ]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for m in mutations.values():
            vafs = m["_vafs"]
            m["n_samples"] = m["_n_samples"]
            m["mean_vaf"]   = f"{sum(vafs)/len(vafs):.4f}" if vafs else "NA"
            sorted_vafs = sorted(vafs)
            n = len(sorted_vafs)
            m["median_vaf"] = f"{sorted_vafs[n//2]:.4f}" if sorted_vafs else "NA"
            writer.writerow(m)

    print(f"Written: {out}", flush=True)


if __name__ == "__main__":
    main()
