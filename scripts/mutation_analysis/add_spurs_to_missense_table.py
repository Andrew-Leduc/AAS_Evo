#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

AA3_TO_1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
    "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
    "Tyr":"Y","Val":"V","Ter":"*",
}
HGVSP_PATTERN = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

def parse_mutation(row):
    """
    Returns (pos_1based:int, mut_aa:str) or None.
    Primary: Amino_acids + Protein_position (present in all_missense_mutations.tsv).
    Backup: HGVSp parse (p.Ile300Val etc).
    """
    aa_swap = (row.get("Amino_acids") or "").strip()
    pos_str = (row.get("Protein_position") or "").strip()

    if "/" in aa_swap and pos_str.isdigit():
        wt, mut = aa_swap.split("/", 1)
        wt = wt.strip()
        mut = mut.strip()
        if len(wt) == 1 and len(mut) == 1:
            return int(pos_str), mut

    hgvsp = (row.get("HGVSp") or "").strip()
    m = HGVSP_PATTERN.search(hgvsp)
    if m:
        alt = AA3_TO_1.get(m.group(3))
        if alt:
            return int(m.group(2)), alt

    return None

def load_ddg_matrix(ddg_path: Path):
    """
    ddg_matrix.tsv format:
      pos_1based wt_aa to_A ... to_Y
    Returns dict: (pos_1based, mut_aa) -> ddg float
    """
    mapping = {}
    with open(ddg_path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        to_aas = [h.replace("to_", "") for h in header[2:]]
        for row in r:
            pos = int(row[0])
            vals = row[2:]
            for aa, v in zip(to_aas, vals):
                try:
                    mapping[(pos, aa)] = float(v)
                except ValueError:
                    mapping[(pos, aa)] = None
    return mapping

def build_gene_to_ddgfile(ddg_dir: Path):
    """
    Files named: ACC.GENE.ddg_matrix.tsv
    Returns dict gene -> Path
    """
    out = {}
    if not ddg_dir.exists():
        return out
    for p in ddg_dir.glob("*.ddg_matrix.tsv"):
        parts = p.name.split(".")
        if len(parts) >= 3:
            gene = parts[1]
            out[gene] = p
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vep-tsv", required=True, help="Input VEP/all_missense_mutations.tsv")
    ap.add_argument("--spurs-dir", required=True, help="Base SPURS dir containing ddg_matrices/")
    ap.add_argument("-o", "--out-tsv", required=True, help="Output TSV with spurs_ddg column")
    args = ap.parse_args()

    vep = Path(args.vep_tsv)
    spurs = Path(args.spurs_dir)
    ddg_dir = spurs / "ddg_matrices"
    out = Path(args.out_tsv)

    gene2file = build_gene_to_ddgfile(ddg_dir)
    ddg_cache = {}

    with open(vep, newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if "spurs_ddg" not in fieldnames:
            fieldnames.append("spurs_ddg")

        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="") as fout:
            writer = csv.DictWriter(fout, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()

            for row in reader:
                gene = (row.get("SYMBOL") or "").strip()
                spurs_ddg = "NA"

                mut = parse_mutation(row)
                if gene and mut and gene in gene2file:
                    pos, mut_aa = mut
                    if gene not in ddg_cache:
                        ddg_cache[gene] = load_ddg_matrix(gene2file[gene])
                    val = ddg_cache[gene].get((pos, mut_aa))
                    if val is not None:
                        spurs_ddg = f"{val:.6f}"

                row["spurs_ddg"] = spurs_ddg
                writer.writerow(row)

    print(f"WROTE: {out}")

if __name__ == "__main__":
    main()
