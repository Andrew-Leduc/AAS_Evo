#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
import urllib.request
from pathlib import Path

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

def read_gene_list(path: str):
    genes = []
    with open(path, "r") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g.split()[0])  # allow extra columns
    return genes

def parse_uniprot_canonical_fasta(fasta_path: str):
    """
    Returns dict: gene_symbol -> uniprot_accession (first occurrence).
    Headers in uniprot_human_canonical.fasta contain GN=<SYMBOL>.
    """
    gene2acc = {}
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.startswith(">"):
                continue
            parts = line.split("|")
            acc = parts[1] if len(parts) >= 3 else line[1:].split()[0]
            m = re.search(r"GN=(\S+)", line)
            gene = m.group(1) if m else None
            if gene and gene not in gene2acc:
                gene2acc[gene] = acc
    return gene2acc

def download_alphafold_pdb(uniprot_acc: str, out_pdb: Path, timeout: int = 60):
    """
    AlphaFold DB (EBI): AF-<ACC>-F1-model_v{4,3,2,1}.pdb
    """
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    if out_pdb.exists() and out_pdb.stat().st_size > 1000:
        return True

    for v in (4, 3, 2, 1):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v{v}.pdb"
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                data = r.read()
            if data and len(data) > 1000:
                out_pdb.write_bytes(data)
                return True
        except Exception:
            continue
    return False

def write_ddg_matrix_tsv(ddg_tensor, seq: str, out_tsv: Path):
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    ddg = ddg_tensor.detach().cpu().numpy()
    with open(out_tsv, "w") as f:
        f.write("pos_1based\twt_aa\t" + "\t".join([f"to_{aa}" for aa in ALPHABET]) + "\n")
        for i, wt in enumerate(seq):
            f.write(f"{i+1}\t{wt}\t" + "\t".join([f"{ddg[i, j]:.6f}" for j in range(20)]) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-list", required=True, help="Path to gene_list.txt (one gene per line)")
    ap.add_argument("--gene-index", type=int, default=None,
                    help="1-based line index for SLURM arrays. If omitted, runs all genes.")
    ap.add_argument("--ref-fasta", required=True, help="uniprot_human_canonical.fasta")
    ap.add_argument("--out-dir", required=True, help="Base SPURS output dir (e.g., /scratch/leduc.an/AAS_Evo/SPURS)")
    ap.add_argument("--chain", default="A", help="Chain to use (AlphaFold PDBs typically chain A)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    pdb_cache = out_dir / "pdb_cache"
    ddg_dir = out_dir / "ddg_matrices"
    hf_cache = out_dir / "hf_cache"
    hf_cache.mkdir(parents=True, exist_ok=True)

    os.environ.setdefault("HF_HOME", str(hf_cache))
    os.environ.setdefault("TRANSFORMERS_CACHE", str(hf_cache))

    genes = read_gene_list(args.gene_list)
    if args.gene_index is not None:
        i = args.gene_index
        if i < 1 or i > len(genes):
            print(f"ERROR: gene-index {i} out of range 1..{len(genes)}", file=sys.stderr)
            sys.exit(1)
        genes = [genes[i - 1]]

    gene2acc = parse_uniprot_canonical_fasta(args.ref_fasta)

    import torch
    from spurs.inference import get_SPURS_from_hub, parse_pdb

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, cfg = get_SPURS_from_hub()
    model = model.to(device).eval()

    for gene in genes:
        acc = gene2acc.get(gene)
        if not acc:
            print(f"[SKIP] {gene}: not found in ref FASTA GN= mapping")
            continue

        out_tsv = ddg_dir / f"{acc}.{gene}.ddg_matrix.tsv"
        if out_tsv.exists() and out_tsv.stat().st_size > 1000:
            print(f"[OK] {gene} ({acc}) already exists")
            continue

        pdb_path = pdb_cache / f"AF-{acc}-F1.pdb"
        if not download_alphafold_pdb(acc, pdb_path):
            print(f"[SKIP] {gene} ({acc}): AlphaFold PDB not found")
            continue

        pdb_name = f"AF-{acc}-F1"
        pdb = parse_pdb(str(pdb_path), pdb_name, args.chain, cfg)
        seq = pdb["seq"]

        with torch.no_grad():
            ddg = model(pdb, return_logist=True)  # (L,20)

        wt1 = seq[0]
        wt_idx = ALPHABET.find(wt1)
        wt_ddg = ddg[0, wt_idx].item() if wt_idx >= 0 else float("nan")
        print(f"[RUN] {gene} ({acc}) L={len(seq)} wt1={wt1} ddg(wt1)={wt_ddg:.6f}")

        write_ddg_matrix_tsv(ddg, seq, out_tsv)
        print(f"[WROTE] {out_tsv}")

if __name__ == "__main__":
    main()
