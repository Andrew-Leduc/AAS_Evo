#!/usr/bin/env python3
"""
Backfill pLDDT column into existing ddg_matrix TSV files.
Reads pLDDT from cached AlphaFold PDB B-factors (no SPURS model needed).
Skips files that already have the plddt column.
"""

import argparse
import os
from pathlib import Path


def extract_plddt(pdb_path: Path) -> list:
    plddt = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            res_num = int(line[22:26].strip())
            if res_num in seen:
                continue
            seen.add(res_num)
            plddt.append(float(line[60:66].strip()))
    return plddt


def backfill(ddg_path: Path, pdb_cache: Path):
    acc = ddg_path.name.split(".")[0]
    pdb_path = pdb_cache / f"AF-{acc}-F1.pdb"

    if not pdb_path.exists():
        print(f"[SKIP] {ddg_path.name}: no cached PDB")
        return False

    # Check if already has plddt column
    with open(ddg_path) as f:
        header = f.readline()
    if "plddt" in header:
        print(f"[OK]   {ddg_path.name}: already has plddt")
        return True

    plddt = extract_plddt(pdb_path)

    # Read existing file
    lines = ddg_path.read_text().splitlines()
    if len(lines) < 2:
        print(f"[SKIP] {ddg_path.name}: empty")
        return False

    n_data = len(lines) - 1  # excluding header
    if len(plddt) != n_data:
        print(f"[WARN] {ddg_path.name}: pLDDT length {len(plddt)} != rows {n_data}, using NA")

    # Rebuild with plddt inserted after wt_aa (col index 1)
    tmp = ddg_path.with_suffix(".tmp")
    with open(tmp, "w") as out:
        # Header: pos_1based, wt_aa, plddt, to_A, to_C, ...
        cols = header.rstrip("\n").split("\t")
        cols.insert(2, "plddt")
        out.write("\t".join(cols) + "\n")

        for i, line in enumerate(lines[1:]):
            fields = line.split("\t")
            pl = f"{plddt[i]:.2f}" if i < len(plddt) else "NA"
            fields.insert(2, pl)
            out.write("\t".join(fields) + "\n")

    tmp.replace(ddg_path)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ddg-dir",   default="/scratch/leduc.an/AAS_Evo/SPURS/ddg_matrices")
    ap.add_argument("--pdb-cache", default="/scratch/leduc.an/AAS_Evo/SPURS/pdb_cache")
    args = ap.parse_args()

    ddg_dir   = Path(args.ddg_dir)
    pdb_cache = Path(args.pdb_cache)

    files = sorted(ddg_dir.glob("*.ddg_matrix.tsv"))
    print(f"Files to process: {len(files)}", flush=True)

    done = skip = warn = 0
    for i, f in enumerate(files):
        if i % 200 == 0:
            print(f"  {i}/{len(files)}...", flush=True)
        ok = backfill(f, pdb_cache)
        if ok:
            done += 1
        else:
            skip += 1

    print(f"\nDone: {done}  Skipped: {skip}")


if __name__ == "__main__":
    main()
