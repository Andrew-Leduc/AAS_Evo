#!/usr/bin/env python3
"""
Build an ESM-MSR (github.com/SKTeamLab/esm-msr) input CSV from our EVE
double-mutant request, so ESM-MSR can score each detected (missense_i, AAS_j)
double mutant's stability ddG + epistasis (dddG) as a fast structural proxy
while the real EVE double-mutant predictions run.

Each row: pdb_file, code, chain, mut_type
  mut_type = "{miss_wt}{miss_pos}{miss_alt}:{aas_wt}{aas_pos}{aas_alt}"  (colon = double)
  pdb_file = the cached AlphaFold structure for that UniProt accession.

Run ESM-MSR on the output, e.g.:
  python src/esm_msr/inference.py --checkpoint_path LoRA_models/esm-msr-small/<ckpt>.ckpt \
      --input_csv esm_msr_input.csv --mode singles+doubles --output_csv esm_msr_scores.csv
(do NOT pass --skip_reverse for doubles — it drops epistasis).
"""

import argparse
import glob
import os
from pathlib import Path
import pandas as pd


def find_pdb(cache, acc):
    """Locate a cached AlphaFold structure for this accession (prefer plain .pdb)."""
    for ext in ("pdb", "cif"):
        for v in ("6", "5", "4"):
            p = Path(cache) / f"AF-{acc}-F1-model_v{v}.{ext}"
            if p.exists() and p.stat().st_size > 0:
                return str(p)
    hits = sorted(glob.glob(os.path.join(cache, f"AF-{acc}-F1-model_v*.pdb")))
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser()
    scratch = "/scratch/leduc.an/AAS_Evo"
    ap.add_argument("--request",
                    default=f"{scratch}/ANALYSIS/contact_saap/eve_double_mutant_request.tsv")
    ap.add_argument("--pdb-cache", default=f"{scratch}/SPURS/pdb_cache")
    ap.add_argument("-o", "--out",
                    default=f"{scratch}/ANALYSIS/contact_saap/esm_msr_input.csv")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--significant-only", action="store_true",
                    help="keep only rows with q < 0.05 (needs a q column in the request)")
    ap.add_argument("--singles", action="store_true",
                    help="emit one row per UNIQUE single mutation (missense and AAS "
                         "singles, mut_type has no ':') instead of the double mutants; "
                         "score with --mode singles to enable double-minus-single analyses")
    args = ap.parse_args()

    req = pd.read_csv(args.request, sep="\t")
    if args.significant_only and "q" in req.columns:
        req = req[pd.to_numeric(req["q"], errors="coerce") < 0.05]

    rows, no_pdb, seen = [], set(), set()
    for r in req.itertuples(index=False):
        acc = r.uniprot
        pdb = find_pdb(args.pdb_cache, acc)
        if not pdb:
            no_pdb.add(acc)
            continue
        miss = f"{r.missense_wt}{int(r.missense_pos)}{r.missense_alt}"
        aas = f"{r.aas_wt}{int(r.aas_pos)}{r.aas_alt}"
        if args.singles:
            for mt in (miss, aas):
                if (acc, mt) not in seen:
                    seen.add((acc, mt))
                    rows.append({"pdb_file": pdb, "code": acc, "chain": args.chain,
                                 "mut_type": mt})
        else:
            rows.append({"pdb_file": pdb, "code": acc, "chain": args.chain,
                         "mut_type": f"{miss}:{aas}"})

    out = pd.DataFrame(rows).drop_duplicates()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    kind = "single" if args.singles else "double-mutant"
    print(f"{len(out):,} {kind} rows -> {args.out}")
    print(f"proteins: {out['code'].nunique() if len(out) else 0} | "
          f"accessions with no cached PDB: {len(no_pdb)}")
    if no_pdb:
        print("  missing PDB (first 10):", sorted(no_pdb)[:10],
              "\n  (re-download via make_contact_maps.py --pdb-dir, or point --pdb-cache elsewhere)")


if __name__ == "__main__":
    main()
