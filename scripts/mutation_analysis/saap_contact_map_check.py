#!/usr/bin/env python3
"""
Check structural contact between SAAP position and co-occurring missense position
using AlphaFold-derived Cα distance maps.

For each (SAAP, missense) co-occurrence record where the missense is both
AM-pathogenic and SPURS-destabilizing, looks up the distance between the
SAAP position and the missense position in the protein's contact map.

Contact threshold: Cα-Cα < 8 Å (standard)

Outputs
-------
saap_contact_annotated.tsv   — co-occurrence records with distance + in_contact columns
saap_contact_summary.txt     — summary statistics
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

REPO_DIR  = Path(__file__).resolve().parents[2]
TSOUR_DIR = REPO_DIR / "metadata" / "Tsour_et_al"

AM_THRESHOLD    = 0.564
SPURS_THRESHOLD = 0.5
CONTACT_ANG     = 8.0   # Cα-Cα distance threshold in Ångstroms

COOC_DEFAULT    = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence/saap_missense_cooccurrence.tsv"
CONTACT_DIR     = "/scratch/leduc.an/AAS_Evo/SPURS/contact_maps"
DDG_DIR         = "/scratch/leduc.an/AAS_Evo/SPURS/ddg_matrices"
OUT_DEFAULT     = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cooccurrence", default=COOC_DEFAULT)
    ap.add_argument("--contact-dir",  default=CONTACT_DIR)
    ap.add_argument("--ddg-dir",      default=DDG_DIR)
    ap.add_argument("--threshold",    type=float, default=CONTACT_ANG)
    ap.add_argument("-o", "--out-dir", default=OUT_DEFAULT)
    return ap.parse_args()


def build_gene_to_acc(ddg_dir: Path) -> dict:
    """Build gene -> UniProt accession from ddg_matrix filenames: {ACC}.{GENE}.ddg_matrix.tsv"""
    g2a = {}
    for f in ddg_dir.glob("*.ddg_matrix.tsv"):
        parts = f.name.split(".")
        if len(parts) >= 3:
            g2a[parts[1]] = parts[0]
    return g2a


def load_contact_map(contact_dir: Path, acc: str):
    """
    Load distance map for a UniProt accession.
    Returns (pos_to_idx dict, distance_matrix) or (None, None) if not found.
    pos_to_idx maps 1-based protein position -> 0-based matrix index.
    """
    prefix = contact_dir / f"AF-{acc}-F1"
    csv_path = Path(str(prefix) + ".csv")
    npy_path = Path(str(prefix) + ".npy")
    if not csv_path.exists() or not npy_path.exists():
        return None, None
    try:
        meta = pd.read_csv(csv_path, index_col=0)
        dm   = np.load(npy_path)
        # id column is 1-based residue position, row index is 0-based matrix index
        pos_to_idx = {int(row["id"]): idx for idx, row in meta.iterrows()
                      if pd.notna(row["id"])}
        return pos_to_idx, dm
    except Exception as e:
        print(f"  [WARN] Failed to load {acc}: {e}")
        return None, None


def main():
    args = parse_args()
    contact_dir = Path(args.contact_dir)
    ddg_dir     = Path(args.ddg_dir)
    out_dir     = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Build gene -> accession map ───────────────────────────────────────────
    print("Building gene -> accession map...", flush=True)
    gene_to_acc = build_gene_to_acc(ddg_dir)
    print(f"  {len(gene_to_acc):,} genes with ddg matrices", flush=True)

    # ── Count available contact maps ──────────────────────────────────────────
    n_maps = len(list(contact_dir.glob("*.npy")))
    print(f"  {n_maps:,} contact maps available in {contact_dir}", flush=True)

    # ── Load co-occurrence, filter to both-destab ─────────────────────────────
    print("Loading co-occurrence table...", flush=True)
    cooc = pd.read_csv(args.cooccurrence, sep="\t")
    cooc["missense_am"]        = pd.to_numeric(cooc["missense_am"],        errors="coerce")
    cooc["missense_spurs_ddg"] = pd.to_numeric(cooc["missense_spurs_ddg"], errors="coerce")

    both = cooc[(cooc["missense_am"]        >= AM_THRESHOLD) |
                (cooc["missense_spurs_ddg"] >= SPURS_THRESHOLD)].copy()
    print(f"  AM or SPURS destabilizing pairs: {len(both):,}", flush=True)

    # Extract missense position from missense_swap (e.g. R273H -> 273)
    def parse_pos(swap):
        import re
        if not isinstance(swap, str): return None
        m = re.match(r"[A-Z](\d+)[A-Z]", swap)
        return int(m.group(1)) if m else None

    both["missense_pos"] = both["missense_swap"].apply(parse_pos)
    both["saap_pos"]     = pd.to_numeric(both["saap_pos"], errors="coerce")

    # ── Check contact for each pair ───────────────────────────────────────────
    print("Checking structural contacts...", flush=True)
    cm_cache = {}  # acc -> (pos_to_idx, dm)

    distances   = []
    in_contacts = []
    has_map     = []

    n_no_acc = n_no_map = n_no_pos = n_checked = 0

    for i, (_, row) in enumerate(both.iterrows()):
        if i % 1000 == 0:
            print(f"  {i:,} / {len(both):,} pairs checked", flush=True)

        gene        = row["gene"]
        saap_pos    = row["saap_pos"]
        miss_pos    = row["missense_pos"]

        if pd.isna(saap_pos) or pd.isna(miss_pos):
            distances.append(np.nan); in_contacts.append(False); has_map.append(False)
            n_no_pos += 1; continue

        saap_pos = int(saap_pos)
        miss_pos = int(miss_pos)

        acc = gene_to_acc.get(gene)
        if not acc:
            distances.append(np.nan); in_contacts.append(False); has_map.append(False)
            n_no_acc += 1; continue

        if acc not in cm_cache:
            cm_cache[acc] = load_contact_map(contact_dir, acc)
        pos_to_idx, dm = cm_cache[acc]

        if dm is None:
            distances.append(np.nan); in_contacts.append(False); has_map.append(False)
            n_no_map += 1; continue

        idx_saap = pos_to_idx.get(saap_pos)
        idx_miss = pos_to_idx.get(miss_pos)

        if idx_saap is None or idx_miss is None:
            distances.append(np.nan); in_contacts.append(False); has_map.append(True)
            n_no_pos += 1; continue

        dist = float(dm[idx_saap, idx_miss])
        distances.append(dist)
        in_contacts.append(dist < args.threshold)
        has_map.append(True)
        n_checked += 1

    both["structural_distance_ang"] = distances
    both["in_structural_contact"]   = in_contacts
    both["has_contact_map"]         = has_map

    # ── Summary ───────────────────────────────────────────────────────────────
    testable   = both["has_contact_map"] & both["structural_distance_ang"].notna()
    n_contact  = both["in_structural_contact"].sum()

    lines = [
        "── STRUCTURAL CONTACT: SAAP × AM-OR-SPURS-DESTAB MISSENSE ────────────",
        f"AM or SPURS destabilizing pairs:  {len(both):,}",
        f"  No gene->accession mapping:     {n_no_acc:,}",
        f"  No contact map file:            {n_no_map:,}",
        f"  Position not in map:            {n_no_pos:,}",
        f"  Distance computed:              {n_checked:,}",
        f"",
        f"In structural contact (<{args.threshold}Å Cα-Cα): "
        f"{n_contact:,} / {testable.sum():,} testable "
        f"({100*n_contact/max(testable.sum(),1):.1f}%)",
        f"",
        f"Distance distribution (testable pairs):",
        f"  mean   = {both.loc[testable,'structural_distance_ang'].mean():.2f} Å",
        f"  median = {both.loc[testable,'structural_distance_ang'].median():.2f} Å",
        f"  <8Å    = {(both.loc[testable,'structural_distance_ang'] < 8).sum():,}",
        f"  <12Å   = {(both.loc[testable,'structural_distance_ang'] < 12).sum():,}",
        f"  <20Å   = {(both.loc[testable,'structural_distance_ang'] < 20).sum():,}",
        "",
        "Top contact pairs (shortest distance):",
    ]
    top = (both[testable & both["in_structural_contact"]]
           .sort_values("structural_distance_ang")
           [["gene","saap_swap","missense_swap","structural_distance_ang",
             "missense_am","missense_spurs_ddg"]]
           .drop_duplicates(["gene","saap_swap","missense_swap"])
           .head(20))
    for _, r in top.iterrows():
        lines.append(f"  {r['gene']:<12} SAAP={r['saap_swap']:<10} "
                     f"miss={r['missense_swap']:<10} "
                     f"dist={r['structural_distance_ang']:.2f}Å  "
                     f"AM={r['missense_am']:.3f}  SPURS={r['missense_spurs_ddg']:.3f}")

    summary = "\n".join(lines)
    print(f"\n{summary}", flush=True)
    (out_dir / "saap_contact_summary.txt").write_text(summary)
    both.to_csv(out_dir / "saap_contact_annotated.tsv", sep="\t", index=False)
    print(f"\nOutputs saved to {out_dir}/", flush=True)


if __name__ == "__main__":
    main()
