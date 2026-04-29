#!/usr/bin/env python3
"""
SAAP × EVCouplings × missense co-occurrence analysis
=====================================================
1. Reports which SAAPs in our Tsour list have EVChits files (coverage)
2. For SAAPs that co-occur with a missense in the same TMT set, checks
   whether the missense position is evolutionarily coupled to the SAAP
   position according to the EVChits files.

EVChits file format (per ENSP):
  i, A_i  — one coupled position
  j, A_j  — other coupled position
  aas_pos — SAAP position (1-based)
  aas_wt, aas_mut — wildtype / mutant AA at SAAP site

A missense at position P is "evolutionarily linked" to a SAAP at position Q if
there is a row in the EVChits file where aas_pos == Q and (i == P or j == P).

Outputs
-------
evc_coverage.tsv        — per-SAAP: has_evc_file, n_coupling_rows in EVChits
evc_cooccurrence.tsv    — co-occurrence pairs annotated with EVC hit info
evc_summary.txt         — console-style summary printed to file
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd

REPO_DIR  = Path(__file__).resolve().parents[2]
TSOUR_DIR = REPO_DIR / "metadata" / "Tsour_et_al"

EVC_DIR_DEFAULT = "/projects/marubi/collabs/slavov_rna/evcouplings_files/EVChits"
COOC_DEFAULT    = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence/saap_missense_cooccurrence.tsv"
OUT_DEFAULT     = "/scratch/leduc.an/AAS_Evo/ANALYSIS/saap_cooccurrence"


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evc-dir",       default=EVC_DIR_DEFAULT)
    ap.add_argument("--cooccurrence",  default=COOC_DEFAULT)
    ap.add_argument("--pep-to-protein",default=str(TSOUR_DIR / "pep_to_protein.csv"))
    ap.add_argument("-o", "--out-dir", default=OUT_DEFAULT)
    return ap.parse_args()


def extract_pos(swap: str):
    """Extract 1-based position from swap string like 'R273H'."""
    if not isinstance(swap, str):
        return None
    m = re.match(r"[A-Z](\d+)[A-Z]", swap)
    return int(m.group(1)) if m else None


def load_evc_index(evc_dir: Path) -> dict[str, pd.DataFrame]:
    """
    Load all EVChits CSVs. Returns dict: ensp_id -> DataFrame.
    Only loads files lazily — we store the path and load on demand.
    """
    files = list(evc_dir.glob("ENSP*_EVChits.csv"))
    return {f.name.split("_")[0]: f for f in files}


def get_evc_hits(evc_path: Path) -> pd.DataFrame:
    df = pd.read_csv(evc_path, index_col=0)
    df["i"] = pd.to_numeric(df["i"], errors="coerce")
    df["j"] = pd.to_numeric(df["j"], errors="coerce")
    df["aas_pos"] = pd.to_numeric(df["aas_pos"], errors="coerce")
    return df


def main():
    args = parse_args()
    evc_dir = Path(args.evc_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Load Tsour pep_to_protein: SAAP → ENSP ───────────────────────────────
    print("Loading pep_to_protein...", flush=True)
    p2prot = pd.read_csv(args.pep_to_protein, index_col=0,
                         usecols=[0, "SAAP", "name", "protein", "protein.position"])
    p2prot = p2prot[["SAAP", "name", "protein", "protein.position"]].drop_duplicates("SAAP")
    p2prot = p2prot.rename(columns={"name": "gene_symbol",
                                    "protein": "ensp_id",
                                    "protein.position": "saap_pos_p2p"})
    saap_to_ensp = dict(zip(p2prot["SAAP"], p2prot["ensp_id"]))
    print(f"  Unique SAAPs in pep_to_protein: {len(saap_to_ensp):,}", flush=True)

    # ── Index available EVChits files ─────────────────────────────────────────
    print(f"Scanning EVChits folder: {evc_dir}", flush=True)
    evc_index = load_evc_index(evc_dir)
    print(f"  EVChits files found: {len(evc_index):,}", flush=True)

    # ── Coverage: which SAAPs have an EVChits file? ───────────────────────────
    print("Computing SAAP coverage...", flush=True)
    coverage_rows = []
    for _, row in p2prot.iterrows():
        ensp = row["ensp_id"]
        has_file = ensp in evc_index
        coverage_rows.append({
            "SAAP":        row["SAAP"],
            "gene":        row["gene_symbol"],
            "ensp_id":     ensp,
            "has_evc_file": has_file,
        })
    coverage = pd.DataFrame(coverage_rows)
    n_total   = len(coverage)
    n_covered = coverage["has_evc_file"].sum()
    print(f"  SAAPs with EVChits file: {n_covered:,} / {n_total:,} "
          f"({100*n_covered/n_total:.1f}%)", flush=True)
    coverage.to_csv(out_dir / "evc_coverage.tsv", sep="\t", index=False)

    # ── Load co-occurrence table ───────────────────────────────────────────────
    print("Loading co-occurrence table...", flush=True)
    cooc = pd.read_csv(args.cooccurrence, sep="\t")
    print(f"  Co-occurrence records: {len(cooc):,}", flush=True)

    # Add ENSP to co-occurrence table
    cooc["ensp_id"] = cooc["SAAP"].map(saap_to_ensp)
    has_ensp = cooc["ensp_id"].notna()
    has_file  = cooc["ensp_id"].isin(evc_index)
    print(f"  Co-occurrence records with ENSP mapping: {has_ensp.sum():,}", flush=True)
    print(f"  Co-occurrence records with EVChits file: {has_file.sum():,}", flush=True)

    # Extract missense position from missense_swap
    cooc["missense_pos"] = cooc["missense_swap"].apply(extract_pos)

    # ── Check EVC coupling for each co-occurrence pair ────────────────────────
    print("Checking evolutionary coupling for co-occurrence pairs...", flush=True)

    # Cache loaded EVChits DataFrames
    evc_cache: dict[str, pd.DataFrame] = {}

    results = []
    testable = cooc[has_file & cooc["missense_pos"].notna()].copy()
    print(f"  Testable pairs (have EVChits + parsed missense pos): {len(testable):,}", flush=True)

    for i, (_, row) in enumerate(testable.iterrows()):
        if i % 1000 == 0:
            print(f"  {i:,} / {len(testable):,} tested", flush=True)

        ensp = row["ensp_id"]
        saap_pos = row["saap_pos"]  # from cooccurrence table (saap_pos column)
        miss_pos = row["missense_pos"]

        # Load EVChits (cached)
        if ensp not in evc_cache:
            evc_cache[ensp] = get_evc_hits(evc_index[ensp])
        evc = evc_cache[ensp]

        # Rows for this SAAP position
        saap_rows = evc[evc["aas_pos"] == saap_pos]
        n_coupling_rows = len(saap_rows)

        # Is the missense position coupled to the SAAP position?
        if n_coupling_rows > 0:
            coupled = saap_rows[(saap_rows["i"] == miss_pos) |
                                (saap_rows["j"] == miss_pos)]
            is_coupled = len(coupled) > 0
            if is_coupled:
                # Take best-scoring coupling row
                best = coupled.sort_values("mad_score", ascending=False).iloc[0]
                best_mad   = best["mad_score"]
                best_prob  = best["probability"]
                best_score = best["score"]
                best_bitscore = best["bitscore"]
            else:
                best_mad = best_prob = best_score = best_bitscore = float("nan")
        else:
            is_coupled = False
            best_mad = best_prob = best_score = best_bitscore = float("nan")

        results.append({
            **row.to_dict(),
            "n_evc_rows_for_saap_pos": n_coupling_rows,
            "is_evc_coupled":          is_coupled,
            "evc_mad_score":           best_mad,
            "evc_probability":         best_prob,
            "evc_score":               best_score,
            "evc_bitscore":            best_bitscore,
        })

    out = pd.DataFrame(results)
    out.to_csv(out_dir / "evc_cooccurrence.tsv", sep="\t", index=False)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_testable  = len(testable)
    n_coupled   = out["is_evc_coupled"].sum()
    n_saap_has_evc_rows = (out["n_evc_rows_for_saap_pos"] > 0).sum()

    lines = [
        "── EVC COUPLING SUMMARY ─────────────────────────────────────────────",
        f"Total co-occurrence records:               {len(cooc):,}",
        f"  With ENSP mapping:                       {has_ensp.sum():,}",
        f"  With EVChits file:                       {has_file.sum():,}",
        f"  Testable (EVChits + parsed miss pos):    {n_testable:,}",
        f"  SAAP position has ≥1 EVC row:            {n_saap_has_evc_rows:,}",
        f"  Evolutionarily coupled (miss pos in EVC):{n_coupled:,} "
        f"({100*n_coupled/n_testable:.1f}% of testable)",
        "",
        "Top coupled genes:",
    ]

    if n_coupled > 0:
        top = (out[out["is_evc_coupled"]]
               .groupby("gene")["SAAP"]
               .count()
               .sort_values(ascending=False)
               .head(15))
        for gene, cnt in top.items():
            lines.append(f"  {gene:<20} {cnt:>5} coupled records")

        lines += [
            "",
            "EVC score distribution (coupled pairs):",
            f"  mad_score:   mean={out.loc[out['is_evc_coupled'],'evc_mad_score'].mean():.3f}  "
            f"median={out.loc[out['is_evc_coupled'],'evc_mad_score'].median():.3f}",
            f"  probability: mean={out.loc[out['is_evc_coupled'],'evc_probability'].mean():.3f}  "
            f"median={out.loc[out['is_evc_coupled'],'evc_probability'].median():.3f}",
        ]

    summary_text = "\n".join(lines)
    print(f"\n{summary_text}", flush=True)
    (out_dir / "evc_summary.txt").write_text(summary_text)
    print(f"\nOutputs saved to {out_dir}/", flush=True)


if __name__ == "__main__":
    main()
