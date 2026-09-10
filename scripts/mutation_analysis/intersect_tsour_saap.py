#!/usr/bin/env python3
"""
Intersect the AAS we detected (contact-SAAP) with the Tsour et al. SAAP set
(metadata/Tsour_et_al/pep_to_protein.csv) as an independent cross-validation.

Match key: (gene, protein position, wt->alt). Tsour encodes this as
name / protein.position / fromto ("P:N" = Pro->Asn). Ours comes from
per_swap_carrier_tests.tsv (gene, contact_pos, swap="W123Y").

Reports exact-key overlap, a position-tolerant (gene, wt->alt) fallback (in case
UniProt vs ENSP numbering differs), and whether the significant hits replicate.
"""
import argparse
import os
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", default=f"{B}/per_swap_carrier_tests.tsv")
    ap.add_argument("--tsour", default=os.path.join(REPO, "metadata/Tsour_et_al/pep_to_protein.csv"))
    ap.add_argument("--out", default=f"{B}/tsour_intersection.tsv")
    a = ap.parse_args()

    ours = pd.read_csv(a.ours, sep="\t")
    sw = ours["swap"].astype(str).str.extract(r"^([A-Z])(\d+)([A-Z])$")
    ours["wt"], ours["alt"] = sw[0], sw[2]
    ours["pos"] = pd.to_numeric(sw[1], errors="coerce")
    ours = ours.dropna(subset=["wt", "alt", "pos"]); ours["pos"] = ours["pos"].astype(int)

    t = pd.read_csv(a.tsour)
    ft = t["fromto"].astype(str).str.split(":", expand=True)
    t["wt"], t["alt"] = ft[0].str.strip(), ft[1].str.strip()
    t["pos"] = pd.to_numeric(t["protein.position"], errors="coerce")
    t = t.dropna(subset=["pos", "wt", "alt"]); t["pos"] = t["pos"].astype(int)

    exact = set(zip(t["name"], t["pos"], t["wt"], t["alt"]))
    loose = set(zip(t["name"], t["wt"], t["alt"]))          # ignore position (numbering offsets)
    ours["in_tsour"] = [(g, p, w, al) in exact
                        for g, p, w, al in zip(ours["gene"], ours["pos"], ours["wt"], ours["alt"])]
    ours["in_tsour_loose"] = [(g, w, al) in loose
                              for g, w, al in zip(ours["gene"], ours["wt"], ours["alt"])]

    n = len(ours)
    print(f"our tested AAS: {n:,}")
    print(f"  in Tsour (exact gene+pos+wt+alt): {ours['in_tsour'].sum():,} "
          f"({100*ours['in_tsour'].mean():.1f}%)")
    print(f"  in Tsour (gene+wt->alt, any pos): {ours['in_tsour_loose'].sum():,} "
          f"({100*ours['in_tsour_loose'].mean():.1f}%)")
    print(f"  (Tsour set: {len(exact):,} unique gene+pos+subs across {t['name'].nunique():,} genes)")

    if "q" in ours.columns:
        sig = ours[ours["q"] < 0.05]
        print(f"\nsignificant hits (q<0.05): {len(sig)} | in Tsour exact: {int(sig['in_tsour'].sum())} "
              f"| loose: {int(sig['in_tsour_loose'].sum())}")
        cols = [c for c in ["gene", "swap", "delta", "q", "in_tsour", "in_tsour_loose"] if c in sig.columns]
        print(sig[cols].to_string(index=False))

    print("\nall AAS that intersect Tsour (exact):")
    hit = ours[ours["in_tsour"]]
    cols = [c for c in ["gene", "swap", "delta", "p", "q"] if c in hit.columns]
    print(hit[cols].to_string(index=False) if len(hit) else "  (none)")
    ours.to_csv(a.out, sep="\t", index=False)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
