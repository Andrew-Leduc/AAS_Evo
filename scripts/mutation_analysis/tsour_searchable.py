#!/usr/bin/env python3
"""
Fair denominator for the Tsour cross-check: of Tsour's SAAPs, how many were even
IN our search space (could we have detected them), and of those how many we did.

Our search space = the contact swap map: for each searched contact site
(gene, contact_pos) we generated ALL allowed substitutions (excluding K/R, the
isobaric N<->D / Q<->E / I<->L, and M-oxidation confounds). A Tsour SAAP is
"searchable" if its (gene, pos) is a searched contact site, its wt matches our
reference residue there, and its substitution is an allowed one.

Funnel:  Tsour total -> in a searched gene -> at a searched contact pos
         -> searchable (wt match + allowed swap) -> detected by us
"""
import argparse
import os
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"

_AA_MASS = {'A':71.03711,'C':103.00919,'D':115.02694,'E':129.04259,'F':147.06841,
            'G':57.02146,'H':137.05891,'I':113.08406,'K':128.09496,'L':113.08406,
            'M':131.04049,'N':114.04293,'P':97.05276,'Q':128.05858,'R':156.10111,
            'S':87.03203,'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06333}
EXCLUDED = {('N','D'),('D','N'),('Q','E'),('E','Q'),('I','L'),('L','I')}
KR = {'K', 'R'}


def allowed(wt, alt):
    if alt == wt or (wt, alt) in EXCLUDED or wt in KR or alt in KR:
        return False
    if wt == 'M' or alt == 'M':
        d = abs(_AA_MASS.get(alt, 0) - _AA_MASS.get(wt, 0))
        if abs(d - 15.9949) < 0.05 or abs(d - 2.01565) < 0.05:   # M oxidation / -2H
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--swap-map", default="/scratch/leduc.an/AAS_Evo/FASTA/contact_saap_swap_map.tsv")
    ap.add_argument("--detected", default=f"{B}/per_swap_raas_summary.tsv")
    ap.add_argument("--tsour", default=os.path.join(REPO, "metadata/Tsour_et_al/pep_to_protein.csv"))
    a = ap.parse_args()

    smap = pd.read_csv(a.swap_map, sep="\t")
    site_wt = {(g, int(p)): w for g, p, w in
               zip(smap["gene"], smap["contact_pos"], smap["wt_aa"])}
    searched_genes = set(smap["gene"])

    det = pd.read_csv(a.detected, sep="\t")
    ds = det["swap"].astype(str).str.extract(r"^([A-Z])(\d+)([A-Z])$")
    det_keys = set(zip(det["gene"], pd.to_numeric(ds[1], errors="coerce"), ds[0], ds[2]))

    t = pd.read_csv(a.tsour)
    ft = t["fromto"].astype(str).str.split(":", expand=True)
    t["wt"], t["alt"] = ft[0].str.strip(), ft[1].str.strip()
    t["pos"] = pd.to_numeric(t["protein.position"], errors="coerce")
    t = t.dropna(subset=["pos", "wt", "alt"]).copy(); t["pos"] = t["pos"].astype(int)
    t = t.drop_duplicates(["name", "pos", "wt", "alt"])

    in_gene = t["name"].isin(searched_genes)
    at_site = [(g, p) in site_wt for g, p in zip(t["name"], t["pos"])]
    t["at_site"] = at_site
    wt_ok = [site_wt.get((g, p)) == w for g, p, w in zip(t["name"], t["pos"], t["wt"])]
    searchable = [s and wok and allowed(w, al)
                  for s, wok, w, al in zip(t["at_site"], wt_ok, t["wt"], t["alt"])]
    t["searchable"] = searchable
    t["detected"] = [s and (g, p, w, al) in det_keys
                     for s, g, p, w, al in
                     zip(t["searchable"], t["name"], t["pos"], t["wt"], t["alt"])]

    N = len(t)
    print(f"Tsour unique SAAPs (gene+pos+subs): {N:,}")
    print(f"  in a gene we searched ..........: {int(in_gene.sum()):,}")
    print(f"  at a searched contact position .: {int(t['at_site'].sum()):,}")
    print(f"  SEARCHABLE (wt match+allowed sub): {int(t['searchable'].sum()):,}   <- fair denominator")
    ns, nd = int(t["searchable"].sum()), int(t["detected"].sum())
    print(f"  DETECTED by us .................: {nd:,}"
          + (f"   = {100*nd/ns:.0f}% of searchable" if ns else ""))
    if nd:
        print("\nTsour SAAPs we searched AND detected:")
        print(t.loc[t["detected"], ["name", "pos", "wt", "alt", "RAAS", "Tissues"]].to_string(index=False))
    if ns:
        miss = t[t["searchable"] & ~t["detected"]]
        print(f"\nsearchable but NOT detected by us: {len(miss):,} "
              f"(we searched the swap but didn't detect it in the MS)")


if __name__ == "__main__":
    main()
