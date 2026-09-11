#!/usr/bin/env python3
"""
Per-plex/per-sample-aware Tsour cross-check.

Our FASTAs were per-plex: a swap was in plex P's database only if a patient in P
carried the contact-driving missense. So a Tsour SAAP was only truly "searchable"
if it was detected (by Tsour) in a sample that maps to one of our plexes AND that
plex's FASTA actually contained the swap. This computes that fair denominator.

Funnel (unique gene+pos+wt->alt):
  Tsour SAAPs with a sample
   -> sample maps to our cohort (case_submitter_id in TMT map)
   -> the swap was searched in that sample's plex (manifest) + wt match + allowed sub
   -> detected by us
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
        if abs(d - 15.9949) < 0.05 or abs(d - 2.01565) < 0.05:
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pep2prot", default=os.path.join(REPO, "metadata/Tsour_et_al/pep_to_protein.csv"))
    ap.add_argument("--pep2pat", default=os.path.join(REPO, "metadata/Tsour_et_al/peptide_to_patient.csv"))
    ap.add_argument("--tmt-map", default=os.path.join(REPO, "metadata/PDC_meta/pdc_file_tmt_map.tsv"))
    ap.add_argument("--swap-map", default="/scratch/leduc.an/AAS_Evo/FASTA/contact_saap_swap_map.tsv")
    ap.add_argument("--manifest", default="/scratch/leduc.an/AAS_Evo/FASTA/contact_saap_manifest.tsv")
    ap.add_argument("--detected", default=f"{B}/per_swap_raas_summary.tsv")
    a = ap.parse_args()

    # Tsour SAAP -> (gene, pos, wt, alt)
    p2p = pd.read_csv(a.pep2prot)
    ft = p2p["fromto"].astype(str).str.split(":", expand=True)
    p2p["wt"], p2p["alt"] = ft[0].str.strip(), ft[1].str.strip()
    p2p["pos"] = pd.to_numeric(p2p["protein.position"], errors="coerce")
    p2p = p2p.dropna(subset=["pos", "wt", "alt"])[["SAAP", "name", "pos", "wt", "alt"]]
    p2p["pos"] = p2p["pos"].astype(int)
    p2p = p2p.rename(columns={"name": "gene"}).drop_duplicates()

    # Tsour SAAP -> samples
    pat = pd.read_csv(a.pep2pat, usecols=["SAAP", "Sample name", "Dataset"])
    pat["case"] = pat["Sample name"].astype(str).str.rsplit("_", n=1).str[0]
    det_rows = pat.merge(p2p, on="SAAP", how="inner")   # SAAP detection x gene mapping

    # our per-plex search space + contact wt
    tmt = pd.read_csv(a.tmt_map, sep="\t")
    case2plex = tmt.groupby("case_submitter_id")["run_metadata_id"].agg(set).to_dict()
    man = pd.read_csv(a.manifest, sep="\t")
    plex_sites = man.groupby("plex_id").apply(
        lambda d: set(zip(d["gene"], d["contact_pos"]))).to_dict()
    smap = pd.read_csv(a.swap_map, sep="\t")
    site_wt = {(g, int(p)): w for g, p, w in zip(smap["gene"], smap["contact_pos"], smap["wt_aa"])}

    det = pd.read_csv(a.detected, sep="\t")
    ds = det["swap"].astype(str).str.extract(r"^([A-Z])(\d+)([A-Z])$")
    det_keys = set(zip(det["gene"], pd.to_numeric(ds[1], errors="coerce"), ds[0], ds[2]))

    def searched_in_sample_plex(gene, pos, wt, alt, case):
        if site_wt.get((gene, pos)) != wt or not allowed(wt, alt):
            return False
        return any((gene, pos) in plex_sites.get(P, set()) for P in case2plex.get(case, set()))

    det_rows["in_cohort"] = det_rows["case"].isin(case2plex)
    det_rows["searched_here"] = [
        searched_in_sample_plex(g, p, w, al, c)
        for g, p, w, al, c in zip(det_rows["gene"], det_rows["pos"], det_rows["wt"],
                                  det_rows["alt"], det_rows["case"])]

    # aggregate to unique SAAP (gene,pos,wt,alt)
    key = ["gene", "pos", "wt", "alt"]
    g = det_rows.groupby(key).agg(in_cohort=("in_cohort", "any"),
                                  searched_here=("searched_here", "any")).reset_index()
    g["detected"] = [(gg, pp, ww, aa) in det_keys
                     for gg, pp, ww, aa in zip(g["gene"], g["pos"], g["wt"], g["alt"])]

    print(f"Tsour SAAPs with a sample (unique gene+pos+subs): {len(g):,}")
    print(f"  sample maps to our cohort .............: {int(g['in_cohort'].sum()):,}")
    ns = int(g["searched_here"].sum())
    print(f"  searched in that sample's plex (fair) .: {ns:,}   <- per-plex denominator")
    nd = int((g["searched_here"] & g["detected"]).sum())
    print(f"  detected by us ........................: {nd:,}"
          + (f"   = {100*nd/ns:.0f}% of per-plex-searchable" if ns else ""))
    if ns:
        print("\nper-plex-searchable AND detected:")
        print(g[g["searched_here"] & g["detected"]][key].to_string(index=False))
        print("\nper-plex-searchable but NOT detected:")
        print(g[g["searched_here"] & ~g["detected"]][key].to_string(index=False))


if __name__ == "__main__":
    main()
