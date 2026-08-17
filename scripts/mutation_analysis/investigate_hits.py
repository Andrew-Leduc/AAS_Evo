#!/usr/bin/env python3
"""
Artifact QC for significant contact-SAAP hits. For each significant swap:
  #3  print the AAS swap (pos, wt>alt) alongside its driving missense(s) (pos, wt>alt)
  #2  coinciding-missense check: any genomic missense near the contact site
      (esp. one whose alt == the swap's alt -> isobaric/localization risk)
  #1  peptide-ID confidence of the hit swap PSMs vs all AAS PSMs vs all PSMs
  #4  PSM/spectral count per hit; swap-peptide uniqueness across genes; ambiguous
      attribution (n_driving_missense > 1)

Usage:
  python3 investigate_hits.py                 # p<0.05 hits from the within-set test
  python3 investigate_hits.py --use-q         # q<0.05 instead
"""
import argparse
import glob
import os
import re
import numpy as np
import pandas as pd

B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"
CONF_CANDIDATES = ["PEP", "Probability", "PeptideProphet Probability"]
LOWER_BETTER = {"PEP"}   # PEP: lower better; Probability: higher better


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sig-file", default=f"{B}/per_swap_carrier_tests.tsv")
    ap.add_argument("--swap-map", default="/scratch/leduc.an/AAS_Evo/FASTA/contact_saap_swap_map.tsv")
    ap.add_argument("--missense", default="/projects/slavov/andrew/AAS_EVO/all_missense_mutations.tsv")
    ap.add_argument("--results-base", default="/scratch/leduc.an/AAS_Evo/MS_SEARCH/results_contact")
    ap.add_argument("--raas-summary", default=f"{B}/per_swap_raas_summary.tsv",
                    help="per-swap raw RAAS summary from the notebook (gene, swap, median ...)")
    ap.add_argument("--pmax", type=float, default=0.05)
    ap.add_argument("--use-q", action="store_true")
    ap.add_argument("--coincide-window", type=int, default=15)
    args = ap.parse_args()

    sig = pd.read_csv(args.sig_file, sep="\t")
    hits = (sig[sig["q"] < 0.05] if args.use_q else sig[sig["p"] < args.pmax]).copy()
    hits = hits.sort_values("p")
    print(f"{len(hits)} hits (from {len(sig)} tested; "
          f"{'q<0.05' if args.use_q else f'p<{args.pmax}'})")
    if not len(hits):
        return
    hit_keys = set(zip(hits["gene"], hits["swap"]))

    # ── swap map: driving + nearby (coinciding) missense ─────────────────────
    smap = pd.read_csv(args.swap_map, sep="\t")
    hits = hits.merge(smap[["gene", "contact_pos", "driving_missense", "n_driving_missense",
                            "nearby_genomic_missense", "nearby_alt_residues"]],
                      on=["gene", "contact_pos"], how="left")

    # ── #2 deeper: coinciding genomic missense from the full missense table ───
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["SYMBOL", "Protein_position", "Amino_acids"])
    miss["pos"] = pd.to_numeric(miss["Protein_position"].astype(str).str.split("-").str[0],
                                errors="coerce")
    aa = miss["Amino_acids"].astype(str).str.split("/", expand=True)
    miss["mwt"], miss["malt"] = aa[0], (aa[1] if aa.shape[1] > 1 else None)
    miss = miss.dropna(subset=["pos", "malt"])
    gm = {g: d for g, d in miss.groupby("SYMBOL")}

    def coinciding(gene, cpos, swap_alt):
        d = gm.get(gene)
        if d is None:
            return "", ""
        near = d[(d["pos"] - cpos).abs() <= args.coincide_window]
        allm = ",".join(sorted({f"{int(p)}:{w}>{a}"
                                for p, w, a in zip(near["pos"], near["mwt"], near["malt"])}))
        same = ",".join(sorted({f"{int(p)}:{w}>{a}"
                                for p, w, a in zip(near["pos"], near["mwt"], near["malt"])
                                if a == swap_alt}))
        return allm, same

    hits[["coincide_near", "coincide_same_residue"]] = hits.apply(
        lambda r: pd.Series(coinciding(r["gene"], r["contact_pos"], str(r["swap"])[-1])), axis=1)

    # ── #1 + #4: scan psm.tsv for ID confidence, PSM count, uniqueness ────────
    files = sorted(glob.glob(os.path.join(args.results_base, "*", "*_1", "psm.tsv"))) or \
            sorted(glob.glob(os.path.join(args.results_base, "*", "psm.tsv")))
    conf_used = None
    hit_conf, aas_conf, all_conf = [], [], []
    hit_psm_count = {}                         # (gene,swap) -> #PSMs
    pep_to_genes = {}                          # swap peptide seq -> set(genes)
    print(f"scanning {len(files)} psm.tsv files for ID confidence...", flush=True)
    for i, f in enumerate(files):
        if i % 40 == 0:
            print(f"  {i}/{len(files)}", flush=True)
        try:
            psm = pd.read_csv(f, sep="\t", low_memory=False)
        except Exception:
            continue
        conf = next((c for c in CONF_CANDIDATES if c in psm.columns), None)
        if conf is None:
            continue
        conf_used = conf
        all_conf.append(pd.to_numeric(psm[conf], errors="coerce").dropna().values)
        cp = psm[psm["Entry Name"].astype(str).str.endswith("-contact")]
        for _, r in cp.iterrows():
            mg = re.match(r"^([A-Z0-9]+)-contact$", str(r.get("Entry Name", "")))
            ms = re.search(r"-([A-Z]\d+[A-Z])-[A-Z0-9]{4}$", str(r.get("Protein ID", "")))
            if not mg or not ms:
                continue
            gene, swap, pep = mg.group(1), ms.group(1), str(r.get("Peptide", ""))
            val = pd.to_numeric(r.get(conf), errors="coerce")
            if pd.isna(val):
                continue
            aas_conf.append(val)
            pep_to_genes.setdefault(pep, set()).add(gene)
            if (gene, swap) in hit_keys:
                hit_conf.append(val)
                hit_psm_count[(gene, swap)] = hit_psm_count.get((gene, swap), 0) + 1

    hits["n_psm"] = [hit_psm_count.get((g, s), 0) for g, s in zip(hits["gene"], hits["swap"])]

    # peptide uniqueness: does any hit swap peptide appear under >1 gene?
    shared = {pep: gs for pep, gs in pep_to_genes.items() if len(gs) > 1}

    # ── raw RAAS: is the hit's RAAS abnormally high vs the overall dist? ──────
    hits["raas_median"] = np.nan
    hits["raas_pctile"] = np.nan
    try:
        rs = pd.read_csv(args.raas_summary, sep="\t")
        overall = rs["median"].dropna().values
        rmap = {(g, s): m for g, s, m in zip(rs["gene"], rs["swap"], rs["median"])}
        hits["raas_median"] = [rmap.get((g, s), np.nan) for g, s in zip(hits["gene"], hits["swap"])]
        hits["raas_pctile"] = [
            round(100 * (overall < v).mean(), 1) if pd.notna(v) else np.nan
            for v in hits["raas_median"]]
        print(f"\n(raw RAAS loaded: overall median log2 RAAS = {np.median(overall):.2f}, "
              f"95th pct = {np.percentile(overall, 95):.2f})")
    except FileNotFoundError:
        print(f"\n(no RAAS summary at {args.raas_summary}; run the notebook RAW RAAS SANITY cell)")

    # ── report ───────────────────────────────────────────────────────────────
    pd.set_option("display.width", 200); pd.set_option("display.max_colwidth", 40)
    cols = ["gene", "contact_pos", "swap", "driving_missense", "n_driving_missense",
            "coincide_same_residue", "n_carrier", "n_noncarrier", "n_sets",
            "delta", "p", "q", "n_psm", "raas_median", "raas_pctile"]
    print("\n=== HITS: AAS swap + driving missense + coinciding same-residue missense ===")
    print(hits[cols].to_string(index=False))

    print("\n=== all coinciding missense within +-{} aa (context) ===".format(args.coincide_window))
    for _, r in hits.iterrows():
        print(f"  {r['gene']} {r['swap']} @ {int(r['contact_pos'])}: {r['coincide_near'] or '(none)'}")

    if conf_used:
        h = np.array(hit_conf); a = np.array(aas_conf)
        al = np.concatenate(all_conf) if all_conf else np.array([])
        d = "lower=better" if conf_used in LOWER_BETTER else "higher=better"
        print(f"\n=== #1 peptide ID confidence ({conf_used}, {d}) ===")
        for lab, v in [("hit swap PSMs", h), ("all AAS/contact PSMs", a), ("all PSMs", al)]:
            if len(v):
                print(f"  {lab:22s} n={len(v):8d}  median={np.median(v):.4f}  "
                      f"q25={np.percentile(v,25):.4f}  q75={np.percentile(v,75):.4f}")

    print("\n=== #4 flags ===")
    weak = hits[hits["n_psm"] <= 1]
    print(f"  hits with <=1 supporting PSM: {len(weak)}"
          + (" -> " + ", ".join(f"{g} {s}" for g, s in zip(weak['gene'], weak['swap'])) if len(weak) else ""))
    amb = hits[pd.to_numeric(hits["n_driving_missense"], errors="coerce") > 1]
    print(f"  hits with >1 driving missense (ambiguous): {len(amb)}"
          + (" -> " + ", ".join(f"{g} {s}" for g, s in zip(amb['gene'], amb['swap'])) if len(amb) else ""))
    coin = hits[hits["coincide_same_residue"].astype(bool) & (hits["coincide_same_residue"] != "")]
    print(f"  hits with a same-residue missense within +-{args.coincide_window} aa: {len(coin)}"
          + (" -> " + ", ".join(f"{g} {s}" for g, s in zip(coin['gene'], coin['swap'])) if len(coin) else ""))
    hit_shared = {p: gs for p, gs in shared.items()}
    print(f"  swap peptides shared across >1 gene (paralog risk): {len(hit_shared)}")
    for p, gs in list(hit_shared.items())[:10]:
        print(f"      {p}: {sorted(gs)}")
    hi_raas = hits[pd.to_numeric(hits["raas_pctile"], errors="coerce") >= 95]
    print(f"  hits with RAAS above the 95th pct of all swaps (abnormally high -> suspect): "
          f"{len(hi_raas)}"
          + (" -> " + ", ".join(f"{g} {s}(p{int(pc)})"
                                for g, s, pc in zip(hi_raas['gene'], hi_raas['swap'],
                                                    hi_raas['raas_pctile'])) if len(hi_raas) else ""))


if __name__ == "__main__":
    main()
