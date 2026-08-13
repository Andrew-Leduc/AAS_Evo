#!/usr/bin/env python3
"""
Estimate how much relaxing the carrier/non-carrier threshold inflates the
contact-SAAP search space — WITHOUT writing any FASTAs.

For each plex and each threshold (e.g. 3v3 vs 2v2) it counts the unique swap
peptides that generate_contact_saap_fastas.py would add. Reuses that script's
exact selection + swap machinery (imported), so the numbers match the real run.

Usage:
    python3 scripts/fasta_gen/count_swap_space.py --thresholds 3:3,2:2
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_contact_saap_fastas as G   # reuse exact machinery

TMT_CH_MAP = {
    "tmt_126":"126","tmt_127n":"127N","tmt_127c":"127C","tmt_128n":"128N",
    "tmt_128c":"128C","tmt_129n":"129N","tmt_129c":"129C","tmt_130n":"130N",
    "tmt_130c":"130C","tmt_131":"131N","tmt_131c":"131C","tmt_126c":"126C","tmt_134n":"134N",
}


def build_canonical(seqs):
    canon = set()
    for seq in seqs.values():
        cuts = [-1]
        for i, aa in enumerate(seq):
            if aa in "KR" and (i + 1 >= len(seq) or seq[i + 1] != "P"):
                cuts.append(i)
        cuts.append(len(seq) - 1)
        for i in range(len(cuts) - 1):
            pep = seq[cuts[i] + 1: cuts[i + 1] + 1]
            if len(pep) >= 6:
                canon.add(pep)
    return canon


def main():
    ap = argparse.ArgumentParser()
    for k, v in G.DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_','-')}", default=v)
    ap.add_argument("--dist", type=float, default=G.DIST_THRESHOLD)
    ap.add_argument("--vaf-threshold", type=float, default=G.VAF_THRESHOLD)
    ap.add_argument("--thresholds", default="1:0",
                    help="comma list of in-plex carrier:noncarrier, e.g. 1:0,3:3")
    ap.add_argument("--min-patients", type=int, default=G.MIN_PATIENTS,
                    help="min patients dataset-wide carrying the driving missense")
    ap.add_argument("--max-patient-frac", type=float, default=G.MAX_PATIENT_FRAC,
                    help="exclude missense carried by >= this fraction of all patients")
    args = ap.parse_args()
    THS = [tuple(int(x) for x in t.split(":")) for t in args.thresholds.split(",")]

    print("Loading reference + canonical peptides...", flush=True)
    seqs, gene2acc, _ = G.load_ref_fasta(args.ref_fasta)
    canon = build_canonical(seqs)

    tmt = pd.read_csv(args.tmt_map, sep="\t")
    gdc = pd.read_csv(args.gdc_meta, sep="\t")
    if "file_id" in gdc.columns and "gdc_file_id" not in gdc.columns:
        gdc = gdc.rename(columns={"file_id": "gdc_file_id"})
    case_sample_to_uuid = (gdc.set_index(["case_submitter_id", "sample_type"])
                           ["gdc_file_id"].to_dict())

    if Path(args.plex_list).exists():
        plex_ids = [l.strip() for l in open(args.plex_list) if l.strip()]
    else:
        plex_ids = sorted(tmt["run_metadata_id"].dropna().astype(str).unique())

    print("Loading missense...", flush=True)
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["sample_id", "SYMBOL", "Protein_position", "Amino_acids", "VAF"])
    miss["VAF"] = pd.to_numeric(miss["VAF"], errors="coerce")
    miss = miss[miss["VAF"] >= args.vaf_threshold]
    miss["pos"] = pd.to_numeric(
        miss["Protein_position"].astype(str).str.split("-").str[0], errors="coerce")
    aa = miss["Amino_acids"].astype(str).str.split("/", expand=True)
    miss["wt"] = aa[0].str.strip()
    miss["alt"] = aa[1].str.strip() if aa.shape[1] > 1 else None
    miss = miss.dropna(subset=["SYMBOL", "pos", "wt", "alt"])
    miss = miss[(miss["wt"].str.len() == 1) & (miss["alt"].str.len() == 1)]
    miss["pos"] = miss["pos"].astype(int)
    processed = set(miss["sample_id"].unique())
    all_plex_uuids = set(gdc.loc[gdc["case_submitter_id"].isin(
        tmt[tmt["run_metadata_id"].isin(plex_ids)]["case_submitter_id"]), "gdc_file_id"])
    miss = miss[miss["sample_id"].isin(all_plex_uuids)]

    # dataset-wide missense eligibility (mirrors generate_contact_saap_fastas)
    total_patients = len(all_plex_uuids & processed)
    _wide = miss.groupby(["SYMBOL", "pos", "wt", "alt"])["sample_id"].nunique()
    eligible = set(_wide[(_wide >= args.min_patients) &
                         (_wide < args.max_patient_frac * total_patients)].index)
    print(f"dataset-wide eligible missense: {len(eligible):,} of {len(_wide):,} "
          f"(>= {args.min_patients} patients & < {args.max_patient_frac:.0%} "
          f"of {total_patients:,} samples)", flush=True)

    gene_to_acc = G.build_gene_to_acc(args.ddg_dir)
    for g, a in gene2acc.items():
        gene_to_acc.setdefault(g, a)
    miss_by_uuid = miss.groupby("sample_id")
    cm_cache, nearby_cache = {}, {}

    totals = {th: [] for th in THS}     # th -> per-plex unique swap-peptide counts
    n_miss = {th: [] for th in THS}     # th -> per-plex testable-missense counts
    uniq   = {th: set() for th in THS}  # th -> unique (gene,pos,wt,alt) testable missense

    print(f"Scanning {len(plex_ids)} plexes for thresholds {THS} ...", flush=True)
    for pi, plex_id in enumerate(plex_ids):
        if pi % 20 == 0:
            print(f"  {pi}/{len(plex_ids)}", flush=True)
        pt = tmt[tmt["run_metadata_id"] == plex_id].copy()
        pt["channel"] = pt["tmt_channel"].map(TMT_CH_MAP)
        pt = pt.dropna(subset=["channel"])
        pt = pt[~pt["case_submitter_id"].astype(str).str.lower().isin(
            ["ref", "reference", "pooled", "pool", "nan"])]
        uuids = set()
        for _, r in pt.iterrows():
            u = case_sample_to_uuid.get((r["case_submitter_id"], r["sample_type"]))
            if u and u in processed:
                uuids.add(u)
        n_gen = len(uuids)
        frames = [miss_by_uuid.get_group(u) for u in uuids if u in miss_by_uuid.groups]
        if not frames:
            for th in THS:
                totals[th].append(0); n_miss[th].append(0)
            continue
        pm = pd.concat(frames)
        cc = pm.groupby(["SYMBOL", "pos", "wt", "alt"])["sample_id"].nunique()

        for th in THS:
            mc, mnc = th
            ent = set()
            testable = cc[(cc >= mc) & ((n_gen - cc) >= mnc)] if n_gen >= mc + mnc else cc.iloc[:0]
            testable = testable[testable.index.isin(eligible)]
            for (gene, pos, wt, alt), _ in testable.items():
                acc = gene_to_acc.get(gene)
                if not acc or acc not in seqs:
                    continue
                seq = seqs[acc]
                if pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
                    continue
                if acc not in cm_cache:
                    cm_cache[acc] = G.load_contact_map(args.contact_dir, acc, seq_len=len(seq))
                p2i, dm = cm_cache[acc]
                if dm is None:
                    continue
                key = (acc, pos)
                if key not in nearby_cache:
                    nearby_cache[key] = G.nearby_positions(p2i, dm, pos, args.dist)
                for jpos in nearby_cache[key]:
                    if jpos < 1 or jpos > len(seq) or abs(jpos - pos) < G.MIN_SEQ_SEP:
                        continue
                    for _, pep in G.make_swap_peptides(seq, acc, gene, jpos, "tumor", "plex",
                                                       source_tag="contact",
                                                       canonical_peptides=canon):
                        ent.add(pep)
            totals[th].append(len(ent))
            n_miss[th].append(len(testable))
            uniq[th].update(testable.index)

    print("\n" + "=" * 64)
    print("SEARCH-SPACE ESTIMATE (unique swap peptides per plex, no FASTA written)")
    print("=" * 64)
    for th in THS:
        v = pd.Series(totals[th]); m = pd.Series(n_miss[th])
        print(f"{th[0]}c/{th[1]}nc : testable missense/plex mean={m.mean():.0f} | "
              f"swap peptides/plex mean={v.mean():.0f} median={v.median():.0f} "
              f"max={int(v.max())} | total across plexes={int(v.sum()):,}")
    if len(THS) >= 2:
        base, alt = THS[0], THS[1]
        sb, sa = sum(totals[base]), sum(totals[alt])
        mb, ma = sum(n_miss[base]), sum(n_miss[alt])
        print(f"\n{alt[0]}c/{alt[1]}nc vs {base[0]}c/{base[1]}nc:")
        print(f"  testable missense : {ma/max(mb,1):.2f}x  ({mb:,} -> {ma:,})")
        print(f"  swap peptides     : {sa/max(sb,1):.2f}x  ({sb:,} -> {sa:,})  "
              f"(+{100*(sa-sb)/max(sb,1):.0f}%)")

    from collections import Counter
    print("\n" + "=" * 64)
    print("PER-MISSENSE CONTACT/SWAP BREAKDOWN (unique testable missense)")
    print("=" * 64)
    for th in THS:
        n_tot = len(uniq[th]); no_acc = no_wt = no_map = 0
        ncontact = Counter(); nswap = []
        for (gene, pos, wt, alt) in uniq[th]:
            acc = gene_to_acc.get(gene)
            if not acc or acc not in seqs:
                no_acc += 1; continue
            seq = seqs[acc]
            if pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
                no_wt += 1; continue
            if acc not in cm_cache:
                cm_cache[acc] = G.load_contact_map(args.contact_dir, acc, seq_len=len(seq))
            p2i, dm = cm_cache[acc]
            if dm is None:
                no_map += 1; continue
            key = (acc, pos)
            if key not in nearby_cache:
                nearby_cache[key] = G.nearby_positions(p2i, dm, pos, args.dist)
            far = [j for j in nearby_cache[key]
                   if 1 <= j <= len(seq) and abs(j - pos) >= G.MIN_SEQ_SEP]
            ncontact[len(far)] += 1
            ps = set()
            for j in far:
                for _, pep in G.make_swap_peptides(seq, acc, gene, j, "tumor", "plex",
                                                   source_tag="contact", canonical_peptides=canon):
                    ps.add(pep)
            nswap.append(len(ps))
        with_c = sum(v for k, v in ncontact.items() if k > 0)
        print(f"\n{th[0]}c/{th[1]}nc: {n_tot:,} unique testable missense")
        print(f"  dropped: no acc/seq={no_acc:,}  wt-mismatch={no_wt:,}  no contact map={no_map:,}")
        print(f"  with >= 1 contact (>= {G.MIN_SEQ_SEP} aa away): {with_c:,} "
              f"({100*with_c/max(n_tot,1):.1f}% of testable)")
        print("  #contacts -> #missense:")
        for k in sorted(ncontact):
            lab = f"{k}" if k < 6 else "6+"
        agg = Counter()
        for k, v in ncontact.items():
            agg["0" if k == 0 else ("1" if k == 1 else ("2" if k == 2 else ("3-5" if k <= 5 else "6+")))] += v
        for lab in ["0", "1", "2", "3-5", "6+"]:
            if agg[lab]:
                print(f"     {lab:>4} contacts: {agg[lab]:,}")
        sw = np.array([x for x in nswap if x > 0])
        if len(sw):
            print(f"  swaps/missense (those with contacts): "
                  f"mean={sw.mean():.1f} median={np.median(sw):.0f} max={int(sw.max())}")


if __name__ == "__main__":
    main()
