#!/usr/bin/env python3
"""
extract_plex_targets.py

For a single TMT plex, identify low-ratio mutant PSMs (ratio 1–2), determine
which TMT channels are "not-have" for each mutation, and output:

  1. targets.tsv  — one row per (not-have UUID, mutation, genomic position)
                    ready for BAM interrogation
  2. bam_manifest.tsv — GDC-format manifest for re-downloading the not-have BAMs

Usage:
    python3 extract_plex_targets.py [options]

Options:
    --plex-id      TMT plex run_metadata_id (default: first plex in plex_list.txt)
    --ratio-lo     Lower ratio bound for "low ratio" filter (default: 1.0)
    --ratio-hi     Upper ratio bound for "low ratio" filter (default: 2.0)
    --out-dir      Output directory for targets.tsv and bam_manifest.tsv
                   (default: MS_SEARCH/bam_interrogate/)
    --data-dir     Root data directory on cluster (default: /scratch/leduc.an/AAS_Evo)
    --repo-dir     Repo root (default: /home/leduc.an/AAS_Evo_project/AAS_Evo)
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

import pandas as pd
import numpy as np


# ── CONSTANTS ─────────────────────────────────────────────────────────────────
CHANNEL_ORDER = ["126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C"]
TMT_CHANNEL_MAP = {
    "tmt_126":"126",  "tmt_127n":"127N", "tmt_127c":"127C",
    "tmt_128n":"128N","tmt_128c":"128C", "tmt_129n":"129N",
    "tmt_129c":"129C","tmt_130n":"130N", "tmt_130c":"130C",
    "tmt_131":"131N", "tmt_131c":"131C",
}
AA3TO1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
}


def parse_hgvsp_to_swap(hgvsp):
    m = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', str(hgvsp))
    if m:
        ref = AA3TO1.get(m.group(1))
        alt = AA3TO1.get(m.group(3))
        if ref and alt:
            return f"{ref}{m.group(2)}{alt}"
    return None


def parse_mapped_proteins(mapped_str):
    """'sp|P01889-S35A-B4B8|HLA-B-mut, ...' → {('P01889','S35A'), ...}"""
    if pd.isna(mapped_str):
        return set()
    pairs = set()
    for entry in str(mapped_str).split(","):
        parts = entry.strip().split("|")
        if len(parts) >= 2:
            pid_parts = parts[1].split("-")
            if len(pid_parts) >= 3:
                pairs.add((pid_parts[0], pid_parts[1]))
    return pairs


def parse_protein_id(protein_id):
    """'P11047-L888P-F262' → {('P11047','L888P')}."""
    if pd.isna(protein_id):
        return set()
    parts = str(protein_id).split("-")
    if len(parts) >= 3:
        return {(parts[0], parts[1])}
    return set()


def find_ri_cols(df):
    found = {}
    for ch in CHANNEL_ORDER:
        if ch in df.columns:
            found[ch] = ch; continue
        cands = [c for c in df.columns if c.startswith("Intensity") and c.endswith(f"_{ch}")]
        if cands:
            found[ch] = cands[0]; continue
        cands = [c for c in df.columns if c.startswith(ch) and "intensity" in c.lower()]
        if cands:
            found[ch] = cands[0]
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plex-id",   default=None,
                        help="TMT plex run_metadata_id. Defaults to first in plex_list.txt.")
    parser.add_argument("--ratio-lo",  type=float, default=1.0)
    parser.add_argument("--ratio-hi",  type=float, default=2.0)
    parser.add_argument("--data-dir",  default="/scratch/leduc.an/AAS_Evo")
    parser.add_argument("--repo-dir",  default="/home/leduc.an/AAS_Evo_project/AAS_Evo")
    parser.add_argument("--out-dir",   default=None)
    args = parser.parse_args()

    DATA_DIR   = args.data_dir
    REPO_DIR   = args.repo_dir
    SEARCH_DIR = os.path.join(DATA_DIR, "MS_SEARCH")
    OUT_DIR    = args.out_dir or os.path.join(SEARCH_DIR, "bam_interrogate")
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Resolve plex ID ───────────────────────────────────────────────────────
    plex_list_path = os.path.join(SEARCH_DIR, "plex_list.txt")
    with open(plex_list_path) as f:
        all_plexes = [l.strip() for l in f if l.strip()]
    PLEX_ID = args.plex_id or all_plexes[0]
    print(f"Plex: {PLEX_ID}")

    # ── Paths ─────────────────────────────────────────────────────────────────
    TMT_MAP    = os.path.join(REPO_DIR, "metadata/PDC_meta/pdc_file_tmt_map.tsv")
    GDC_META   = os.path.join(REPO_DIR, "metadata/GDC_meta/gdc_meta_matched.tsv")
    FASTA_DIR  = os.path.join(DATA_DIR, "FASTA/per_sample")
    MISSENSE   = os.path.join(DATA_DIR, "VEP/all_missense_mutations.tsv")
    REF_FASTA  = os.path.join(DATA_DIR, "SEQ_FILES/uniprot_human_canonical.fasta")
    RESULTS_DIR = os.path.join(SEARCH_DIR, "results", PLEX_ID)

    psm_matches = sorted(glob.glob(os.path.join(RESULTS_DIR, "*_1/psm.tsv")))
    if not psm_matches:
        sys.exit(f"ERROR: No psm.tsv found under {RESULTS_DIR}")
    PSM_FILE = psm_matches[0]
    print(f"PSM file: {PSM_FILE}")

    # ── Load shared resources ─────────────────────────────────────────────────
    print("Loading metadata...")
    tmt = pd.read_csv(TMT_MAP, sep="\t")
    gdc = pd.read_csv(GDC_META, sep="\t")
    if "file_id" in gdc.columns and "gdc_file_id" not in gdc.columns:
        gdc = gdc.rename(columns={"file_id": "gdc_file_id"})

    # gene symbol → UniProt accession
    gene_to_acc = {}
    with open(REF_FASTA) as f:
        for line in f:
            if line.startswith(">"):
                m = re.search(r'GN=(\S+)', line)
                if m:
                    gene_to_acc[m.group(1)] = line.split("|")[1]
    acc_to_gene = {v: k for k, v in gene_to_acc.items()}
    print(f"  Gene→acc entries: {len(gene_to_acc):,}")

    # ── Build plex_meta: channel → uuid ──────────────────────────────────────
    print("Building plex metadata...")
    plex_tmt = (tmt[tmt["run_metadata_id"] == PLEX_ID]
                [["tmt_channel","case_submitter_id","sample_type"]].drop_duplicates())
    plex_tmt = plex_tmt[~plex_tmt["case_submitter_id"].str.lower()
                        .isin(["ref","reference","pooled","pool","nan"])]
    plex_tmt["channel"] = plex_tmt["tmt_channel"].map(TMT_CHANNEL_MAP)
    plex_meta = plex_tmt.merge(
        gdc[["gdc_file_id","case_submitter_id","sample_type"]],
        on=["case_submitter_id","sample_type"], how="left")

    uuid_to_channel = (plex_meta.dropna(subset=["gdc_file_id","channel"])
                                .set_index("gdc_file_id")["channel"].to_dict())
    channel_to_uuid = {v: k for k, v in uuid_to_channel.items()}
    channel_to_case = (plex_meta.dropna(subset=["channel"])
                                .set_index("channel")["case_submitter_id"].to_dict())

    print(f"  Channels in plex: {len(plex_tmt)},  matched to GDC: {len(uuid_to_channel)}")

    # ── Build mutation_to_channels and channels_with_fastas ──────────────────
    print("Reading per-sample FASTAs...")
    mutation_to_channels = defaultdict(set)
    channels_with_fastas = set()
    for _, row in plex_meta.iterrows():
        uuid, channel = row["gdc_file_id"], row["channel"]
        if pd.isna(uuid) or pd.isna(channel):
            continue
        fasta_path = os.path.join(FASTA_DIR, f"{uuid}_mutant.fasta")
        if not os.path.exists(fasta_path):
            continue
        channels_with_fastas.add(channel)
        with open(fasta_path) as f:
            for line in f:
                if not line.startswith(">"): continue
                parts = line[1:].strip().split("|")
                if len(parts) >= 4 and parts[0] == "mut":
                    mutation_to_channels[(parts[1], parts[3])].add(channel)

    print(f"  FASTA channels: {len(channels_with_fastas)}")
    print(f"  Unique mutations mapped: {len(mutation_to_channels):,}")

    # ── Build genomic position lookup: (uuid, gene, swap) → (CHROM,POS,REF,ALT,VAF)
    print("Loading missense table...")
    missense = pd.read_csv(MISSENSE, sep="\t", low_memory=False)
    plex_missense = missense[missense["sample_id"].isin(uuid_to_channel)]
    print(f"  Missense rows for this plex: {len(plex_missense):,}")

    # Build: (have_uuid, gene, swap) → (CHROM, POS, REF, ALT, VAF)
    genomic_lookup = {}
    for _, vrow in plex_missense.iterrows():
        uuid  = vrow["sample_id"]
        gene  = str(vrow.get("SYMBOL", ""))
        swap  = parse_hgvsp_to_swap(str(vrow.get("HGVSp", "")))
        if not swap:
            continue
        chrom = str(vrow.get("CHROM", ""))
        pos   = vrow.get("POS")
        ref   = str(vrow.get("REF", ""))
        alt   = str(vrow.get("ALT", ""))
        vaf   = float(vrow.get("VAF", np.nan)) if not pd.isna(vrow.get("VAF")) else np.nan
        if chrom and pos and ref and alt:
            genomic_lookup[(uuid, gene, swap)] = {
                "CHROM": chrom, "POS": int(pos), "REF": ref, "ALT": alt, "have_vaf": vaf
            }

    print(f"  Genomic position entries: {len(genomic_lookup):,}")

    # ── Load PSM and compute ratios ───────────────────────────────────────────
    print("Loading PSM table...")
    psm = pd.read_csv(PSM_FILE, sep="\t", low_memory=False)
    ri_col_map = find_ri_cols(psm)
    ri_cols    = list(ri_col_map.values())

    mut_all = psm[psm["Entry Name"].str.endswith("-mut", na=False)].copy()
    mut_psm = mut_all[(mut_all[ri_cols].fillna(0) > 0).any(axis=1)].copy() \
              if ri_cols else mut_all.copy()
    print(f"  Mutant PSMs with signal: {len(mut_psm):,}")

    # ── Extract targets ───────────────────────────────────────────────────────
    print(f"Extracting targets (ratio {args.ratio_lo}–{args.ratio_hi})...")
    target_rows = []
    seen = set()   # deduplicate on (not_have_uuid, CHROM, POS)

    for _, row in mut_psm.iterrows():
        mutations = parse_mapped_proteins(row.get("Mapped Proteins", float("nan")))
        if not mutations:
            mutations = parse_protein_id(row.get("Protein ID", float("nan")))
        if not mutations:
            continue

        have_ch = set()
        for mk in mutations:
            have_ch |= mutation_to_channels.get(mk, set())
        have_ch &= channels_with_fastas
        not_have_ch = channels_with_fastas - have_ch
        if not have_ch or not not_have_ch:
            continue

        have_ri  = [row[ri_col_map[ch]] for ch in have_ch
                    if ch in ri_col_map and pd.notna(row[ri_col_map[ch]])]
        not_ri   = [row[ri_col_map[ch]] for ch in not_have_ch
                    if ch in ri_col_map and pd.notna(row[ri_col_map[ch]])]
        if not have_ri or not not_ri:
            continue
        mean_not  = np.mean(not_ri)
        mean_have = np.mean(have_ri)
        if mean_not == 0:
            continue
        ratio = mean_have / mean_not
        if not (args.ratio_lo <= ratio <= args.ratio_hi):
            continue

        pep = str(row.get("Peptide", ""))

        # For each (accession, swap) in this PSM
        for accession, swap in mutations:
            gene = acc_to_gene.get(accession, accession)
            # Find genomic position from any have-channel UUID
            geo = None
            have_uuid_used, have_case_used = None, None
            for have_channel in have_ch:
                have_uuid = channel_to_uuid.get(have_channel)
                if have_uuid:
                    geo = genomic_lookup.get((have_uuid, gene, swap))
                    if geo:
                        have_uuid_used = have_uuid
                        have_case_used = channel_to_case.get(have_channel, "")
                        break

            if not geo:
                continue  # can't find genomic coordinates for this mutation

            # Emit one row per not-have UUID
            for not_have_channel in not_have_ch:
                not_have_uuid = channel_to_uuid.get(not_have_channel)
                if not not_have_uuid:
                    continue
                dedup_key = (not_have_uuid, geo["CHROM"], geo["POS"])
                if dedup_key in seen:
                    continue
                seen.add(dedup_key)
                target_rows.append({
                    "not_have_uuid":    not_have_uuid,
                    "not_have_case":    channel_to_case.get(not_have_channel, ""),
                    "not_have_channel": not_have_channel,
                    "have_uuid":        have_uuid_used,
                    "have_case":        have_case_used,
                    "have_channel":     next(iter(have_ch)),
                    "gene":             gene,
                    "accession":        accession,
                    "swap":             swap,
                    "peptide":          pep,
                    "ratio":            round(ratio, 3),
                    "have_vaf":         round(geo["have_vaf"], 4) if not np.isnan(geo["have_vaf"]) else "",
                    "CHROM":            geo["CHROM"],
                    "POS":              geo["POS"],
                    "REF":              geo["REF"],
                    "ALT":              geo["ALT"],
                })

    targets = pd.DataFrame(target_rows)
    n_unique_uuids = targets["not_have_uuid"].nunique() if len(targets) else 0
    print(f"  Target rows:        {len(targets):,}")
    print(f"  Unique not-have UUIDs: {n_unique_uuids}")

    if len(targets) == 0:
        print("WARNING: No targets found. Check plex ID and ratio range.")
        sys.exit(0)

    targets_path = os.path.join(OUT_DIR, "targets.tsv")
    targets.to_csv(targets_path, sep="\t", index=False)
    print(f"  → Wrote {targets_path}")

    # ── Build GDC manifest for not-have UUIDs ────────────────────────────────
    print("Building BAM manifest for not-have UUIDs...")
    unique_not_have = set(targets["not_have_uuid"])
    manifest_rows = gdc[gdc["gdc_file_id"].isin(unique_not_have)][
        ["gdc_file_id", "file_name", "md5", "size", "state"]
    ].rename(columns={"gdc_file_id": "id", "file_name": "filename"})

    print(f"  UUIDs requested: {len(unique_not_have)},  matched in GDC meta: {len(manifest_rows)}")
    missing = unique_not_have - set(manifest_rows["id"])
    if missing:
        print(f"  WARNING: {len(missing)} UUIDs not found in gdc_meta_matched.tsv:")
        for u in sorted(missing):
            print(f"    {u}")

    manifest_path = os.path.join(OUT_DIR, "bam_manifest.tsv")
    manifest_rows.to_csv(manifest_path, sep="\t", index=False)
    print(f"  → Wrote {manifest_path}  ({len(manifest_rows)} BAMs to download)")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\nDone.")
    print(f"  targets.tsv:      {targets_path}")
    print(f"  bam_manifest.tsv: {manifest_path}")
    print(f"\nNext steps:")
    print(f"  1. Download BAMs:")
    print(f"       python3 scripts/download/gdc/download.py --manifest {manifest_path}")
    print(f"  2. Run BAM interrogation:")
    print(f"       N=$(wc -l < {targets_path})")
    print(f"       sbatch --array=1-${{N}} scripts/ms_search/fragpipe/validation/interrogate_bams.sh")
    print(f"  3. Visualize results:")
    print(f"       python3 scripts/ms_search/fragpipe/validation/plot_bam_support.py --out-dir {OUT_DIR}")


if __name__ == "__main__":
    main()
