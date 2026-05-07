#!/usr/bin/env python3
"""
Generate per-plex FASTAs containing predicted SAAP candidates at positions
structurally proximal (<4Å Cα-Cα) to destabilizing missense mutations.

For each TMT plex:
  1. Find destabilizing missense mutations in patients in that plex
     (AM>=0.564 OR SPURS>=0.5, rare gnomAD<0.01, VAF>=0.3)
  2. For each missense at position i, find all positions j with CA-CA < DIST_THRESHOLD
  3. Add all 19 possible AA swaps at position j as tryptic peptides
  4. Add equal-sized random neutral missense sample (AM<0.1, gnomAD>0.1) with
     their proximal positions as negative control
  5. Write per-plex FASTA for FragPipe search

Output headers: >contact_pred|{ACC}|{GENE}|{SWAP}|{source}|{patient}|{sample_type}
"""

import argparse
import re
import random
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

REPO_DIR = Path(__file__).resolve().parents[2]

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

AM_THRESHOLD    = 0.564
SPURS_THRESHOLD = 0.5
AM_BENIGN_MAX   = 0.1
GNOMAD_NEUTRAL  = 0.1
GNOMAD_MAX      = 0.01
VAF_THRESHOLD   = 0.3
DIST_THRESHOLD  = 4.0   # Å Cα-Cα

DEFAULTS = dict(
    missense   = "/scratch/leduc.an/AAS_Evo/ANALYSIS/all_missense_with_spurs.tsv",
    ref_fasta  = "/scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta",
    contact_dir= "/scratch/leduc.an/AAS_Evo/SPURS/contact_maps",
    ddg_dir    = "/scratch/leduc.an/AAS_Evo/SPURS/ddg_matrices",
    tmt_map    = str(REPO_DIR / "metadata/PDC_meta/pdc_file_tmt_map.tsv"),
    gdc_meta   = str(REPO_DIR / "metadata/GDC_meta/gdc_meta_matched.tsv"),
    out_dir    = "/scratch/leduc.an/AAS_Evo/FASTA/contact_saap",
    plex_list  = "/scratch/leduc.an/AAS_Evo/MS_SEARCH/plex_list.txt",
)


def parse_args():
    ap = argparse.ArgumentParser()
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_','-')}", default=v)
    ap.add_argument("--dist", type=float, default=DIST_THRESHOLD)
    ap.add_argument("--seed", type=int, default=42)
    return ap.parse_args()


# ── Tryptic digestion ────────────────────────────────────────────────────────
def tryptic_peptides(seq, pos_1based, max_missed=1):
    """Return tryptic peptides (0-missed + 1-missed) that cover pos_1based."""
    cuts = [-1]
    for i, aa in enumerate(seq):
        if aa in "KR" and (i + 1 >= len(seq) or seq[i + 1] != "P"):
            cuts.append(i)
    cuts.append(len(seq) - 1)

    idx = pos_1based - 1
    results = []
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 2 + max_missed, len(cuts))):
            start = cuts[i] + 1
            end   = cuts[j] + 1
            if start <= idx < end:
                results.append((start, end))
    return list(set(results))


def make_swap_peptides(seq, acc, gene, pos_1based, sample_type, patient, source_tag):
    """Generate tryptic peptides with all 19 AA swaps at pos_1based."""
    entries = []
    wt = seq[pos_1based - 1]
    for alt in ALPHABET:
        if alt == wt:
            continue
        mut_seq = seq[:pos_1based - 1] + alt + seq[pos_1based:]
        for start, end in tryptic_peptides(seq, pos_1based):
            pep = mut_seq[start:end]
            if len(pep) < 6:
                continue
            swap = f"{wt}{pos_1based}{alt}"
            header = f">contact_pred|{acc}|{gene}|{swap}|{source_tag}|{patient}|{sample_type}"
            entries.append((header, pep))
    return entries


# ── Reference FASTA ──────────────────────────────────────────────────────────
def load_ref_fasta(path):
    seqs, accs, gene2acc = {}, {}, {}
    cur_acc = cur_gene = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                parts = line.split("|")
                cur_acc  = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                m = re.search(r"GN=(\S+)", line)
                cur_gene = m.group(1) if m else None
                seqs[cur_acc] = []
                if cur_gene:
                    gene2acc[cur_gene] = cur_acc
                    accs[cur_acc] = cur_gene
            elif cur_acc:
                seqs[cur_acc].append(line)
    return {acc: "".join(s) for acc, s in seqs.items()}, gene2acc, accs


# ── Contact map helpers ───────────────────────────────────────────────────────
def build_gene_to_acc(ddg_dir):
    return {f.name.split(".")[1]: f.name.split(".")[0]
            for f in Path(ddg_dir).glob("*.ddg_matrix.tsv")}


def load_contact_map(contact_dir, acc):
    cdir = Path(contact_dir)
    candidates = sorted(cdir.glob(f"AF-{acc}-*F1.npy"))
    if not candidates:
        return None, None
    for npy_path in candidates:
        csv_path = npy_path.with_suffix(".csv")
        if not csv_path.exists():
            continue
        try:
            meta = pd.read_csv(csv_path, index_col=0)
            dm   = np.load(npy_path)
            pos_to_idx = {int(row["id"]): i for i, row in meta.iterrows()
                          if pd.notna(row["id"])}
            return pos_to_idx, dm
        except Exception:
            continue
    return None, None


def nearby_positions(pos_to_idx, dm, pos_1based, threshold):
    """Return list of 1-based positions within threshold Å of pos_1based."""
    idx = pos_to_idx.get(pos_1based)
    if idx is None:
        return []
    dists = dm[idx]
    nearby = []
    for p, j in pos_to_idx.items():
        d = dists[j]
        if 0 < d < threshold and p != pos_1based:
            nearby.append(p)
    return nearby


def main():
    args  = parse_args()
    random.seed(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading reference FASTA...", flush=True)
    seqs, gene2acc, acc2gene = load_ref_fasta(args.ref_fasta)
    print(f"  {len(seqs):,} sequences loaded", flush=True)

    print("Loading metadata...", flush=True)
    tmt = pd.read_csv(args.tmt_map, sep="\t")
    gdc = pd.read_csv(args.gdc_meta, sep="\t")
    if "file_id" in gdc.columns and "gdc_file_id" not in gdc.columns:
        gdc = gdc.rename(columns={"file_id": "gdc_file_id"})

    case_sample_to_uuid = (gdc.set_index(["case_submitter_id","sample_type"])
                           ["gdc_file_id"].to_dict())

    with open(args.plex_list) as f:
        plex_ids = [l.strip() for l in f if l.strip()]
    print(f"  {len(plex_ids)} plexes", flush=True)

    print("Loading missense table...", flush=True)
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["sample_id","SYMBOL","Protein_position",
                                "Amino_acids","VAF","gnomADe_AF",
                                "am_pathogenicity","spurs_ddg"])
    miss["VAF"]            = pd.to_numeric(miss["VAF"],            errors="coerce")
    miss["gnomADe_AF"]     = pd.to_numeric(miss["gnomADe_AF"],     errors="coerce").fillna(0)
    miss["am_pathogenicity"]= pd.to_numeric(miss["am_pathogenicity"], errors="coerce")
    miss["spurs_ddg"]      = pd.to_numeric(miss["spurs_ddg"],      errors="coerce")
    miss["pos"]            = pd.to_numeric(
        miss["Protein_position"].astype(str).str.split("-").str[0], errors="coerce")

    base_ok   = (miss["VAF"] >= VAF_THRESHOLD) & (miss["gnomADe_AF"] < GNOMAD_MAX)
    am_pos    = miss["am_pathogenicity"] >= AM_THRESHOLD
    spurs_pos = miss["spurs_ddg"]        >= SPURS_THRESHOLD
    destab    = miss[(am_pos | spurs_pos) & base_ok].copy()
    neutral   = miss[(miss["am_pathogenicity"] <= AM_BENIGN_MAX) &
                     (miss["gnomADe_AF"] >= GNOMAD_NEUTRAL)].copy()

    print(f"  Destab: {len(destab):,} | Neutral: {len(neutral):,}", flush=True)

    gene_to_acc = build_gene_to_acc(args.ddg_dir)
    # supplement from ref FASTA
    for gene, acc in gene2acc.items():
        if gene not in gene_to_acc:
            gene_to_acc[gene] = acc

    cm_cache = {}

    TMT_CH_MAP = {
        "tmt_126":"126","tmt_127n":"127N","tmt_127c":"127C",
        "tmt_128n":"128N","tmt_128c":"128C","tmt_129n":"129N",
        "tmt_129c":"129C","tmt_130n":"130N","tmt_130c":"130C",
        "tmt_131":"131N","tmt_131c":"131C","tmt_126c":"126C","tmt_134n":"134N",
    }

    print(f"\nGenerating FASTAs for {len(plex_ids)} plexes...", flush=True)

    for pi, plex_id in enumerate(plex_ids):
        if pi % 20 == 0:
            print(f"  {pi}/{len(plex_ids)} plexes...", flush=True)

        pt = tmt[tmt["run_metadata_id"] == plex_id].copy()
        pt["channel"] = pt["tmt_channel"].map(TMT_CH_MAP)
        pt = pt.dropna(subset=["channel"])
        pt = pt[~pt["case_submitter_id"].str.lower().isin(
            ["ref","reference","pooled","pool","nan"])]

        # Map to GDC UUIDs
        uuids = set()
        for _, r in pt.iterrows():
            uuid = case_sample_to_uuid.get((r["case_submitter_id"], r["sample_type"]))
            if uuid:
                uuids.add(uuid)
        if not uuids:
            continue

        plex_destab  = destab[destab["sample_id"].isin(uuids)]
        plex_neutral = neutral[neutral["sample_id"].isin(uuids)]

        # Sample neutral to match destab size
        n_destab = len(plex_destab)
        if len(plex_neutral) > n_destab:
            plex_neutral = plex_neutral.sample(n_destab, random_state=args.seed)

        entries = {}  # (header, pep) deduplicated by peptide seq

        def process_mutations(df, source_tag):
            for _, row in df.iterrows():
                gene = str(row["SYMBOL"])
                pos  = row["pos"]
                if pd.isna(pos):
                    continue
                pos = int(pos)
                acc = gene_to_acc.get(gene)
                if not acc or acc not in seqs:
                    continue
                seq = seqs[acc]
                if pos < 1 or pos > len(seq):
                    continue

                # Get contact map
                if acc not in cm_cache:
                    cm_cache[acc] = load_contact_map(args.contact_dir, acc)
                pos_to_idx, dm = cm_cache[acc]
                if dm is None:
                    continue

                nearby = nearby_positions(pos_to_idx, dm, pos, args.dist)
                if not nearby:
                    continue

                patient    = row["sample_id"]
                sample_type = gdc.loc[gdc["gdc_file_id"] == patient,
                                      "sample_type"].iloc[0] \
                              if patient in gdc["gdc_file_id"].values else "unknown"

                for jpos in nearby:
                    if jpos < 1 or jpos > len(seq):
                        continue
                    for header, pep in make_swap_peptides(
                            seq, acc, gene, jpos,
                            sample_type, patient, source_tag):
                        if pep not in entries:
                            entries[pep] = (header, pep)

        process_mutations(plex_destab,  "destab_contact")
        process_mutations(plex_neutral, "neutral_contact")

        if not entries:
            continue

        out_fasta = out_dir / f"{plex_id}.fasta"
        with open(out_fasta, "w") as f:
            for header, pep in entries.values():
                f.write(f"{header}\n{pep}\n")

    print(f"\nDone. FASTAs written to {out_dir}/", flush=True)


if __name__ == "__main__":
    main()
