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
import hashlib
import re
import random
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

REPO_DIR = Path(__file__).resolve().parents[2]

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Monoisotopic residue masses (Da)
_AA_MASS = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G':  57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P':  97.05276, 'Q': 128.05858, 'R': 156.10111, 'S':  87.03203,
    'T': 101.04768, 'V':  99.06841, 'W': 186.07931, 'Y': 163.06333,
}
_OXIDATION = 15.9949   # mass of one oxygen atom

_DEHYDROGENATION = 2.01565  # loss of 2H (methionine → dehydromethionine)

def _is_oxidation_confound(wt, alt):
    """True if one side is M and the mass difference overlaps with a known
    methionine modification artifact:
      M↔F (+16.03 Da): oxidation/sulfoxide
      M↔D (−16.01 Da): oxidation artifact
      M↔E (−2.00 Da):  dehydromethionine (−2H) ≈ E within typical search tolerance
    """
    if wt != 'M' and alt != 'M':
        return False
    delta = abs(_AA_MASS.get(alt, 0) - _AA_MASS.get(wt, 0))
    return (abs(delta - _OXIDATION) < 0.05 or
            abs(delta - _DEHYDROGENATION) < 0.05)

# Residues that are tryptic cleavage sites or carry TMT label —
# swapping to/from these changes peptide boundaries or labeling state.
_KR = {'K', 'R'}

# Swaps excluded because they are indistinguishable from common modifications
# or isobaric with the wildtype peptide at typical Orbitrap resolution:
#   N->D / Q->E  : deamidation artifact (+0.984 Da), routinely searched as PTM
#   D->N / E->Q  : reverse deamidation (same mass shift)
#   I<->L        : isobaric — identical mass, cannot be distinguished by MS/MS
EXCLUDED_SWAPS = {
    ('N', 'D'), ('D', 'N'),
    ('Q', 'E'), ('E', 'Q'),
    ('I', 'L'), ('L', 'I'),
}

AM_THRESHOLD    = 0.564
SPURS_THRESHOLD = 0.5
AM_BENIGN_MAX   = 0.1
GNOMAD_NEUTRAL  = 0.1
GNOMAD_MAX      = 0.01
VAF_THRESHOLD   = 0.3
DIST_THRESHOLD  = 3.0   # Å Cα-Cα
MIN_SEQ_SEP     = 21    # min AA separation between contact pos and missense pos

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
def tryptic_peptides(seq, pos_1based, max_missed=0):
    """Return tryptic peptides covering pos_1based (no missed cleavages by default).
    FragPipe handles missed cleavages during search so the FASTA only needs the
    fully-tryptic peptide."""
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


def make_swap_peptides(seq, acc, gene, pos_1based, sample_type, patient, source_tag,
                       canonical_peptides=None):
    """Generate tryptic peptides with all 19 AA swaps at pos_1based.

    Header format matches Philosopher-compatible mock-UniProt format:
      >sp|{ACC}-{SWAP}-{HASH}|{GENE}-contact {description} GN={GENE} ...
    Accession uses real UniProt accession as base so Philosopher classifies
    it as a UniProt variant (not generic), stores it in db.bin, and doesn't
    crash during filter. No underscores in accession field.
    """
    entries = []
    wt = seq[pos_1based - 1]
    for alt in ALPHABET:
        if alt == wt:
            continue
        if (wt, alt) in EXCLUDED_SWAPS:
            continue
        if wt in _KR or alt in _KR:
            continue
        if _is_oxidation_confound(wt, alt):
            continue
        mut_seq = seq[:pos_1based - 1] + alt + seq[pos_1based:]
        for start, end in tryptic_peptides(seq, pos_1based):
            pep = mut_seq[start:end]
            if len(pep) < 6:
                continue
            # Skip if this swap peptide sequence exists as a canonical tryptic
            # peptide in the human proteome — it would be indistinguishable
            # from a real protein and create false positives.
            if canonical_peptides is not None and pep in canonical_peptides:
                continue
            swap     = f"{wt}{pos_1based}{alt}"
            seq_hash = hashlib.md5(pep.encode()).hexdigest()[:4].upper()
            accession = f"{acc}-{swap}-{seq_hash}"
            header = (f">sp|{accession}|{gene}-contact "
                      f"{gene} contact prediction {swap} {source_tag} "
                      f"OS=Homo sapiens OX=9606 GN={gene} PE=1 SV=1")
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


def load_contact_map(contact_dir, acc, seq_len=None):
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
            if not pos_to_idx:
                continue
            # Sanity check: contact map positions must be consistent with the
            # reference sequence. A large overshoot means a numbering mismatch
            # (e.g. wrong isoform) and would silently corrupt position lookups.
            if seq_len is not None:
                map_max = max(pos_to_idx.keys())
                if map_max > seq_len + 5:
                    print(f"  [SKIP] {acc}: contact map max pos {map_max} > "
                          f"seq len {seq_len} — likely isoform mismatch")
                    return None, None
            return pos_to_idx, dm
        except Exception:
            continue
    return None, None


def nearby_positions(pos_to_idx, dm, pos_1based, threshold):
    """Return list of 1-based positions within threshold Å of pos_1based (vectorized)."""
    idx = pos_to_idx.get(pos_1based)
    if idx is None:
        return []
    dists    = dm[idx]
    pos_arr  = np.array(list(pos_to_idx.keys()),   dtype=np.int32)
    idx_arr  = np.array(list(pos_to_idx.values()),  dtype=np.int32)
    d_vals   = dists[idx_arr]
    mask     = (d_vals > 0) & (d_vals < threshold) & (pos_arr != pos_1based)
    return pos_arr[mask].tolist()


def main():
    args  = parse_args()
    random.seed(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading reference FASTA...", flush=True)
    seqs, gene2acc, acc2gene = load_ref_fasta(args.ref_fasta)
    print(f"  {len(seqs):,} sequences loaded", flush=True)

    # Build set of all canonical tryptic peptides (0 missed cleavage, len>=6).
    # Any swap peptide matching one of these is indistinguishable from a real
    # canonical peptide in another protein and must be excluded from the FASTA.
    print("Building canonical tryptic peptide set...", flush=True)
    canonical_peptides = set()
    for seq in seqs.values():
        cuts = [-1]
        for i, aa in enumerate(seq):
            if aa in "KR" and (i + 1 >= len(seq) or seq[i + 1] != "P"):
                cuts.append(i)
        cuts.append(len(seq) - 1)
        for i in range(len(cuts) - 1):
            pep = seq[cuts[i] + 1 : cuts[i + 1] + 1]
            if len(pep) >= 6:
                canonical_peptides.add(pep)
    print(f"  {len(canonical_peptides):,} canonical tryptic peptides indexed", flush=True)

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

    # Pre-filter to only sample_ids that appear in any plex (avoids 10M groupby)
    all_plex_uuids = set(gdc.loc[gdc["case_submitter_id"].isin(
        tmt[tmt["run_metadata_id"].isin(plex_ids)]["case_submitter_id"]), "gdc_file_id"])
    destab  = destab[destab["sample_id"].isin(all_plex_uuids)]
    neutral = neutral[neutral["sample_id"].isin(all_plex_uuids)]
    print(f"  After plex filter — Destab: {len(destab):,} | Neutral: {len(neutral):,}", flush=True)

    gene_to_acc = build_gene_to_acc(args.ddg_dir)
    for gene, acc in gene2acc.items():
        if gene not in gene_to_acc:
            gene_to_acc[gene] = acc

    # Pre-build UUID -> sample_type lookup
    uuid_to_stype = gdc.set_index("gdc_file_id")["sample_type"].to_dict()

    # Pre-index missense by sample_id for fast plex subsetting
    destab_by_uuid  = destab.groupby("sample_id")
    neutral_by_uuid = neutral.groupby("sample_id")

    cm_cache    = {}   # acc -> (pos_to_idx, dm)
    nearby_cache = {}  # (acc, pos) -> list of nearby positions

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

        # Collect unique (gene, pos) pairs from plex patients — faster than per-row
        def get_unique_gene_pos(uuid_set, groupby_obj):
            frames = []
            for uid in uuid_set:
                if uid in groupby_obj.groups:
                    frames.append(groupby_obj.get_group(uid))
            if not frames:
                return pd.DataFrame()
            return pd.concat(frames)[["SYMBOL","pos"]].dropna().drop_duplicates()

        destab_gp  = get_unique_gene_pos(uuids, destab_by_uuid)
        neutral_gp = get_unique_gene_pos(uuids, neutral_by_uuid)

        # Sample neutral to match destab count
        if len(neutral_gp) > len(destab_gp):
            neutral_gp = neutral_gp.sample(len(destab_gp), random_state=args.seed)

        entries = {}  # pep_seq -> (header, pep)

        def process_gene_pos(df, source_tag):
            for _, row in df.iterrows():
                gene = str(row["SYMBOL"])
                pos  = int(row["pos"])
                acc  = gene_to_acc.get(gene)
                if not acc or acc not in seqs:
                    continue
                seq = seqs[acc]
                if pos < 1 or pos > len(seq):
                    continue

                if acc not in cm_cache:
                    cm_cache[acc] = load_contact_map(args.contact_dir, acc, seq_len=len(seq))
                pos_to_idx, dm = cm_cache[acc]
                if dm is None:
                    continue

                cache_key = (acc, pos)
                if cache_key not in nearby_cache:
                    nearby_cache[cache_key] = nearby_positions(
                        pos_to_idx, dm, pos, args.dist)
                nearby = nearby_cache[cache_key]
                if not nearby:
                    continue

                for jpos in nearby:
                    if jpos < 1 or jpos > len(seq):
                        continue
                    # Skip if contact position is within max peptide length of the
                    # missense — they would land in the same tryptic peptide.
                    if abs(jpos - pos) < MIN_SEQ_SEP:
                        continue
                    for header, pep in make_swap_peptides(
                            seq, acc, gene, jpos,
                            "tumor", "plex_patient", source_tag,
                            canonical_peptides=canonical_peptides):
                        if pep not in entries:
                            entries[pep] = (header, pep)

        process_gene_pos(destab_gp,  "destab_contact")
        process_gene_pos(neutral_gp, "neutral_contact")

        if not entries:
            continue

        out_fasta = out_dir / f"{plex_id}.fasta"
        with open(out_fasta, "w") as f:
            for header, pep in entries.values():
                f.write(f"{header}\n{pep}\n")

    print(f"\nDone. Contact SAAP FASTAs written to {out_dir}/", flush=True)

    # ── Combine with existing per-plex FASTAs ─────────────────────────────────
    # Appends contact SAAP entries to per_plex/ FASTAs → per_plex_contact/
    # Then add_decoys.py should be run on per_plex_contact/ before FragPipe.
    per_plex_dir     = Path(args.out_dir).parent / "per_plex"
    combined_dir     = Path(args.out_dir).parent / "per_plex_contact"
    combined_dir.mkdir(parents=True, exist_ok=True)

    n_combined = 0
    for plex_id in plex_ids:
        base_fasta    = per_plex_dir    / f"{plex_id}.fasta"
        contact_fasta = out_dir         / f"{plex_id}.fasta"
        combined_fasta= combined_dir    / f"{plex_id}.fasta"

        if not base_fasta.exists():
            continue
        with open(combined_fasta, "w") as out_f:
            out_f.write(base_fasta.read_text())
            if contact_fasta.exists():
                out_f.write("\n")
                out_f.write(contact_fasta.read_text())
        n_combined += 1

    print(f"Combined FASTAs written to {combined_dir}/ ({n_combined} plexes)", flush=True)
    print("Next step: run add_decoys.py on per_plex_contact/ then resubmit FragPipe.", flush=True)


if __name__ == "__main__":
    main()
