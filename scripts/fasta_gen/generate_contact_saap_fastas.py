#!/usr/bin/env python3
"""
Generate per-plex FASTAs containing predicted SAAP candidates at positions
structurally proximal to missense mutations that are testable in that plex.

Eligibility (dataset-wide, computed once):
  A missense is eligible if it is carried by >= --min-patients patients across the
  whole on-plex cohort AND by < --max-patient-frac of all patients (the upper bound
  drops near-fixed common variants that have essentially no non-carriers).

For each TMT plex:
  1. Add contact swaps for any dataset-wide-eligible missense that has >= --min-carrier
     carrier channels in that plex (VAF>=0.3). The old within-set 3v3 requirement is
     dropped, so a 1-carrier-vs-many-non-carrier plex still contributes.
     No stability (AM/SPURS) filter — EVE is applied later at analysis time.
  2. For each missense at position i, find positions j whose MINIMUM HEAVY-ATOM
     distance to i is < DIST_THRESHOLD (Angstrom) and at least --min-seq-sep
     residues away in sequence.
  3. Add ALL allowed AA swaps at position j as tryptic peptides (excluding K/R,
     isobaric, and M-modification confounds; dropping canonical peptides).
  4. Write per-plex search FASTA = reference proteome + those swap peptides, plus
     a manifest of (missense -> contact site, carrier/non-carrier counts).

Output headers: >sp|{ACC}-{SWAP}-{HASH}|{GENE}-contact ... GN={GENE} ...
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
DIST_THRESHOLD  = 5.0   # Å, minimum heavy-atom distance (maps are min-heavy-atom, not Cα-Cα)
MIN_SEQ_SEP     = 21    # min AA separation between contact pos and missense pos
MIN_PATIENTS    = 5     # min patients dataset-wide carrying the driving missense
MAX_PATIENT_FRAC= 0.75  # exclude missense carried by >= this fraction of all patients

DEFAULTS = dict(
    missense   = "/projects/slavov/andrew/AAS_EVO/all_missense_mutations.tsv",
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
    ap.add_argument("--min-carrier", type=int, default=1,
                    help="min carrier channels of the missense IN a plex to add its "
                         "contact swaps to that plex's FASTA (within-set 3v3 dropped)")
    ap.add_argument("--min-noncarrier", type=int, default=0,
                    help="min non-carrier channels in the plex (0 = no within-plex "
                         "non-carrier requirement)")
    ap.add_argument("--min-patients", type=int, default=MIN_PATIENTS,
                    help="min patients dataset-wide carrying the driving missense")
    ap.add_argument("--max-patient-frac", type=float, default=MAX_PATIENT_FRAC,
                    help="exclude missense carried by >= this fraction of all patients "
                         "(kills near-fixed common variants)")
    ap.add_argument("--min-seq-sep", type=int, default=MIN_SEQ_SEP,
                    help="min sequence separation (aa) between contact and missense")
    ap.add_argument("--vaf-threshold", type=float, default=VAF_THRESHOLD,
                    help="min VAF for a confident carrier; set 0 to disable")
    ap.add_argument("--gnomad-max", type=float, default=None,
                    help="drop variants with gnomADe_AF above this (e.g. 0.1) to "
                         "exclude common polymorphisms; default: no rarity filter")
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

    if Path(args.plex_list).exists():
        with open(args.plex_list) as f:
            plex_ids = [l.strip() for l in f if l.strip()]
    else:
        plex_ids = sorted(tmt["run_metadata_id"].dropna().astype(str).unique())
        print(f"  (plex_list not found; using all {len(plex_ids)} plexes from TMT map)",
              flush=True)
    print(f"  {len(plex_ids)} plexes", flush=True)

    print("Loading missense table...", flush=True)
    miss = pd.read_csv(args.missense, sep="\t", low_memory=False,
                       usecols=["sample_id", "SYMBOL", "Protein_position",
                                "Amino_acids", "VAF", "gnomADe_AF"])
    miss["VAF"] = pd.to_numeric(miss["VAF"], errors="coerce")
    miss = miss[miss["VAF"] >= args.vaf_threshold]
    if args.gnomad_max is not None:
        af = pd.to_numeric(miss["gnomADe_AF"], errors="coerce").fillna(0)
        miss = miss[af <= args.gnomad_max]
        print(f"  gnomADe_AF <= {args.gnomad_max}: {len(miss):,} rows", flush=True)
    miss["pos"] = pd.to_numeric(
        miss["Protein_position"].astype(str).str.split("-").str[0], errors="coerce")
    aa = miss["Amino_acids"].astype(str).str.split("/", expand=True)
    miss["wt"]  = aa[0].str.strip()
    miss["alt"] = aa[1].str.strip() if aa.shape[1] > 1 else None
    miss = miss.dropna(subset=["SYMBOL", "pos", "wt", "alt"])
    miss = miss[(miss["wt"].str.len() == 1) & (miss["alt"].str.len() == 1)]
    miss["pos"] = miss["pos"].astype(int)
    processed_uuids = set(miss["sample_id"].unique())
    print(f"  {len(miss):,} VAF>={args.vaf_threshold} missense rows "
          f"({len(processed_uuids):,} genomics samples)", flush=True)

    # restrict to sample_ids that appear in any plex (avoids scanning all patients)
    all_plex_uuids = set(gdc.loc[gdc["case_submitter_id"].isin(
        tmt[tmt["run_metadata_id"].isin(plex_ids)]["case_submitter_id"]),
        "gdc_file_id"])
    miss = miss[miss["sample_id"].isin(all_plex_uuids)]
    print(f"  {len(miss):,} on plex patients", flush=True)

    # ── dataset-wide missense eligibility (replaces within-set 3v3) ──────────
    # denominator = every genomics sample on a plex (whether or not it has this
    # missense); numerator = distinct samples carrying the exact (gene,pos,wt,alt).
    total_patients = len(all_plex_uuids & processed_uuids)
    wide = miss.groupby(["SYMBOL", "pos", "wt", "alt"])["sample_id"].nunique()
    hi = args.max_patient_frac * total_patients
    eligible = set(wide[(wide >= args.min_patients) & (wide < hi)].index)
    print(f"  dataset-wide eligible missense: {len(eligible):,} of {len(wide):,} "
          f"(>= {args.min_patients} patients & < {args.max_patient_frac:.0%} "
          f"of {total_patients:,} samples)", flush=True)

    gene_to_acc = build_gene_to_acc(args.ddg_dir)
    for gene, acc in gene2acc.items():
        gene_to_acc.setdefault(gene, acc)

    miss_by_uuid = miss.groupby("sample_id")

    cm_cache     = {}   # acc -> (pos_to_idx, dm)
    nearby_cache = {}   # (acc, pos) -> nearby positions

    TMT_CH_MAP = {
        "tmt_126":"126","tmt_127n":"127N","tmt_127c":"127C",
        "tmt_128n":"128N","tmt_128c":"128C","tmt_129n":"129N",
        "tmt_129c":"129C","tmt_130n":"130N","tmt_130c":"130C",
        "tmt_131":"131N","tmt_131c":"131C","tmt_126c":"126C","tmt_134n":"134N",
    }

    ref_text = Path(args.ref_fasta).read_text()
    if not ref_text.endswith("\n"):
        ref_text += "\n"
    combined_dir = Path(args.out_dir).parent / "per_plex_contact"
    combined_dir.mkdir(parents=True, exist_ok=True)

    manifest = []       # rows describing what was searched
    n_with_fasta = 0

    print(f"\nGenerating FASTAs for {len(plex_ids)} plexes "
          f"(eligible missense with >= {args.min_carrier} in-plex carrier(s), "
          f"contacts < {args.dist}A, >= {args.min_seq_sep} aa away)...", flush=True)

    for pi, plex_id in enumerate(plex_ids):
        if pi % 20 == 0:
            print(f"  {pi}/{len(plex_ids)} plexes, {n_with_fasta} FASTAs...", flush=True)

        pt = tmt[tmt["run_metadata_id"] == plex_id].copy()
        pt["channel"] = pt["tmt_channel"].map(TMT_CH_MAP)
        pt = pt.dropna(subset=["channel"])
        pt = pt[~pt["case_submitter_id"].astype(str).str.lower().isin(
            ["ref","reference","pooled","pool","nan"])]

        # genomics-processed UUIDs on distinct channels in this plex
        uuids = set()
        for _, r in pt.iterrows():
            u = case_sample_to_uuid.get((r["case_submitter_id"], r["sample_type"]))
            if u and u in processed_uuids:
                uuids.add(u)
        n_gen = len(uuids)
        if n_gen < args.min_carrier + args.min_noncarrier:
            continue

        frames = [miss_by_uuid.get_group(u) for u in uuids
                  if u in miss_by_uuid.groups]
        if not frames:
            continue
        pm = pd.concat(frames)

        # carriers of each variant within this plex; keep dataset-wide-eligible ones
        cc = pm.groupby(["SYMBOL","pos","wt","alt"])["sample_id"].nunique()
        testable = cc[(cc >= args.min_carrier) &
                      ((n_gen - cc) >= args.min_noncarrier)]
        testable = testable[testable.index.isin(eligible)]
        if testable.empty:
            continue

        entries = {}   # pep -> (header, pep)
        for (gene, pos, wt, alt), ncarr in testable.items():
            acc = gene_to_acc.get(gene)
            if not acc or acc not in seqs:
                continue
            seq = seqs[acc]
            if pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
                continue
            if acc not in cm_cache:
                cm_cache[acc] = load_contact_map(args.contact_dir, acc, seq_len=len(seq))
            pos_to_idx, dm = cm_cache[acc]
            if dm is None:
                continue
            key = (acc, pos)
            if key not in nearby_cache:
                nearby_cache[key] = nearby_positions(pos_to_idx, dm, pos, args.dist)
            for jpos in nearby_cache[key]:
                if jpos < 1 or jpos > len(seq) or abs(jpos - pos) < args.min_seq_sep:
                    continue
                n_added = 0
                for header, pep in make_swap_peptides(
                        seq, acc, gene, jpos, "tumor", "plex",
                        source_tag="contact",
                        canonical_peptides=canonical_peptides):
                    if pep not in entries:
                        entries[pep] = (header, pep)
                    n_added += 1
                if n_added:
                    manifest.append((plex_id, gene, acc, pos, wt, alt,
                                     int(ncarr), int(n_gen - ncarr), jpos, n_added))

        if not entries:
            continue

        # per-plex swap-only FASTA (for inspection)
        with open(out_dir / f"{plex_id}.fasta", "w") as f:
            for header, pep in entries.values():
                f.write(f"{header}\n{pep}\n")
        # combined search FASTA = reference proteome + contact swap peptides
        with open(combined_dir / f"{plex_id}.fasta", "w") as f:
            f.write(ref_text)
            for header, pep in entries.values():
                f.write(f"{header}\n{pep}\n")
        n_with_fasta += 1

    man = pd.DataFrame(manifest, columns=[
        "plex_id","gene","acc","miss_pos","miss_wt","miss_alt",
        "n_carrier","n_noncarrier","contact_pos","n_swap_peptides"])
    man_path = Path(args.out_dir).parent / "contact_saap_manifest.tsv"
    man.to_csv(man_path, sep="\t", index=False)

    print(f"\nDone. {n_with_fasta} plex FASTAs -> {combined_dir}/", flush=True)
    print(f"  swap-only FASTAs -> {out_dir}/", flush=True)
    print(f"  manifest         -> {man_path}", flush=True)
    if len(man):
        n_uniq = man[["gene","miss_pos","miss_wt","miss_alt"]].drop_duplicates().shape[0]
        print(f"  missense-contact rows: {len(man):,} | unique missense: {n_uniq:,} "
              f"| genes: {man['gene'].nunique():,}", flush=True)
    print("Next: run add_decoys.py on per_plex_contact/ then FragPipe.", flush=True)


if __name__ == "__main__":
    main()
