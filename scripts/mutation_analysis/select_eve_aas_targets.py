#!/usr/bin/env python3
"""
Select AAS (translation-error) targets that EVE predicts will interact strongly
— or not at all — with an observed missense mutation.

Reframed pipeline (contacts narrow candidates, EVE scores the interaction):

  1. Take missense mutations observed in the CPTAC cohort (any consequence for
     stability — AM/SPURS are annotation only, NOT filters).
  2. Keep missense that are "powered": there is >= MIN_CARRIER carrier channels
     and >= MIN_NONCARRIER non-carrier channels in at least one TMT set, so a
     within-set carrier vs non-carrier RAAS contrast is possible.
  3. For each such missense at position i, take the structurally proximal
     positions j (Ca-Ca < DIST_THRESHOLD from existing AlphaFold contact maps),
     at least MIN_SEQ_SEP residues away in sequence (so the AAS peptide is a
     different tryptic peptide than the missense).
  4. For each candidate AAS (position j -> alt_j):
        epistasis = dmm_score - (missense_single + swap_single)
     computed from the protein's EVE double-mutant matrix.
        epistasis > 0  => double LESS deleterious than additive => synergistic
        epistasis < 0  => double MORE deleterious than additive => aggravating
        epistasis ~ 0  => additive / no interaction (null control)
  5. Two-arm design:
        --n-syn  most synergistic (top +epistasis)
        --n-del  most deleterious  (top -epistasis)
        --n-null nearest-zero |epistasis|
     tagged so each arm can be powered separately in the RAAS readout.

Stage 1 (this run): select arms from the PRECOMPUTED EVE matrices.
Stage 2 (this run): emit the list of proteins that have eligible missense +
     contact candidates but NO precomputed EVE matrix, to run EVE on next.

Outputs (OUT_DIR):
  eve_aas_targets.tsv        one row per selected (missense, AAS) pair, arm-labeled
  eve_aas_candidates.tsv     ALL scored candidate pairs (pre-selection)
  eve_proteins_to_compute.tsv proteins needing EVE (eligible missense, no matrix)
  selection_summary.txt      counts / distributions
"""

import argparse
import re
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

REPO_DIR = Path(__file__).resolve().parents[2]
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# ── Monoisotopic residue masses (for MS confound exclusions) ─────────────────
_AA_MASS = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G':  57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P':  97.05276, 'Q': 128.05858, 'R': 156.10111, 'S':  87.03203,
    'T': 101.04768, 'V':  99.06841, 'W': 186.07931, 'Y': 163.06333,
}
_OXIDATION       = 15.9949
_DEHYDROGENATION = 2.01565
_KR = {'K', 'R'}
# isobaric / deamidation-confounded swaps (same as generate_contact_saap_fastas.py)
EXCLUDED_SWAPS = {('N','D'),('D','N'),('Q','E'),('E','Q'),('I','L'),('L','I')}


def _is_oxidation_confound(wt, alt):
    if wt != 'M' and alt != 'M':
        return False
    delta = abs(_AA_MASS.get(alt, 0) - _AA_MASS.get(wt, 0))
    return abs(delta - _OXIDATION) < 0.05 or abs(delta - _DEHYDROGENATION) < 0.05


def swap_allowed(wt, alt):
    """Cheap swap-level detectability filter (independent of position)."""
    if wt == alt:
        return False
    if wt in _KR or alt in _KR:            # changes tryptic cleavage / TMT labeling
        return False
    if (wt, alt) in EXCLUDED_SWAPS:        # isobaric / deamidation artifact
        return False
    if _is_oxidation_confound(wt, alt):    # M +/-16 / -2 Da artifacts
        return False
    return True


# ── Reference proteome + gene<->accession ────────────────────────────────────
def load_reference(ref_fasta, ddg_dir):
    gene_to_acc = {}
    if ddg_dir and Path(ddg_dir).exists():
        for f in Path(ddg_dir).glob('*.ddg_matrix.tsv'):
            parts = f.name.split('.')
            if len(parts) >= 2:
                gene_to_acc[parts[1]] = parts[0]
    seqs = {}
    cur = None
    with open(ref_fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                p = line.split('|')
                cur = p[1] if len(p) >= 3 else line[1:].split()[0]
                m = re.search(r'GN=(\S+)', line)
                seqs[cur] = []
                if m and m.group(1) not in gene_to_acc:
                    gene_to_acc[m.group(1)] = cur
            elif cur:
                seqs[cur].append(line)
    seqs = {a: ''.join(s) for a, s in seqs.items()}
    return seqs, gene_to_acc


def build_canonical_peptides(seqs):
    canon = set()
    for seq in seqs.values():
        cuts = [-1]
        for i, aa in enumerate(seq):
            if aa in 'KR' and (i + 1 >= len(seq) or seq[i + 1] != 'P'):
                cuts.append(i)
        cuts.append(len(seq) - 1)
        for i in range(len(cuts) - 1):
            pep = seq[cuts[i] + 1: cuts[i + 1] + 1]
            if len(pep) >= 6:
                canon.add(pep)
    return canon


def mutant_peptides(seq, pos_1based, alt_aa):
    """Tryptic peptides (0 and 1 missed cleavage, >=6 aa) covering the swap,
    computed on the MUTATED sequence."""
    if pos_1based < 1 or pos_1based > len(seq):
        return []
    mseq = seq[:pos_1based - 1] + alt_aa + seq[pos_1based:]
    cuts = [-1]
    for i, aa in enumerate(mseq):
        if aa in 'KR' and (i + 1 >= len(mseq) or mseq[i + 1] != 'P'):
            cuts.append(i)
    cuts.append(len(mseq) - 1)
    idx = pos_1based - 1
    out = set()
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 3, len(cuts))):
            start, end = cuts[i] + 1, cuts[j] + 1
            if start <= idx < end:
                pep = mseq[start:end]
                if len(pep) >= 6:
                    out.add(pep)
    return list(out)


# ── Contact maps (AlphaFold Ca-Ca) ───────────────────────────────────────────
class ContactMaps:
    def __init__(self, contact_dir, dist_threshold):
        self.dir = Path(contact_dir)
        self.thr = dist_threshold
        self._cache = {}
        # index available accessions once (AF-{acc}-*F1.npy with matching .csv)
        self._avail = set()
        for npy in self.dir.glob('AF-*F1.npy'):
            if npy.with_suffix('.csv').exists():
                m = re.match(r'AF-([A-Z0-9]+)-', npy.name)
                if m:
                    self._avail.add(m.group(1))

    def _load(self, acc):
        for npy in sorted(self.dir.glob(f'AF-{acc}-*F1.npy')):
            csv = npy.with_suffix('.csv')
            if not csv.exists():
                continue
            try:
                meta = pd.read_csv(csv, index_col=0)
                dm = np.load(npy)
                p2i = {int(r['id']): i for i, r in meta.iterrows() if pd.notna(r['id'])}
                return p2i, dm
            except Exception:
                continue
        return None, None

    def has(self, acc):
        return acc in self._avail

    def _get(self, acc):
        if acc not in self._cache:
            self._cache[acc] = self._load(acc)
        return self._cache[acc]

    def nearby(self, acc, pos):
        p2i, dm = self._get(acc)
        if dm is None:
            return []
        idx = p2i.get(pos)
        if idx is None:
            return []
        pos_arr = np.array(list(p2i.keys()), dtype=np.int32)
        idx_arr = np.array(list(p2i.values()), dtype=np.int32)
        d = dm[idx][idx_arr]
        mask = (d > 0) & (d < self.thr) & (pos_arr != pos)
        return pos_arr[mask].tolist()


# ── EVE double-mutant matrix -> endpoint index ───────────────────────────────
def load_eve_files(eve_dir, beta):
    files = {}
    ensp = {}
    pat = re.compile(r'destab_(ENSP\d+)_([A-Z0-9]+)_b' + re.escape(beta)
                     + r'_double_mutant_matrix\.csv$')
    for fp in Path(eve_dir).glob(f'destab_*_b{beta}_double_mutant_matrix.csv'):
        m = pat.match(fp.name)
        if m:
            files[m.group(2)] = str(fp)
            ensp[m.group(2)] = m.group(1)
    return files, ensp


def build_eve_index(fp):
    """(pos, aa) -> list of (partner_pos, partner_aa, dmm, self_single, partner_single)."""
    d = pd.read_csv(fp, usecols=['i', 'aa_i', 'j', 'aa_j', 'i_score', 'j_score', 'dmm_score'])
    idx = defaultdict(list)
    for i, ai, j, aj, si, sj, dm in zip(d['i'], d['aa_i'], d['j'], d['aa_j'],
                                        d['i_score'], d['j_score'], d['dmm_score']):
        i, j = int(i), int(j)
        idx[(i, str(ai))].append((j, str(aj), float(dm), float(si), float(sj)))
        idx[(j, str(aj))].append((i, str(ai), float(dm), float(sj), float(si)))
    return idx


# ── Metadata: patients, TMT sets, carrier/non-carrier accounting ─────────────
TMT_CHANNEL_MAP = {
    'tmt_126':'126','tmt_127n':'127N','tmt_127c':'127C','tmt_128n':'128N',
    'tmt_128c':'128C','tmt_129n':'129N','tmt_129c':'129C','tmt_130n':'130N',
    'tmt_130c':'130C','tmt_131':'131N','tmt_131c':'131C','tmt_126c':'126C',
    'tmt_134n':'134N',
}


def build_plex_uuidsets(tmt, gdc, processed_uuids):
    """run_metadata_id -> set of genomics-processed GDC UUIDs in that plex."""
    gdc_map = gdc[['gdc_file_id', 'case_submitter_id', 'sample_type']].drop_duplicates()
    plex_uuids = {}
    for plex_id, grp in tmt.groupby('run_metadata_id'):
        pt = grp[['tmt_channel', 'case_submitter_id', 'sample_type']].drop_duplicates()
        pt = pt[~pt['case_submitter_id'].astype(str).str.lower()
                .isin(['ref', 'reference', 'pooled', 'pool', 'nan'])]
        pt = pt.assign(channel=pt['tmt_channel'].map(TMT_CHANNEL_MAP)).dropna(subset=['channel'])
        merged = pt.merge(gdc_map, on=['case_submitter_id', 'sample_type'], how='left')
        uuids = {u for u in merged['gdc_file_id'].dropna().unique() if u in processed_uuids}
        if uuids:
            plex_uuids[plex_id] = uuids
    return plex_uuids


def powered(carrier_uuids, plex_uuids, min_c, min_nc):
    """Return (is_eligible, n_usable_sets, best_min_balance)."""
    n_sets = 0
    best = 0
    for uuids in plex_uuids.values():
        c = len(carrier_uuids & uuids)
        nc = len(uuids) - c
        if c >= min_c and nc >= min_nc:
            n_sets += 1
            best = max(best, min(c, nc))
    return (n_sets > 0), n_sets, best


def main():
    ap = argparse.ArgumentParser()
    scratch = '/scratch/leduc.an/AAS_Evo'
    ap.add_argument('--missense', default=f'{scratch}/ANALYSIS/all_missense_with_spurs.tsv')
    ap.add_argument('--ref-fasta', default=f'{scratch}/SEQ_FILES/uniprot_human_canonical.fasta')
    ap.add_argument('--contact-dir', default=f'{scratch}/SPURS/contact_maps')
    ap.add_argument('--ddg-dir', default=f'{scratch}/SPURS/ddg_matrices')
    ap.add_argument('--eve-dir',
                    default='/projects/marubi/collabs/slavov_rna/aa_mut_preds/eve_alns_double_muts/destab')
    ap.add_argument('--tmt-map', default=str(REPO_DIR / 'metadata/PDC_meta/pdc_file_tmt_map.tsv'))
    ap.add_argument('--gdc-meta', default=str(REPO_DIR / 'metadata/GDC_meta/gdc_meta_matched.tsv'))
    ap.add_argument('-o', '--out-dir', default=f'{scratch}/ANALYSIS/eve_aas_targets')
    ap.add_argument('--eve-beta', default='0.1')
    ap.add_argument('--dist-threshold', type=float, default=3.0)
    ap.add_argument('--min-seq-sep', type=int, default=21)
    ap.add_argument('--vaf-threshold', type=float, default=0.3)
    ap.add_argument('--min-carrier', type=int, default=2)
    ap.add_argument('--min-noncarrier', type=int, default=2)
    ap.add_argument('--n-syn', type=int, default=2500)
    ap.add_argument('--n-del', type=int, default=2500)
    ap.add_argument('--n-null', type=int, default=5000)
    ap.add_argument('--null-abs-max', type=float, default=None,
                    help='optional cap on |epistasis| to qualify as null')
    ap.add_argument('--max-per-missense', type=int, default=None,
                    help='cap candidate pairs contributed by a single missense per arm')
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print('Loading reference proteome...', flush=True)
    seqs, gene_to_acc = load_reference(args.ref_fasta, args.ddg_dir)
    canon = build_canonical_peptides(seqs)
    print(f'  {len(seqs):,} seqs | {len(gene_to_acc):,} gene->acc | {len(canon):,} canonical peptides')

    print('Loading EVE matrices...', flush=True)
    eve_files, eve_ensp = load_eve_files(args.eve_dir, args.eve_beta)
    print(f'  {len(eve_files)} precomputed accessions (beta={args.eve_beta})')

    contacts = ContactMaps(args.contact_dir, args.dist_threshold)

    # ── missense: full sample set (genomics-processed) + restricted rows ──────
    print('Loading missense table...', flush=True)
    full = pd.read_csv(args.missense, sep='\t', usecols=['sample_id'], low_memory=False)
    processed_uuids = set(full['sample_id'].unique())
    print(f'  {len(processed_uuids):,} genomics-processed samples')

    miss = pd.read_csv(args.missense, sep='\t', low_memory=False,
                       usecols=['sample_id', 'SYMBOL', 'Protein_position', 'Amino_acids',
                                'VAF', 'gnomADe_AF', 'am_pathogenicity', 'spurs_ddg'])
    miss['VAF'] = pd.to_numeric(miss['VAF'], errors='coerce')
    miss = miss[miss['VAF'] >= args.vaf_threshold]
    miss['acc'] = miss['SYMBOL'].map(gene_to_acc)
    miss = miss.dropna(subset=['acc'])
    # only proteins with a contact map are usable (contacts narrow candidates)
    miss = miss[miss['acc'].apply(contacts.has)]
    miss['pos'] = pd.to_numeric(
        miss['Protein_position'].astype(str).str.split('-').str[0], errors='coerce')
    aa = miss['Amino_acids'].astype(str).str.split('/', expand=True)
    miss['wt_i'] = aa[0].str.strip()
    miss['alt_i'] = aa[1].str.strip() if aa.shape[1] > 1 else None
    miss = miss.dropna(subset=['pos', 'wt_i', 'alt_i'])
    miss = miss[(miss['wt_i'].str.len() == 1) & (miss['alt_i'].str.len() == 1)]
    miss['pos'] = miss['pos'].astype(int)
    print(f'  {len(miss):,} missense rows on contact-mapped proteins '
          f'({miss["SYMBOL"].nunique():,} genes)')

    # unique missense -> carrier uuids + annotations
    grp = miss.groupby(['SYMBOL', 'acc', 'pos', 'wt_i', 'alt_i'])
    uniq = grp.agg(carrier_uuids=('sample_id', lambda s: frozenset(s.unique())),
                   gnomADe_AF=('gnomADe_AF', 'first'),
                   am_pathogenicity=('am_pathogenicity', 'first'),
                   spurs_ddg=('spurs_ddg', 'first')).reset_index()
    print(f'  {len(uniq):,} unique missense', flush=True)

    print('Loading TMT / GDC metadata...', flush=True)
    tmt = pd.read_csv(args.tmt_map, sep='\t')
    gdc = pd.read_csv(args.gdc_meta, sep='\t')
    if 'file_id' in gdc.columns and 'gdc_file_id' not in gdc.columns:
        gdc = gdc.rename(columns={'file_id': 'gdc_file_id'})
    plex_uuids = build_plex_uuidsets(tmt, gdc, processed_uuids)
    print(f'  {len(plex_uuids)} TMT sets with genomics')

    # ── power filter ─────────────────────────────────────────────────────────
    elig_rows = []
    for r in uniq.itertuples(index=False):
        ok, n_sets, best = powered(set(r.carrier_uuids), plex_uuids,
                                   args.min_carrier, args.min_noncarrier)
        if ok:
            elig_rows.append((r.SYMBOL, r.acc, r.pos, r.wt_i, r.alt_i,
                              n_sets, best, len(r.carrier_uuids),
                              r.gnomADe_AF, r.am_pathogenicity, r.spurs_ddg))
    elig = pd.DataFrame(elig_rows, columns=[
        'gene', 'acc', 'miss_pos', 'miss_wt', 'miss_alt', 'n_usable_sets',
        'best_balance', 'n_carrier_samples', 'gnomADe_AF', 'am_pathogenicity', 'spurs_ddg'])
    print(f'\nPowered missense (>= {args.min_carrier}c / {args.min_noncarrier}nc in a set): '
          f'{len(elig):,}  ({elig["gene"].nunique():,} genes)')
    has_eve = elig['acc'].isin(set(eve_files))
    print(f'  on precomputed EVE proteins : {has_eve.sum():,}  '
          f'({elig.loc[has_eve, "gene"].nunique()} genes)')
    print(f'  need EVE computed           : {(~has_eve).sum():,}  '
          f'({elig.loc[~has_eve, "gene"].nunique()} genes)')

    # ── Stage 1: scored candidate pairs from precomputed EVE ─────────────────
    eve_idx_cache = {}
    cand = []
    for r in elig[has_eve].itertuples(index=False):
        acc, mp, mwt, malt = r.acc, r.miss_pos, r.miss_wt, r.miss_alt
        if acc not in eve_idx_cache:
            eve_idx_cache[acc] = build_eve_index(eve_files[acc])
        partners = eve_idx_cache[acc].get((mp, malt))
        if not partners:
            continue
        nearby = set(contacts.nearby(acc, mp))
        if not nearby:
            continue
        seq = seqs.get(acc, '')
        for (jpos, jaa, dmm, ms, ss) in partners:
            if jpos not in nearby:
                continue
            if abs(jpos - mp) < args.min_seq_sep:
                continue
            if jpos < 1 or jpos > len(seq):
                continue
            jwt = seq[jpos - 1]
            if not swap_allowed(jwt, jaa):
                continue
            # detectable + specific tryptic peptide?
            peps = [p for p in mutant_peptides(seq, jpos, jaa) if p not in canon]
            if not peps:
                continue
            cand.append({
                'gene': r.gene, 'acc': acc,
                'miss_pos': mp, 'miss_wt': mwt, 'miss_alt': malt,
                'aas_pos': jpos, 'aas_wt': jwt, 'aas_alt': jaa,
                'epistasis': dmm - ms - ss, 'dmm': dmm,
                'miss_single': ms, 'swap_single': ss,
                'n_usable_sets': r.n_usable_sets, 'best_balance': r.best_balance,
                'n_carrier_samples': r.n_carrier_samples,
                'gnomADe_AF': r.gnomADe_AF, 'am_pathogenicity': r.am_pathogenicity,
                'spurs_ddg': r.spurs_ddg,
                'swap_peptide': max(peps, key=len),
            })
    cand = pd.DataFrame(cand)
    if len(cand):
        cand = cand.drop_duplicates(['acc', 'miss_pos', 'miss_alt', 'aas_pos', 'aas_alt'])
    print(f'\nScored candidate (missense, AAS) pairs: {len(cand):,}')
    cand.to_csv(out / 'eve_aas_candidates.tsv', sep='\t', index=False)

    # ── arm selection ────────────────────────────────────────────────────────
    def _cap(df):
        if args.max_per_missense is None or len(df) == 0:
            return df
        return (df.groupby(['acc', 'miss_pos', 'miss_alt'], group_keys=False)
                  .head(args.max_per_missense))

    targets = pd.DataFrame()
    if len(cand):
        syn = _cap(cand[cand['epistasis'] > 0].sort_values('epistasis', ascending=False)) \
            .head(args.n_syn).assign(arm='synergistic')
        dele = _cap(cand[cand['epistasis'] < 0].sort_values('epistasis', ascending=True)) \
            .head(args.n_del).assign(arm='deleterious')
        nullp = cand.copy()
        nullp['abs_ep'] = nullp['epistasis'].abs()
        if args.null_abs_max is not None:
            nullp = nullp[nullp['abs_ep'] <= args.null_abs_max]
        nullp = _cap(nullp.sort_values('abs_ep')).head(args.n_null).assign(arm='null') \
            .drop(columns='abs_ep')
        targets = pd.concat([syn, dele, nullp], ignore_index=True)
        targets.to_csv(out / 'eve_aas_targets.tsv', sep='\t', index=False)

    # ── Stage 2: proteins needing EVE (eligible missense, no matrix) ──────────
    need = (elig[~has_eve]
            .groupby(['gene', 'acc'])
            .agg(n_missense=('miss_pos', 'nunique'),
                 n_usable_sets=('n_usable_sets', 'max'),
                 best_balance=('best_balance', 'max'))
            .reset_index()
            .sort_values(['n_missense', 'best_balance'], ascending=False))
    # annotate contact-candidate count (how many scored pairs EVE could yield)
    cc = []
    for r in need.itertuples(index=False):
        seq = seqs.get(r.acc, '')
        n_pairs = 0
        for m in elig[(elig['acc'] == r.acc) & (~has_eve)].itertuples(index=False):
            nearby = [j for j in contacts.nearby(r.acc, m.miss_pos)
                      if abs(j - m.miss_pos) >= args.min_seq_sep and 1 <= j <= len(seq)]
            n_pairs += len(nearby)
        cc.append(n_pairs)
    need['n_contact_positions'] = cc
    need.to_csv(out / 'eve_proteins_to_compute.tsv', sep='\t', index=False)

    # ── summary ──────────────────────────────────────────────────────────────
    lines = []
    def _p(s=''):
        print(s); lines.append(s)

    _p('\n' + '=' * 62)
    _p('SELECTION SUMMARY')
    _p('=' * 62)
    _p(f'Powered missense           : {len(elig):,} ({elig["gene"].nunique()} genes)')
    _p(f'  on precomputed EVE       : {has_eve.sum():,}')
    _p(f'  need EVE (stage 2)       : {(~has_eve).sum():,} '
       f'over {need.shape[0]} proteins, {need["n_contact_positions"].sum():,} contact positions')
    _p(f'Scored candidate pairs     : {len(cand):,}')
    if len(targets):
        _p('\nArm sizes (target vs actual):')
        for arm, tgt in [('synergistic', args.n_syn), ('deleterious', args.n_del),
                         ('null', args.n_null)]:
            a = targets[targets['arm'] == arm]
            _p(f'  {arm:<12} {len(a):>6,} / {tgt:<6}  '
               f'unique peptides={a["swap_peptide"].nunique():>5}  '
               f'genes={a["gene"].nunique():>3}  '
               f'epistasis[min/med/max]='
               f'{a["epistasis"].min():+.2f}/{a["epistasis"].median():+.2f}/{a["epistasis"].max():+.2f}'
               if len(a) else f'  {arm:<12} {0:>6} / {tgt}')
        _p(f'\nTotal targets              : {len(targets):,}')
        _p(f'Unique AAS peptides        : {targets["swap_peptide"].nunique():,}')
        _p(f'Unique proteins            : {targets["acc"].nunique():,}')
        _p(f'Unique seeding missense    : '
           f'{targets[["acc","miss_pos","miss_alt"]].drop_duplicates().shape[0]:,}')
    _p('\nTop proteins needing EVE (stage 2):')
    _p(need.head(25).to_string(index=False))
    _p(f'\nWrote -> {out}/')

    (out / 'selection_summary.txt').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
