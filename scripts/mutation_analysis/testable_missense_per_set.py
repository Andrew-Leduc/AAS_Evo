#!/usr/bin/env python3
"""
Feasibility diagnostic for the carrier-vs-non-carrier RAAS design.

Panel A: for each TMT set (plex), how many distinct missense mutations have at
    least --min-patients CARRIER channels AND at least --min-patients NON-CARRIER
    channels within that set (so a per-event carrier-vs-non-carrier test is
    possible). Histogram over TMT sets.

Panel B: for the SUBSET of testable missense (>=k carriers & >=k non-carriers in
    >=1 set), how many contact-proximal positions (Ca-Ca < --dist-threshold) sit
    at least --min-seq-sep residues away in sequence — i.e. how many AAS
    candidate sites each testable missense actually has. Histogram over that
    subset. Shows how tied our hands are once contacts are required.

A "missense" is a (gene, protein position, wt, alt) protein-level variant with
VAF >= --vaf-threshold.  No EVE filtering — this is the ceiling before EVE.
"""

import argparse
import re
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_DIR = Path(__file__).resolve().parents[2]

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


def build_gene_to_acc(ref_fasta, ddg_dir):
    g2a = {}
    if ddg_dir and Path(ddg_dir).exists():
        for f in Path(ddg_dir).glob('*.ddg_matrix.tsv'):
            parts = f.name.split('.')
            if len(parts) >= 2:
                g2a[parts[1]] = parts[0]
    if not (ref_fasta and Path(ref_fasta).exists()):
        return g2a
    with open(ref_fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                p = line.split('|')
                acc = p[1] if len(p) >= 3 else line[1:].split()[0]
                m = re.search(r'GN=(\S+)', line)
                if m and m.group(1) not in g2a:
                    g2a[m.group(1)] = acc
    return g2a


class ContactMaps:
    def __init__(self, contact_dir, dist_threshold):
        self.dir = Path(contact_dir)
        self.thr = dist_threshold
        self._cache = {}
        self._avail = set()
        for npy in self.dir.glob('AF-*F1.npy'):
            if npy.with_suffix('.csv').exists():
                m = re.match(r'AF-([A-Z0-9]+)-', npy.name)
                if m:
                    self._avail.add(m.group(1))

    def has(self, acc):
        return acc in self._avail

    def _get(self, acc):
        if acc not in self._cache:
            res = (None, None)
            for npy in sorted(self.dir.glob(f'AF-{acc}-*F1.npy')):
                csv = npy.with_suffix('.csv')
                if not csv.exists():
                    continue
                try:
                    meta = pd.read_csv(csv, index_col=0)
                    dm = np.load(npy)
                    p2i = {int(r['id']): i for i, r in meta.iterrows() if pd.notna(r['id'])}
                    res = (p2i, dm)
                    break
                except Exception:
                    continue
            self._cache[acc] = res
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


def main():
    ap = argparse.ArgumentParser()
    scratch = '/scratch/leduc.an/AAS_Evo'
    ap.add_argument('--missense', default=f'{scratch}/ANALYSIS/all_missense_with_spurs.tsv')
    ap.add_argument('--tmt-map', default=str(REPO_DIR / 'metadata/PDC_meta/pdc_file_tmt_map.tsv'))
    ap.add_argument('--gdc-meta', default=str(REPO_DIR / 'metadata/GDC_meta/gdc_meta_matched.tsv'))
    ap.add_argument('--ref-fasta', default=f'{scratch}/SEQ_FILES/uniprot_human_canonical.fasta')
    ap.add_argument('--contact-dir', default=f'{scratch}/SPURS/contact_maps')
    ap.add_argument('--ddg-dir', default=f'{scratch}/SPURS/ddg_matrices')
    ap.add_argument('-o', '--out-dir', default=f'{scratch}/ANALYSIS/eve_aas_targets')
    ap.add_argument('--min-patients', type=int, default=3)
    ap.add_argument('--vaf-threshold', type=float, default=0.3)
    ap.add_argument('--dist-threshold', type=float, default=3.0)
    ap.add_argument('--min-seq-sep', type=int, default=30)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    k = args.min_patients

    print('Loading missense...', flush=True)
    miss = pd.read_csv(args.missense, sep='\t', low_memory=False,
                       usecols=['sample_id', 'SYMBOL', 'Protein_position', 'Amino_acids', 'VAF'])
    processed_uuids = set(miss['sample_id'].unique())
    miss['VAF'] = pd.to_numeric(miss['VAF'], errors='coerce')
    miss = miss[miss['VAF'] >= args.vaf_threshold]
    miss['pos'] = pd.to_numeric(
        miss['Protein_position'].astype(str).str.split('-').str[0], errors='coerce')
    aa = miss['Amino_acids'].astype(str).str.split('/', expand=True)
    miss['wt'] = aa[0].str.strip()
    miss['alt'] = aa[1].str.strip() if aa.shape[1] > 1 else None
    miss = miss.dropna(subset=['pos', 'wt', 'alt'])
    miss = miss[(miss['wt'].str.len() == 1) & (miss['alt'].str.len() == 1)]
    miss['pos'] = miss['pos'].astype(int)
    miss['M'] = list(zip(miss['SYMBOL'], miss['pos'], miss['wt'], miss['alt']))
    print(f'  {len(processed_uuids):,} genomics samples | '
          f'{len(miss):,} missense rows | {miss["M"].nunique():,} unique missense')

    uuid_missense = defaultdict(set)
    for uid, M in zip(miss['sample_id'], miss['M']):
        uuid_missense[uid].add(M)

    print('Loading TMT / GDC metadata...', flush=True)
    tmt = pd.read_csv(args.tmt_map, sep='\t')
    gdc = pd.read_csv(args.gdc_meta, sep='\t')
    if 'file_id' in gdc.columns and 'gdc_file_id' not in gdc.columns:
        gdc = gdc.rename(columns={'file_id': 'gdc_file_id'})
    plex_uuids = build_plex_uuidsets(tmt, gdc, processed_uuids)
    print(f'  {len(plex_uuids)} TMT sets with genomics')

    # ── Panel A: per-set testable missense ──────────────────────────────────
    per_set = []
    testable_pairs = 0
    testable_missense = set()
    obs_carr, obs_nonc = [], []
    for plex_id, uuids in plex_uuids.items():
        n = len(uuids)
        carrier_count = Counter()
        for uid in uuids:
            for M in uuid_missense.get(uid, ()):
                carrier_count[M] += 1
        cnt = 0
        for M, c in carrier_count.items():
            nc = n - c
            obs_carr.append(c)
            obs_nonc.append(nc)
            if c >= k and nc >= k:
                cnt += 1
                testable_pairs += 1
                testable_missense.add(M)
        per_set.append({'plex_id': plex_id, 'n_channels': n, 'n_testable': cnt})
    obs_carr = np.array(obs_carr)
    obs_nonc = np.array(obs_nonc)

    ps = pd.DataFrame(per_set).sort_values('n_testable', ascending=False)
    ps.to_csv(out / 'testable_missense_per_set.tsv', sep='\t', index=False)
    counts = ps['n_testable'].values

    print('\n' + '=' * 60)
    print(f'PANEL A — testable missense per TMT set (>= {k}c AND >= {k}nc)')
    print('=' * 60)
    print(f'TMT sets                          : {len(ps)}')
    print(f'  with >= 1 testable missense     : {(counts > 0).sum()}')
    print(f'per set  min/median/mean/max      : '
          f'{counts.min()} / {np.median(counts):.1f} / {counts.mean():.1f} / {counts.max()}')
    print(f'testable (missense, set) pairs    : {testable_pairs:,}')
    print(f'unique missense testable in >=1 set : {len(testable_missense):,}')

    # ── sanity: genomics-processed samples per TMT set ──────────────────────
    n_plex_total = tmt['run_metadata_id'].nunique()
    nch = ps['n_channels'].values
    print('\n' + '-' * 60)
    print('GENOMICS SAMPLES PER TMT SET (coverage sanity check)')
    print('-' * 60)
    print(f'TMT sets in map                   : {n_plex_total}')
    print(f'  with >= 1 genomics sample       : {len(ps)}')
    print(f'  with >= 6 (min for a 3v3 split) : {(nch >= 6).sum()}')
    print(f'genomics samples/set min/median/mean/max : '
          f'{nch.min()} / {np.median(nch):.1f} / {nch.mean():.1f} / {nch.max()}')
    print('top sets by genomics samples:')
    print(ps.sort_values('n_channels', ascending=False)
            [['plex_id', 'n_channels', 'n_testable']].head(8).to_string(index=False))
    figc, axc = plt.subplots(figsize=(8, 5))
    axc.hist(nch, bins=np.arange(0, int(nch.max()) + 2) - 0.5,
             color='#9467bd', edgecolor='white')
    axc.axvline(6, color='r', ls='--', lw=1.2, label='3v3 needs >= 6')
    axc.set_xlabel('# genomics-processed samples (channels) in the set')
    axc.set_ylabel('# TMT sets')
    axc.set_title(f'Genomics samples per TMT set\n'
                  f'({len(ps)}/{n_plex_total} sets have >=1; '
                  f'{(nch >= 6).sum()} have >=6)')
    axc.legend()
    plt.tight_layout()
    figc.savefig(out / 'samples_per_set.png', dpi=150, bbox_inches='tight')
    plt.close(figc)
    print(f'Wrote -> {out}/samples_per_set.png')

    # ── breakdown: carrier / non-carrier counts across ALL (missense, set) ──
    print('\n' + '-' * 60)
    print('BREAKDOWN over ALL (missense, set) observations (missense present)')
    print('-' * 60)
    tot = len(obs_carr)
    print(f'(missense, set) observations       : {tot:,}')
    print(f'{"threshold":>9} | {">= carriers":>12} | {">= non-carriers":>15}')
    for t in [1, 2, 3, 4, 5]:
        print(f'{t:>9} | {int((obs_carr >= t).sum()):>12,} | '
              f'{int((obs_nonc >= t).sum()):>15,}')
    print(f'>= {k} carriers AND >= {k} non-carriers : {testable_pairs:,}')

    figd, axd = plt.subplots(1, 2, figsize=(14, 5))
    for ax, data, lab, col in [
            (axd[0], obs_carr, 'carriers', '#d62728'),
            (axd[1], obs_nonc, 'non-carriers', '#1f77b4')]:
        cap = min(int(data.max()) if len(data) else 1, 12)
        ax.hist(np.clip(data, 0, cap), bins=np.arange(0, cap + 2) - 0.5,
                color=col, edgecolor='white')
        ax.axvline(k, color='k', ls='--', lw=1.2, label=f'threshold = {k}')
        ax.set_xlabel(f'# {lab} in the set (clipped at {cap})')
        ax.set_ylabel('# (missense, set) observations (log)')
        ax.set_title(f'{lab} per (missense, set)')
        ax.set_yscale('log')
        ax.legend()
    plt.tight_layout()
    figd.savefig(out / 'carrier_breakdown.png', dpi=150, bbox_inches='tight')
    plt.close(figd)
    print(f'Wrote -> {out}/carrier_breakdown.png')

    # ── Panel B: contact candidates for the testable subset ─────────────────
    print('\nLoading contact maps / gene->acc ...', flush=True)
    gene_to_acc = build_gene_to_acc(args.ref_fasta, args.ddg_dir)
    contacts = ContactMaps(args.contact_dir, args.dist_threshold)

    rows = []
    n_no_acc = n_no_map = 0
    for (gene, pos, wt, alt) in testable_missense:
        acc = gene_to_acc.get(gene)
        if not acc:
            n_no_acc += 1
            rows.append((gene, pos, wt, alt, '', 0, False, False))
            continue
        if not contacts.has(acc):
            n_no_map += 1
            rows.append((gene, pos, wt, alt, acc, 0, False, False))
            continue
        nearby = contacts.nearby(acc, pos)
        n_far = sum(1 for j in nearby if abs(j - pos) >= args.min_seq_sep)
        rows.append((gene, pos, wt, alt, acc, n_far, True, n_far > 0))

    cm = pd.DataFrame(rows, columns=['gene', 'pos', 'wt', 'alt', 'acc',
                                     'n_contacts_ge_sep', 'has_map', 'has_candidate'])
    cm.to_csv(out / 'testable_missense_contacts.tsv', sep='\t', index=False)

    n_tot = len(cm)
    n_map = int(cm['has_map'].sum())
    n_cand = int(cm['has_candidate'].sum())
    far_counts = cm.loc[cm['has_map'], 'n_contacts_ge_sep'].values
    print('\n' + '=' * 60)
    print(f'PANEL B — contact candidates (Ca-Ca < {args.dist_threshold}A, '
          f'>= {args.min_seq_sep} aa away) for testable missense')
    print('=' * 60)
    print(f'testable missense                 : {n_tot:,}')
    print(f'  no gene->acc mapping            : {n_no_acc:,}')
    print(f'  no AlphaFold contact map        : {n_no_map:,}')
    print(f'  with contact map                : {n_map:,}')
    print(f'  with >= 1 long-range candidate  : {n_cand:,}  '
          f'({100 * n_cand / n_tot:.1f}% of testable)')
    if len(far_counts):
        print(f'candidates/missense (mapped) min/median/mean/max : '
              f'{far_counts.min()} / {np.median(far_counts):.1f} / '
              f'{far_counts.mean():.1f} / {far_counts.max()}')

    # ── figure: Panel A (+ Panel B if contacts available) ───────────────────
    two = len(far_counts) > 0
    fig, axes = plt.subplots(1, 2 if two else 1,
                             figsize=(14, 5) if two else (8, 5), squeeze=False)
    axA = axes[0][0]
    hi = int(counts.max()) if len(counts) and counts.max() > 0 else 1
    axA.hist(counts, bins=np.arange(0, hi + 2) - 0.5,
             color='#1f77b4', edgecolor='white')
    axA.axvline(counts.mean(), color='k', ls='--', lw=1,
                label=f'mean = {counts.mean():.1f}')
    axA.set_xlabel(f'# missense with >= {k} carriers AND >= {k} non-carriers')
    axA.set_ylabel('# TMT sets')
    axA.set_title(f'A. {k}v{k}-testable missense per TMT set\n'
                  f'({len(ps)} sets; {testable_pairs:,} testable (missense,set) pairs)')
    axA.legend()

    if two:
        axB = axes[0][1]
        cap = min(int(far_counts.max()), 40)
        clipped = np.clip(far_counts, 0, cap)
        axB.hist(clipped, bins=np.arange(0, cap + 2) - 0.5,
                 color='#2ca02c', edgecolor='white')
        axB.set_xlabel(f'# contact positions >= {args.min_seq_sep} aa away '
                       f'(clipped at {cap})')
        axB.set_ylabel('# testable missense (with contact map)')
        axB.set_title(f'B. AAS candidate sites per testable missense\n'
                      f'{n_cand:,}/{n_tot:,} have >=1 candidate; '
                      f'{n_no_map:,} lack a contact map')
    plt.tight_layout()
    png = out / 'testable_missense_feasibility.png'
    plt.savefig(png, dpi=150, bbox_inches='tight')
    plt.savefig(png.with_suffix('.pdf'), bbox_inches='tight')
    print(f'\nWrote -> {png}\n        {out}/testable_missense_per_set.tsv'
          f'\n        {out}/testable_missense_contacts.tsv')


if __name__ == '__main__':
    main()
