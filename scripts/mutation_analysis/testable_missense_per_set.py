#!/usr/bin/env python3
"""
Feasibility diagnostic for the carrier-vs-non-carrier RAAS design.

For each TMT set (plex), count how many distinct missense mutations have at
least --min-patients CARRIER channels AND at least --min-patients NON-CARRIER
channels within that set (so a per-event carrier-vs-non-carrier test is
possible). Plots a histogram over TMT sets of that count.

A "missense" is a (gene, protein position, wt, alt) protein-level variant with
VAF >= --vaf-threshold. Carriers = genomics-processed channels in the set whose
GDC sample carries it; non-carriers = the remaining genomics channels in the set.

No EVE / contact filtering — this is the upper bound before those narrow it.
"""

import argparse
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


def main():
    ap = argparse.ArgumentParser()
    scratch = '/scratch/leduc.an/AAS_Evo'
    ap.add_argument('--missense', default=f'{scratch}/ANALYSIS/all_missense_with_spurs.tsv')
    ap.add_argument('--tmt-map', default=str(REPO_DIR / 'metadata/PDC_meta/pdc_file_tmt_map.tsv'))
    ap.add_argument('--gdc-meta', default=str(REPO_DIR / 'metadata/GDC_meta/gdc_meta_matched.tsv'))
    ap.add_argument('-o', '--out-dir', default=f'{scratch}/ANALYSIS/eve_aas_targets')
    ap.add_argument('--min-patients', type=int, default=3)
    ap.add_argument('--vaf-threshold', type=float, default=0.3)
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

    # uuid -> set of missense
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

    # per-set count of testable missense (>=k carriers AND >=k non-carriers)
    per_set = []
    testable_pairs = 0
    testable_missense = set()
    for plex_id, uuids in plex_uuids.items():
        n = len(uuids)
        carrier_count = Counter()
        for uid in uuids:
            for M in uuid_missense.get(uid, ()):
                carrier_count[M] += 1
        cnt = 0
        for M, c in carrier_count.items():
            if c >= k and (n - c) >= k:
                cnt += 1
                testable_pairs += 1
                testable_missense.add(M)
        per_set.append({'plex_id': plex_id, 'n_channels': n, 'n_testable': cnt})

    ps = pd.DataFrame(per_set).sort_values('n_testable', ascending=False)
    ps.to_csv(out / 'testable_missense_per_set.tsv', sep='\t', index=False)

    counts = ps['n_testable'].values
    print('\n' + '=' * 56)
    print(f'TESTABLE MISSENSE per TMT set  (>= {k} carriers AND >= {k} non-carriers)')
    print('=' * 56)
    print(f'TMT sets                         : {len(ps)}')
    print(f'  with >= 1 testable missense    : {(counts > 0).sum()}')
    print(f'testable missense per set  min/median/mean/max : '
          f'{counts.min()} / {np.median(counts):.1f} / {counts.mean():.1f} / {counts.max()}')
    print(f'testable (missense, set) pairs   : {testable_pairs:,}')
    print(f'unique missense testable in >=1 set : {len(testable_missense):,}')
    print(f'\nTop sets:\n{ps.head(10).to_string(index=False)}')

    # histogram over TMT sets
    fig, ax = plt.subplots(figsize=(8, 5))
    hi = int(counts.max()) if len(counts) and counts.max() > 0 else 1
    bins = np.arange(0, hi + 2) - 0.5
    ax.hist(counts, bins=bins, color='#1f77b4', edgecolor='white')
    ax.set_xlabel(f'# missense with >= {k} carriers AND >= {k} non-carriers in the set')
    ax.set_ylabel('# TMT sets')
    ax.set_title(f'Per-TMT-set count of {k}v{k}-testable missense mutations\n'
                 f'({len(ps)} sets; {testable_pairs:,} testable (missense,set) pairs)')
    ax.axvline(counts.mean(), color='k', ls='--', lw=1,
               label=f'mean = {counts.mean():.1f}')
    ax.legend()
    plt.tight_layout()
    png = out / 'testable_missense_per_set.png'
    plt.savefig(png, dpi=150, bbox_inches='tight')
    plt.savefig(png.with_suffix('.pdf'), bbox_inches='tight')
    print(f'\nWrote -> {png}\n        {out}/testable_missense_per_set.tsv')


if __name__ == '__main__':
    main()
