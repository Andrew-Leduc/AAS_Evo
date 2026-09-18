#!/usr/bin/env python3
"""
Count contact-swap peptides IDENTIFIED across all FragPipe results.

A swap PSM is one whose protein accession matches the mock-UniProt swap format
emitted by generate_contact_saap_fastas.py:  {ACC}-{SWAP}-{HASH}  (e.g.
P04637-R273H-9AF2), with a "{gene}-contact" description. We scan every plex's
psm.tsv, pull those PSMs, and report totals + uniques through a small funnel:

  1. Identified (raw)        - all swap PSMs passing FragPipe 1% FDR
  2. After PTM-mass filter   - drop swaps whose delta-mass ~= a common PTM
                               (same _PTM_MASSES / 0.05 Da rule as the notebook)

For each stage we print:
  total PSMs | unique swap peptide sequences | unique (gene, swap) events
"""
import glob
import os
import re
import pandas as pd

RESULTS_BASE = "/scratch/leduc.an/AAS_Evo/MS_SEARCH/results_contact"
OUT_DIR = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"

_AA_MASS = {'A':71.03711,'C':103.00919,'D':115.02694,'E':129.04259,'F':147.06841,
            'G':57.02146,'H':137.05891,'I':113.08406,'K':128.09496,'L':113.08406,
            'M':131.04049,'N':114.04293,'P':97.05276,'Q':128.05858,'R':156.10111,
            'S':87.03203,'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06333}
_PTM_MASSES = [14.01565, 15.99491, 2.01565, 0.98402, 42.01057, 79.96633, 18.01056,
               14.99967, 28.03130, 16.97893, 31.98983]
_PTM_TOL = 0.05

# swap accession signature: ...-{WT}{pos}{ALT}-{4 hex}
SWAP_RE = re.compile(r'\b([A-Z0-9]+)-([A-Z]\d+[A-Z])-[0-9A-F]{4}\b')


def is_suspicious(swap):
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', swap)
    if not m:
        return False
    wt, alt = m.group(1), m.group(3)
    if wt not in _AA_MASS or alt not in _AA_MASS:
        return False
    d = abs(_AA_MASS[alt] - _AA_MASS[wt])
    return any(abs(d - p) < _PTM_TOL for p in _PTM_MASSES)


def main():
    plex_dirs = sorted(d for d in glob.glob(os.path.join(RESULTS_BASE, '*'))
                       if os.path.isdir(d))
    parts = []          # per-file DataFrames of swap PSMs
    n_files = 0
    for pd_dir in plex_dirs:
        psm_files = sorted(glob.glob(os.path.join(pd_dir, '*_1', 'psm.tsv'))) \
                    or sorted(glob.glob(os.path.join(pd_dir, 'psm.tsv')))
        for pf in psm_files:
            n_files += 1
            if n_files % 20 == 0:
                print(f'  {n_files} files, {sum(len(p) for p in parts):,} swap PSMs so far', flush=True)
            try:
                head = pd.read_csv(pf, sep='\t', nrows=0)
            except Exception as e:
                print(f'  WARN cannot read {pf}: {e}')
                continue
            cols = {c.lower(): c for c in head.columns}
            c_prot = cols.get('protein')
            c_pid = cols.get('protein id')
            c_map = cols.get('mapped proteins')
            c_pep = cols.get('peptide') or cols.get('modified peptide')
            use = [c for c in (c_pid, c_prot, c_map, c_pep) if c]
            t = pd.read_csv(pf, sep='\t', usecols=use, dtype=str,
                            low_memory=False).fillna('')
            # combined protein string, one regex-extract for the whole column
            prot_all = t[c_pid] if c_pid else ''
            for extra in (c_prot, c_map):
                if extra:
                    prot_all = prot_all + ',' + t[extra]
            ex = prot_all.str.extract(SWAP_RE.pattern)   # cols 0=acc, 1=swap
            mask = ex[1].notna()
            if not mask.any():
                continue
            sub = pd.DataFrame({
                'acc': ex.loc[mask, 0].values,
                'swap': ex.loc[mask, 1].values,
                'peptide': t.loc[mask, c_pep].values if c_pep else '',
            })
            parts.append(sub)
    df = pd.concat(parts, ignore_index=True) if parts else \
        pd.DataFrame(columns=['acc', 'swap', 'peptide'])
    print(f'scanned {n_files} psm.tsv files across {len(plex_dirs)} plex dirs\n')

    def report(name, sub):
        tot = len(sub)
        uniq_pep = sub['peptide'].nunique()
        uniq_ev = sub.drop_duplicates(['acc', 'swap']).shape[0]
        print(f'{name}')
        print(f'  total swap PSMs .............: {tot:,}')
        print(f'  unique swap peptide seqs ....: {uniq_pep:,}')
        print(f'  unique (acc, swap) events ...: {uniq_ev:,}\n')

    report('1. Identified (raw, FragPipe 1% FDR)', df)
    clean = df[~df['swap'].map(is_suspicious)]
    report('2. After PTM-mass suspicious filter', clean)
    n_susp = df.drop_duplicates(['acc', 'swap'])['swap'].map(is_suspicious).sum()
    print(f'(dropped {n_susp:,} unique swaps flagged as PTM-mass confounds)')

    write_heatmap(clean, 'clean')
    write_heatmap(df, 'raw')


# 20 aa in a chemically-grouped order (nonpolar / polar / acidic / basic / aromatic)
AA_ORDER = list('GAVLIPMCSTNQDEKRHFYW')


def write_heatmap(sub, tag):
    """20x20 WT(X)->ALT(Z) substitution matrix of UNIQUE swaps + PNG heatmap."""
    u = sub.drop_duplicates(['acc', 'swap']).copy()
    u['wt'] = u['swap'].str[0]
    u['alt'] = u['swap'].str[-1]
    mat = (pd.crosstab(u['wt'], u['alt'])
             .reindex(index=AA_ORDER, columns=AA_ORDER, fill_value=0))
    tsv = os.path.join(OUT_DIR, f'swap_heatmap_{tag}.tsv')
    mat.to_csv(tsv, sep='\t')
    print(f'wrote {tsv}  (n unique swaps = {int(mat.values.sum()):,})')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        fig, ax = plt.subplots(figsize=(9, 8))
        M = mat.values.astype(float)
        im = ax.imshow(np.ma.masked_equal(M, 0), cmap='viridis', aspect='equal')
        ax.set_xticks(range(20)); ax.set_xticklabels(AA_ORDER, fontsize=9)
        ax.set_yticks(range(20)); ax.set_yticklabels(AA_ORDER, fontsize=9)
        ax.set_xlabel('ALT (Z)', fontweight='bold')
        ax.set_ylabel('WT (X)', fontweight='bold')
        ax.set_title(f'Identified X→Z swaps ({tag}): '
                     f'{int(M.sum()):,} unique', fontweight='bold')
        for i in range(20):
            for j in range(20):
                v = int(M[i, j])
                if v:
                    ax.text(j, i, v, ha='center', va='center', fontsize=6,
                            color='white' if v < M.max() * 0.6 else 'black')
        fig.colorbar(im, ax=ax, shrink=0.8, label='unique swaps')
        plt.tight_layout()
        png = os.path.join(OUT_DIR, f'swap_heatmap_{tag}.png')
        plt.savefig(png, dpi=200, bbox_inches='tight')
        print(f'wrote {png}')
    except Exception as e:
        print(f'  (heatmap PNG skipped: {e})')


if __name__ == '__main__':
    main()
