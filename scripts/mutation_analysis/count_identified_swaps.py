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


def swap_from_row(prot, prot_id, mapped, gene_desc):
    """Return (gene, swap, acc) if any protein field is a swap accession, else None."""
    for field in (prot_id, prot, mapped):
        if not isinstance(field, str):
            continue
        for chunk in field.split(','):
            m = SWAP_RE.search(chunk.strip())
            if m:
                return m.group(1), m.group(2)      # acc, swap
    return None


def main():
    plex_dirs = sorted(d for d in glob.glob(os.path.join(RESULTS_BASE, '*'))
                       if os.path.isdir(d))
    rows = []           # (gene/acc, swap, peptide, spectrum)
    n_files = 0
    for pd_dir in plex_dirs:
        psm_files = sorted(glob.glob(os.path.join(pd_dir, '*_1', 'psm.tsv'))) \
                    or sorted(glob.glob(os.path.join(pd_dir, 'psm.tsv')))
        for pf in psm_files:
            n_files += 1
            try:
                t = pd.read_csv(pf, sep='\t', low_memory=False)
            except Exception as e:
                print(f'  WARN cannot read {pf}: {e}')
                continue
            cols = {c.lower(): c for c in t.columns}
            c_prot = cols.get('protein')
            c_pid = cols.get('protein id')
            c_map = cols.get('mapped proteins')
            c_pep = cols.get('peptide') or cols.get('modified peptide')
            c_spec = cols.get('spectrum')
            for r in t.itertuples(index=False):
                d = r._asdict()
                hit = swap_from_row(d.get(c_prot), d.get(c_pid), d.get(c_map), None)
                if hit is None:
                    continue
                acc, swap = hit
                rows.append((acc, swap,
                             d.get(c_pep, ''), d.get(c_spec, '')))
    df = pd.DataFrame(rows, columns=['acc', 'swap', 'peptide', 'spectrum'])
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


if __name__ == '__main__':
    main()
