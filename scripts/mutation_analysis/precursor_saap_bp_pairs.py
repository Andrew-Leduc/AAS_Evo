#!/usr/bin/env python3
"""
Foundation for the SAAP validation analyses (#1, #2, and the RAAS axis for #3-#6).

For every identified contact-swap (SAAP) peptide, find its base peptide (BP =
the WT tryptic peptide covering the same position) and pair them PER MS RUN
(RAW file / fraction). Precursor-level, MS1 apex intensity (psm.tsv `Intensity`),
NOT reporter ions.

Per (SAAP, BP, run) we record:
  - saap_intensity, bp_intensity  -> raas_precursor = log2(saap/bp)
  - saap_rt, bp_rt               -> rt_shift = saap_rt - bp_rt
  - run (fraction), plex

Outputs:
  contact_saap/precursor_pairs.tsv     one row per (acc,pos,wt,alt,run) with a BP in the SAME run
  contact_saap/coid_summary.tsv        per SAAP: #runs seen, #runs BP also seen (answers #1)
  + prints the #1 co-identification summary

Downstream annotation (codon mismatch, codon freq, phyloP, gnomAD) joins on
(acc,pos,wt,alt); see build_position_codon_map.py.
"""
import glob
import os
import re
import numpy as np
import pandas as pd

RESULTS_BASE = "/scratch/leduc.an/AAS_Evo/MS_SEARCH/results_contact"
REF_FASTA = "/scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta"
OUT_DIR = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"

SWAP_RE = re.compile(r'\b([A-Z0-9]+)-([A-Z]\d+[A-Z])-[0-9A-F]{4}\b')
MIN_PEP_LEN = 6


# ── reference sequences + tryptic BP peptides ────────────────────────────────
def load_ref(path):
    seqs, cur = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                parts = line.split('|')
                cur = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return {a: ''.join(s) for a, s in seqs.items()}


def tryptic_cuts(seq):
    cuts = [-1]
    for i, aa in enumerate(seq):
        if aa in 'KR' and (i + 1 >= len(seq) or seq[i + 1] != 'P'):
            cuts.append(i)
    cuts.append(len(seq) - 1)
    return cuts


def wt_peptides_covering(seq, pos1):
    """WT tryptic peptides (0 and 1 missed cleavage) covering 1-based pos1."""
    cuts = tryptic_cuts(seq)
    idx = pos1 - 1
    out = set()
    for mc in (0, 1):
        for k in range(len(cuts) - 1 - mc):
            start = cuts[k] + 1
            end = cuts[k + 1 + mc] + 1
            if start <= idx < end:
                pep = seq[start:end]
                if len(pep) >= MIN_PEP_LEN:
                    out.add(pep)
    return out


def parse_run(spectrum):
    """FragPipe Spectrum = '<rawbasename>.<scan>.<scan>.<charge>' -> rawbasename."""
    if not isinstance(spectrum, str):
        return ''
    return spectrum.rsplit('.', 3)[0]


def main():
    seqs = load_ref(REF_FASTA)
    print(f'loaded {len(seqs):,} reference sequences')

    plex_dirs = sorted(d for d in glob.glob(os.path.join(RESULTS_BASE, '*'))
                       if os.path.isdir(d))

    # First pass: collect identified swaps -> the (acc,pos) sites we need BPs for.
    # Build BP-sequence -> (acc,pos) map lazily as we discover sites.
    bp_map = {}            # bp_seq -> set of (acc,pos)
    site_wt = {}           # (acc,pos) -> wt_aa
    saap_rows = []         # (acc,pos,wt,alt,run,plex,intensity,rt)
    bp_hits = []           # (acc,pos,run,plex,intensity,rt) — filled 2nd pass
    n_files = 0

    def ensure_site(acc, pos, wt):
        key = (acc, pos)
        if key in site_wt:
            return
        site_wt[key] = wt
        seq = seqs.get(acc, '')
        if not seq or pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
            return
        for pep in wt_peptides_covering(seq, pos):
            bp_map.setdefault(pep, set()).add(key)

    # ---- pass 1: swaps + register BP sites ----
    per_file_cols = {}
    for pd_dir in plex_dirs:
        plex = os.path.basename(pd_dir)
        psm_files = sorted(glob.glob(os.path.join(pd_dir, '*_1', 'psm.tsv'))) \
                    or sorted(glob.glob(os.path.join(pd_dir, 'psm.tsv')))
        for pf in psm_files:
            n_files += 1
            head = pd.read_csv(pf, sep='\t', nrows=0)
            cols = {c.lower(): c for c in head.columns}
            c_pid = cols.get('protein id'); c_prot = cols.get('protein')
            c_map = cols.get('mapped proteins'); c_pep = cols.get('peptide')
            c_int = cols.get('intensity'); c_rt = cols.get('retention')
            c_spec = cols.get('spectrum')
            per_file_cols[pf] = (plex, c_pid, c_prot, c_map, c_pep, c_int, c_rt, c_spec)
            use = [c for c in (c_pid, c_prot, c_map, c_pep) if c]
            t = pd.read_csv(pf, sep='\t', usecols=use, dtype=str).fillna('')
            prot_all = t[c_pid] if c_pid else pd.Series([''] * len(t))
            for extra in (c_prot, c_map):
                if extra:
                    prot_all = prot_all + ',' + t[extra]
            ex = prot_all.str.extract(SWAP_RE.pattern)
            for acc, swap in ex.dropna().itertuples(index=False):
                m = re.match(r'^([A-Z])(\d+)([A-Z])$', swap)
                if m:
                    ensure_site(acc, int(m.group(2)), m.group(1))
    print(f'registered {len(site_wt):,} contact sites; {len(bp_map):,} BP peptide seqs')

    # ---- pass 2: pull intensities/RT for swaps and BPs ----
    for pf, (plex, c_pid, c_prot, c_map, c_pep, c_int, c_rt, c_spec) in per_file_cols.items():
        use = [c for c in (c_pid, c_prot, c_map, c_pep, c_int, c_rt, c_spec) if c]
        t = pd.read_csv(pf, sep='\t', usecols=use, dtype=str).fillna('')
        inten = pd.to_numeric(t[c_int], errors='coerce') if c_int else np.nan
        rt = pd.to_numeric(t[c_rt], errors='coerce') if c_rt else np.nan
        run = t[c_spec].map(parse_run) if c_spec else ''
        prot_all = t[c_pid] if c_pid else pd.Series([''] * len(t))
        for extra in (c_prot, c_map):
            if extra:
                prot_all = prot_all + ',' + t[extra]
        ex = prot_all.str.extract(SWAP_RE.pattern)
        for i in range(len(t)):
            sw = ex.iat[i, 1]
            if isinstance(sw, str):
                m = re.match(r'^([A-Z])(\d+)([A-Z])$', sw)
                if m:
                    saap_rows.append((ex.iat[i, 0], int(m.group(2)), m.group(1),
                                      m.group(3), run.iat[i], plex,
                                      inten.iat[i], rt.iat[i]))
                    continue
            sites = bp_map.get(t[c_pep].iat[i])
            if sites:
                for acc, pos in sites:
                    bp_hits.append((acc, pos, run.iat[i], plex,
                                    inten.iat[i], rt.iat[i]))

    saap = pd.DataFrame(saap_rows, columns=['acc', 'pos', 'wt', 'alt', 'run',
                                            'plex', 'intensity', 'rt'])
    bp = pd.DataFrame(bp_hits, columns=['acc', 'pos', 'run', 'plex',
                                        'intensity', 'rt'])
    print(f'swap PSMs: {len(saap):,} | BP PSMs: {len(bp):,}')

    # aggregate to one value per (site[,alt], run): sum intensity, intensity-wt RT
    def agg(df, keys):
        df = df.dropna(subset=['intensity'])
        g = df.groupby(keys)
        out = g['intensity'].sum().rename('intensity').reset_index()
        rtw = (g.apply(lambda x: np.average(x['rt'], weights=x['intensity'])
                       if x['intensity'].sum() > 0 else x['rt'].median())
               .rename('rt').reset_index())
        return out.merge(rtw, on=keys)

    saap_a = agg(saap, ['acc', 'pos', 'wt', 'alt', 'run', 'plex'])
    bp_a = agg(bp, ['acc', 'pos', 'run', 'plex']).rename(
        columns={'intensity': 'bp_intensity', 'rt': 'bp_rt'})

    pairs = saap_a.merge(bp_a, on=['acc', 'pos', 'run', 'plex'], how='left')
    pairs = pairs.rename(columns={'intensity': 'saap_intensity', 'rt': 'saap_rt'})
    same = pairs['bp_intensity'].notna()
    pairs.loc[same, 'raas_precursor'] = np.log2(
        pairs.loc[same, 'saap_intensity'] / pairs.loc[same, 'bp_intensity'])
    pairs.loc[same, 'rt_shift'] = pairs.loc[same, 'saap_rt'] - pairs.loc[same, 'bp_rt']

    os.makedirs(OUT_DIR, exist_ok=True)
    out_pairs = pairs[same].copy()
    out_pairs['swap'] = out_pairs['wt'] + out_pairs['pos'].astype(str) + out_pairs['alt']
    out_pairs.to_csv(os.path.join(OUT_DIR, 'precursor_pairs.tsv'), sep='\t', index=False)

    # ── #1 co-identification summary ──
    runs_seen = saap_a.groupby(['acc', 'pos', 'wt', 'alt'])['run'].nunique()
    runs_bpsame = (pairs[same].groupby(['acc', 'pos', 'wt', 'alt'])['run']
                   .nunique())
    coid = pd.DataFrame({'runs_saap_seen': runs_seen}).join(
        runs_bpsame.rename('runs_bp_same_run')).fillna(0)
    coid['frac_runs_with_bp'] = coid['runs_bp_same_run'] / coid['runs_saap_seen']
    coid.reset_index().to_csv(os.path.join(OUT_DIR, 'coid_summary.tsv'),
                              sep='\t', index=False)

    tot_obs = int(runs_seen.sum())
    tot_bp = int(pairs[same].shape[0])
    print('\n#1  SAME-RUN CO-IDENTIFICATION')
    print(f'  unique SAAPs (acc,pos,wt,alt) ....: {len(coid):,}')
    print(f'  SAAP x run observations ..........: {tot_obs:,}')
    print(f'  of those, BP in the SAME run .....: {tot_bp:,} '
          f'({100*tot_bp/max(tot_obs,1):.1f}%)')
    print(f'  SAAPs with >=1 same-run BP .......: {int((coid["runs_bp_same_run"]>0).sum()):,} '
          f'({100*(coid["runs_bp_same_run"]>0).mean():.1f}%)')
    print(f'  median frac of a SAAP\'s runs w/ BP: {coid["frac_runs_with_bp"].median():.2f}')
    print(f'\nwrote {OUT_DIR}/precursor_pairs.tsv ({len(out_pairs):,} same-run pairs)')
    print(f'wrote {OUT_DIR}/coid_summary.tsv')


if __name__ == '__main__':
    main()
