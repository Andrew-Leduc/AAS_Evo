#!/usr/bin/env python3
"""
Regenerate AlphaFold contact maps (reconstruction of the lost contact_make.py).

For each UniProt accession, downloads the AlphaFold structure and writes a
residue-residue distance map in the exact format the pipeline consumes
(select_eve_aas_targets.py, contact_saap_analysis.ipynb, generate_contact_saap_fastas.py):

  AF-{acc}-F1.npy : L x L float32 matrix of MINIMUM HEAVY-ATOM distance (Angstrom)
                    between residues (diagonal = 0).
  AF-{acc}-F1.csv : first column = matrix row index 0..L-1 (this is what the
                    readers use as index_col=0); columns id (1-based residue
                    number), aa (1-letter), plddt (AlphaFold confidence).

Distance metric note: the pipeline's contact threshold is ~3 A, which is far
below the ~3.8 A minimum for Ca-Ca — so contacts are defined on the MINIMUM
HEAVY-ATOM inter-residue distance (a real atomic contact / H-bond), which is
what reproduces the original behaviour. The readers do `d > 0 & d < threshold`
and filter sequence-adjacent pairs separately (MIN_SEQ_SEP), so trivial
backbone-neighbour contacts are handled downstream.

Accessions are derived from the missense table (genes -> UniProt accession via
the reference proteome GN= field), or supplied via --acc-list.

Resumable (skips accessions whose .npy already exists) and shardable for SLURM
arrays via --nshards N (+ SLURM_ARRAY_TASK_ID) or --shard i --nshards N.
"""

import argparse
import os
import re
import sys
import urllib.request
import urllib.error
import numpy as np
import pandas as pd
from pathlib import Path

AF_URL = 'https://alphafold.ebi.ac.uk/files/AF-{acc}-F{frag}-model_v4.pdb'

_THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}


def load_gene_to_acc(ref_fasta):
    g2a = {}
    with open(ref_fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                p = line.split('|')
                acc = p[1] if len(p) >= 3 else line[1:].split()[0]
                m = re.search(r'GN=(\S+)', line)
                if m and m.group(1) not in g2a:
                    g2a[m.group(1)] = acc
    return g2a


def accs_from_missense(missense, gene_to_acc, vaf_threshold):
    miss = pd.read_csv(missense, sep='\t', low_memory=False, usecols=['SYMBOL', 'VAF'])
    miss['VAF'] = pd.to_numeric(miss['VAF'], errors='coerce')
    genes = miss.loc[miss['VAF'] >= vaf_threshold, 'SYMBOL'].dropna().unique()
    accs = {gene_to_acc[g] for g in genes if g in gene_to_acc}
    return sorted(accs)


def download_pdb(acc, cache, frag=1):
    dst = cache / f'AF-{acc}-F{frag}-model_v4.pdb'
    if dst.exists() and dst.stat().st_size > 0:
        return dst
    url = AF_URL.format(acc=acc, frag=frag)
    try:
        urllib.request.urlretrieve(url, dst)
    except (urllib.error.HTTPError, urllib.error.URLError, OSError):
        if dst.exists():
            dst.unlink()
        return None
    return dst


def parse_pdb(pdb):
    """Return (resnums sorted, aa1 list, list-of-atom-coord-arrays, plddt array)."""
    atoms = {}    # resnum -> list of (x,y,z)
    resname = {}  # resnum -> 3-letter
    plddt = {}    # resnum -> B-factor (AlphaFold pLDDT)
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith('ATOM'):
                continue
            element = line[76:78].strip()
            atom = line[12:16].strip()
            if element == 'H' or (not element and atom.startswith('H')):
                continue
            resnum = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            atoms.setdefault(resnum, []).append((x, y, z))
            if resnum not in resname:
                resname[resnum] = line[17:20].strip()
                try:
                    plddt[resnum] = float(line[60:66])
                except ValueError:
                    plddt[resnum] = np.nan
    resnums = sorted(atoms)
    aa1 = [_THREE2ONE.get(resname[r], 'X') for r in resnums]
    coords = [np.asarray(atoms[r], dtype=np.float64) for r in resnums]
    pl = np.array([plddt[r] for r in resnums], dtype=np.float32)
    return resnums, aa1, coords, pl


def min_heavy_atom_matrix(coords):
    """L x L minimum heavy-atom distance between residues."""
    L = len(coords)
    dm = np.zeros((L, L), dtype=np.float32)
    for i in range(L):
        ci = coords[i]
        for j in range(i + 1, L):
            cj = coords[j]
            d = np.sqrt(((ci[:, None, :] - cj[None, :, :]) ** 2).sum(-1)).min()
            dm[i, j] = dm[j, i] = d
    return dm


def make_one(acc, cache, out_dir):
    npy = out_dir / f'AF-{acc}-F1.npy'
    if npy.exists():
        return 'skip'
    pdb = download_pdb(acc, cache)
    if pdb is None:
        return 'no_structure'
    resnums, aa1, coords, pl = parse_pdb(pdb)
    if len(resnums) < 2:
        return 'empty'
    dm = min_heavy_atom_matrix(coords)
    meta = pd.DataFrame({'id': resnums, 'aa': aa1, 'plddt': pl})
    meta.to_csv(out_dir / f'AF-{acc}-F1.csv')   # index 0..L-1 becomes col 0
    np.save(npy, dm)
    return 'ok'


def main():
    ap = argparse.ArgumentParser()
    scratch = '/scratch/leduc.an/AAS_Evo'
    ap.add_argument('--missense', default='/projects/slavov/andrew/AAS_EVO/all_missense_mutations.tsv')
    ap.add_argument('--ref-fasta', default=f'{scratch}/SEQ_FILES/uniprot_human_canonical.fasta')
    ap.add_argument('--acc-list', default=None, help='file with one UniProt acc per line (overrides missense)')
    ap.add_argument('--out-dir', default=f'{scratch}/SPURS/contact_maps')
    ap.add_argument('--pdb-cache', default=f'{scratch}/SPURS/pdb_cache')
    ap.add_argument('--vaf-threshold', type=float, default=0.3)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--shard', type=int, default=None,
                    help='0-based shard index; defaults to SLURM_ARRAY_TASK_ID-1')
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    cache = Path(args.pdb_cache); cache.mkdir(parents=True, exist_ok=True)

    if args.acc_list:
        accs = sorted({l.strip() for l in open(args.acc_list) if l.strip()})
    else:
        g2a = load_gene_to_acc(args.ref_fasta)
        accs = accs_from_missense(args.missense, g2a, args.vaf_threshold)
    print(f'{len(accs)} target accessions', flush=True)

    # shard for SLURM array
    shard = args.shard
    if shard is None and 'SLURM_ARRAY_TASK_ID' in os.environ:
        shard = int(os.environ['SLURM_ARRAY_TASK_ID']) - 1
    if shard is not None and args.nshards > 1:
        accs = accs[shard::args.nshards]
        print(f'shard {shard}/{args.nshards}: {len(accs)} accessions', flush=True)

    tally = {}
    for i, acc in enumerate(accs):
        try:
            r = make_one(acc, cache, out_dir)
        except Exception as e:
            r = 'error'
            print(f'  {acc}: ERROR {e}', file=sys.stderr, flush=True)
        tally[r] = tally.get(r, 0) + 1
        if (i + 1) % 25 == 0:
            print(f'  {i + 1}/{len(accs)}  {tally}', flush=True)
    print(f'DONE  {tally}', flush=True)


if __name__ == '__main__':
    main()
