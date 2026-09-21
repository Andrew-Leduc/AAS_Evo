#!/usr/bin/env python3
"""
Map each SAAP contact site (UniProt acc, protein position) -> reference codon and
genomic coordinates, via Ensembl canonical CDS (pyensembl).

Strategy that avoids UniProt<->Ensembl id headaches: for each gene, pick the
Ensembl protein-coding transcript whose translated protein SEQUENCE matches our
UniProt canonical sequence (exact, else best prefix match). That guarantees the
protein-position numbering is identical, so codon = coding_sequence[3*(p-1):3*p].
We verify translate(codon) == wt and drop mismatches.

Genomic coordinates (for phyloP #5 / gnomAD #6) are a best-effort walk of the
transcript exons; if pyensembl can't give the CDS start offset the codon columns
are still valid and coords are left blank.

Setup (compute node has proxy internet):
    module load anaconda3/2024.06
    pip install --user pyensembl
    export https_proxy=http://10.99.0.130:3128 http_proxy=http://10.99.0.130:3128
    pyensembl install --release 111 --species homo_sapiens

Run:
    python3 build_codon_map.py            # sites from precursor_pairs.tsv
Output: contact_saap/codon_map.tsv
"""
import argparse
import os
import re
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"
REF_FASTA = "/scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta"

CODON_TABLE = {}  # filled below
_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for _i, _a in enumerate(_AAS):
    CODON_TABLE[_BASES[_i >> 4] + _BASES[(_i >> 2) & 3] + _BASES[_i & 3]] = _a


def load_ref(path):
    seqs, gene = {}, {}
    cur = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                parts = line.split("|")
                cur = parts[1] if len(parts) >= 3 else line[1:].split()[0]
                m = re.search(r"GN=(\S+)", line)
                gene[cur] = m.group(1) if m else None
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return {a: "".join(s) for a, s in seqs.items()}, gene


def coding_positions(tx):
    """Genomic position of each mRNA base in 5'->3' translation order (incl UTR)."""
    exons = sorted(tx.exons, key=lambda e: e.start)
    if tx.strand == "-":
        exons = exons[::-1]
    pos = []
    for e in exons:
        rng = range(e.start, e.end + 1) if tx.strand == "+" else range(e.end, e.start - 1, -1)
        pos.extend(rng)
    return pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=f"{B}/precursor_pairs.tsv")
    ap.add_argument("--release", type=int, default=111)
    ap.add_argument("--out", default=f"{B}/codon_map.tsv")
    a = ap.parse_args()

    from pyensembl import EnsemblRelease
    data = EnsemblRelease(a.release, species="homo_sapiens")

    seqs, acc2gene = load_ref(REF_FASTA)
    pairs = pd.read_csv(a.pairs, sep="\t", usecols=["acc", "pos", "wt"]).drop_duplicates()
    sites = pairs.groupby("acc")["pos"].agg(lambda s: sorted(set(s))).to_dict()
    wt_of = {(r.acc, int(r.pos)): r.wt for r in pairs.itertuples()}
    print(f"{len(pairs):,} sites across {len(sites):,} proteins")

    tx_cache = {}

    def pick_transcript(acc):
        if acc in tx_cache:
            return tx_cache[acc]
        gene = acc2gene.get(acc)
        up = seqs.get(acc, "")
        chosen = None
        if gene and up:
            try:
                cands = [t for g in data.genes_by_name(gene) for t in g.transcripts
                         if t.is_protein_coding and t.complete]
            except Exception:
                cands = []
            best = None
            for t in cands:
                try:
                    prot = t.protein_sequence or ""
                except Exception:
                    continue
                if prot == up:
                    chosen = t; break
                # best partial: longest common prefix length
                k = 0
                for x, y in zip(prot, up):
                    if x != y:
                        break
                    k += 1
                if best is None or k > best[0]:
                    best = (k, t)
            if chosen is None and best and best[0] >= 0.9 * len(up):
                chosen = best[1]
        tx_cache[acc] = chosen
        return chosen

    rows = []
    n_ok = n_notx = n_mismatch = 0
    for acc, positions in sites.items():
        tx = pick_transcript(acc)
        if tx is None:
            n_notx += 1
            continue
        try:
            cds = tx.coding_sequence or ""
        except Exception:
            cds = ""
        gpos = None
        try:
            off = tx.first_start_codon_spliced_offset
            allpos = coding_positions(tx)
            gpos = allpos[off:]           # coding bases in translation order
        except Exception:
            gpos = None
        for p in positions:
            wt = wt_of[(acc, p)]
            codon = cds[3 * (p - 1): 3 * p]
            if len(codon) != 3 or CODON_TABLE.get(codon) != wt:
                n_mismatch += 1
                continue
            row = dict(acc=acc, gene=acc2gene.get(acc), pos=p, wt=wt,
                       transcript_id=tx.transcript_id, codon=codon,
                       chrom=tx.contig, strand=tx.strand,
                       g1="", g2="", g3="")
            if gpos is not None and 3 * p <= len(gpos):
                g = gpos[3 * (p - 1): 3 * p]
                row["g1"], row["g2"], row["g3"] = g[0], g[1], g[2]
            rows.append(row)
            n_ok += 1

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    print(f"mapped codons: {n_ok:,} | proteins w/o transcript: {n_notx:,} | "
          f"aa mismatches dropped: {n_mismatch:,}")
    print(f"  with genomic coords: {int((out['g1']!='').sum()):,}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
