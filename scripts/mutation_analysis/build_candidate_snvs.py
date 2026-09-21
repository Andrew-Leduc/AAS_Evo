#!/usr/bin/env python3
"""
#6 (part 1)  Enumerate the single-nucleotide variants that could produce each
detected SAAP swap, as a VCF to run through VEP for gnomAD allele frequencies.

For each unique SAAP (acc,pos,wt,alt) whose site is in codon_map.tsv, look at the
reference codon and emit every single-nt change (at any of the 3 codon positions)
that yields a codon translating to `alt`. Genomic ref/alt are strand-corrected
from the codon (codon is the mRNA codon, verified to translate to wt).

VCF uses chr-prefixed coordinates to match the UCSC hg38 FASTA used by VEP.
The ID field encodes the swap:  acc|pos|wt|alt  (so max-AF-per-swap is a groupby).

Then run VEP (see header of analysis_raas_gnomad.py) and parse the output there.
"""
import argparse
import os
import pandas as pd

B = "/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap"

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
for _i, _a in enumerate(_AAS):
    CODON_TABLE[_BASES[_i >> 4] + _BASES[(_i >> 2) & 3] + _BASES[_i & 3]] = _a
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def ucsc_chrom(c):
    c = str(c)
    return "chrM" if c == "MT" else (c if c.startswith("chr") else "chr" + c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=f"{B}/precursor_pairs.tsv")
    ap.add_argument("--codons", default=f"{B}/codon_map.tsv")
    ap.add_argument("--out", default=f"{B}/candidate_snvs.vcf")
    a = ap.parse_args()

    cm = pd.read_csv(a.codons, sep="\t").dropna(subset=["g1", "g2", "g3", "codon"])
    for c in ("g1", "g2", "g3"):
        cm[c] = cm[c].astype(int)
    cm = cm[cm["codon"].str.len() == 3]
    site = cm.set_index(["acc", "pos", "wt"])

    swaps = (pd.read_csv(a.pairs, sep="\t", usecols=["acc", "pos", "wt", "alt"])
             .drop_duplicates())
    print(f"unique swaps: {len(swaps):,} | codon sites: {len(cm):,}")

    recs = set()
    n_swaps_with_snv = 0
    for r in swaps.itertuples(index=False):
        key = (r.acc, r.pos, r.wt)
        if key not in site.index:
            continue
        row = site.loc[key]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        codon, strand = row["codon"], row["strand"]
        gpos = [row["g1"], row["g2"], row["g3"]]
        hit = False
        for i in range(3):
            for b in "ACGT":
                if b == codon[i]:
                    continue
                new = codon[:i] + b + codon[i + 1:]
                if CODON_TABLE.get(new) != r.alt:
                    continue
                ref_g = codon[i] if strand == "+" else COMP[codon[i]]
                alt_g = b if strand == "+" else COMP[b]
                recs.add((ucsc_chrom(row["chrom"]), int(gpos[i]), ref_g, alt_g,
                          f"{r.acc}|{r.pos}|{r.wt}|{r.alt}"))
                hit = True
        n_swaps_with_snv += hit

    def chrom_key(c):
        s = c[3:]
        return (0, int(s)) if s.isdigit() else (1, s)
    rows = sorted(recs, key=lambda x: (chrom_key(x[0]), x[1]))

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt, sid in rows:
            f.write(f"{chrom}\t{pos}\t{sid}\t{ref}\t{alt}\t.\t.\t.\n")
    print(f"wrote {a.out}: {len(rows):,} SNVs for "
          f"{n_swaps_with_snv:,} swaps reachable by a single nt change "
          f"({100*n_swaps_with_snv/len(swaps):.0f}% of swaps)")


if __name__ == "__main__":
    main()
