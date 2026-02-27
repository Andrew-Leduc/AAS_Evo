# -*- coding: utf-8 -*-
"""
count_plex_peptides.py

Count unique peptides in a per-plex FASTA, broken down by entry type:
  - ref:  standard reference proteome entries (sp|accession|...)
  - mut:  mutant tryptic peptides (sp|accession-SWAP-HASH|gene-mut ...)
  - comp: compensatory peptides (sp|accession-comp-...|gene-comp ...)

Usage:
    python3 count_plex_peptides.py /path/to/FASTA/per_plex/PLEX_ID.fasta
    python3 count_plex_peptides.py /path/to/FASTA/per_plex/  # first FASTA found
"""

import os
import sys


def classify_header(header):
    """Return 'mut', 'comp', or 'ref' based on the FASTA header."""
    # Mutant: third pipe field starts with gene-mut
    # Comp:   third pipe field starts with gene-comp
    parts = header.split("|")
    if len(parts) >= 3:
        desc = parts[2]
        if "-mut " in desc or desc.strip().endswith("-mut"):
            return "mut"
        if "-comp " in desc or desc.strip().endswith("-comp"):
            return "comp"
    return "ref"


def count_fasta(path):
    ref_seqs  = set()
    mut_seqs  = set()
    comp_seqs = set()

    mut_headers  = set()   # unique (accession, swap) combos
    comp_headers = set()

    header, seq_parts = None, []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts)
                    kind = classify_header(header)
                    if kind == "mut":
                        mut_seqs.add(seq)
                        # accession field: P04637-R273H-A3F2 → (P04637, R273H)
                        acc_field = header.split("|")[1] if "|" in header else ""
                        parts = acc_field.split("-")
                        if len(parts) >= 2:
                            mut_headers.add((parts[0], parts[1]))
                    elif kind == "comp":
                        comp_seqs.add(seq)
                        acc_field = header.split("|")[1] if "|" in header else ""
                        comp_headers.add(acc_field)
                    else:
                        ref_seqs.add(seq)
                header, seq_parts = line, []
            else:
                seq_parts.append(line)

        # flush last entry
        if header is not None:
            seq = "".join(seq_parts)
            kind = classify_header(header)
            if kind == "mut":
                mut_seqs.add(seq)
                acc_field = header.split("|")[1] if "|" in header else ""
                parts = acc_field.split("-")
                if len(parts) >= 2:
                    mut_headers.add((parts[0], parts[1]))
            elif kind == "comp":
                comp_seqs.add(seq)
            else:
                ref_seqs.add(seq)

    return {
        "ref_proteins":       len(ref_seqs),
        "mut_unique_seqs":    len(mut_seqs),
        "mut_unique_mutations": len(mut_headers),
        "comp_unique_seqs":   len(comp_seqs),
        "comp_unique_entries": len(comp_headers),
        "total_entries":      len(ref_seqs) + len(mut_seqs) + len(comp_seqs),
    }


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]

    if os.path.isdir(path):
        fastas = sorted(f for f in os.listdir(path) if f.endswith(".fasta"))
        if not fastas:
            sys.exit(f"No .fasta files found in {path}")
        path = os.path.join(path, fastas[0])
        print(f"Using first FASTA: {fastas[0]}\n")

    if not os.path.isfile(path):
        sys.exit(f"File not found: {path}")

    size_mb = os.path.getsize(path) / 1e6
    print(f"FASTA: {os.path.basename(path)}  ({size_mb:.1f} MB)")
    print(f"{'─' * 45}")

    stats = count_fasta(path)

    print(f"Reference proteins:          {stats['ref_proteins']:>8,}")
    print(f"Mutant unique peptide seqs:  {stats['mut_unique_seqs']:>8,}")
    print(f"Mutant unique (acc, swap):   {stats['mut_unique_mutations']:>8,}")
    print(f"Compensatory unique seqs:    {stats['comp_unique_seqs']:>8,}")
    print(f"Compensatory unique entries: {stats['comp_unique_entries']:>8,}")
    print(f"{'─' * 45}")
    print(f"Total unique sequences:      {stats['total_entries']:>8,}")


if __name__ == "__main__":
    main()
