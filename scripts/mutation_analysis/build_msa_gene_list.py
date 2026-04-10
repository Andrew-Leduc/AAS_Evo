#!/usr/bin/env python3
"""
Map MSA_have_list.txt (UniProt entry names like 1433G_HUMAN_b03.a2m)
to accessions and gene symbols via uniprot_human_canonical.fasta.

Outputs:
  - metadata/MSA_have_list_mapped.tsv  (entry_name, accession, gene_symbol)
  - ANALYSIS/gene_list_for_spurs.txt   (gene symbols, one per line)
"""

import argparse
import re
from pathlib import Path


def parse_uniprot_fasta(fasta_path):
    """Returns dicts: entry_name -> accession, entry_name -> gene_symbol"""
    entry2acc = {}
    entry2gene = {}
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            # >sp|P61981|1433G_HUMAN Gene=YWHAG ...
            parts = line.split("|")
            if len(parts) < 3:
                continue
            acc = parts[1].strip()
            rest = parts[2]
            entry_name = rest.split()[0].strip()  # e.g. 1433G_HUMAN
            m = re.search(r"GN=(\S+)", line)
            gene = m.group(1) if m else None
            entry2acc[entry_name] = acc
            entry2gene[entry_name] = gene
    return entry2acc, entry2gene


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--msa-list", required=True, help="metadata/MSA_have_list.txt")
    ap.add_argument("--ref-fasta", required=True, help="uniprot_human_canonical.fasta")
    ap.add_argument("--mapped-out", required=True, help="Output TSV with accession+gene columns")
    ap.add_argument("--gene-list-out", required=True, help="Output gene list for SPURS")
    args = ap.parse_args()

    entry2acc, entry2gene = parse_uniprot_fasta(args.ref_fasta)

    mapped_out = Path(args.mapped_out)
    gene_list_out = Path(args.gene_list_out)
    mapped_out.parent.mkdir(parents=True, exist_ok=True)
    gene_list_out.parent.mkdir(parents=True, exist_ok=True)

    missing = 0
    written = 0
    genes_written = set()

    with open(args.msa_list) as fin, \
         open(mapped_out, "w") as fmap, \
         open(gene_list_out, "w") as fgene:

        fmap.write("entry_name\taccession\tgene_symbol\tmsa_file\n")

        for line in fin:
            msa_file = line.strip()
            if not msa_file:
                continue

            # Strip batch suffix and extension: 1433G_HUMAN_b03.a2m -> 1433G_HUMAN
            entry_name = re.sub(r"_b\d+\.a2m$", "", msa_file)

            acc = entry2acc.get(entry_name)
            gene = entry2gene.get(entry_name)

            if not acc:
                missing += 1
                fmap.write(f"{entry_name}\tNA\tNA\t{msa_file}\n")
                continue

            fmap.write(f"{entry_name}\t{acc}\t{gene or 'NA'}\t{msa_file}\n")
            written += 1

            if gene and gene not in genes_written:
                fgene.write(f"{gene}\n")
                genes_written.add(gene)

    print(f"Mapped:   {written}")
    print(f"Missing:  {missing}")
    print(f"Genes:    {len(genes_written)}")
    print(f"Mapped TSV: {mapped_out}")
    print(f"Gene list:  {gene_list_out}")


if __name__ == "__main__":
    main()
