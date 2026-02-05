# -*- coding: utf-8 -*-
"""
generate_msas.py

Generate multiple sequence alignments via MMseqs2 for coevolution analysis.

Two modes:

  1) --make-gene-list: Identify genes that need MSAs, write list to file.

     python3 generate_msas.py --make-gene-list \
         --vep-tsv /path/to/all_missense_mutations.tsv \
         --ref-fasta /path/to/uniprot_human_canonical.fasta \
         --msa-dir /path/to/MSA/ \
         -o gene_list.txt

  2) Single-gene MSA generation (called by SLURM array job):

     python3 generate_msas.py \
         --gene-list gene_list.txt \
         --index 1 \
         --ref-fasta /path/to/uniprot_human_canonical.fasta \
         --target-db /path/to/uniref30_2302 \
         --msa-dir /path/to/MSA/ \
         --tmp-dir /path/to/tmp/ \
         --threads 4
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile


# Extensions recognized as existing MSA files
MSA_EXTENSIONS = {".a3m", ".sto", ".stockholm", ".fasta", ".fa", ".aln"}


# ---------------------------------------------------------------------------
# Reference FASTA parsing (duplicated from coevolution_analysis.py)
# ---------------------------------------------------------------------------

def parse_reference_fasta(fasta_path):
    """Parse UniProt FASTA. Returns gene_to_protein and accession/entry maps."""
    gene_to_protein = {}
    accession_to_gene = {}
    entry_to_gene = {}

    current_acc = None
    current_entry = None
    current_gene = None
    lines = []

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_acc is not None:
                    seq = "".join(lines)
                    if current_gene and current_gene not in gene_to_protein:
                        gene_to_protein[current_gene] = (
                            current_acc, current_entry, seq
                        )
                    if current_gene:
                        accession_to_gene[current_acc] = current_gene
                        entry_to_gene[current_entry] = current_gene

                parts = line.split("|")
                if len(parts) >= 3:
                    current_acc = parts[1]
                    current_entry = parts[2].split()[0]
                else:
                    current_acc = line[1:].split()[0]
                    current_entry = current_acc

                gn_match = re.search(r"GN=(\S+)", line)
                current_gene = gn_match.group(1) if gn_match else None
                lines = []
            else:
                lines.append(line.strip())

    if current_acc is not None:
        seq = "".join(lines)
        if current_gene and current_gene not in gene_to_protein:
            gene_to_protein[current_gene] = (current_acc, current_entry, seq)
        if current_gene:
            accession_to_gene[current_acc] = current_gene
            entry_to_gene[current_entry] = current_gene

    return gene_to_protein, accession_to_gene, entry_to_gene


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_genes_with_mutations(vep_tsv_path):
    """Read VEP TSV and return set of unique gene symbols."""
    genes = set()
    with open(vep_tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            symbol = row.get("SYMBOL", "").strip()
            if symbol:
                genes.add(symbol)
    return genes


def get_existing_msas(msa_dir):
    """Scan MSA directory and return set of filename stems with MSA extensions."""
    stems = set()
    if not os.path.isdir(msa_dir):
        return stems
    for fname in os.listdir(msa_dir):
        stem, ext = os.path.splitext(fname)
        if ext.lower() in MSA_EXTENSIONS:
            stems.add(stem)
    return stems


def get_genes_with_existing_msas(msa_dir, accession_to_gene, entry_to_gene):
    """
    Return set of gene symbols that already have MSA files.

    MSA files may be named by gene symbol, UniProt accession, or entry name.
    """
    genes = set()
    stems = get_existing_msas(msa_dir)
    all_gene_symbols = set(accession_to_gene.values()) | set(entry_to_gene.values())

    for stem in stems:
        if stem in all_gene_symbols:
            genes.add(stem)
        elif stem in accession_to_gene:
            genes.add(accession_to_gene[stem])
        elif stem in entry_to_gene:
            genes.add(entry_to_gene[stem])
    return genes


def write_query_fasta(accession, entry_name, gene, sequence, output_path):
    """Write a single-protein FASTA file for MMseqs2 query."""
    with open(output_path, "w") as f:
        f.write(f">sp|{accession}|{entry_name} GN={gene}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i + 60] + "\n")


def count_msa_sequences(a3m_path):
    """Count number of sequences in an A3M/FASTA file."""
    count = 0
    with open(a3m_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def run_command(cmd, description):
    """Run a shell command and check for errors."""
    print(f"    {description}")
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"    ERROR: {description} failed (exit {result.returncode})")
        if result.stderr:
            for line in result.stderr.strip().split("\n")[:20]:
                print(f"      {line}")
        return False
    return True


# ---------------------------------------------------------------------------
# Mode 1: Generate gene list
# ---------------------------------------------------------------------------

def make_gene_list(args):
    """Identify genes needing MSAs and write them to a file."""
    print(f"Loading mutations: {args.vep_tsv}")
    mutated_genes = get_genes_with_mutations(args.vep_tsv)
    print(f"  {len(mutated_genes)} genes with mutations")

    print(f"Loading reference proteome: {args.ref_fasta}")
    gene_to_protein, accession_to_gene, entry_to_gene = parse_reference_fasta(
        args.ref_fasta
    )
    print(f"  {len(gene_to_protein)} genes in reference")

    # Keep only genes that have both mutations AND a UniProt entry
    valid_genes = mutated_genes & set(gene_to_protein.keys())
    no_uniprot = mutated_genes - set(gene_to_protein.keys())
    if no_uniprot:
        print(f"  {len(no_uniprot)} mutated genes not in reference (skipped)")

    # Check for existing MSAs (handles files named by gene, accession, or entry)
    genes_with_msas = get_genes_with_existing_msas(
        args.msa_dir, accession_to_gene, entry_to_gene
    )
    already_done = valid_genes & genes_with_msas
    need_msa = sorted(valid_genes - genes_with_msas)

    print(f"\nSummary:")
    print(f"  Genes with mutations + UniProt entry: {len(valid_genes)}")
    print(f"  Already have MSA:                     {len(already_done)}")
    print(f"  Need MSA generation:                  {len(need_msa)}")

    # Write output
    out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)

    with open(args.output, "w") as f:
        for gene in need_msa:
            f.write(gene + "\n")

    print(f"\nGene list written: {args.output} ({len(need_msa)} genes)")


# ---------------------------------------------------------------------------
# Mode 2: Single-gene MSA generation
# ---------------------------------------------------------------------------

def run_mmseqs_pipeline(query_fasta, target_db, output_a3m, tmp_dir,
                        threads, num_iterations, sensitivity):
    """Run the MMseqs2 search + MSA pipeline for a single query protein."""
    query_db = os.path.join(tmp_dir, "queryDB")
    result_db = os.path.join(tmp_dir, "resultDB")
    msa_db = os.path.join(tmp_dir, "msaDB")
    msa_out = os.path.join(tmp_dir, "msa_out")
    mmseqs_tmp = os.path.join(tmp_dir, "mmseqs_tmp")

    os.makedirs(msa_out, exist_ok=True)
    os.makedirs(mmseqs_tmp, exist_ok=True)

    # Step 1: Create query database
    if not run_command(
        f"mmseqs createdb {query_fasta} {query_db}",
        "Creating query database"
    ):
        return False

    # Step 2: Search against target database
    if not run_command(
        f"mmseqs search {query_db} {target_db} {result_db} {mmseqs_tmp} "
        f"--num-iterations {num_iterations} -s {sensitivity} "
        f"--threads {threads}",
        f"Searching target DB (iters={num_iterations}, sens={sensitivity})"
    ):
        return False

    # Step 3: Convert results to MSA (A3M format)
    if not run_command(
        f"mmseqs result2msa {query_db} {target_db} {result_db} {msa_db} "
        f"--msa-format-mode 6",
        "Converting results to A3M"
    ):
        return False

    # Step 4: Unpack MSA database to individual files
    if not run_command(
        f"mmseqs unpackdb {msa_db} {msa_out} --unpack-name-mode 0",
        "Unpacking MSA database"
    ):
        return False

    # Find the output file (single query = one file, usually named "0")
    out_files = [f for f in os.listdir(msa_out) if not f.endswith(".index")]
    if not out_files:
        print("    ERROR: No MSA output produced")
        return False

    # Copy to final destination
    src = os.path.join(msa_out, out_files[0])
    shutil.copy2(src, output_a3m)
    return True


def generate_single_msa(args):
    """Generate MSA for a single gene (called by SLURM array task)."""
    # Check mmseqs is available
    if shutil.which("mmseqs") is None:
        sys.exit("Error: mmseqs not found on PATH. "
                 "Load the module: module load mmseqs2")

    # Read gene from list
    with open(args.gene_list) as f:
        lines = f.read().strip().split("\n")

    if args.index < 1 or args.index > len(lines):
        sys.exit(f"Error: --index {args.index} out of range "
                 f"(gene list has {len(lines)} entries)")

    gene = lines[args.index - 1].strip()
    print(f"Gene: {gene} (index {args.index}/{len(lines)})")

    # Look up protein sequence to get accession for output naming
    print(f"Loading reference proteome: {args.ref_fasta}")
    gene_to_protein, accession_to_gene, entry_to_gene = parse_reference_fasta(
        args.ref_fasta
    )

    if gene not in gene_to_protein:
        print(f"  WARNING: {gene} not found in reference FASTA. Skipping.")
        return

    accession, entry_name, sequence = gene_to_protein[gene]
    print(f"  {accession} ({entry_name}), {len(sequence)} aa")

    # Skip-if-done: check for accession-named or gene-named MSA
    output_a3m = os.path.join(args.msa_dir, f"{accession}.a3m")
    gene_a3m = os.path.join(args.msa_dir, f"{gene}.a3m")
    for existing_path in (output_a3m, gene_a3m):
        if os.path.isfile(existing_path):
            n_seqs = count_msa_sequences(existing_path)
            print(f"  Already exists: {existing_path} ({n_seqs} sequences)")
            print("  Skipping.")
            return

    # Create temp directory
    os.makedirs(args.tmp_dir, exist_ok=True)
    tmp_dir = tempfile.mkdtemp(prefix=f"msa_{gene}_", dir=args.tmp_dir)

    try:
        # Write query FASTA
        query_fasta = os.path.join(tmp_dir, "query.fasta")
        write_query_fasta(accession, entry_name, gene, sequence, query_fasta)

        # Run MMseqs2
        os.makedirs(args.msa_dir, exist_ok=True)
        success = run_mmseqs_pipeline(
            query_fasta, args.target_db, output_a3m, tmp_dir,
            args.threads, args.num_iterations, args.sensitivity
        )

        if not success:
            # Clean up partial output
            if os.path.isfile(output_a3m):
                os.remove(output_a3m)
            sys.exit(f"Error: MMseqs2 pipeline failed for {gene}")

        # Report result
        n_seqs = count_msa_sequences(output_a3m)
        print(f"\n  Output: {output_a3m}")
        print(f"  Sequences: {n_seqs}")

    finally:
        # Clean up temp directory
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate MSAs via MMseqs2 for coevolution analysis."
    )

    # Mode selection
    parser.add_argument("--make-gene-list", action="store_true",
                        help="Mode 1: generate gene list file instead of MSA")

    # Shared arguments
    parser.add_argument("--ref-fasta", required=True,
                        help="UniProt reference proteome FASTA")
    parser.add_argument("--msa-dir", required=True,
                        help="Directory for MSA files (input for checking, "
                             "output for generation)")

    # Mode 1 arguments
    parser.add_argument("--vep-tsv",
                        help="Consolidated missense mutations TSV "
                             "(required for --make-gene-list)")
    parser.add_argument("-o", "--output",
                        help="Output gene list file "
                             "(required for --make-gene-list)")

    # Mode 2 arguments
    parser.add_argument("--gene-list",
                        help="Gene list file (one gene per line)")
    parser.add_argument("--index", type=int,
                        help="1-based index into gene list "
                             "(from SLURM_ARRAY_TASK_ID)")
    parser.add_argument("--target-db",
                        help="Path to MMseqs2 target database prefix "
                             "(e.g. UniRef30 or UniRef90)")
    parser.add_argument("--uniref90-db", dest="target_db",
                        help="(deprecated, use --target-db) "
                             "Path to MMseqs2 database prefix")
    parser.add_argument("--tmp-dir", default="/tmp",
                        help="Temporary directory for MMseqs2 working files")
    parser.add_argument("--threads", type=int, default=4,
                        help="CPU threads for MMseqs2 (default: 4)")
    parser.add_argument("--num-iterations", type=int, default=3,
                        help="MMseqs2 search iterations (default: 3)")
    parser.add_argument("--sensitivity", type=float, default=7.5,
                        help="MMseqs2 sensitivity (default: 7.5)")

    args = parser.parse_args()

    if args.make_gene_list:
        # Mode 1: generate gene list
        if not args.vep_tsv:
            sys.exit("Error: --vep-tsv required for --make-gene-list mode")
        if not args.output:
            sys.exit("Error: -o/--output required for --make-gene-list mode")
        if not os.path.isfile(args.vep_tsv):
            sys.exit(f"Error: VEP TSV not found: {args.vep_tsv}")
        if not os.path.isfile(args.ref_fasta):
            sys.exit(f"Error: Reference FASTA not found: {args.ref_fasta}")

        make_gene_list(args)
    else:
        # Mode 2: generate single MSA
        if not args.gene_list:
            sys.exit("Error: --gene-list required for MSA generation mode")
        if args.index is None:
            sys.exit("Error: --index required for MSA generation mode")
        if not args.target_db:
            sys.exit("Error: --target-db required for MSA generation mode")
        if not os.path.isfile(args.gene_list):
            sys.exit(f"Error: Gene list not found: {args.gene_list}")
        if not os.path.isfile(args.ref_fasta):
            sys.exit(f"Error: Reference FASTA not found: {args.ref_fasta}")

        generate_single_msa(args)


if __name__ == "__main__":
    main()
