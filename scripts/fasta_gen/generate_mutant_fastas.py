# -*- coding: utf-8 -*-
"""
generate_mutant_fastas.py

Generate per-sample mutant protein FASTA files from VEP missense output.

For each sample's missense mutations:
  - Look up the gene in the reference proteome (UniProt canonical)
  - Validate the reference amino acid at the mutation position
  - Apply the substitution and write a mutant FASTA entry

Usage:
    python3 generate_mutant_fastas.py \
        --ref-fasta /path/to/uniprot_human_canonical.fasta \
        --vep-dir /path/to/VEP/ \
        -o /path/to/FASTA/per_sample/

Optional:
    --min-vaf 0.1           Minimum variant allele frequency
    --max-gnomad-af 0.01    Maximum gnomAD population frequency
"""

import argparse
import csv
import os
import re
import sys

# Three-letter to one-letter amino acid mapping
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}

# HGVSp pattern: e.g. p.Arg273His or ENSP00000269305.4:p.Arg273His
HGVSP_PATTERN = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


def parse_reference_fasta(fasta_path):
    """
    Parse UniProt FASTA file.

    Returns:
        gene_to_protein: dict mapping gene symbol -> (accession, entry_name, sequence)
            If multiple proteins share a gene name, the first is used.
        accession_to_seq: dict mapping accession -> sequence
    """
    gene_to_protein = {}
    accession_to_seq = {}
    duplicates = []

    current_acc = None
    current_entry = None
    current_gene = None
    current_header = None
    lines = []

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # Save previous entry
                if current_acc is not None:
                    seq = "".join(lines)
                    accession_to_seq[current_acc] = seq
                    if current_gene and current_gene not in gene_to_protein:
                        gene_to_protein[current_gene] = (current_acc, current_entry, seq)
                    elif current_gene and current_gene in gene_to_protein:
                        duplicates.append(current_gene)

                # Parse new header
                # Format: >sp|P04637|P53_HUMAN Cellular tumor ... GN=TP53 PE=1 SV=4
                # or:     >tr|A0A0K0K1L2|A0A0K0K1L2_HUMAN ... GN=GENE ...
                current_header = line
                parts = line.split("|")
                if len(parts) >= 3:
                    current_acc = parts[1]
                    current_entry = parts[2].split()[0]
                else:
                    current_acc = line[1:].split()[0]
                    current_entry = current_acc

                # Extract gene name from GN= field
                gn_match = re.search(r"GN=(\S+)", line)
                current_gene = gn_match.group(1) if gn_match else None
                lines = []
            else:
                lines.append(line.strip())

    # Save last entry
    if current_acc is not None:
        seq = "".join(lines)
        accession_to_seq[current_acc] = seq
        if current_gene and current_gene not in gene_to_protein:
            gene_to_protein[current_gene] = (current_acc, current_entry, seq)

    if duplicates:
        n_dup = len(set(duplicates))
        print(f"  Note: {n_dup} gene symbols mapped to multiple proteins (first used)")

    return gene_to_protein, accession_to_seq


def parse_mutation(hgvsp, amino_acids, protein_position):
    """
    Parse a missense mutation from VEP output.

    Returns (ref_aa_1letter, position_1based, alt_aa_1letter) or None if unparseable.
    """
    # Try HGVSp first
    if hgvsp and hgvsp not in ("", "-", "NA"):
        m = HGVSP_PATTERN.search(hgvsp)
        if m:
            ref_3 = m.group(1)
            pos = int(m.group(2))
            alt_3 = m.group(3)
            ref_1 = AA3_TO_1.get(ref_3)
            alt_1 = AA3_TO_1.get(alt_3)
            if ref_1 and alt_1:
                return ref_1, pos, alt_1

    # Fallback: Amino_acids (e.g. "R/H") + Protein_position (e.g. "273")
    if amino_acids and "/" in amino_acids and protein_position:
        parts = amino_acids.split("/")
        if len(parts) == 2 and len(parts[0]) == 1 and len(parts[1]) == 1:
            try:
                pos = int(protein_position)
                return parts[0], pos, parts[1]
            except ValueError:
                pass

    return None


def process_sample(tsv_path, gene_to_protein, out_dir, min_vaf, max_gnomad_af):
    """
    Process one VEP TSV file and generate a mutant FASTA.

    Returns dict with counts: applied, skipped_gene, skipped_mismatch,
        skipped_vaf, skipped_gnomad, skipped_parse, skipped_dup
    """
    counts = {
        "applied": 0,
        "skipped_gene": 0,
        "skipped_mismatch": 0,
        "skipped_vaf": 0,
        "skipped_gnomad": 0,
        "skipped_parse": 0,
        "skipped_dup": 0,
    }
    log_entries = []  # (status, gene, detail)
    seen_mutations = set()  # (gene, ref, pos, alt) to deduplicate

    sample_id = os.path.basename(tsv_path).replace(".vep.tsv", "")
    out_path = os.path.join(out_dir, f"{sample_id}_mutant.fasta")

    entries = []

    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            symbol = row.get("SYMBOL", "").strip()
            hgvsp = row.get("HGVSp", "").strip()
            amino_acids = row.get("Amino_acids", "").strip()
            protein_position = row.get("Protein_position", "").strip()
            vaf_str = row.get("VAF", "NA").strip()
            gnomad_str = row.get("gnomADe_AF", "").strip()

            # Parse VAF filter
            if min_vaf is not None and vaf_str not in ("", "NA", "-"):
                try:
                    if float(vaf_str) < min_vaf:
                        counts["skipped_vaf"] += 1
                        continue
                except ValueError:
                    pass

            # Parse gnomAD filter
            if max_gnomad_af is not None and gnomad_str not in ("", "NA", "-"):
                try:
                    if float(gnomad_str) > max_gnomad_af:
                        counts["skipped_gnomad"] += 1
                        continue
                except ValueError:
                    pass

            # Parse mutation
            mutation = parse_mutation(hgvsp, amino_acids, protein_position)
            if mutation is None:
                counts["skipped_parse"] += 1
                log_entries.append(("PARSE_FAIL", symbol,
                                    f"HGVSp={hgvsp} AA={amino_acids} pos={protein_position}"))
                continue

            ref_aa, pos, alt_aa = mutation

            # Deduplicate within sample
            mut_key = (symbol, ref_aa, pos, alt_aa)
            if mut_key in seen_mutations:
                counts["skipped_dup"] += 1
                continue
            seen_mutations.add(mut_key)

            # Look up gene in reference
            if symbol not in gene_to_protein:
                counts["skipped_gene"] += 1
                log_entries.append(("GENE_NOT_FOUND", symbol, ""))
                continue

            accession, entry_name, ref_seq = gene_to_protein[symbol]

            # Validate position (1-based)
            if pos < 1 or pos > len(ref_seq):
                counts["skipped_mismatch"] += 1
                log_entries.append(("POS_OUT_OF_RANGE", symbol,
                                    f"pos={pos} len={len(ref_seq)}"))
                continue

            # Validate reference amino acid
            if ref_seq[pos - 1] != ref_aa:
                counts["skipped_mismatch"] += 1
                log_entries.append(("AA_MISMATCH", symbol,
                                    f"expected={ref_aa} found={ref_seq[pos - 1]} pos={pos}"))
                continue

            # Apply substitution
            mutant_seq = ref_seq[:pos - 1] + alt_aa + ref_seq[pos:]

            # Build header: >mut|P04637|TP53_R273H OS=Homo sapiens GN=TP53
            mut_label = f"{symbol}_{ref_aa}{pos}{alt_aa}"
            header = f">mut|{accession}|{mut_label} OS=Homo sapiens GN={symbol}"

            entries.append((header, mutant_seq))
            counts["applied"] += 1

    # Write FASTA (only if there are entries)
    if entries:
        with open(out_path, "w") as out:
            for header, seq in entries:
                out.write(header + "\n")
                # Write sequence in 60-char lines
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")

    return counts, log_entries, sample_id


def main():
    parser = argparse.ArgumentParser(
        description="Generate per-sample mutant protein FASTAs from VEP output."
    )
    parser.add_argument("--ref-fasta", required=True,
                        help="UniProt reference proteome FASTA")
    parser.add_argument("--vep-dir", required=True,
                        help="Directory containing *.vep.tsv files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for per-sample mutant FASTAs")
    parser.add_argument("--min-vaf", type=float, default=None,
                        help="Minimum variant allele frequency (default: no filter)")
    parser.add_argument("--max-gnomad-af", type=float, default=None,
                        help="Maximum gnomAD population frequency (default: no filter)")

    args = parser.parse_args()

    if not os.path.isfile(args.ref_fasta):
        sys.exit(f"Error: Reference FASTA not found: {args.ref_fasta}")
    if not os.path.isdir(args.vep_dir):
        sys.exit(f"Error: VEP directory not found: {args.vep_dir}")

    os.makedirs(args.output, exist_ok=True)

    # Load reference proteome
    print(f"Loading reference proteome: {args.ref_fasta}")
    gene_to_protein, _ = parse_reference_fasta(args.ref_fasta)
    print(f"  Loaded {len(gene_to_protein)} gene->protein mappings")

    # Find all VEP TSV files
    tsv_files = sorted([
        os.path.join(args.vep_dir, f)
        for f in os.listdir(args.vep_dir)
        if f.endswith(".vep.tsv")
    ])

    if not tsv_files:
        sys.exit(f"Error: No .vep.tsv files found in {args.vep_dir}")

    print(f"Processing {len(tsv_files)} samples...")

    # Summary log
    summary_path = os.path.join(args.output, "generation_summary.tsv")
    all_logs = []
    totals = {
        "applied": 0, "skipped_gene": 0, "skipped_mismatch": 0,
        "skipped_vaf": 0, "skipped_gnomad": 0, "skipped_parse": 0,
        "skipped_dup": 0,
    }

    with open(summary_path, "w") as summary_f:
        summary_f.write("sample_id\tapplied\tskipped_gene\tskipped_mismatch\t"
                        "skipped_vaf\tskipped_gnomad\tskipped_parse\tskipped_dup\n")

        for tsv in tsv_files:
            counts, log_entries, sample_id = process_sample(
                tsv, gene_to_protein, args.output,
                args.min_vaf, args.max_gnomad_af
            )

            summary_f.write(f"{sample_id}\t{counts['applied']}\t"
                            f"{counts['skipped_gene']}\t{counts['skipped_mismatch']}\t"
                            f"{counts['skipped_vaf']}\t{counts['skipped_gnomad']}\t"
                            f"{counts['skipped_parse']}\t{counts['skipped_dup']}\n")

            for key in totals:
                totals[key] += counts[key]

            all_logs.extend(log_entries)

    # Write detailed log of issues
    if all_logs:
        log_path = os.path.join(args.output, "generation_issues.tsv")
        with open(log_path, "w") as log_f:
            log_f.write("status\tgene\tdetail\n")
            for status, gene, detail in all_logs:
                log_f.write(f"{status}\t{gene}\t{detail}\n")
        print(f"Issue log: {log_path}")

    # Print summary
    print(f"\n{'='*50}")
    print(f"Mutant FASTA Generation Summary")
    print(f"{'='*50}")
    print(f"Samples processed:  {len(tsv_files)}")
    print(f"Mutations applied:  {totals['applied']}")
    print(f"Skipped (no gene):  {totals['skipped_gene']}")
    print(f"Skipped (mismatch): {totals['skipped_mismatch']}")
    print(f"Skipped (low VAF):  {totals['skipped_vaf']}")
    print(f"Skipped (gnomAD):   {totals['skipped_gnomad']}")
    print(f"Skipped (parse):    {totals['skipped_parse']}")
    print(f"Skipped (dup):      {totals['skipped_dup']}")
    print(f"Summary: {summary_path}")
    print(f"Output:  {args.output}")


if __name__ == "__main__":
    main()
