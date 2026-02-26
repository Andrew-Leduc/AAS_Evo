# -*- coding: utf-8 -*-
"""
generate_mutant_fastas.py

Generate per-sample mutant tryptic peptide FASTA files from VEP missense output.

For each sample's missense mutations:
  - Look up the gene in the reference proteome (UniProt canonical)
  - Validate the reference amino acid at the mutation position
  - Extract the tryptic peptide(s) containing the mutation (with 1 missed cleavage)
  - Write mutant peptide FASTA entries

This outputs tryptic peptides only (not full proteins) to minimize FASTA size
and improve FDR statistics during MS database search.

Header format:
    >mut|{accession}|{gene}|{swap}|genetic

Example:
    >mut|P04637|TP53|R273H|genetic
    SVTCTYSPALNKMFCQLAK

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


def find_tryptic_cleavage_sites(sequence):
    """
    Find tryptic cleavage sites in a protein sequence.

    Trypsin cleaves after K or R, unless followed by P.

    Returns list of cleavage positions (0-indexed, position AFTER which cleavage occurs).
    Includes 0 (start) and len(sequence) (end) as boundaries.
    """
    sites = [0]  # Start of protein
    for i, aa in enumerate(sequence):
        # Cleave after K or R, unless followed by P
        if aa in ('K', 'R'):
            if i + 1 < len(sequence) and sequence[i + 1] == 'P':
                continue  # No cleavage before proline
            sites.append(i + 1)  # Position after this residue
    sites.append(len(sequence))  # End of protein
    return sites


def get_reference_peptides(sequence, max_missed_cleavages=0):
    """
    Get all tryptic peptides from a reference protein sequence.

    Used to filter mutant entries that would generate peptides identical
    to reference peptides (e.g. when a K/R cleavage site is mutated,
    the missed-cleavage peptide logic can still emit the pre-mutation
    fragment under certain conditions).

    Returns a set of peptide sequences.
    """
    sites = find_tryptic_cleavage_sites(sequence)
    peptides = set()
    for i in range(len(sites) - 1):
        for missed in range(max_missed_cleavages + 1):
            if i + missed + 1 >= len(sites):
                break
            start = sites[i]
            end = sites[i + missed + 1]
            pep = sequence[start:end]
            if 6 <= len(pep) <= 50:
                peptides.add(pep)
    return peptides


def extract_tryptic_peptides(sequence, mutation_pos, max_missed_cleavages=0):
    """
    Extract tryptic peptide(s) containing a mutation position.

    Args:
        sequence: Full protein sequence (with mutation already applied)
        mutation_pos: 1-based position of the mutation
        max_missed_cleavages: Include peptides with up to this many missed cleavages

    Returns:
        List of unique peptide sequences containing the mutation
    """
    sites = find_tryptic_cleavage_sites(sequence)
    mut_idx = mutation_pos - 1  # Convert to 0-based

    peptides = set()

    # Find all peptide windows that contain the mutation position
    for i in range(len(sites) - 1):
        for missed in range(max_missed_cleavages + 1):
            if i + missed + 1 >= len(sites):
                break

            start = sites[i]
            end = sites[i + missed + 1]

            # Check if this peptide contains the mutation
            if start <= mut_idx < end:
                peptide = sequence[start:end]
                # Filter out very short or very long peptides
                if 6 <= len(peptide) <= 50:
                    peptides.add(peptide)

    return list(peptides)


def parse_mutation(hgvsp, amino_acids, protein_position):
    """
    Parse a missense mutation from VEP output.

    Returns (ref_aa_1letter, position_1based, alt_aa_1letter, vep_protein_length)
    or None if unparseable. vep_protein_length may be None if not available.
    """
    vep_length = None

    # Parse VEP protein length from Protein_position (e.g., "273/393" -> 393)
    if protein_position and "/" in protein_position:
        try:
            vep_length = int(protein_position.split("/")[1])
        except (ValueError, IndexError):
            pass

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
                return ref_1, pos, alt_1, vep_length

    # Fallback: Amino_acids (e.g. "R/H") + Protein_position (e.g. "273" or "273/393")
    if amino_acids and "/" in amino_acids and protein_position:
        parts = amino_acids.split("/")
        if len(parts) == 2 and len(parts[0]) == 1 and len(parts[1]) == 1:
            try:
                # Handle both "273" and "273/393" formats
                pos_str = protein_position.split("/")[0] if "/" in protein_position else protein_position
                pos = int(pos_str)
                return parts[0], pos, parts[1], vep_length
            except ValueError:
                pass

    return None


def process_sample(tsv_path, gene_to_protein, out_dir, min_vaf, max_gnomad_af,
                   ref_pep_cache=None):
    """
    Process one VEP TSV file and generate a mutant tryptic peptide FASTA.

    ref_pep_cache: shared dict mapping accession -> set of reference tryptic
        peptide sequences. Lazily populated; pass the same dict across all
        samples to avoid redundant computation.

    Returns dict with counts: applied, skipped_gene, skipped_mismatch,
        skipped_vaf, skipped_gnomad, skipped_parse, skipped_dup
    """
    if ref_pep_cache is None:
        ref_pep_cache = {}

    counts = {
        "applied": 0,
        "peptides": 0,
        "skipped_gene": 0,
        "skipped_mismatch": 0,
        "skipped_length_mismatch": 0,
        "skipped_vaf": 0,
        "skipped_gnomad": 0,
        "skipped_parse": 0,
        "skipped_dup": 0,
        "skipped_ref_peptide": 0,
    }
    log_entries = []  # (status, gene, detail)
    seen_mutations = set()  # (gene, ref, pos, alt) to deduplicate
    seen_peptides = set()  # Deduplicate peptide sequences

    sample_id = os.path.basename(tsv_path).replace(".vep.tsv", "")
    out_path = os.path.join(out_dir, f"{sample_id}_mutant.fasta")

    if os.path.exists(out_path):
        return {"applied": 0, "peptides": 0, "skipped_gene": 0, "skipped_mismatch": 0,
                "skipped_length_mismatch": 0, "skipped_vaf": 0, "skipped_gnomad": 0,
                "skipped_parse": 0, "skipped_dup": 0,
                "skipped_ref_peptide": 0}, [], sample_id

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

            ref_aa, pos, alt_aa, vep_length = mutation

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
            uniprot_length = len(ref_seq)

            # Cache reference tryptic peptides for this protein (lazy)
            if accession not in ref_pep_cache:
                ref_pep_cache[accession] = get_reference_peptides(ref_seq)

            # Validate protein length matches (VEP uses Ensembl, we use UniProt)
            # Allow small differences (±5 AA) for signal peptide/initiator Met variations
            if vep_length is not None and abs(vep_length - uniprot_length) > 5:
                counts["skipped_length_mismatch"] += 1
                log_entries.append(("LENGTH_MISMATCH", symbol,
                                    f"VEP_len={vep_length} UniProt_len={uniprot_length} "
                                    f"diff={vep_length - uniprot_length} pos={pos}"))
                continue

            # Validate position (1-based)
            if pos < 1 or pos > uniprot_length:
                counts["skipped_mismatch"] += 1
                log_entries.append(("POS_OUT_OF_RANGE", symbol,
                                    f"pos={pos} len={uniprot_length}"))
                continue

            # Validate reference amino acid
            if ref_seq[pos - 1] != ref_aa:
                counts["skipped_mismatch"] += 1
                log_entries.append(("AA_MISMATCH", symbol,
                                    f"expected={ref_aa} found={ref_seq[pos - 1]} pos={pos}"))
                continue

            # Apply substitution to get mutant sequence
            mutant_seq = ref_seq[:pos - 1] + alt_aa + ref_seq[pos:]

            # Extract tryptic peptide(s) containing the mutation
            peptides = extract_tryptic_peptides(mutant_seq, pos, max_missed_cleavages=0)

            if not peptides:
                log_entries.append(("NO_PEPTIDE", symbol,
                                    f"pos={pos} - no valid tryptic peptide"))
                continue

            # Build header: >mut|P04637|TP53|R273H|genetic
            swap = f"{ref_aa}{pos}{alt_aa}"

            ref_peptides = ref_pep_cache[accession]
            for peptide in peptides:
                # Skip peptides identical to any reference tryptic peptide.
                # This can happen when a mutation removes a K/R cleavage site:
                # the missed-cleavage extraction may still emit a pre-mutation
                # fragment that is unchanged from the reference sequence.
                # Such peptides cannot distinguish mutant from reference in MS.
                if peptide in ref_peptides:
                    counts["skipped_ref_peptide"] += 1
                    continue

                # Deduplicate peptides across mutations
                pep_key = (accession, swap, peptide)
                if pep_key in seen_peptides:
                    continue
                seen_peptides.add(pep_key)

                header = f">mut|{accession}|{symbol}|{swap}|genetic"
                entries.append((header, peptide))
                counts["peptides"] += 1

            counts["applied"] += 1

    # Write FASTA atomically via temp file to avoid partial writes on job kill
    if entries:
        tmp_path = out_path + ".tmp"
        with open(tmp_path, "w") as out:
            for header, seq in entries:
                out.write(header + "\n")
                out.write(seq + "\n")
        os.rename(tmp_path, out_path)

    return counts, log_entries, sample_id


def main():
    parser = argparse.ArgumentParser(
        description="Generate per-sample mutant tryptic peptide FASTAs from VEP output."
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

    # Find all VEP TSV files (walks chunk subdirectories)
    tsv_files = sorted([
        os.path.join(root, f)
        for root, _, files in os.walk(args.vep_dir)
        for f in files
        if f.endswith(".vep.tsv")
    ])

    if not tsv_files:
        sys.exit(f"Error: No .vep.tsv files found in {args.vep_dir}")

    print(f"Processing {len(tsv_files)} samples...")

    # Summary log
    summary_path = os.path.join(args.output, "generation_summary.tsv")
    all_logs = []
    totals = {
        "applied": 0, "peptides": 0, "skipped_gene": 0, "skipped_mismatch": 0,
        "skipped_length_mismatch": 0, "skipped_vaf": 0, "skipped_gnomad": 0,
        "skipped_parse": 0, "skipped_dup": 0, "skipped_ref_peptide": 0,
    }

    # Shared cache of reference tryptic peptides per protein accession.
    # Avoids recomputing the same protein's digest for every sample.
    ref_pep_cache = {}

    with open(summary_path, "w") as summary_f:
        summary_f.write("sample_id\tmutations_applied\tpeptides_written\tskipped_gene\t"
                        "skipped_mismatch\tskipped_length_mismatch\tskipped_vaf\t"
                        "skipped_gnomad\tskipped_parse\tskipped_dup\tskipped_ref_peptide\n")

        for tsv in tsv_files:
            counts, log_entries, sample_id = process_sample(
                tsv, gene_to_protein, args.output,
                args.min_vaf, args.max_gnomad_af,
                ref_pep_cache=ref_pep_cache,
            )

            summary_f.write(f"{sample_id}\t{counts['applied']}\t{counts['peptides']}\t"
                            f"{counts['skipped_gene']}\t{counts['skipped_mismatch']}\t"
                            f"{counts['skipped_length_mismatch']}\t{counts['skipped_vaf']}\t"
                            f"{counts['skipped_gnomad']}\t{counts['skipped_parse']}\t"
                            f"{counts['skipped_dup']}\t{counts['skipped_ref_peptide']}\n")
            summary_f.flush()

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
    print(f"Mutant Tryptic Peptide FASTA Generation Summary")
    print(f"{'='*50}")
    print(f"Samples processed:    {len(tsv_files)}")
    print(f"Mutations applied:    {totals['applied']}")
    print(f"Peptides written:     {totals['peptides']}")
    print(f"Skipped (no gene):    {totals['skipped_gene']}")
    print(f"Skipped (AA mismatch):{totals['skipped_mismatch']}")
    print(f"Skipped (len mismatch):{totals['skipped_length_mismatch']}")
    print(f"Skipped (low VAF):    {totals['skipped_vaf']}")
    print(f"Skipped (gnomAD):     {totals['skipped_gnomad']}")
    print(f"Skipped (parse):      {totals['skipped_parse']}")
    print(f"Skipped (dup):        {totals['skipped_dup']}")
    print(f"Skipped (ref peptide):{totals['skipped_ref_peptide']}")
    print(f"Summary: {summary_path}")
    print(f"Output:  {args.output}")


if __name__ == "__main__":
    main()
