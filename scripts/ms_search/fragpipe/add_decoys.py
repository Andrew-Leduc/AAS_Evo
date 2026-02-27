# -*- coding: utf-8 -*-
"""
add_decoys.py

Copy per-plex FASTA files to an output directory with reversed decoy sequences
appended. FragPipe requires decoy sequences (prefix 'rev_') in the search
database for target-decoy FDR estimation.

By default writes to --output-dir (per_plex_fragpipe/), leaving the source
per_plex/ FASTAs untouched so MaxQuant (which handles decoys internally) can
use them without modification.

Usage:
    python3 add_decoys.py \
        --fasta-dir /path/to/FASTA/per_plex/ \
        --output-dir /path/to/FASTA/per_plex_fragpipe/

Options:
    --fasta-dir    Source directory with per-plex .fasta files (required)
    --output-dir   Destination directory for FASTAs with decoys (required)
    --prefix       Decoy header prefix (default: rev_)
    --force        Overwrite output files even if already present
"""

import argparse
import os
import shutil
import sys


def write_with_decoys(src_path, dst_path, prefix, force=False):
    """
    Copy src_path to dst_path and append reversed decoy entries.

    If dst_path already exists and force=False, skips the file.
    Returns True if written, False if skipped.
    """
    if not force and os.path.exists(dst_path):
        return False

    # Stream through source, writing targets then decoys to destination.
    with open(src_path, "r") as f_in, open(dst_path, "w") as f_out:
        entries = []
        header, seq_parts = None, []
        for line in f_in:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_parts)))
                header, seq_parts = line, []
            else:
                seq_parts.append(line)
        if header is not None:
            entries.append((header, "".join(seq_parts)))

        # Write all target entries first
        for header, seq in entries:
            f_out.write(f"{header}\n{seq}\n")

        # Append reversed decoy entries
        for header, seq in entries:
            f_out.write(f">{prefix}{header[1:]}\n{seq[::-1]}\n")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Copy per-plex FASTAs to a new directory with reversed decoys appended."
    )
    parser.add_argument("--fasta-dir", required=True,
                        help="Source directory containing per-plex .fasta files")
    parser.add_argument("--output-dir", required=True,
                        help="Destination directory for FASTAs with decoys")
    parser.add_argument("--prefix", default="rev_",
                        help="Decoy header prefix (default: rev_)")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite output files even if already present")

    args = parser.parse_args()

    if not os.path.isdir(args.fasta_dir):
        sys.exit(f"Error: FASTA directory not found: {args.fasta_dir}")

    os.makedirs(args.output_dir, exist_ok=True)

    fastas = sorted(f for f in os.listdir(args.fasta_dir) if f.endswith(".fasta"))
    if not fastas:
        sys.exit(f"Error: No .fasta files found in {args.fasta_dir}")

    print(f"Writing FASTAs with decoys (prefix='{args.prefix}')...")
    print(f"  Source:      {args.fasta_dir}")
    print(f"  Destination: {args.output_dir}")

    written, skipped = 0, 0
    for fname in fastas:
        src = os.path.join(args.fasta_dir, fname)
        dst = os.path.join(args.output_dir, fname)
        if write_with_decoys(src, dst, args.prefix, force=args.force):
            written += 1
        else:
            skipped += 1

    print(f"  Written: {written}")
    print(f"  Skipped (already exist, use --force to overwrite): {skipped}")


if __name__ == "__main__":
    main()
