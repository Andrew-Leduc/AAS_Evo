# -*- coding: utf-8 -*-
"""
add_decoys.py

Append reversed decoy sequences to per-plex FASTA files.

FragPipe requires decoy sequences (prefix 'rev_') in the search database for
target-decoy FDR estimation. This script appends reversed versions of all
entries to any per-plex FASTA that does not already contain decoys.

Usage:
    python3 add_decoys.py --fasta-dir /path/to/FASTA/per_plex/

Options:
    --fasta-dir   Directory containing per-plex .fasta files (required)
    --prefix      Decoy header prefix (default: rev_)
    --force       Re-add decoys even if already present
"""

import argparse
import os
import sys


def parse_fasta(path):
    """Read a FASTA file. Returns list of (header, sequence) tuples."""
    entries = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq)))
                header, seq = line, []
            else:
                seq.append(line)
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


def has_decoys(path, prefix):
    """Return True if the file already contains decoy entries."""
    decoy_header = f">{prefix}"
    with open(path) as f:
        for line in f:
            if line.startswith(decoy_header):
                return True
    return False


def add_decoys(path, prefix, force=False):
    """
    Append reversed decoy entries to a FASTA file in place.

    Returns True if decoys were added, False if skipped.
    """
    if not force and has_decoys(path, prefix):
        return False

    entries = parse_fasta(path)

    with open(path, "a") as f:
        for header, seq in entries:
            decoy_header = f">{prefix}{header[1:]}"
            f.write(f"{decoy_header}\n{seq[::-1]}\n")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Append reversed decoy sequences to per-plex FASTA files."
    )
    parser.add_argument("--fasta-dir", required=True,
                        help="Directory containing per-plex .fasta files")
    parser.add_argument("--prefix", default="rev_",
                        help="Decoy header prefix (default: rev_)")
    parser.add_argument("--force", action="store_true",
                        help="Re-add decoys even if already present")

    args = parser.parse_args()

    if not os.path.isdir(args.fasta_dir):
        sys.exit(f"Error: FASTA directory not found: {args.fasta_dir}")

    fastas = sorted(f for f in os.listdir(args.fasta_dir)
                    if f.endswith(".fasta"))
    if not fastas:
        sys.exit(f"Error: No .fasta files found in {args.fasta_dir}")

    print(f"Adding decoys (prefix='{args.prefix}') to {len(fastas)} FASTAs...")

    added, skipped = 0, 0
    for fname in fastas:
        path = os.path.join(args.fasta_dir, fname)
        if add_decoys(path, args.prefix, force=args.force):
            added += 1
        else:
            skipped += 1

    print(f"  Added decoys: {added}")
    print(f"  Already had decoys (skipped): {skipped}")


if __name__ == "__main__":
    main()
