#!/usr/bin/env python3
"""
PDC Data Download Script

Downloads proteomics RAW files using pdc-client.

Usage:
    python download.py                           # Uses default matched manifest
    python download.py -m path/to/manifest.tsv   # Custom manifest
    python download.py --dry-run                 # Show command without running
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from config import META_DIR, RAW_DIR, ensure_dirs

DEFAULT_MANIFEST = META_DIR / 'PDC_meta' / 'pdc_all_files_matched.tsv'


def download_with_pdc_client(manifest_path, output_dir, n_processes=4, dry_run=False):
    """Download files using pdc-client."""
    cmd = [
        'pdc-client', 'download',
        '-m', str(manifest_path),
        '-d', str(output_dir),
        '-n', str(n_processes),
    ]

    print(f"Manifest: {manifest_path}")
    print(f"Output directory: {output_dir}")
    print(f"Parallel downloads: {n_processes}")
    print(f"\nCommand: {' '.join(cmd)}")

    if dry_run:
        print("\n[DRY RUN - command not executed]")
        return

    print("\nStarting download...\n")
    subprocess.run(cmd)


def main():
    parser = argparse.ArgumentParser(
        description='Download PDC RAW files using pdc-client',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python download.py                              # Download matched files
    python download.py -m custom_manifest.tsv       # Custom manifest
    python download.py -o /path/to/output           # Custom output directory
    python download.py -n 8                         # 8 parallel downloads
    python download.py --dry-run                    # Preview command
        """
    )
    parser.add_argument(
        '-m', '--manifest',
        default=str(DEFAULT_MANIFEST),
        help=f'Path to manifest file (default: {DEFAULT_MANIFEST})'
    )
    parser.add_argument(
        '-o', '--output',
        default=str(RAW_DIR),
        help=f'Output directory (default: {RAW_DIR})'
    )
    parser.add_argument(
        '-n', '--n-processes',
        type=int,
        default=4,
        help='Number of parallel download connections (default: 4)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show command without executing'
    )

    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"Error: Manifest not found: {manifest_path}")
        sys.exit(1)

    ensure_dirs()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    download_with_pdc_client(
        manifest_path,
        output_dir,
        n_processes=args.n_processes,
        dry_run=args.dry_run
    )


if __name__ == '__main__':
    main()
