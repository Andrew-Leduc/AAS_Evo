#!/usr/bin/env python3
"""
GDC Data Download Script

Downloads genomics data (BAM files) from the Genomic Data Commons portal.
Requires gdc-client to be installed: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool

Usage:
    python download.py --manifest manifest.txt
    python download.py --manifest manifest.txt --token token.txt
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from config import BAMS_DIR, ensure_dirs

def download_from_manifest(manifest_path: str, token_path: str = None):
    """Download files from GDC using a manifest file."""
    ensure_dirs()

    cmd = [
        'gdc-client', 'download',
        '-m', manifest_path,
        '-d', str(BAMS_DIR),
    ]

    if token_path:
        cmd.extend(['-t', token_path])

    print(f"Downloading to: {BAMS_DIR}")
    print(f"Command: {' '.join(cmd)}")

    # Uncomment to execute
    # subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description='Download GDC data')
    parser.add_argument('-m', '--manifest', required=True, help='Path to manifest file')
    parser.add_argument('-t', '--token', help='Path to GDC token file (for controlled access data)')
    parser.add_argument('--dry-run', action='store_true', help='Print command without executing')

    args = parser.parse_args()

    if not Path(args.manifest).exists():
        print(f"Error: Manifest file not found: {args.manifest}")
        sys.exit(1)

    download_from_manifest(args.manifest, args.token)

if __name__ == '__main__':
    main()
