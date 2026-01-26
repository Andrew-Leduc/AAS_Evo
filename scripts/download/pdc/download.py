#!/usr/bin/env python3
"""
PDC Data Download Script

Downloads proteomics RAW files from signed URLs in a manifest file.

Usage:
    python download.py                           # Uses default matched manifest
    python download.py -m path/to/manifest.tsv   # Custom manifest
    python download.py --dry-run                 # Show what would be downloaded
"""

import argparse
import csv
import os
import sys
import time
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from config import META_DIR, RAW_DIR, ensure_dirs

DEFAULT_MANIFEST = META_DIR / 'PDC_meta' / 'pdc_all_files_matched.tsv'


def load_manifest(manifest_path):
    """Load file manifest and return list of files to download."""
    files = []
    with open(manifest_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            files.append({
                'file_name': row.get('file_name', ''),
                'file_size': int(row.get('file_size', 0) or 0),
                'md5sum': row.get('md5sum', ''),
                'download_url': row.get('download_url', ''),
            })
    return files


def file_exists_with_size(filepath, expected_size):
    """Check if file exists and has the expected size."""
    if not filepath.exists():
        return False
    return filepath.stat().st_size == expected_size


def download_file(file_info, output_dir, dry_run=False):
    """Download a single file."""
    file_name = file_info['file_name']
    url = file_info['download_url']
    expected_size = file_info['file_size']
    output_path = output_dir / file_name

    # Skip if no URL
    if not url:
        return {'file': file_name, 'status': 'skipped', 'reason': 'no URL'}

    # Skip if already exists with correct size
    if file_exists_with_size(output_path, expected_size):
        return {'file': file_name, 'status': 'exists', 'size': expected_size}

    if dry_run:
        return {'file': file_name, 'status': 'would_download', 'size': expected_size}

    # Download file
    try:
        urllib.request.urlretrieve(url, output_path)
        actual_size = output_path.stat().st_size
        return {'file': file_name, 'status': 'downloaded', 'size': actual_size}
    except Exception as e:
        # Clean up partial download
        if output_path.exists():
            output_path.unlink()
        return {'file': file_name, 'status': 'failed', 'error': str(e)}


def download_files(manifest_path, output_dir, dry_run=False, workers=4):
    """Download all files from manifest."""
    print(f"Loading manifest: {manifest_path}")
    files = load_manifest(manifest_path)
    print(f"Found {len(files)} files in manifest")

    if not files:
        print("No files to download")
        return

    # Calculate total size
    total_size = sum(f['file_size'] for f in files)
    print(f"Total size: {total_size / 1e9:.1f} GB")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    if dry_run:
        print("\n[DRY RUN - no files will be downloaded]\n")

    # Track progress
    downloaded = 0
    skipped = 0
    failed = 0
    downloaded_size = 0

    start_time = time.time()

    # Download files (sequential for now - parallel can cause issues with signed URLs)
    for i, file_info in enumerate(files, 1):
        file_name = file_info['file_name']
        print(f"[{i}/{len(files)}] {file_name}...", end=' ', flush=True)

        result = download_file(file_info, output_dir, dry_run)

        if result['status'] == 'downloaded':
            downloaded += 1
            downloaded_size += result.get('size', 0)
            print(f"OK ({result.get('size', 0) / 1e6:.1f} MB)")
        elif result['status'] == 'exists':
            skipped += 1
            print("exists")
        elif result['status'] == 'would_download':
            print(f"would download ({result.get('size', 0) / 1e6:.1f} MB)")
        elif result['status'] == 'skipped':
            skipped += 1
            print(f"skipped ({result.get('reason', '')})")
        else:
            failed += 1
            print(f"FAILED: {result.get('error', 'unknown error')}")

    elapsed = time.time() - start_time

    # Summary
    print("\n" + "=" * 50)
    print("DOWNLOAD SUMMARY")
    print("=" * 50)
    print(f"Downloaded: {downloaded} files ({downloaded_size / 1e9:.2f} GB)")
    print(f"Skipped:    {skipped} files (already exist)")
    print(f"Failed:     {failed} files")
    print(f"Time:       {elapsed:.1f} seconds")
    if downloaded > 0 and elapsed > 0:
        print(f"Speed:      {downloaded_size / 1e6 / elapsed:.1f} MB/s")


def main():
    parser = argparse.ArgumentParser(
        description='Download PDC RAW files from manifest',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python download.py                              # Download matched files
    python download.py -m custom_manifest.tsv       # Custom manifest
    python download.py -o /path/to/output           # Custom output directory
    python download.py --dry-run                    # Preview without downloading
        """
    )
    parser.add_argument(
        '-m', '--manifest',
        default=str(DEFAULT_MANIFEST),
        help=f'Path to manifest TSV file (default: {DEFAULT_MANIFEST})'
    )
    parser.add_argument(
        '-o', '--output',
        default=str(RAW_DIR),
        help=f'Output directory (default: {RAW_DIR})'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be downloaded without actually downloading'
    )

    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"Error: Manifest not found: {manifest_path}")
        sys.exit(1)

    ensure_dirs()
    download_files(manifest_path, args.output, dry_run=args.dry_run)


if __name__ == '__main__':
    main()
