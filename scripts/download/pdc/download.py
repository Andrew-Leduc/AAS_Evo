#!/usr/bin/env python3
"""
PDC Data Download Script

Downloads proteomics data (RAW files) from the Proteomics Data Commons portal.
Uses the PDC API: https://pdc.cancer.gov/pdc/

Usage:
    python download.py --study-id <study_id>
    python download.py --file-ids file1,file2,file3
"""

import argparse
import requests
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from config import RAW_DIR, ensure_dirs

PDC_API_URL = "https://pdc.cancer.gov/graphql"

def get_study_files(study_id: str):
    """Query PDC API for files in a study."""
    query = """
    query StudyFiles($study_id: String!) {
        filesPerStudy(study_id: $study_id) {
            file_id
            file_name
            file_size
            md5sum
        }
    }
    """

    response = requests.post(
        PDC_API_URL,
        json={'query': query, 'variables': {'study_id': study_id}}
    )
    response.raise_for_status()
    return response.json()

def download_file(file_id: str, output_dir: Path):
    """Download a single file from PDC."""
    # PDC download URL pattern
    url = f"https://pdc.cancer.gov/data-download/{file_id}"

    print(f"Downloading: {file_id}")
    print(f"URL: {url}")
    print(f"Output: {output_dir}")

    # Uncomment to execute actual download
    # response = requests.get(url, stream=True)
    # response.raise_for_status()
    # with open(output_dir / f"{file_id}.raw", 'wb') as f:
    #     for chunk in response.iter_content(chunk_size=8192):
    #         f.write(chunk)

def main():
    parser = argparse.ArgumentParser(description='Download PDC data')
    parser.add_argument('--study-id', help='PDC Study ID')
    parser.add_argument('--file-ids', help='Comma-separated list of file IDs')
    parser.add_argument('--dry-run', action='store_true', help='Print info without downloading')

    args = parser.parse_args()

    if not args.study_id and not args.file_ids:
        parser.error("Either --study-id or --file-ids required")

    ensure_dirs()

    if args.study_id:
        print(f"Querying files for study: {args.study_id}")
        files = get_study_files(args.study_id)
        print(f"Found files: {files}")

    if args.file_ids:
        for file_id in args.file_ids.split(','):
            download_file(file_id.strip(), RAW_DIR)

if __name__ == '__main__':
    main()
