#!/usr/bin/env python3
"""
refresh_urls.py

Refresh expired PDC signed URLs by querying the GraphQL API.

This script reads an existing manifest file, extracts the unique study IDs,
queries the PDC GraphQL API for fresh signed URLs, and writes an updated
manifest with the new URLs.

Usage:
    python3 refresh_urls.py /path/to/pdc_all_files.tsv

    # Or specify output path:
    python3 refresh_urls.py /path/to/pdc_all_files.tsv -o /path/to/refreshed.tsv

    # Dry run (check without writing):
    python3 refresh_urls.py /path/to/pdc_all_files.tsv --dry-run
"""

import argparse
import csv
import json
import requests
import sys
import time
from datetime import datetime
from collections import defaultdict


PDC_GRAPHQL_URL = "https://pdc.cancer.gov/graphql"

# Rate limiting to be nice to the PDC API
REQUEST_DELAY = 1.0  # seconds between requests


def query_pdc(query, variables=None, retries=3):
    """Send a GraphQL query to PDC API with retries."""
    payload = {"query": query}
    if variables:
        payload["variables"] = variables

    headers = {"Content-Type": "application/json"}

    for attempt in range(retries):
        try:
            resp = requests.post(PDC_GRAPHQL_URL, json=payload, headers=headers, timeout=120)
            resp.raise_for_status()
            result = resp.json()

            if "errors" in result:
                print(f"  GraphQL errors: {result['errors']}")
                return None

            return result

        except requests.RequestException as e:
            if attempt < retries - 1:
                print(f"  Retry {attempt + 1}/{retries} after error: {e}")
                time.sleep(2 ** attempt)
            else:
                print(f"  ERROR: API request failed after {retries} attempts: {e}")
                return None

    return None


def get_files_for_study(pdc_study_id):
    """Query all files and their fresh signed URLs for a study."""
    query = """
    query FilesPerStudy($pdc_study_id: String!) {
        filesPerStudy(pdc_study_id: $pdc_study_id) {
            file_id
            file_name
            file_type
            data_category
            md5sum
            file_size
            signedUrl {
                url
            }
        }
    }
    """
    return query_pdc(query, {"pdc_study_id": pdc_study_id})


def get_unique_studies(manifest_path, delimiter='\t'):
    """Extract unique PDC study IDs from manifest."""
    studies = set()
    file_count = 0

    with open(manifest_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            file_count += 1
            study_id = row.get('PDC Study ID', '')
            if study_id:
                studies.add(study_id)

    return sorted(studies), file_count


def build_url_map(studies, progress_callback=None):
    """Query PDC API for all studies and build file_id -> signedUrl map."""
    url_map = {}  # file_id -> signedUrl
    file_info = {}  # file_id -> {file_name, file_type, etc.}

    for i, study_id in enumerate(studies):
        if progress_callback:
            progress_callback(i + 1, len(studies), study_id)

        result = get_files_for_study(study_id)

        if result and "data" in result:
            files = result["data"].get("filesPerStudy", [])

            for f in files:
                file_id = f.get("file_id")
                signed_url = f.get("signedUrl", {})
                url = signed_url.get("url") if signed_url else None

                if file_id and url:
                    url_map[file_id] = url
                    file_info[file_id] = {
                        "file_name": f.get("file_name"),
                        "file_type": f.get("file_type"),
                        "data_category": f.get("data_category"),
                        "md5sum": f.get("md5sum"),
                        "file_size": f.get("file_size"),
                    }

        # Rate limiting
        if i < len(studies) - 1:
            time.sleep(REQUEST_DELAY)

    return url_map, file_info


def refresh_manifest(manifest_path, output_path, delimiter='\t', dry_run=False):
    """Refresh all URLs in a manifest file."""

    print("=" * 60)
    print("PDC URL Refresh Tool")
    print("=" * 60)
    print(f"Input:  {manifest_path}")
    print(f"Output: {output_path}")
    print(f"Time:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Step 1: Get unique study IDs
    print("Step 1: Scanning manifest for study IDs...")
    studies, total_files = get_unique_studies(manifest_path, delimiter)
    print(f"  Found {len(studies)} unique studies, {total_files} files")
    print()

    # Step 2: Query API for fresh URLs
    print("Step 2: Querying PDC API for fresh signed URLs...")

    def progress(current, total, study_id):
        print(f"  [{current}/{total}] Fetching {study_id}...")

    url_map, file_info = build_url_map(studies, progress_callback=progress)
    print(f"  Retrieved {len(url_map)} signed URLs")
    print()

    if dry_run:
        print("DRY RUN - not writing output file")
        return len(url_map)

    # Step 3: Update manifest with new URLs
    print("Step 3: Writing updated manifest...")

    updated_count = 0
    missing_count = 0
    rows_written = 0

    with open(manifest_path, newline='', encoding='utf-8') as infile:
        reader = csv.DictReader(infile, delimiter=delimiter)
        fieldnames = reader.fieldnames

        with open(output_path, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()

            for row in reader:
                file_id = row.get('File ID', '')

                if file_id in url_map:
                    old_url = row.get('File Download Link', '')
                    new_url = url_map[file_id]

                    if old_url != new_url:
                        row['File Download Link'] = new_url
                        updated_count += 1
                else:
                    missing_count += 1

                writer.writerow(row)
                rows_written += 1

    print(f"  Rows written: {rows_written}")
    print(f"  URLs updated: {updated_count}")
    print(f"  Files not found in API: {missing_count}")
    print()

    print("=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Output file: {output_path}")
    print()
    print("Next steps:")
    print("  1. Resume downloads:")
    print("     sbatch scripts/download/pdc/submit_download.sh")
    print("     # or locally:")
    print(f"     python3 scripts/download/pdc/download.py {output_path} -o /path/to/RAW")

    return updated_count


def main():
    parser = argparse.ArgumentParser(
        description="Refresh expired PDC signed URLs via GraphQL API"
    )
    parser.add_argument("manifest", help="Path to existing manifest file (TSV or CSV)")
    parser.add_argument("-o", "--output", help="Output path (default: overwrite input)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Check API without writing output")

    args = parser.parse_args()

    # Determine delimiter
    if args.manifest.endswith('.csv'):
        delimiter = ','
    else:
        delimiter = '\t'

    # Output path
    output_path = args.output or args.manifest

    if not args.output and not args.dry_run:
        print(f"WARNING: Will overwrite {args.manifest}")
        print("Press Ctrl+C to cancel, or wait 3 seconds to continue...")
        try:
            time.sleep(3)
        except KeyboardInterrupt:
            print("\nCancelled.")
            sys.exit(0)

    refresh_manifest(args.manifest, output_path, delimiter, args.dry_run)


if __name__ == "__main__":
    main()
