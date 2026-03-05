#!/usr/bin/env python3
"""
Fetch metadata for GDC WXS BAM files from the GDC API.

Creates a metadata file with:
- gdc_file_id, file_name, patient_id, sample_type, tissue_type

Usage:
    python fetch_metadata.py
"""

import json
import requests
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from config import META_DIR

# Request settings
REQUEST_TIMEOUT = 30  # seconds
RETRY_ATTEMPTS = 3
DELAY_BETWEEN_BATCHES = 1  # seconds

WXS_MANIFEST = META_DIR / 'GDC_meta' / 'manifest_wxs_bams.tsv'
OUTPUT_META = META_DIR / 'GDC_meta' / 'gdc_meta.tsv'

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"

def get_file_ids_from_manifest():
    """Extract file IDs from WXS manifest."""
    file_ids = []
    with open(WXS_MANIFEST, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                file_ids.append(parts[0])
    return file_ids

def fetch_metadata_batch(file_ids, batch_size=100):
    """Fetch metadata from GDC API in batches."""
    all_results = []
    total_batches = (len(file_ids) + batch_size - 1) // batch_size

    for i in range(0, len(file_ids), batch_size):
        batch = file_ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        print(f"Fetching batch {batch_num}/{total_batches} ({len(batch)} files)...")

        filters = {
            "op": "in",
            "content": {
                "field": "file_id",
                "value": batch
            }
        }

        params = {
            "filters": json.dumps(filters),
            "fields": ",".join([
                "file_id",
                "file_name",
                "md5sum",
                "file_size",
                "state",
                "cases.submitter_id",
                "cases.samples.submitter_id",
                "cases.samples.sample_type",
                "cases.samples.tissue_type",
                "cases.samples.portions.analytes.aliquots.submitter_id",
                "cases.project.project_id",
                "cases.disease_type",
                "cases.primary_site",
            ]),
            "size": batch_size
        }

        # Retry logic
        for attempt in range(RETRY_ATTEMPTS):
            try:
                response = requests.get(GDC_FILES_ENDPOINT, params=params, timeout=REQUEST_TIMEOUT)
                response.raise_for_status()
                data = response.json()
                all_results.extend(data.get('data', {}).get('hits', []))
                break
            except (requests.Timeout, requests.RequestException) as e:
                if attempt < RETRY_ATTEMPTS - 1:
                    wait_time = 2 ** attempt  # exponential backoff
                    print(f"  Retry {attempt + 1}/{RETRY_ATTEMPTS} after {wait_time}s: {e}")
                    time.sleep(wait_time)
                else:
                    print(f"  Failed after {RETRY_ATTEMPTS} attempts: {e}")
                    raise

        # Delay between batches to avoid rate limiting
        if i + batch_size < len(file_ids):
            time.sleep(DELAY_BETWEEN_BATCHES)

    return all_results

def parse_results(results):
    """Parse API results into metadata rows."""
    rows = []

    for hit in results:
        file_id = hit.get('file_id', '')
        file_name = hit.get('file_name', '')
        md5 = hit.get('md5sum', '')
        size = hit.get('file_size', '')
        state = hit.get('state', '')

        cases = hit.get('cases', [])
        if not cases:
            continue

        case = cases[0]
        case_id = case.get('submitter_id', '')
        project_id = case.get('project', {}).get('project_id', '')
        disease_type = case.get('disease_type', '')
        primary_site = case.get('primary_site', '')

        samples = case.get('samples', [])
        if not samples:
            continue

        sample = samples[0]
        sample_id = sample.get('submitter_id', '')
        sample_type = sample.get('sample_type', '')  # Primary Tumor, Solid Tissue Normal, Blood Derived Normal
        tissue_type = sample.get('tissue_type', '')  # Tumor, Normal

        # Get aliquot ID
        aliquot_id = ''
        portions = sample.get('portions', [])
        if portions:
            analytes = portions[0].get('analytes', [])
            if analytes:
                aliquots = analytes[0].get('aliquots', [])
                if aliquots:
                    aliquot_id = aliquots[0].get('submitter_id', '')

        rows.append({
            'gdc_file_id': file_id,
            'file_name': file_name,
            'md5': md5,
            'size': size,
            'state': state,
            'case_submitter_id': case_id,
            'sample_submitter_id': sample_id,
            'aliquot_submitter_id': aliquot_id,
            'sample_type': sample_type,
            'tissue_type': tissue_type,
            'tissue_or_organ_of_origin': '',
            'project_id': project_id,
            'disease_type': disease_type,
            'primary_site': primary_site,
            'unifier_id': aliquot_id,
        })

    return rows

def write_metadata(rows):
    """Write metadata to TSV file."""
    headers = [
        'gdc_file_id', 'file_name', 'md5', 'size', 'state',
        'case_submitter_id', 'sample_submitter_id', 'aliquot_submitter_id',
        'sample_type', 'tissue_type', 'tissue_or_organ_of_origin',
        'project_id', 'disease_type', 'primary_site', 'unifier_id'
    ]

    with open(OUTPUT_META, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write('\t'.join(str(row.get(h, '')) for h in headers) + '\n')

    print(f"Wrote {len(rows)} rows to {OUTPUT_META}")

def main():
    print("Reading WXS manifest...")
    file_ids = get_file_ids_from_manifest()
    print(f"Found {len(file_ids)} WXS BAM files")

    print("\nFetching metadata from GDC API...")
    results = fetch_metadata_batch(file_ids)
    print(f"Retrieved {len(results)} results")

    print("\nParsing results...")
    rows = parse_results(results)

    print("\nWriting metadata...")
    write_metadata(rows)

if __name__ == '__main__':
    main()
