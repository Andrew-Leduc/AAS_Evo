#!/usr/bin/env python3
"""
fetch_unmatched_bams.py

Find and fetch metadata for PDC patients that are missing from the GDC
matched manifest. Queries the GDC API by case_submitter_id to find WXS
BAM files, then outputs metadata in the same format as gdc_meta.tsv.

The resulting file can be appended to gdc_meta.tsv, then mapping_report.py
re-run to produce an updated gdc_meta_matched.tsv with better coverage.

Usage:
    python3 fetch_unmatched_bams.py \
        --gdc-meta /path/to/gdc_meta.tsv \
        --tmt-map /path/to/pdc_file_tmt_map.tsv \
        -o /path/to/gdc_meta_unmatched.tsv

    # Then merge and re-run matching:
    cat gdc_meta.tsv > gdc_meta_combined.tsv
    tail -n +2 gdc_meta_unmatched.tsv >> gdc_meta_combined.tsv
    cp gdc_meta_combined.tsv gdc_meta.tsv
    python3 scripts/download/mapping_report.py

    # Re-chunk for downloading:
    bash scripts/download/gdc/setup_chunks.sh
"""

import argparse
import csv
import json
import os
import ssl
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict


GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
BATCH_SIZE = 20  # Case IDs per API request
REQUEST_TIMEOUT = 30
RETRY_ATTEMPTS = 3
DELAY_BETWEEN_REQUESTS = 0.5


# Output columns (same as fetch_metadata.py / gdc_meta.tsv)
HEADERS = [
    'gdc_file_id', 'file_name', 'md5', 'size', 'state',
    'case_submitter_id', 'sample_submitter_id', 'aliquot_submitter_id',
    'sample_type', 'tissue_type', 'tissue_or_organ_of_origin',
    'project_id', 'disease_type', 'primary_site', 'unifier_id'
]


def load_gdc_case_ids(gdc_meta_path):
    """Load existing GDC case_submitter_ids."""
    case_ids = set()
    with open(gdc_meta_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            case_ids.add(row.get('case_submitter_id', '').strip())
    return case_ids


def load_pdc_case_ids(tmt_map_path):
    """Load unique PDC case_submitter_ids (excluding reference channels)."""
    case_ids = set()
    skip = {'ref', 'reference', 'pooled sample', 'pool', ''}
    with open(tmt_map_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            case_id = row.get('case_submitter_id', '').strip()
            if case_id.lower() not in skip:
                case_ids.add(case_id)
    return case_ids


def query_gdc_api(url, timeout=REQUEST_TIMEOUT):
    """Query GDC API with retry logic. Uses curl fallback for SSL issues."""
    for attempt in range(RETRY_ATTEMPTS):
        try:
            # Try urllib first
            ctx = ssl.create_default_context()
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
                return json.loads(resp.read())
        except ssl.SSLError:
            # Fall back to curl (common on macOS)
            try:
                import subprocess
                result = subprocess.run(
                    ['curl', '-sk', '--max-time', str(timeout), url],
                    capture_output=True, text=True
                )
                if result.returncode == 0 and result.stdout.strip():
                    return json.loads(result.stdout)
            except Exception:
                pass
        except Exception as e:
            if attempt < RETRY_ATTEMPTS - 1:
                time.sleep(2 ** attempt)
            else:
                raise

    return None


def find_wxs_bams_for_cases(case_ids):
    """
    Query GDC API to find WXS BAM files for a list of case_submitter_ids.

    Returns list of file metadata dicts.
    """
    all_files = []

    # Process in batches
    case_list = sorted(case_ids)
    total_batches = (len(case_list) + BATCH_SIZE - 1) // BATCH_SIZE

    for batch_idx in range(0, len(case_list), BATCH_SIZE):
        batch = case_list[batch_idx:batch_idx + BATCH_SIZE]
        batch_num = batch_idx // BATCH_SIZE + 1
        print(f"  Querying batch {batch_num}/{total_batches} "
              f"({len(batch)} cases)...")

        filters = json.dumps({
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "cases.submitter_id",
                        "value": batch
                    }
                },
                {
                    "op": "=",
                    "content": {
                        "field": "experimental_strategy",
                        "value": "WXS"
                    }
                },
                {
                    "op": "=",
                    "content": {
                        "field": "data_format",
                        "value": "BAM"
                    }
                },
            ]
        })

        fields = ",".join([
            "file_id", "file_name", "md5sum", "file_size", "state",
            "cases.submitter_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
            "cases.samples.tissue_type",
            "cases.samples.portions.analytes.aliquots.submitter_id",
            "cases.project.project_id",
            "cases.disease_type",
            "cases.primary_site",
        ])

        params = urllib.parse.urlencode({
            "filters": filters,
            "fields": fields,
            "size": 500,  # Max per request
        })

        url = f"{GDC_FILES_ENDPOINT}?{params}"

        try:
            data = query_gdc_api(url)
            if data:
                hits = data.get('data', {}).get('hits', [])
                total = data.get('data', {}).get('pagination', {}).get('total', 0)
                all_files.extend(hits)
                print(f"    Found {len(hits)} files (total available: {total})")

                # Handle pagination if more than 500 results
                if total > 500:
                    for offset in range(500, total, 500):
                        page_params = urllib.parse.urlencode({
                            "filters": filters,
                            "fields": fields,
                            "size": 500,
                            "from": offset,
                        })
                        page_url = f"{GDC_FILES_ENDPOINT}?{page_params}"
                        page_data = query_gdc_api(page_url)
                        if page_data:
                            page_hits = page_data.get('data', {}).get('hits', [])
                            all_files.extend(page_hits)
                        time.sleep(DELAY_BETWEEN_REQUESTS)
            else:
                print(f"    No response")
        except Exception as e:
            print(f"    ERROR: {e}")

        time.sleep(DELAY_BETWEEN_REQUESTS)

    return all_files


def parse_hits(hits):
    """Parse GDC API hits into metadata rows (same format as gdc_meta.tsv)."""
    rows = []

    for hit in hits:
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
        sample_type = sample.get('sample_type', '')
        tissue_type = sample.get('tissue_type', '')

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
            'size': str(size),
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


def main():
    parser = argparse.ArgumentParser(
        description="Find WXS BAMs for PDC patients missing from GDC metadata"
    )
    parser.add_argument("--gdc-meta", required=True,
                        help="Existing gdc_meta.tsv (or gdc_meta_matched.tsv)")
    parser.add_argument("--tmt-map", required=True,
                        help="PDC TMT mapping file (pdc_file_tmt_map.tsv)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV for newly found BAM metadata")

    args = parser.parse_args()

    print("=" * 60)
    print("Find Unmatched PDC Patients in GDC")
    print("=" * 60)

    # Load existing GDC case IDs
    print(f"\nLoading GDC metadata: {args.gdc_meta}")
    gdc_cases = load_gdc_case_ids(args.gdc_meta)
    print(f"  {len(gdc_cases)} unique GDC cases")

    # Load PDC case IDs
    print(f"Loading PDC TMT map: {args.tmt_map}")
    pdc_cases = load_pdc_case_ids(args.tmt_map)
    print(f"  {len(pdc_cases)} unique PDC cases")

    # Find unmatched
    unmatched = pdc_cases - gdc_cases
    print(f"\nUnmatched PDC cases (have proteomics, no genomics): {len(unmatched)}")

    if not unmatched:
        print("All PDC cases are already in GDC metadata. Nothing to do.")
        return

    # Query GDC API
    print(f"\nQuerying GDC API for WXS BAMs...")
    hits = find_wxs_bams_for_cases(unmatched)
    print(f"\nTotal hits from API: {len(hits)}")

    # Parse results
    rows = parse_hits(hits)

    # Filter to only the unmatched cases (API might return extras from batch queries)
    rows = [r for r in rows if r['case_submitter_id'] in unmatched]

    # Deduplicate by file_id
    seen = set()
    unique_rows = []
    for r in rows:
        if r['gdc_file_id'] not in seen:
            seen.add(r['gdc_file_id'])
            unique_rows.append(r)
    rows = unique_rows

    # Summary
    found_cases = set(r['case_submitter_id'] for r in rows)
    still_missing = unmatched - found_cases

    print(f"\n{'=' * 60}")
    print("RESULTS")
    print(f"{'=' * 60}")
    print(f"Unmatched PDC cases queried:   {len(unmatched)}")
    print(f"Cases with WXS BAMs found:     {len(found_cases)}")
    print(f"Cases with NO WXS BAMs at GDC: {len(still_missing)}")
    print(f"New BAM files found:           {len(rows)}")

    # Sample type breakdown
    type_counts = defaultdict(int)
    for r in rows:
        type_counts[r['sample_type']] += 1
    print(f"\nBy sample type:")
    for st, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {st}: {count}")

    # Write output
    if rows:
        with open(args.output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=HEADERS, delimiter='\t')
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

        print(f"\nOutput: {args.output}")
        print(f"\nTo merge with existing metadata and re-match:")
        print(f"  # Append new rows to gdc_meta.tsv")
        gdc_dir = os.path.dirname(args.gdc_meta)
        print(f"  tail -n +2 {args.output} >> {os.path.join(gdc_dir, 'gdc_meta.tsv')}")
        print(f"  # Re-run matching")
        print(f"  python3 scripts/download/mapping_report.py")
        print(f"  # Filter out blood, re-chunk")
        print(f"  head -1 gdc_meta_matched.tsv > gdc_meta_matched_tissue.tsv")
        print(f"  awk -F'\\t' '$9 != \"Blood Derived Normal\"' gdc_meta_matched.tsv | tail -n +2 >> gdc_meta_matched_tissue.tsv")
        print(f"  bash scripts/download/gdc/setup_chunks.sh")
    else:
        print("\nNo new BAMs found.")

    if still_missing:
        print(f"\n{len(still_missing)} cases have no WXS BAMs at GDC:")
        for case_id in sorted(still_missing)[:20]:
            print(f"  {case_id}")
        if len(still_missing) > 20:
            print(f"  ... and {len(still_missing) - 20} more")


if __name__ == "__main__":
    main()
