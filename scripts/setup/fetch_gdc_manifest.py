#!/usr/bin/env python3
"""
Fetch GDC WXS BAM file manifest via the GDC REST API.

Queries for all CPTAC Whole Exome Sequencing BAM files and writes a manifest
in standard GDC format (id, filename, md5, size, state). This replaces the
manual "download manifest from portal" step.

Usage:
    python3 scripts/setup/fetch_gdc_manifest.py
"""

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import META_DIR

OUTPUT = META_DIR / 'GDC_meta' / 'manifests' / 'manifest_wxs_bams.tsv'
STUDIES_FILE = META_DIR / 'studies.tsv'
GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"

PAGE_SIZE = 500
REQUEST_TIMEOUT = 30
RETRY_ATTEMPTS = 3


def query_gdc(url, timeout=REQUEST_TIMEOUT):
    """Query GDC API with retry logic."""
    for attempt in range(RETRY_ATTEMPTS):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read())
        except Exception as e:
            if attempt < RETRY_ATTEMPTS - 1:
                time.sleep(2 ** attempt)
            else:
                raise


def load_gdc_projects():
    """Read unique GDC project IDs from studies.tsv."""
    projects = set()
    with open(STUDIES_FILE, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pid = row.get('gdc_project_id', '').strip()
            if pid:
                projects.add(pid)
    return sorted(projects)


def fetch_wxs_bam_manifest():
    """Fetch all CPTAC WXS BAM files from GDC API."""
    projects = load_gdc_projects()
    if not projects:
        print("  WARNING: No gdc_project_id values found in studies.tsv")
        return []

    filters = json.dumps({
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": projects,
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "experimental_strategy",
                    "value": "WXS",
                }
            },
            {
                "op": "=",
                "content": {
                    "field": "data_format",
                    "value": "BAM",
                }
            },
        ]
    })

    fields = "file_id,file_name,md5sum,file_size,state"

    all_files = []
    offset = 0

    # First request to get total count
    params = urllib.parse.urlencode({
        "filters": filters,
        "fields": fields,
        "size": PAGE_SIZE,
        "from": offset,
    })
    url = f"{GDC_FILES_ENDPOINT}?{params}"
    data = query_gdc(url)

    total = data.get("data", {}).get("pagination", {}).get("total", 0)
    hits = data.get("data", {}).get("hits", [])
    all_files.extend(hits)
    print(f"  Total WXS BAMs at GDC: {total}")
    print(f"  Fetched {len(all_files)}/{total}...")

    # Page through remaining results
    while len(all_files) < total:
        offset += PAGE_SIZE
        params = urllib.parse.urlencode({
            "filters": filters,
            "fields": fields,
            "size": PAGE_SIZE,
            "from": offset,
        })
        url = f"{GDC_FILES_ENDPOINT}?{params}"
        data = query_gdc(url)
        hits = data.get("data", {}).get("hits", [])
        if not hits:
            break
        all_files.extend(hits)
        print(f"  Fetched {len(all_files)}/{total}...")
        time.sleep(0.3)

    return all_files


def write_manifest(files, output_path):
    """Write files to GDC manifest format."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w") as f:
        f.write("id\tfilename\tmd5\tsize\tstate\n")
        for hit in files:
            f.write("\t".join([
                hit.get("file_id", ""),
                hit.get("file_name", ""),
                hit.get("md5sum", ""),
                str(hit.get("file_size", "")),
                hit.get("state", ""),
            ]) + "\n")

    print(f"\nManifest written: {output_path}")
    print(f"  {len(files)} WXS BAM files")


def main():
    print("=" * 60)
    print("Fetch GDC WXS BAM Manifest")
    print("=" * 60)
    print(f"\nQuerying GDC API for CPTAC WXS BAMs...")
    projects = load_gdc_projects()
    print(f"  Projects (from studies.tsv): {', '.join(projects)}")

    files = fetch_wxs_bam_manifest()
    write_manifest(files, OUTPUT)


if __name__ == "__main__":
    main()
