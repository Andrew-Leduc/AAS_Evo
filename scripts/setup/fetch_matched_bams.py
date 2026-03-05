#!/usr/bin/env python3
"""
fetch_matched_bams.py

PDC-first GDC metadata fetch.

Reads pdc_file_tmt_map.tsv (source of truth for all searched patients),
queries the GDC API for WXS BAMs belonging to those cases, and produces
gdc_meta_matched.tsv + the tissue download manifest directly.

This replaces the old GDC-first pipeline (fetch_gdc_manifest.py →
fetch_metadata.py → mapping_report.py → fetch_unmatched_bams.py) with a
single pass that has no missing-patient blind spots.

Usage:
    python3 scripts/setup/fetch_matched_bams.py

Prerequisites:
    - metadata/PDC_meta/pdc_file_tmt_map.tsv  (run fetch_pdc_metadata.py +
      consolidate_metadata.py first)

Outputs:
    - metadata/GDC_meta/gdc_meta_matched.tsv            (15-column metadata)
    - metadata/GDC_meta/manifests/manifest_wxs_bams_tissue.tsv  (5-col, no blood)
"""

import csv
import json
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import META_DIR


# ── Constants ──────────────────────────────────────────────────────────────────

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
BATCH_SIZE = 20      # case IDs per API request
PAGE_SIZE = 500      # results per page (GDC max)
REQUEST_TIMEOUT = 30
RETRY_ATTEMPTS = 3
DELAY_BETWEEN_REQUESTS = 0.5  # seconds

# 15-column schema — exact match to existing gdc_meta_matched.tsv
HEADERS = [
    'gdc_file_id', 'file_name', 'md5', 'size', 'state',
    'case_submitter_id', 'sample_submitter_id', 'aliquot_submitter_id',
    'sample_type', 'tissue_type', 'tissue_or_organ_of_origin',
    'project_id', 'disease_type', 'primary_site', 'unifier_id',
]

MANIFEST_HEADERS = ['id', 'filename', 'md5', 'size', 'state']

# Patient-channel filtering (reused from find_missing_bams.py)
CPTAC_CASE_RE = re.compile(r'^C3[LN]-\d{5}$')
SKIP_CASE_LOWER = {
    'ref', 'reference', 'pooled', 'pooled sample', 'pool',
    'nan', 'not reported', '',
}
SKIP_SAMPLE_TYPES = {'Not Reported', 'Cell Lines'}  # Blood kept in metadata; excluded from manifest
BLOOD_SAMPLE_TYPE = 'Blood Derived Normal'

# Output paths
GDC_META_OUT = META_DIR / 'GDC_meta' / 'gdc_meta_matched.tsv'
MANIFEST_OUT  = META_DIR / 'GDC_meta' / 'manifests' / 'manifest_wxs_bams_tissue.tsv'
TMT_MAP       = META_DIR / 'PDC_meta' / 'pdc_file_tmt_map.tsv'


# ── PDC input ──────────────────────────────────────────────────────────────────

def load_pdc_cases(tmt_map_path):
    """
    Return dict: case_submitter_id -> frozenset of sample_types.

    Only includes real CPTAC patient channels (C3L/C3N IDs), skipping
    reference/pooled channels and excluded sample types.
    """
    cases = defaultdict(set)
    with open(tmt_map_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            case_id     = (row.get('case_submitter_id') or '').strip()
            sample_type = (row.get('sample_type') or '').strip()

            if case_id.lower() in SKIP_CASE_LOWER:
                continue
            if not CPTAC_CASE_RE.match(case_id):
                continue
            if sample_type in SKIP_SAMPLE_TYPES:
                continue

            cases[case_id].add(sample_type)

    return {k: frozenset(v) for k, v in cases.items()}


# ── GDC API ────────────────────────────────────────────────────────────────────

def _query_url(url, timeout=REQUEST_TIMEOUT):
    """Single URL fetch with SSL fallback to curl."""
    for attempt in range(RETRY_ATTEMPTS):
        try:
            ctx = ssl.create_default_context()
            with urllib.request.urlopen(
                urllib.request.Request(url), timeout=timeout, context=ctx
            ) as resp:
                return json.loads(resp.read())
        except ssl.SSLError:
            try:
                import subprocess
                res = subprocess.run(
                    ['curl', '-sk', '--max-time', str(timeout), url],
                    capture_output=True, text=True,
                )
                if res.returncode == 0 and res.stdout.strip():
                    return json.loads(res.stdout)
            except Exception:
                pass
            return None
        except Exception:
            if attempt < RETRY_ATTEMPTS - 1:
                time.sleep(2 ** attempt)
            else:
                raise
    return None


def _build_query_params(batch, fields, size, from_=0):
    filters = json.dumps({
        'op': 'and',
        'content': [
            {'op': 'in',  'content': {'field': 'cases.submitter_id',    'value': batch}},
            {'op': '=',   'content': {'field': 'experimental_strategy', 'value': 'WXS'}},
            {'op': '=',   'content': {'field': 'data_format',           'value': 'BAM'}},
        ],
    })
    p = {'filters': filters, 'fields': fields, 'size': size}
    if from_:
        p['from'] = from_
    return urllib.parse.urlencode(p)


def fetch_wxs_bams_for_cases(case_ids):
    """
    Batch-query the GDC API for all WXS BAMs belonging to case_ids.
    Returns a list of raw GDC hit dicts.
    """
    fields = ','.join([
        'file_id', 'file_name', 'md5sum', 'file_size', 'state',
        'cases.submitter_id',
        'cases.samples.submitter_id',
        'cases.samples.sample_type',
        'cases.samples.tissue_type',
        'cases.samples.tissue_or_organ_of_origin',
        'cases.samples.portions.analytes.aliquots.submitter_id',
        'cases.project.project_id',
        'cases.project.disease_type',
        'cases.project.primary_site',
    ])

    case_list     = sorted(case_ids)
    total_batches = (len(case_list) + BATCH_SIZE - 1) // BATCH_SIZE
    all_hits      = []

    for batch_idx in range(0, len(case_list), BATCH_SIZE):
        batch     = case_list[batch_idx : batch_idx + BATCH_SIZE]
        batch_num = batch_idx // BATCH_SIZE + 1
        print(f'  Batch {batch_num}/{total_batches} ({len(batch)} cases)... ',
              end='', flush=True)

        url  = f'{GDC_FILES_ENDPOINT}?{_build_query_params(batch, fields, PAGE_SIZE)}'
        data = _query_url(url)
        if not data:
            print('no response')
            time.sleep(DELAY_BETWEEN_REQUESTS)
            continue

        hits  = data.get('data', {}).get('hits', [])
        total = data.get('data', {}).get('pagination', {}).get('total', 0)
        all_hits.extend(hits)
        print(f'{len(hits)}/{total} files')

        # Paginate if total > PAGE_SIZE
        for offset in range(PAGE_SIZE, total, PAGE_SIZE):
            page_url  = f'{GDC_FILES_ENDPOINT}?{_build_query_params(batch, fields, PAGE_SIZE, offset)}'
            page_data = _query_url(page_url)
            if page_data:
                all_hits.extend(page_data.get('data', {}).get('hits', []))
            time.sleep(DELAY_BETWEEN_REQUESTS)

        time.sleep(DELAY_BETWEEN_REQUESTS)

    return all_hits


# ── Hit parsing ────────────────────────────────────────────────────────────────

def parse_hit(hit, allowed_sample_types):
    """
    Convert a GDC API hit to zero or more 15-column metadata rows.

    Only samples whose sample_type is in allowed_sample_types are kept.
    Returns a list so that multi-sample hits produce multiple rows (rare
    for WXS, but handled correctly).
    """
    rows = []

    gdc_cases = hit.get('cases', [])
    if not gdc_cases:
        return rows

    case       = gdc_cases[0]
    case_id    = case.get('submitter_id', '')
    proj       = case.get('project', {})
    project_id = proj.get('project_id', '')
    disease    = proj.get('disease_type', '')
    site       = proj.get('primary_site', '')

    for sample in case.get('samples', []):
        sample_type = sample.get('sample_type', '')
        if sample_type not in allowed_sample_types:
            continue

        # Drill down to aliquot ID
        aliquot_id = ''
        portions   = sample.get('portions', [])
        if portions:
            analytes = portions[0].get('analytes', [])
            if analytes:
                aliquots = analytes[0].get('aliquots', [])
                if aliquots:
                    aliquot_id = aliquots[0].get('submitter_id', '')

        rows.append({
            'gdc_file_id':             hit.get('file_id', ''),
            'file_name':               hit.get('file_name', ''),
            'md5':                     hit.get('md5sum', ''),
            'size':                    str(hit.get('file_size', '')),
            'state':                   hit.get('state', ''),
            'case_submitter_id':       case_id,
            'sample_submitter_id':     sample.get('submitter_id', ''),
            'aliquot_submitter_id':    aliquot_id,
            'sample_type':             sample_type,
            'tissue_type':             sample.get('tissue_type', ''),
            'tissue_or_organ_of_origin': sample.get('tissue_or_organ_of_origin', ''),
            'project_id':              project_id,
            'disease_type':            disease,
            'primary_site':            site,
            'unifier_id':              aliquot_id,
        })

    return rows


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print('=' * 60)
    print('PDC-First GDC Metadata Fetch')
    print('=' * 60)

    if not TMT_MAP.exists():
        print(f'ERROR: TMT map not found: {TMT_MAP}')
        print('Run fetch_pdc_metadata.py + consolidate_metadata.py first.')
        sys.exit(1)

    # ── Step 1: Load PDC cases ─────────────────────────────────────────────────
    print(f'\n[1] Loading PDC TMT map: {TMT_MAP}')
    pdc_cases = load_pdc_cases(str(TMT_MAP))
    print(f'    {len(pdc_cases):,} unique CPTAC cases')

    stype_counts = defaultdict(int)
    for stypes in pdc_cases.values():
        for st in stypes:
            stype_counts[st] += 1
    print('    Sample types:')
    for st, n in sorted(stype_counts.items(), key=lambda x: -x[1]):
        print(f'      {st}: {n}')

    # ── Step 2: Query GDC API ──────────────────────────────────────────────────
    n_batches = (len(pdc_cases) + BATCH_SIZE - 1) // BATCH_SIZE
    print(f'\n[2] Querying GDC API ({len(pdc_cases):,} cases, {n_batches} batches)...')
    hits = fetch_wxs_bams_for_cases(list(pdc_cases.keys()))
    print(f'    {len(hits):,} raw hits returned')

    # ── Step 3: Parse and filter hits ─────────────────────────────────────────
    print('\n[3] Parsing hits and filtering to PDC sample types...')
    all_rows      = []
    seen_file_ids = set()

    for hit in hits:
        gdc_cases_list = hit.get('cases', [])
        if not gdc_cases_list:
            continue
        case_id = gdc_cases_list[0].get('submitter_id', '')
        if case_id not in pdc_cases:
            continue  # unexpected case from API — skip

        for row in parse_hit(hit, pdc_cases[case_id]):
            fid = row['gdc_file_id']
            if fid and fid not in seen_file_ids:
                seen_file_ids.add(fid)
                all_rows.append(row)

    # ── Step 4: Summary ────────────────────────────────────────────────────────
    found_cases = {r['case_submitter_id'] for r in all_rows}
    not_found   = set(pdc_cases.keys()) - found_cases

    type_counts = defaultdict(int)
    for r in all_rows:
        type_counts[r['sample_type']] += 1

    print(f'\n{"=" * 60}')
    print('RESULTS')
    print(f'{"=" * 60}')
    print(f'  PDC cases queried:       {len(pdc_cases):>6,}')
    print(f'  Cases found in GDC:      {len(found_cases):>6,}')
    print(f'  Cases NOT found in GDC:  {len(not_found):>6,}')
    print(f'  Total BAM rows:          {len(all_rows):>6,}')
    print(f'  BAMs by sample type:')
    for st, n in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f'    {st}: {n}')

    if not_found:
        print(f'\n  Cases with no WXS BAMs at GDC ({len(not_found)}):')
        for cid in sorted(not_found)[:20]:
            print(f'    {cid}')
        if len(not_found) > 20:
            print(f'    ... and {len(not_found) - 20} more')

    # ── Step 5: Write gdc_meta_matched.tsv ────────────────────────────────────
    GDC_META_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(GDC_META_OUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=HEADERS, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_rows)
    print(f'\n  gdc_meta_matched.tsv  → {GDC_META_OUT}  ({len(all_rows)} rows)')

    # ── Step 6: Write tissue manifest (blood excluded) ────────────────────────
    tissue_rows = [r for r in all_rows if r['sample_type'] != BLOOD_SAMPLE_TYPE]
    MANIFEST_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST_OUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_HEADERS, delimiter='\t')
        writer.writeheader()
        for r in tissue_rows:
            writer.writerow({
                'id':       r['gdc_file_id'],
                'filename': r['file_name'],
                'md5':      r['md5'],
                'size':     r['size'],
                'state':    r['state'],
            })
    print(f'  manifest_wxs_bams_tissue.tsv → {MANIFEST_OUT}  ({len(tissue_rows)} rows)')
    print()


if __name__ == '__main__':
    main()
