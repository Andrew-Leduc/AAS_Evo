#!/usr/bin/env python3
"""
Generate a report on GDC (genomics) to PDC (proteomics) sample mapping.

Compares BAM files from GDC with TMT-labeled samples from PDC to find
matching case_submitter_id + sample_type combinations.

Usage:
    python mapping_report.py
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import META_DIR

GDC_META = META_DIR / 'GDC_meta' / 'gdc_meta.tsv'
PDC_TMT_MAP = META_DIR / 'PDC_meta' / 'pdc_file_tmt_map.tsv'
OUTPUT_REPORT = META_DIR / 'mapping_report.tsv'

def load_gdc_samples():
    """Load GDC samples keyed by (case_id, sample_type)."""
    samples = defaultdict(list)

    with open(GDC_META, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            case_id = row.get('case_submitter_id', '')
            sample_type = row.get('sample_type', '')
            file_name = row.get('file_name', '')
            tissue_type = row.get('tissue_type', '')
            primary_site = row.get('primary_site', '')

            if case_id and sample_type:
                key = (case_id, sample_type)
                samples[key].append({
                    'file_name': file_name,
                    'tissue_type': tissue_type,
                    'primary_site': primary_site,
                })

    return samples

def load_pdc_samples():
    """Load PDC samples keyed by (case_id, sample_type)."""
    samples = defaultdict(list)

    with open(PDC_TMT_MAP, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            case_id = row.get('case_submitter_id', '')
            sample_type = row.get('sample_type', '')
            file_name = row.get('file_name', '')
            tmt_channel = row.get('tmt_channel', '')
            tissue_folder = row.get('tissue_folder', '')

            if case_id and sample_type:
                key = (case_id, sample_type)
                samples[key].append({
                    'file_name': file_name,
                    'tmt_channel': tmt_channel,
                    'tissue_folder': tissue_folder,
                })

    return samples

def normalize_sample_type(sample_type):
    """Normalize sample type for matching."""
    st = sample_type.lower().strip()
    if 'tumor' in st and 'normal' not in st:
        return 'tumor'
    elif 'normal' in st:
        return 'normal'
    elif 'blood' in st:
        return 'blood'
    return st

def is_blood_sample(sample_type):
    """Check if sample type is blood-derived (no proteomics expected)."""
    return 'blood' in sample_type.lower()

def generate_report(gdc_samples, pdc_samples):
    """Generate mapping report."""

    # Build normalized lookup for PDC
    pdc_normalized = defaultdict(list)
    for (case_id, sample_type), entries in pdc_samples.items():
        norm_type = normalize_sample_type(sample_type)
        pdc_normalized[(case_id, norm_type)].extend(entries)

    # Statistics
    gdc_total = len(gdc_samples)
    gdc_blood = 0  # Blood samples (no proteomics expected)
    gdc_tissue = 0  # Tissue samples (proteomics expected)
    matched = 0
    unmatched_gdc = []
    matched_details = []

    for (case_id, sample_type), gdc_entries in gdc_samples.items():
        norm_type = normalize_sample_type(sample_type)
        pdc_key = (case_id, norm_type)

        # Track blood vs tissue samples
        if is_blood_sample(sample_type):
            gdc_blood += 1
        else:
            gdc_tissue += 1

        if pdc_key in pdc_normalized:
            matched += 1
            pdc_entries = pdc_normalized[pdc_key]

            # Count unique raw files for this sample
            unique_raw_files = len(set(e['file_name'] for e in pdc_entries))

            matched_details.append({
                'case_id': case_id,
                'sample_type': sample_type,
                'gdc_bam_count': len(gdc_entries),
                'pdc_raw_count': unique_raw_files,
                'pdc_tmt_entries': len(pdc_entries),
                'gdc_bam': gdc_entries[0]['file_name'],
                'primary_site': gdc_entries[0]['primary_site'],
                'is_blood': is_blood_sample(sample_type),
            })
        else:
            unmatched_gdc.append({
                'case_id': case_id,
                'sample_type': sample_type,
                'gdc_bam_count': len(gdc_entries),
                'primary_site': gdc_entries[0]['primary_site'],
                'is_blood': is_blood_sample(sample_type),
            })

    # PDC samples not in GDC
    gdc_normalized = set()
    for (case_id, sample_type) in gdc_samples.keys():
        gdc_normalized.add((case_id, normalize_sample_type(sample_type)))

    unmatched_pdc = []
    for (case_id, sample_type), entries in pdc_samples.items():
        norm_type = normalize_sample_type(sample_type)
        if (case_id, norm_type) not in gdc_normalized:
            unique_raw = len(set(e['file_name'] for e in entries))
            unmatched_pdc.append({
                'case_id': case_id,
                'sample_type': sample_type,
                'pdc_raw_count': unique_raw,
                'tissue_folder': entries[0]['tissue_folder'] if entries else '',
            })

    # Count unique raw files for matched samples (avoiding double-counting)
    matched_raw_files = set()
    for m in matched_details:
        case_id = m['case_id']
        sample_type = m['sample_type']
        norm_type = normalize_sample_type(sample_type)
        pdc_key = (case_id, norm_type)
        if pdc_key in pdc_normalized:
            for entry in pdc_normalized[pdc_key]:
                matched_raw_files.add(entry['file_name'])

    return {
        'gdc_total': gdc_total,
        'gdc_blood': gdc_blood,
        'gdc_tissue': gdc_tissue,
        'pdc_total': len(pdc_samples),
        'matched': matched,
        'unmatched_gdc': unmatched_gdc,
        'unmatched_pdc': unmatched_pdc,
        'matched_details': matched_details,
        'matched_raw_files': len(matched_raw_files),
    }

def print_report(report):
    """Print report to console."""
    print("=" * 60)
    print("GDC to PDC SAMPLE MAPPING REPORT")
    print("=" * 60)

    # Count blood in unmatched (expected) vs tissue unmatched (unexpected)
    unmatched_blood = sum(1 for u in report['unmatched_gdc'] if u['is_blood'])
    unmatched_tissue = sum(1 for u in report['unmatched_gdc'] if not u['is_blood'])

    # Matched tissue samples (excluding any blood that somehow matched)
    matched_tissue = sum(1 for m in report['matched_details'] if not m['is_blood'])

    print(f"\nSUMMARY:")
    print(f"  GDC samples (BAM files):     {report['gdc_total']:,}")
    print(f"    - Tissue samples:          {report['gdc_tissue']:,}")
    print(f"    - Blood samples:           {report['gdc_blood']:,} (no proteomics expected)")
    print(f"  PDC samples (TMT channels):  {report['pdc_total']:,}")
    print(f"  Matched samples:             {report['matched']:,}")
    print(f"  GDC-only (no proteomics):    {len(report['unmatched_gdc']):,}")
    print(f"    - Blood (expected):        {unmatched_blood:,}")
    print(f"    - Tissue (missing):        {unmatched_tissue:,}")
    print(f"  PDC-only (no genomics):      {len(report['unmatched_pdc']):,}")

    # Overall match rate
    match_rate = (report['matched'] / report['gdc_total'] * 100) if report['gdc_total'] > 0 else 0
    print(f"\n  Overall match rate:          {match_rate:.1f}%")

    # Tissue-only match rate (expected match rate)
    tissue_match_rate = (matched_tissue / report['gdc_tissue'] * 100) if report['gdc_tissue'] > 0 else 0
    print(f"  Tissue match rate:           {tissue_match_rate:.1f}% (excludes blood samples)")

    # Sample type breakdown
    print(f"\nMATCHED BY SAMPLE TYPE:")
    type_counts = defaultdict(int)
    for m in report['matched_details']:
        type_counts[m['sample_type']] += 1
    for st, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {st}: {count}")

    # Tissue breakdown for matched
    print(f"\nMATCHED BY TISSUE (primary_site):")
    tissue_counts = defaultdict(int)
    for m in report['matched_details']:
        tissue_counts[m['primary_site']] += 1
    for tissue, count in sorted(tissue_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {tissue}: {count}")

    # File counts
    if report['matched_details']:
        total_bam = sum(m['gdc_bam_count'] for m in report['matched_details'])
        unique_raw = report.get('matched_raw_files', 0)

        print(f"\nFILE COUNTS:")
        print(f"  GDC BAM files (matched):     {total_bam:,}")
        print(f"  PDC RAW files (matched):     {unique_raw:,}")

        print(f"\n  Structure: Each TMT plex has ~11 samples × ~25 fractions")
        print(f"             {unique_raw:,} RAW files contain data for {report['matched']:,} matched samples")

def write_matched_report(report):
    """Write detailed matched samples to TSV."""
    with open(OUTPUT_REPORT, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['case_id', 'sample_type', 'primary_site',
                        'gdc_bam_count', 'pdc_raw_count', 'pdc_tmt_entries', 'gdc_bam'])

        for m in sorted(report['matched_details'], key=lambda x: x['case_id']):
            writer.writerow([
                m['case_id'],
                m['sample_type'],
                m['primary_site'],
                m['gdc_bam_count'],
                m['pdc_raw_count'],
                m['pdc_tmt_entries'],
                m['gdc_bam'],
            ])

    print(f"\nDetailed report written to: {OUTPUT_REPORT}")

def main():
    print("Loading GDC metadata...")
    gdc_samples = load_gdc_samples()
    print(f"  Loaded {len(gdc_samples)} unique GDC samples")

    print("Loading PDC TMT mapping...")
    pdc_samples = load_pdc_samples()
    print(f"  Loaded {len(pdc_samples)} unique PDC samples")

    print("\nGenerating report...")
    report = generate_report(gdc_samples, pdc_samples)

    print_report(report)
    write_matched_report(report)

if __name__ == '__main__':
    main()
