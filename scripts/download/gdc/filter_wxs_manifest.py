#!/usr/bin/env python3
"""
Filter GDC manifest to only include WXS BAM files.

Usage:
    python filter_wxs_manifest.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from config import META_DIR

INPUT_MANIFEST = META_DIR / 'GDC_meta' / 'gdc_manifest.2026-01-26.133020.txt'
OUTPUT_MANIFEST = META_DIR / 'GDC_meta' / 'manifest_wxs_bams.tsv'

def filter_wxs_bams():
    """Filter manifest to only WXS BAM files."""
    with open(INPUT_MANIFEST, 'r') as f:
        lines = f.readlines()

    header = lines[0]
    wxs_lines = [header]

    for line in lines[1:]:
        # WXS BAM files have _wxs_gdc_realn.bam in filename
        if '_wxs_gdc_realn.bam' in line:
            wxs_lines.append(line)

    with open(OUTPUT_MANIFEST, 'w') as f:
        f.writelines(wxs_lines)

    print(f"Filtered {len(lines)-1} total files to {len(wxs_lines)-1} WXS BAM files")
    print(f"Output: {OUTPUT_MANIFEST}")

if __name__ == '__main__':
    filter_wxs_manifest()
