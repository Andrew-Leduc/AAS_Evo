#!/usr/bin/env python3
"""
Consolidate PDC metadata from all tissue folders.

Each tissue folder contains:
- PDC_file_manifest_*.csv: Raw file download info
- PDC_study_biospecimen_*.csv: Sample/patient metadata
- PDC_study_experimental_*.csv: TMT plex layout

Creates unified metadata file with:
- file_name, file_id, patient_id (case_submitter_id), sample_type, tissue_type

Usage:
    python consolidate_metadata.py
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from config import META_DIR

PDC_META_DIR = META_DIR / 'PDC_meta'
OUTPUT_FILE_MANIFEST = META_DIR / 'PDC_meta' / 'pdc_all_files.tsv'
OUTPUT_SAMPLE_META = META_DIR / 'PDC_meta' / 'pdc_meta.tsv'

def find_csv_files(folder, prefix):
    """Find CSV files matching prefix in folder."""
    for f in folder.glob(f'{prefix}*.csv'):
        return f
    return None

def read_csv_with_bom(filepath):
    """Read CSV handling BOM marker."""
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def consolidate_all_tissues():
    """Consolidate data from all tissue folders."""
    all_files = []
    all_biospecimens = []
    all_experiments = []

    tissue_folders = [d for d in PDC_META_DIR.iterdir() if d.is_dir()]

    for tissue_dir in tissue_folders:
        tissue_name = tissue_dir.name
        print(f"Processing {tissue_name}...")

        # Find files
        file_manifest = find_csv_files(tissue_dir, 'PDC_file_manifest')
        biospecimen = find_csv_files(tissue_dir, 'PDC_study_biospecimen')
        experimental = find_csv_files(tissue_dir, 'PDC_study_experimental')

        if file_manifest:
            files = read_csv_with_bom(file_manifest)
            for f in files:
                f['tissue_folder'] = tissue_name
            all_files.extend(files)
            print(f"  - {len(files)} files")

        if biospecimen:
            samples = read_csv_with_bom(biospecimen)
            for s in samples:
                s['tissue_folder'] = tissue_name
            all_biospecimens.extend(samples)
            print(f"  - {len(samples)} biospecimens")

        if experimental:
            exps = read_csv_with_bom(experimental)
            for e in exps:
                e['tissue_folder'] = tissue_name
            all_experiments.extend(exps)
            print(f"  - {len(exps)} experimental runs")

    return all_files, all_biospecimens, all_experiments

def build_aliquot_to_sample_map(biospecimens):
    """Map aliquot ID to sample metadata."""
    aliquot_map = {}
    for b in biospecimens:
        aliquot_id = b.get('Aliquot Submitter ID', '')
        if aliquot_id:
            aliquot_map[aliquot_id] = {
                'case_submitter_id': b.get('Case Submitter ID', ''),
                'sample_submitter_id': b.get('Sample Submitter ID', ''),
                'sample_type': b.get('Sample Type', ''),
                'tissue_type': b.get('Tissue Type', ''),
                'primary_site': b.get('Primary Site', ''),
                'disease_type': b.get('Disease Type', ''),
            }
    return aliquot_map

def build_run_to_samples_map(experiments, aliquot_map):
    """Map Run Metadata ID to samples in that plex."""
    tmt_channels = [
        'tmt_126', 'tmt_127n', 'tmt_127c', 'tmt_128n', 'tmt_128c',
        'tmt_129n', 'tmt_129c', 'tmt_130n', 'tmt_130c', 'tmt_131', 'tmt_131c'
    ]

    run_map = {}
    for exp in experiments:
        run_id = exp.get('Run Metadata ID', '')
        if not run_id:
            continue

        samples_in_run = []
        for channel in tmt_channels:
            value = exp.get(channel, '')
            if value:
                # Format is typically "AliquotID\nSampleType" or just aliquot ID
                lines = value.strip().split('\n')
                aliquot_id = lines[0].strip() if lines else ''

                # Look up sample info
                if aliquot_id in aliquot_map:
                    sample_info = aliquot_map[aliquot_id].copy()
                    sample_info['aliquot_id'] = aliquot_id
                    sample_info['tmt_channel'] = channel
                    samples_in_run.append(sample_info)

        run_map[run_id] = samples_in_run

    return run_map

def write_file_manifest(files):
    """Write consolidated file manifest."""
    headers = ['file_id', 'file_name', 'run_metadata_id', 'study_name',
               'pdc_study_id', 'file_size', 'md5sum', 'tissue_folder']

    with open(OUTPUT_FILE_MANIFEST, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)

        for file in files:
            writer.writerow([
                file.get('File ID', ''),
                file.get('File Name', ''),
                file.get('Run Metadata ID', ''),
                file.get('Study Name', ''),
                file.get('PDC Study ID', ''),
                file.get('File Size (in bytes)', ''),
                file.get('Md5sum', ''),
                file.get('tissue_folder', ''),
            ])

    print(f"\nWrote {len(files)} files to {OUTPUT_FILE_MANIFEST}")

def write_sample_metadata(biospecimens):
    """Write consolidated sample metadata."""
    headers = ['aliquot_submitter_id', 'sample_submitter_id', 'case_submitter_id',
               'sample_type', 'tissue_type', 'primary_site', 'disease_type', 'tissue_folder']

    with open(OUTPUT_SAMPLE_META, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)

        for b in biospecimens:
            writer.writerow([
                b.get('Aliquot Submitter ID', ''),
                b.get('Sample Submitter ID', ''),
                b.get('Case Submitter ID', ''),
                b.get('Sample Type', ''),
                b.get('Tissue Type', ''),
                b.get('Primary Site', ''),
                b.get('Disease Type', ''),
                b.get('tissue_folder', ''),
            ])

    print(f"Wrote {len(biospecimens)} samples to {OUTPUT_SAMPLE_META}")

def main():
    print("Consolidating PDC metadata from all tissue folders...\n")

    all_files, all_biospecimens, all_experiments = consolidate_all_tissues()

    print(f"\nTotals:")
    print(f"  Files: {len(all_files)}")
    print(f"  Biospecimens: {len(all_biospecimens)}")
    print(f"  Experimental runs: {len(all_experiments)}")

    # Build mapping
    aliquot_map = build_aliquot_to_sample_map(all_biospecimens)
    run_map = build_run_to_samples_map(all_experiments, aliquot_map)

    # Write outputs
    write_file_manifest(all_files)
    write_sample_metadata(all_biospecimens)

    print("\nDone!")

if __name__ == '__main__':
    main()
