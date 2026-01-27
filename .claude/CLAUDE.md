# AAS_Evo Project Context

## Project Overview

Multi-omics data pipeline for CPTAC3 (Clinical Proteomic Tumor Analysis Consortium) data. Downloads and processes matched genomics (WXS BAM files from GDC) and proteomics (TMT-labeled RAW files from PDC) data for integrated analysis.

## Data Sources

### GDC (Genomic Data Commons)
- **Data type**: Whole Exome Sequencing (WXS) BAM files
- **Portal**: https://portal.gdc.cancer.gov
- **Download tool**: `gdc-client`
- **Files**: ~3,619 BAM files total, ~2,098 matched to proteomics

### PDC (Proteomics Data Commons)
- **Data type**: TMT-labeled mass spectrometry RAW files
- **Portal**: https://pdc.cancer.gov
- **Download**: Signed URLs (expire after 7 days)
- **Files**: ~5,215 RAW files total, ~5,020 matched to genomics (~3.77 TB)

### Sample Matching
- Samples matched by `case_submitter_id` + `sample_type` (normalized to tumor/normal/blood)
- TMT multiplexing: ~11 samples per plex, ~25 fractions per plex
- Blood samples have genomics but no proteomics (expected)
- Match rate: ~72% for tissue samples

## Directory Structure

### Code Repository
```
AAS_Evo/                              # This repo
├── config/
│   ├── __init__.py
│   └── paths.py                      # Environment-aware path config
├── scripts/
│   ├── download/
│   │   ├── gdc/
│   │   │   ├── download.py           # GDC download via gdc-client
│   │   │   ├── fetch_metadata.py     # Fetch sample metadata from GDC API
│   │   │   └── filter_wxs_manifest.py # Filter manifest to WXS BAMs only
│   │   ├── pdc/
│   │   │   ├── download.py           # PDC download (streaming, rate-limited)
│   │   │   ├── consolidate_metadata.py # Merge per-tissue PDC manifests
│   │   │   ├── submit_download.sh    # SLURM job wrapper
│   │   │   └── download_files.sh     # Alternative bash downloader
│   │   └── mapping_report.py         # GDC-PDC sample matching & pruning
│   └── preproc/                      # Preprocessing (WIP)
│       ├── ms_raw/                   # RAW file processing
│       └── wxs/                      # BAM file processing
└── utils/
```

### Metadata Directory (AAS_Evo_meta/)
```
AAS_Evo_meta/
├── GDC_meta/
│   ├── gdc_manifest.*.txt            # Full GDC manifest from portal
│   ├── manifest_wxs_bams.tsv         # Filtered to WXS BAMs
│   ├── gdc_meta.tsv                  # Enriched metadata (from API)
│   └── gdc_meta_matched.tsv          # Pruned to matched samples only
├── PDC_meta/
│   ├── {tissue}/                     # Per-tissue folders (LSCC, LUAD, etc.)
│   │   ├── PDC_file_manifest_*.csv   # Raw file info + download URLs
│   │   ├── PDC_study_biospecimen_*.csv # Sample/patient metadata
│   │   └── PDC_study_experimental_*.csv # TMT plex layouts
│   ├── pdc_all_files.tsv             # Consolidated manifest (all tissues)
│   ├── pdc_all_files_matched.tsv     # Pruned to matched samples only
│   ├── pdc_file_tmt_map.tsv          # File -> TMT channel -> sample mapping
│   └── pdc_meta.tsv                  # Consolidated sample metadata
└── mapping_report.tsv                # Detailed match report
```

### Data Directory (cluster: /scratch/leduc.an/AAS_Evo/)
```
/scratch/leduc.an/AAS_Evo/
├── BAMS/                             # GDC BAM files
├── RAW/                              # PDC RAW files
└── logs/                             # SLURM job logs
```

## Environments

### Local (macOS)
- Scripts: `/Users/andrewleduc/Desktop/Github/AAS_Evo`
- Meta: `/Users/andrewleduc/Desktop/AAS_Evo_meta`
- Auto-detected when `/scratch/leduc.an` doesn't exist

### Cluster (Northeastern Discovery)
- Scripts: `/home/leduc.an/AAS_Evo_project/AAS_Evo`
- Meta: `/home/leduc.an/AAS_Evo_project/AAS_Evo_meta`
- Data: `/scratch/leduc.an/AAS_Evo`
- Partition: `slavov`

## Key Workflows

### 1. Prepare GDC Metadata
```bash
# Download manifest from GDC portal, then:
python scripts/download/gdc/filter_wxs_manifest.py  # Filter to WXS BAMs
python scripts/download/gdc/fetch_metadata.py       # Enrich with API metadata
```

### 2. Prepare PDC Metadata
```bash
# Download manifests from PDC portal for each tissue, then:
python scripts/download/pdc/consolidate_metadata.py # Merge all tissues
```

### 3. Generate Matched Manifests
```bash
python scripts/download/mapping_report.py
# Creates: gdc_meta_matched.tsv, pdc_all_files_matched.tsv
```

### 4. Download Data (Cluster)
```bash
# PDC RAW files
sbatch scripts/download/pdc/submit_download.sh

# GDC BAM files (requires gdc-client + token for controlled access)
python scripts/download/gdc/download.py -m manifest.txt -t token.txt
```

## PDC Download Script Details

The `download.py` script has built-in rate limiting to avoid PDC restrictions:
- 10 downloads per 10-minute window (configurable in lines 48-50)
- 2-second pause between files
- Streaming downloads to avoid memory issues
- Auto-retry on transient errors

**Important**: PDC download URLs expire after 7 days. Re-export manifests from portal if downloads fail.

## Tissue Types in Dataset

PDC tissues: CCRCC, GMB_conf, GMG_disco, HNSCC, LSCC, LUAD, LUAD_disco, PDA, PDAC_enrich, UCEC, UCEC_disco, non-ccRCC

Missing from PDC (no proteomics): Stomach, some Kidney cases

## Column Name Conventions

### PDC Manifest (download.py compatible)
`File ID`, `File Name`, `Run Metadata ID`, `Study Name`, `PDC Study ID`, `PDC Study Version`, `Data Category`, `File Type`, `File Size (in bytes)`, `Md5sum`, `tissue_folder`, `File Download Link`

### TMT Mapping (internal)
`file_name`, `run_metadata_id`, `tmt_channel`, `case_submitter_id`, `sample_type`, `tissue_type`, `aliquot_submitter_id`, `tissue_folder`

### GDC Metadata
`file_id`, `file_name`, `case_submitter_id`, `sample_type`, `tissue_type`, `primary_site`, `disease_type`, `project_id`
