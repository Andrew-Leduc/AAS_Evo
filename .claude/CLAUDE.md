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
│   └── proc_bams/                    # BAM processing pipeline
│       ├── submit_variant_call.sh    # SLURM array job for variant calling
│       ├── submit_vep.sh             # SLURM array job for VEP annotation
│       ├── consolidate_missense.sh   # Merge all missense mutations
│       ├── align_and_variant_call.sh # Single-sample test script
│       └── run_vep.sh                # Single-sample VEP test
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
├── BAMS/                             # GDC BAM files (by UUID subdirectory)
├── RAW/                              # PDC RAW files
├── VCF/                              # Variant calls (from BAM processing)
├── VEP/                              # VEP annotations + missense tables
├── ref/                              # Reference files
│   ├── hg38.fa                       # Human reference genome
│   ├── hg38.fa.fai                   # Reference index
│   └── cds.bed                       # CDS regions (speeds up variant calling)
├── bam_list.txt                      # List of BAM paths for array jobs
├── vcf_list.txt                      # List of VCF paths for VEP
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

### 5. Process BAM Files (Variant Calling)
```bash
# Generate BAM list
find /scratch/leduc.an/AAS_Evo/BAMS -name "*.bam" > /scratch/leduc.an/AAS_Evo/bam_list.txt

# Submit variant calling array job
NUM_BAMS=$(wc -l < /scratch/leduc.an/AAS_Evo/bam_list.txt)
sbatch --array=1-${NUM_BAMS}%10 scripts/proc_bams/submit_variant_call.sh
```

### 6. Annotate Variants with VEP
```bash
# Generate VCF list (after step 5 completes)
ls /scratch/leduc.an/AAS_Evo/VCF/*.vcf.gz > /scratch/leduc.an/AAS_Evo/vcf_list.txt

# Submit VEP annotation array job
NUM_VCFS=$(wc -l < /scratch/leduc.an/AAS_Evo/vcf_list.txt)
sbatch --array=1-${NUM_VCFS}%10 scripts/proc_bams/submit_vep.sh
```

### 7. Consolidate Missense Mutations
```bash
# Merge all per-sample missense TSVs into one table
bash scripts/proc_bams/consolidate_missense.sh
# Output: /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv
```

## BAM Processing Pipeline Details

The BAM processing pipeline extracts missense mutations for multi-omics integration:

1. **Variant Calling** (`submit_variant_call.sh`)
   - Uses bcftools mpileup/call on CDS regions only (via `-R cds.bed`)
   - Filters: QUAL≥20, depth≥10
   - Outputs VCF + TSV with VAF per sample

2. **VEP Annotation** (`submit_vep.sh`)
   - Runs Ensembl VEP via Apptainer container
   - Adds gene symbols, HGVS notation, protein changes
   - Includes gnomAD population frequencies (`--af_gnomad`)
   - Extracts missense variants to per-sample TSV

3. **Final Output Columns** (`all_missense_mutations.tsv`)
   - `sample_id`: GDC UUID
   - `CHROM`, `POS`, `REF`, `ALT`: Variant location
   - `SYMBOL`: Gene symbol (e.g., TP53)
   - `HGVSp`: Protein change (e.g., p.Arg273His)
   - `Amino_acids`: R/H format
   - `gnomAD_AF`: Population frequency
   - `VAF`: Variant allele frequency in sample

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
