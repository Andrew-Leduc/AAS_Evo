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
│   │   │   ├── submit_download.sh    # Split manifest & submit SLURM jobs
│   │   │   ├── download_chunk.sh     # Single-chunk SLURM download job
│   │   │   ├── download.py           # GDC download via gdc-client
│   │   │   ├── fetch_metadata.py     # Fetch sample metadata from GDC API
│   │   │   ├── filter_wxs_manifest.py # Filter manifest to WXS BAMs only
│   │   │   └── setup_chunks.sh       # Split manifest into 500-BAM chunks
│   │   ├── pdc/
│   │   │   ├── download.py           # PDC download (streaming, rate-limited)
│   │   │   ├── consolidate_metadata.py # Merge per-tissue PDC manifests
│   │   │   ├── submit_download.sh    # SLURM job wrapper
│   │   │   └── download_files.sh     # Alternative bash downloader
│   │   └── mapping_report.py         # GDC-PDC sample matching & pruning
│   ├── proc_bams/                    # BAM processing pipeline
│   │   ├── run_pipeline.sh           # Wrapper: auto-generates file lists + submits jobs
│   │   ├── submit_variant_call.sh    # SLURM array job for variant calling
│   │   ├── submit_vep.sh             # SLURM array job for VEP annotation
│   │   ├── consolidate_missense.sh   # Merge all missense mutations
│   │   ├── align_and_variant_call.sh # Single-sample test script
│   │   └── run_vep.sh                # Single-sample VEP test
│   └── proteogenomics/               # Custom FASTA generation
│       ├── generate_mutant_fastas.py # Per-sample mutant FASTAs from VEP
│       ├── combine_plex_fastas.py    # Combine by TMT plex
│       └── submit_proteogenomics.sh  # SLURM wrapper
└── utils/
```

### Metadata Directory (AAS_Evo_meta/)
```
AAS_Evo_meta/
├── GDC_meta/
│   ├── gdc_manifest.*.txt            # Full GDC manifest from portal
│   ├── manifest_wxs_bams.tsv         # Filtered to WXS BAMs
│   ├── gdc_meta.tsv                  # Enriched metadata (from API)
│   ├── gdc_meta_matched.tsv          # Pruned to matched samples only
│   └── gdc-user-token_AL.txt         # GDC access token (controlled data)
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
├── SEQ_FILES/                        # Reference files
│   ├── hg38.fa                       # Human reference genome (GRCh38, UCSC)
│   ├── hg38.fa.fai                   # Reference index
│   ├── cds.chr.bed                   # CDS regions (merged, GENCODE)
│   └── uniprot_human_canonical.fasta # UniProt reviewed proteome
├── FASTA/                            # Custom proteogenomics FASTAs
│   ├── per_sample/                   # {uuid}_mutant.fasta (mutants only)
│   └── per_plex/                     # {run_metadata_id}.fasta (ref + mutants)
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

### 4. Download & Process BAMs (Chunked)

BAMs are processed in chunks of ~500 to stay within scratch storage limits.
VCF/VEP outputs persist across chunks; only BAMs are deleted between chunks.

```bash
# One-time: split manifest into 500-BAM chunks
bash scripts/download/gdc/setup_chunks.sh

# Per chunk (repeat for each chunk manifest):
bash scripts/download/gdc/submit_download.sh path/to/chunk_00.tsv

# Step 1: Variant calling (finds BAMs automatically)
bash scripts/proc_bams/run_pipeline.sh variant-call

# Step 2: VEP annotation (finds VCFs automatically)
bash scripts/proc_bams/run_pipeline.sh vep

# Step 3: Delete BAMs, download next chunk, repeat
rm -rf /scratch/leduc.an/AAS_Evo/BAMS/*
```

PDC RAW files (separate, not chunked):
```bash
sbatch scripts/download/pdc/submit_download.sh
```

### 5. Consolidate Missense Mutations
```bash
# After ALL chunks are processed:
bash scripts/proc_bams/consolidate_missense.sh
# Output: /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv
```

### 6. Generate Custom Proteogenomics FASTAs
```bash
# One-time: download UniProt reference proteome
wget -O /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29"

# Generate per-sample mutant FASTAs + per-plex search databases
sbatch scripts/proteogenomics/submit_proteogenomics.sh
```

## BAM Processing Pipeline Details

The BAM processing pipeline extracts missense mutations for multi-omics integration:

1. **Variant Calling** (`submit_variant_call.sh`)
   - Uses bcftools mpileup/call on CDS regions only (via `-R cds.chr.bed`)
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

## Proteogenomics Pipeline Details

Generates custom protein FASTA search databases with sample-specific missense mutations.

1. **Per-sample mutant FASTAs** (`generate_mutant_fastas.py`)
   - Parses VEP `SYMBOL` column → looks up gene in UniProt via `GN=` header field
   - Parses mutation from `HGVSp` (e.g., `p.Arg273His`) or falls back to `Amino_acids` + `Protein_position`
   - Validates reference AA at position before substituting
   - Headers: `>mut|P04637|TP53_R273H OS=Homo sapiens GN=TP53`
   - Multiple mutations in same gene → separate entries (standard SAV approach)
   - Logs issues to `generation_summary.tsv` and `generation_issues.tsv`

2. **Per-plex FASTAs** (`combine_plex_fastas.py`)
   - Linking: `run_metadata_id` → `case_submitter_id` (via TMT map) → GDC UUID (via GDC meta) → mutant FASTA
   - Each plex FASTA = full reference proteome + deduplicated mutant entries from all plex samples
   - Deduplication by mutation identity (accession + mutation label)

3. **Reference proteome**: UniProt reviewed human canonical (~20,400 proteins, ~25MB)

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

## Reference File Acquisition

Reference files are stored in `/scratch/leduc.an/AAS_Evo/SEQ_FILES/`.

### Human Reference Genome (hg38.fa)
Source: UCSC Genome Browser (uses `chr` prefix matching GDC BAM alignments)
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

### CDS Regions BED (cds.chr.bed)
Source: GENCODE annotation (human gene annotations, coding sequences only).
GENCODE uses `chr` prefixes (chr1, chr2, ...) matching the UCSC hg38 reference and GDC BAM alignments.
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

# Extract CDS regions, filter to standard chromosomes, sort, and merge overlapping intervals
zcat gencode.v46.annotation.gtf.gz \
    | awk '$3 == "CDS" {print $1"\t"$4-1"\t"$5}' \
    | grep -E '^chr([0-9]+|X|Y|M)\b' \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > cds.chr.bed
```
The `grep` step removes alt/patch contigs (chrUn, chrGL, chrKI) not present in GDC BAMs.
The `bedtools merge` step is critical — without it, overlapping CDS intervals cause `bcftools mpileup -R` to emit unsorted positions.

### UniProt Reference Proteome (uniprot_human_canonical.fasta)
Source: UniProt reviewed (Swiss-Prot) human canonical sequences
```bash
wget -O uniprot_human_canonical.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29"
```
~25MB, ~20,400 proteins. Headers contain `GN=GENE_SYMBOL` used for VEP→protein mapping.

### VEP Container and Cache
```bash
# Ensembl VEP Apptainer image
apptainer pull ensembl-vep.sif docker://ensemblorg/ensembl-vep
# VEP cache (offline mode)
mkdir -p cache && cd cache
apptainer exec ensembl-vep.sif vep_install --AUTO c --ASSEMBLY GRCh38 --SPECIES homo_sapiens
```
Stored at `/scratch/leduc.an/tools/vep/`

## Chunk Workflow Details

BAMs are processed in batches of ~500 to stay within scratch storage limits.

- `setup_chunks.sh` (in `download/gdc/`) splits the GDC manifest into persistent chunk files at `META_DIR/GDC_meta/manifests/chunks/chunk_NN.tsv`
- Chunk manifests are valid GDC format (with header), usable directly by `submit_download.sh`
- `run_pipeline.sh` (in `proc_bams/`) auto-generates file lists and submits SLURM array jobs for variant calling and VEP
- Both `submit_variant_call.sh` and `submit_vep.sh` have skip-if-done logic, so re-running is safe
- VCF/VEP outputs accumulate across chunks; `consolidate_missense.sh` and proteogenomics scripts run once after all chunks

## GDC Download Details

The GDC download uses `gdc-client` for controlled-access data (requires dbGaP authorization).

- Manifest is split into ~20 sub-manifests and submitted as separate SLURM jobs
- Each chunk runs on `short` partition with 24-hour time limit
- `gdc-client` auto-skips already-downloaded files (resumable)
- Token file: `gdc-user-token_AL.txt` (expires periodically, re-download from GDC portal)
- Each BAM is downloaded into its own UUID subdirectory under `BAMS/`
