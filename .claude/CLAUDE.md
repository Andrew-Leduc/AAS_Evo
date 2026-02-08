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
- 211 unique TMT plexes, ~1,308 unique PDC patients, ~1,057 matched to GDC
- Blood samples (426 BAMs) have genomics but no proteomics — excluded from processing
- Tissue-only BAMs after filtering: ~1,672
- ~251 PDC patients unmatched; ~135 recoverable via `fetch_unmatched_bams.py`

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
│   │   │   ├── fetch_unmatched_bams.py # Find WXS BAMs for unmatched PDC patients
│   │   │   ├── filter_wxs_manifest.py # Filter manifest to WXS BAMs only
│   │   │   └── setup_chunks.sh       # Split manifest into 500-BAM chunks
│   │   ├── pdc/
│   │   │   ├── download.py           # PDC download (streaming, rate-limited)
│   │   │   ├── consolidate_metadata.py # Merge per-tissue PDC manifests
│   │   │   ├── refresh_urls.py       # Refresh expired PDC signed URLs via API
│   │   │   ├── submit_download.sh    # SLURM job wrapper
│   │   │   └── download_files.sh     # Alternative bash downloader
│   │   └── mapping_report.py         # GDC-PDC sample matching & pruning
│   ├── proc_bams/                    # BAM processing pipeline
│   │   ├── run_pipeline.sh           # Wrapper: auto-generates file lists + submits jobs
│   │   ├── submit_variant_call.sh    # SLURM array job for variant calling
│   │   ├── submit_vep.sh             # SLURM array job for VEP annotation
│   │   └── consolidate_missense.sh   # Merge all missense mutations
│   ├── fasta_gen/                    # Custom proteogenomics FASTAs
│   │   ├── generate_mutant_fastas.py # Per-sample mutant FASTAs from VEP
│   │   ├── combine_plex_fastas.py    # Combine by TMT plex (ref + mut + comp)
│   │   ├── generate_compensatory_fastas.py # Compensatory mutation FASTAs
│   │   ├── submit_compensatory_fastas.sh   # SLURM wrapper
│   │   └── submit_proteogenomics.sh  # SLURM wrapper for FASTA generation
│   ├── mutation_analysis/            # Filtering, MSA generation & coevolution
│   │   ├── filter_and_rank.py        # Rank mutations by composite pathogenicity score
│   │   ├── generate_msas.py          # MSA generation via MMseqs2
│   │   ├── submit_msa_generation.sh  # SLURM array: one gene per task
│   │   ├── coevolution_analysis.py   # MI+APC covariation & compensatory prediction
│   │   └── submit_coevolution.sh     # SLURM wrapper for coevolution analysis
│   ├── ms_search/                    # FragPipe MS database search
│   │   ├── generate_manifests.py     # Per-plex FragPipe manifests + TMT annotations
│   │   ├── submit_fragpipe.sh        # SLURM array: one plex per task
│   │   └── run_ms_search.sh          # Orchestrator: generate manifests + submit
│   └── setup/
│       └── setup_seq_files.sh        # Download all external reference files
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
├── VCF/                              # Variant calls (per-chunk subdirectories)
│   ├── chunk_00/                     # VCFs from first BAM chunk
│   ├── chunk_01/                     # VCFs from second BAM chunk
│   └── ...
├── VEP/                              # VEP annotations (per-chunk subdirectories)
│   ├── chunk_00/                     # VEP output from first chunk
│   ├── chunk_01/                     # VEP output from second chunk
│   ├── ...
│   └── all_missense_mutations.tsv    # Consolidated (top-level, from consolidate_missense.sh)
├── SEQ_FILES/                        # Reference files
│   ├── hg38.fa                       # Human reference genome (GRCh38, UCSC)
│   ├── hg38.fa.fai                   # Reference index
│   ├── cds.chr.bed                   # CDS regions (merged, GENCODE)
│   ├── uniprot_human_canonical.fasta # UniProt reviewed proteome
│   ├── uniref50                      # MMseqs2 UniRef90 database
│   ├── AlphaMissense_hg38.tsv.gz     # AlphaMissense pathogenicity data
│   └── AlphaMissense_hg38.tsv.gz.tbi # tabix index
├── ANALYSIS/                         # Mutation filtering & ranking output
│   ├── ranked_mutations.tsv          # All mutations with composite scores
│   ├── top_5000_mutations.tsv        # Top N subset
│   ├── gene_list_for_msa.txt         # Unique genes from top N (used by MSA/coevolution)
│   ├── mutation_burden.tsv           # Per-protein-per-patient counts
│   └── ranking_summary.txt           # Statistics
├── FASTA/                            # Custom proteogenomics FASTAs
│   ├── per_sample/                   # {uuid}_mutant.fasta (mutants only)
│   ├── per_plex/                     # {run_metadata_id}.fasta (ref + mut + comp)
│   └── compensatory/                 # Compensatory mutation FASTAs
│       ├── {gene}_compensatory.fasta # Per-gene compensatory entries
│       └── all_compensatory.fasta    # Consolidated
├── MSA/                              # Per-gene MSAs ({accession}.a3m)
├── COEVOL/                           # Coevolution analysis output
├── MS_SEARCH/                        # FragPipe MS search setup & results
│   ├── manifests/                    # Per-plex .fp-manifest files
│   ├── annotations/                  # Per-plex TMT channel annotations
│   ├── workflows/                    # Per-plex workflow files (patched FASTA paths)
│   ├── results/                      # Per-plex FragPipe output
│   └── plex_list.txt                 # Plex IDs for SLURM array
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

### 4. Download Reference Files (One-Time Setup)
```bash
# Downloads and indexes all external reference files (~100 GB total):
# hg38.fa, cds.chr.bed, UniProt proteome, AlphaMissense, UniRef90, VEP container
sbatch scripts/setup/setup_seq_files.sh

# Or skip large files (UniRef90, VEP) for initial testing:
bash scripts/setup/setup_seq_files.sh --skip-large
```

### 5. Download & Process BAMs (Chunked)

BAMs are processed in chunks of ~400 to stay within scratch storage limits.
VCF/VEP outputs are stored in per-chunk subdirectories (e.g. `VCF/chunk_00/`,
`VEP/chunk_00/`). The chunk name must be passed to `run_pipeline.sh`.

```bash
# One-time: split manifest into 400-BAM chunks
bash scripts/download/gdc/setup_chunks.sh 400 path/to/manifest_wxs_bams_tissue.tsv

# Per chunk (repeat for each chunk manifest):
bash scripts/download/gdc/submit_download.sh path/to/chunks/chunk_00.tsv

# Step 1: Variant calling (finds BAMs automatically, outputs to VCF/chunk_00/)
bash scripts/proc_bams/run_pipeline.sh variant-call chunk_00

# Step 2: VEP annotation (reads VCF/chunk_00/, outputs to VEP/chunk_00/)
bash scripts/proc_bams/run_pipeline.sh vep chunk_00

# Step 3: Delete BAMs, download next chunk, repeat
rm -rf /scratch/leduc.an/AAS_Evo/BAMS/*
```

PDC RAW files (separate, not chunked):
```bash
sbatch scripts/download/pdc/submit_download.sh
```

### 6. Consolidate & Filter Mutations
```bash
# After ALL BAM chunks are processed:
bash scripts/proc_bams/consolidate_missense.sh
# Output: /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv (18 columns)

# Rank mutations by composite pathogenicity score
python3 scripts/mutation_analysis/filter_and_rank.py \
    --vep-tsv /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv \
    --ref-fasta /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta \
    -o /scratch/leduc.an/AAS_Evo/ANALYSIS
# Output: ANALYSIS/top_5000_mutations.tsv, ANALYSIS/gene_list_for_msa.txt
```

### 7. Generate MSAs for Coevolution Analysis
```bash
# Submit MSA generation (auto-finds gene list from ANALYSIS/)
NUM_GENES=$(wc -l < /scratch/leduc.an/AAS_Evo/ANALYSIS/gene_list_for_msa.txt)
sbatch --array=1-${NUM_GENES}%10 scripts/mutation_analysis/submit_msa_generation.sh
```

### 8. Coevolution Analysis (Compensatory Prediction)
```bash
# After MSA generation completes (auto-finds gene list from ANALYSIS/):
sbatch scripts/mutation_analysis/submit_coevolution.sh

# Or limit to specific genes:
sbatch scripts/mutation_analysis/submit_coevolution.sh TP53 BRCA1
# Output: COEVOL/compensatory_predictions.tsv
```

### 9. Generate Compensatory FASTAs
```bash
# After coevolution analysis completes:
sbatch scripts/fasta_gen/submit_compensatory_fastas.sh
# Output: FASTA/compensatory/all_compensatory.fasta
```

### 10. Generate Custom Proteogenomics FASTAs
```bash
# Prerequisites: reference proteome downloaded (included in setup_seq_files.sh)

# Generate per-sample mutant FASTAs + per-plex search databases
# Automatically includes compensatory entries if FASTA/compensatory/ exists
sbatch scripts/fasta_gen/submit_proteogenomics.sh
# Output: FASTA/per_plex/{run_metadata_id}.fasta (ref + mut + comp per plex)
```

### 11. MS Database Search (FragPipe)
```bash
# Set up per-plex manifests and submit FragPipe searches
# Option A: provide a FragPipe workflow template (auto-patches FASTA path per plex)
bash scripts/ms_search/run_ms_search.sh /path/to/template.workflow

# Option B: generate manifests first, configure workflow separately
bash scripts/ms_search/run_ms_search.sh
# -> Place workflow at MS_SEARCH/fragpipe.workflow, then:
NUM_PLEXES=$(wc -l < /scratch/leduc.an/AAS_Evo/MS_SEARCH/plex_list.txt)
sbatch --array=1-${NUM_PLEXES}%5 scripts/ms_search/submit_fragpipe.sh
# Output: MS_SEARCH/results/{plex_id}/
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
   - Includes gnomAD exome allele frequencies (`--af_gnomade`)
   - Includes AlphaMissense pathogenicity scores (`--plugin AlphaMissense`)
   - Extracts missense variants to per-sample TSV (18 columns)

3. **Final Output Columns** (`all_missense_mutations.tsv`, 18 columns)
   - `sample_id`: GDC UUID
   - `CHROM`, `POS`, `REF`, `ALT`: Variant location
   - `Consequence`: VEP consequence type (missense_variant)
   - `SYMBOL`: Gene symbol (e.g., TP53)
   - `Gene`: Ensembl gene ID
   - `HGVSp`: Protein change (e.g., p.Arg273His)
   - `Amino_acids`: R/H format
   - `Protein_position`: Position in protein
   - `gnomADe_AF`: gnomAD exome allele frequency
   - `am_pathogenicity`: AlphaMissense score (0–1, higher = more pathogenic)
   - `am_class`: AlphaMissense classification (likely_benign/ambiguous/likely_pathogenic)
   - `DP`, `AD_ref`, `AD_alt`: Read depth and allelic depths
   - `VAF`: Variant allele frequency in sample

## Mutation Filtering & Ranking Details

`filter_and_rank.py` prioritizes mutations for downstream coevolution analysis.

- **Mutation burden**: For each (sample, gene), counts mutations with VAF > 0.3
- **Pre-filter**: Excludes common polymorphisms (gnomADe_AF > 0.01)
- **Composite score** (0–1):
  - 50% AlphaMissense pathogenicity (missing = 0.5)
  - 20% gnomAD rarity (1 - AF, missing = 1.0)
  - 15% sample recurrence (min(1, n_samples/20))
  - 15% protein mutation burden context
- **Output**: Top N mutations (default 5000), gene list for MSA generation

## Proteogenomics Pipeline Details

Generates custom FASTA search databases with tryptic peptides for sample-specific missense mutations and predicted compensatory entries.

**Tryptic peptide approach**: Instead of adding full mutant proteins (~500 AA each), the pipeline extracts only tryptic peptides containing mutations (~15 AA avg, with 1 missed cleavage). This minimizes database size and improves FDR statistics.

**Header format**: `>type|accession|gene|swap|source|patient|sample_type`

1. **Per-sample mutant FASTAs** (`generate_mutant_fastas.py`)
   - Parses VEP `SYMBOL` column → looks up gene in UniProt via `GN=` header field
   - Parses mutation from `HGVSp` (e.g., `p.Arg273His`) or falls back to `Amino_acids` + `Protein_position`
   - Validates reference AA at position before substituting
   - Extracts tryptic peptide(s) containing the mutation (K/R cleavage, not before P, 1 missed cleavage)
   - Filters common variants by `gnomADe_AF` threshold
   - Headers: `>mut|P04637|TP53|R273H|genetic`
   - Logs issues to `generation_summary.tsv` and `generation_issues.tsv`

2. **Compensatory FASTAs** (`generate_compensatory_fastas.py`)
   - Reads coevolution predictions → applies both original + compensatory mutation to reference
   - Extracts tryptic peptide(s) containing mutation positions
   - Headers: `>comp|P04637|TP53|R273H_G245S|predicted`
   - Validates both mutation positions against reference, deduplicates within gene
   - Output: per-gene FASTAs + consolidated `all_compensatory.fasta`

3. **Per-plex FASTAs** (`combine_plex_fastas.py`)
   - Linking: `run_metadata_id` → `case_submitter_id` (via TMT map) → GDC UUID (via GDC meta) → mutant FASTA
   - Each plex FASTA = reference proteome + deduplicated mutant peptides + plex-specific compensatory peptides
   - Compensatory entries are **plex-specific**: only included if the original mutation is observed in that plex
   - Adds patient info to headers: `>mut|P04637|TP53|R273H|genetic|C3L-00001|tumor`
   - Deduplication by peptide sequence + mutation identity
   - Three FASTA header prefixes: `sp|`/`tr|` (reference), `mut|` (observed), `comp|` (compensatory)

4. **Reference proteome**: UniProt reviewed human canonical (~20,400 proteins, ~25MB)

## MS Search Pipeline Details

FragPipe MS searches are run per TMT plex, each against its own custom FASTA database.

- `generate_manifests.py` reads the TMT map to group RAW files by plex (`run_metadata_id`), checks for matching per-plex FASTAs, and generates FragPipe manifests + TMT channel annotations
- Each plex gets a separate FragPipe run because each has a different FASTA database
- Workflow template is patched per plex with the correct `database.db-path`
- TMT annotations map channels to `{case_submitter_id}_{sample_type}` (e.g., `C3L-00001_tumor`)
- SLURM array: 16 CPUs, 64G RAM, 24h per plex, 5 concurrent

## MSA Generation & Coevolution Pipeline Details

Predicts compensatory translation errors: given a destabilizing missense mutation at position i, finds covarying positions j via MI+APC and predicts which amino acid substitution at j could compensate.

### MSA Generation (`generate_msas.py`)

Two modes in one script:

1. **`--make-gene-list` mode**: Reads VEP mutations → filters to genes with UniProt entries → excludes genes with existing MSAs → writes gene list file
2. **Default mode**: Called by SLURM array job. Reads gene from list by index, extracts protein sequence from UniProt FASTA, runs MMseqs2 (`createdb` → `search` → `result2msa` → `unpackdb`), outputs `{accession}.a3m`

MSA files are named by UniProt accession (e.g., `P04637.a3m`). The gene list uses gene symbols from VEP; `build_gene_to_msa_map()` in the analysis script handles the mapping. Pre-existing MSAs named by gene symbol, accession, or entry name are all auto-detected.

MMseqs2 search parameters: 3 iterations, sensitivity 7.5, against UniRef90. Resources: 4 CPUs, 32G RAM, 4h per gene.

### Coevolution Analysis (`coevolution_analysis.py`)

For each gene with both mutations and an MSA:
1. Read MSA → numpy int8 array (N sequences x L columns, 0-19 for AAs, -1 for gaps)
2. Compute Neff (effective sequence count, identity threshold 0.8) → skip if < 50
3. Compute MI+APC matrix (L x L coupling scores). Memory-efficient: one (20,20) pair table at a time
4. For each mutation at position i: find top-k covarying positions j, compute conditional P(b at j | mut_aa at i), rank by preference shift vs wildtype baseline
5. Output predictions TSV

Pluggable backend via `CouplingBackend` class — MI+APC is default (pure numpy). Future: EVcouplings/plmDCA, MSA Transformer.

**Output columns** (`compensatory_predictions.tsv`):
`gene`, `uniprot_accession`, `mutation`, `mutation_hgvsp`, `n_samples`, `covarying_pos`, `wildtype_aa`, `predicted_compensatory_aa`, `coupling_score`, `conditional_score`, `preference_shift`, `neff`, `msa_depth`

**Additional outputs**: `_summary.tsv` (per-gene stats), `_skipped.tsv` (genes excluded with reasons)

## PDC Download Script Details

The `download.py` script has built-in rate limiting to avoid PDC restrictions:
- 10 downloads per 10-minute window (configurable in lines 48-50)
- 2-second pause between files
- Streaming downloads to avoid memory issues
- Auto-retry on transient errors

**Important**: PDC download URLs expire after 7 days. Use `refresh_urls.py` to programmatically refresh them (see below).

## PDC URL Refresh Workflow

PDC signed URLs are CloudFront URLs that expire after 7 days. When downloads start returning HTTP 403, refresh them via the PDC GraphQL API:

```bash
# Refresh URLs in place (backs up old file first):
python3 scripts/download/pdc/refresh_urls.py /path/to/pdc_all_files.tsv

# Or write to a new file:
python3 scripts/download/pdc/refresh_urls.py /path/to/pdc_all_files.tsv -o refreshed.tsv

# Dry run (check without writing):
python3 scripts/download/pdc/refresh_urls.py /path/to/pdc_all_files.tsv --dry-run
```

The script queries `https://pdc.cancer.gov/graphql` using the `filesPerStudy` query with `signedUrl { url }` to get fresh URLs for each PDC study ID in the manifest. No portal login required.

## Recovering Unmatched PDC Patients

~251 PDC patients have proteomics but are missing from the GDC matched metadata. ~135 of these DO have WXS BAMs at GDC that were missed by the initial manifest download. To recover them:

```bash
# Step 1: Find WXS BAMs for unmatched PDC patients via GDC API
python3 scripts/download/gdc/fetch_unmatched_bams.py \
    --gdc-meta /path/to/GDC_meta/gdc_meta.tsv \
    --tmt-map /path/to/PDC_meta/pdc_file_tmt_map.tsv \
    -o /path/to/GDC_meta/gdc_meta_unmatched.tsv

# Step 2: Append to existing metadata
tail -n +2 /path/to/GDC_meta/gdc_meta_unmatched.tsv >> /path/to/GDC_meta/gdc_meta.tsv

# Step 3: Re-run sample matching
python3 scripts/download/mapping_report.py

# Step 4: Filter out blood samples and re-chunk
cd /path/to/GDC_meta
head -1 gdc_meta_matched.tsv > gdc_meta_matched_tissue.tsv
awk -F'\t' '$9 != "Blood Derived Normal"' gdc_meta_matched.tsv | tail -n +2 >> gdc_meta_matched_tissue.tsv

# Step 5: Rebuild chunk manifests from tissue-only metadata
bash scripts/download/gdc/setup_chunks.sh
```

This adds ~135 patients (~272 BAMs, ~180 after excluding blood), bringing coverage from ~1,057 to ~1,192 patients (~89% of PDC patients). The remaining ~116 PDC patients genuinely have no WXS data at GDC.

## Blood Sample Exclusion

Blood Derived Normal samples (426 of 2,098 BAMs, ~20%) have no proteomics match and are not used as germline references in the current single-sample bcftools pipeline. They should be excluded before chunking:

```bash
head -1 gdc_meta_matched.tsv > gdc_meta_matched_tissue.tsv
awk -F'\t' '$9 != "Blood Derived Normal"' gdc_meta_matched.tsv | tail -n +2 >> gdc_meta_matched_tissue.tsv
```

This reduces processing from ~2,098 to ~1,672 BAMs. Already-processed blood BAM VEP files cause no harm (they just won't link to any TMT plex).

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
Source: UniProt reference proteome UP000005640 (reviewed/Swiss-Prot, canonical only)
```bash
wget -O uniprot_human_canonical.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29"
```
~25MB, ~20,400 proteins. One canonical protein per gene. Headers contain `GN=GENE_SYMBOL` used for VEP SYMBOL → protein mapping.

### UniRef90 Database (for MSA generation)
Source: UniProt Reference Clusters at 90% identity
```bash
# Download FASTA (~28 GB compressed)
wget -O uniref50.fasta.gz \
    https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz

# Build MMseqs2 database (~60 GB on disk)
mmseqs createdb uniref50.fasta.gz SEQ_FILES/uniref50
```
Used by `generate_msas.py` for MSA generation via MMseqs2 profile search.

### AlphaMissense Data (for VEP plugin)
Source: DeepMind AlphaMissense pathogenicity predictions
```bash
cd /scratch/leduc.an/AAS_Evo/SEQ_FILES
wget -O AlphaMissense_hg38.tsv.gz \
    "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
# Index with tabix (requires htslib)
tabix -s 1 -b 2 -e 2 -S 1 AlphaMissense_hg38.tsv.gz
```
~6 GB data file + `.tbi` index. Used by VEP `--plugin AlphaMissense`.

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

BAMs are processed in batches of ~400 to stay within scratch storage limits.

- `setup_chunks.sh` (in `download/gdc/`) splits the GDC manifest into persistent chunk files at `META_DIR/GDC_meta/manifests/chunks/chunk_NN.tsv`
- Chunk manifests are valid GDC format (with header), usable directly by `submit_download.sh`
- `run_pipeline.sh` (in `proc_bams/`) requires a chunk name argument (e.g. `chunk_00`), auto-generates file lists, and submits SLURM array jobs with `CHUNK_NAME` exported
- VCF and VEP outputs go into per-chunk subdirectories: `VCF/chunk_00/`, `VEP/chunk_00/`, etc.
- Both `submit_variant_call.sh` and `submit_vep.sh` have skip-if-done logic, so re-running is safe
- `consolidate_missense.sh` scans all `VEP/chunk_*/*.vep.tsv` and writes to top-level `VEP/all_missense_mutations.tsv`
- `generate_mutant_fastas.py` walks VEP subdirectories automatically to find all `.vep.tsv` files

## GDC Download Details

The GDC download uses `gdc-client` for controlled-access data (requires dbGaP authorization).

- Manifest is split into ~20 sub-manifests and submitted as separate SLURM jobs
- Each chunk runs on `short` partition with 24-hour time limit
- `gdc-client` auto-skips already-downloaded files (resumable)
- Token file: `gdc-user-token_AL.txt` (expires periodically, re-download from GDC portal)
- Each BAM is downloaded into its own UUID subdirectory under `BAMS/`
