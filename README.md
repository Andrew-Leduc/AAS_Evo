# AAS_Evo

Multi-omics data pipeline for CPTAC3 (Clinical Proteomic Tumor Analysis Consortium). Downloads and processes matched genomics (WXS BAM files from GDC) and proteomics (TMT-labeled RAW files from PDC) data for integrated analysis of amino acid substitutions.

## Repository Structure

```
AAS_Evo/
├── config/
│   └── paths.py                      # Environment-aware path config
├── metadata/                            # Tracked metadata (in-repo)
│   ├── studies.tsv                      # Study registry (tissue → PDC study ID)
│   ├── GDC_meta/                        # GDC manifests & sample metadata
│   │   └── manifests/manifest_wxs_bams.tsv
│   ├── PDC_meta/                        # PDC per-tissue CSVs & consolidated files
│   │   ├── {tissue}/PDC_file_manifest.csv
│   │   ├── {tissue}/PDC_study_biospecimen.csv
│   │   └── {tissue}/PDC_study_experimental.csv
│   └── mapping_report.tsv              # GDC↔PDC matching report
├── scripts/
│   ├── setup/
│   │   ├── setup_metadata.sh            # Orchestrator: fetch all metadata from APIs
│   │   ├── setup_seq_files.sh           # Download all external reference files
│   │   ├── fetch_gdc_manifest.py        # Query GDC REST API for WXS BAM manifest
│   │   └── fetch_pdc_metadata.py        # Query PDC GraphQL API for per-tissue metadata
│   ├── download/
│   │   ├── gdc/
│   │   │   ├── submit_download.sh       # Submit SLURM download jobs
│   │   │   ├── download_chunk.sh        # Single-chunk SLURM download job
│   │   │   ├── download.py              # GDC download via gdc-client
│   │   │   ├── fetch_metadata.py        # Fetch sample metadata from GDC API
│   │   │   ├── fetch_unmatched_bams.py  # Find WXS BAMs for unmatched PDC patients
│   │   │   └── setup_chunks.sh          # Split manifest into chunks
│   │   ├── pdc/
│   │   │   ├── submit_download.sh       # SLURM job wrapper
│   │   │   ├── download.py              # PDC download (rate-limited)
│   │   │   ├── refresh_urls.py          # Refresh expired PDC signed URLs via API
│   │   │   └── consolidate_metadata.py  # Merge per-tissue CSVs
│   │   └── mapping_report.py            # GDC-PDC sample matching
│   ├── proc_bams/
│   │   ├── run_pipeline.sh              # Wrapper: auto-generates file lists + submits jobs
│   │   ├── submit_variant_call.sh       # SLURM array: variant calling
│   │   ├── submit_vep.sh               # SLURM array: VEP annotation + AlphaMissense
│   │   └── consolidate_missense.sh      # Merge missense mutations
│   ├── mutation_analysis/               # Filtering, MSA generation & coevolution
│   │   ├── filter_and_rank.py
│   │   ├── generate_msas.py
│   │   ├── submit_msa_generation.sh
│   │   ├── coevolution_analysis.py
│   │   └── submit_coevolution.sh
│   ├── fasta_gen/                       # Custom proteogenomics FASTAs
│   │   ├── generate_mutant_fastas.py
│   │   ├── combine_plex_fastas.py
│   │   ├── generate_compensatory_fastas.py
│   │   ├── submit_compensatory_fastas.sh
│   │   └── submit_proteogenomics.sh
│   └── ms_search/                       # FragPipe MS database search
│       ├── generate_manifests.py
│       ├── submit_fragpipe.sh
│       └── run_ms_search.sh
└── .claude/
    └── CLAUDE.md                        # Detailed project context
```

## Cluster Data Layout

```
/scratch/leduc.an/AAS_Evo/
├── BAMS/              # GDC BAM files (by UUID subdirectory)
├── RAW/               # PDC RAW files (flattened)
├── VCF/               # Variant calls (per-chunk subdirectories)
│   ├── chunk_00/      # VCFs from first BAM chunk
│   └── chunk_01/ ...
├── VEP/               # VEP annotations (per-chunk subdirectories)
│   ├── chunk_00/      # VEP output from first chunk
│   ├── chunk_01/ ...
│   └── all_missense_mutations.tsv  # Consolidated (top-level)
├── SEQ_FILES/         # Reference files
│   ├── hg38.fa        # Human reference genome (GRCh38)
│   ├── cds.chr.bed    # CDS regions for targeted variant calling
│   ├── uniprot_human_canonical.fasta  # UniProt reviewed proteome
│   ├── uniref50       # MMseqs2 UniRef50 database
│   └── AlphaMissense_hg38.tsv.gz     # AlphaMissense pathogenicity data
├── ANALYSIS/          # Mutation filtering & ranking output
├── FASTA/             # Custom proteogenomics FASTAs
│   ├── per_sample/    # Per-sample mutant entries
│   ├── per_plex/      # Reference + plex-specific mutants + compensatory
│   └── compensatory/  # Predicted compensatory mutation entries
├── MSA/               # Per-gene multiple sequence alignments (A3M)
├── COEVOL/            # Coevolution analysis output
├── MS_SEARCH/         # FragPipe MS search setup & results
│   ├── manifests/     # Per-plex .fp-manifest files
│   ├── annotations/   # Per-plex TMT channel annotations
│   └── results/       # Per-plex FragPipe output
└── logs/              # SLURM job logs
```

## Pipeline Overview

```
Download BAMs → Variant Call → VEP (+ AlphaMissense) → Consolidate
    ↓                                                       ↓
    ↓                                           All Missense Mutations
    ↓                                                       ↓
Filter & Rank (top 5000)                      Per-Sample Mutant FASTAs
    ↓                                                       ↓
Generate MSAs (MMseqs2)                                     ↓
    ↓                                                       ↓
Coevolution Analysis                                        ↓
    ↓                                                       ↓
Compensatory FASTAs ─────────────────────────────────────→  ↓
                                                            ↓
                              Combine Per-Plex FASTAs (ref + observed + compensatory)
                                                            ↓
                                              FragPipe MS Search (per plex)
```

## Workflows

### 0. Metadata Setup (One-Time)

Fetches all metadata from GDC and PDC APIs programmatically. No manual portal downloads needed.

```bash
bash scripts/setup/setup_metadata.sh
```

This runs six steps automatically:
1. **GDC manifest** — queries GDC REST API for all CPTAC WXS BAM files (projects from `studies.tsv`)
2. **GDC metadata** — enriches manifest with case/sample info
3. **PDC metadata** — queries PDC GraphQL API for file manifests, biospecimen data, and TMT plex layouts (all 12 tissues)
4. **Consolidate** — merges per-tissue PDC CSVs into unified files
5. **Matching** — cross-references GDC↔PDC samples, generates pruned manifests
6. **Recover unmatched** — finds WXS BAMs for PDC patients missing from the initial GDC manifest, merges them in, and re-runs matching

To add a new tissue/study, add a row to `metadata/studies.tsv` and re-run.

**GDC token**: BAM downloads require dbGaP authorization. Download your token from the [GDC portal](https://portal.gdc.cancer.gov) (username → Download Token) and pass the file path to `submit_download.sh`. The token is NOT needed for this metadata setup step.

### 1. Reference File Setup (One-Time)

Downloads and indexes all external reference files (~100 GB total):

```bash
# As SLURM job (recommended for large downloads):
sbatch scripts/setup/setup_seq_files.sh

# Or skip large files (UniRef90, VEP) for initial testing:
bash scripts/setup/setup_seq_files.sh --skip-large
```

Downloads: hg38 genome, GENCODE CDS regions, UniProt canonical proteome, AlphaMissense data, UniRef90 database, and VEP container + cache.

### 2. Data Download (Chunked)

BAMs are downloaded in chunks of ~400 to stay within storage limits:

```bash
# One-time: split manifest into chunks (400 BAMs each, tissue-only)
bash scripts/download/gdc/setup_chunks.sh 400 path/to/manifest_wxs_bams_tissue.tsv

# Per chunk (repeat for each chunk manifest):
bash scripts/download/gdc/submit_download.sh path/to/chunks/chunk_00.tsv
```

**PDC RAW files** (open access):
```bash
sbatch scripts/download/pdc/submit_download.sh
```

### 3. BAM Processing Pipeline

Processes all BAMs currently in `BAMS/`, outputs to per-chunk subdirectories (`VCF/chunk_00/`, `VEP/chunk_00/`):

```bash
# Step 1: Variant calling (finds BAMs automatically, outputs to VCF/chunk_00/)
bash scripts/proc_bams/run_pipeline.sh variant-call chunk_00

# Step 2: VEP annotation with AlphaMissense (reads VCF/chunk_00/, outputs to VEP/chunk_00/)
bash scripts/proc_bams/run_pipeline.sh vep chunk_00

# Step 3: Delete BAMs, download next chunk, repeat with chunk_01, etc.
rm -rf /scratch/leduc.an/AAS_Evo/BAMS/*
```

VCF/VEP outputs persist across chunks in their subdirectories. Both scripts skip already-processed samples.

**Final output** (`all_missense_mutations.tsv`, 18 columns): sample_id, genomic position, consequence, gene symbol, protein change (HGVSp), amino acid swap, gnomADe_AF, AlphaMissense pathogenicity + class, read depths, VAF.

### 4. Consolidate & Filter Mutations

```bash
# After ALL chunks are processed: merge per-sample VEP TSVs
bash scripts/proc_bams/consolidate_missense.sh

# Rank mutations by composite pathogenicity score (AlphaMissense, gnomAD, recurrence)
python3 scripts/mutation_analysis/filter_and_rank.py \
    --vep-tsv /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv \
    --ref-fasta /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta \
    -o /scratch/leduc.an/AAS_Evo/ANALYSIS
```

Output: `ANALYSIS/top_5000_mutations.tsv`, `ANALYSIS/gene_list_for_msa.txt` (used automatically by downstream steps).

### 5. MSA Generation

```bash
# Submit MSA generation (auto-finds gene list from ANALYSIS/)
NUM_GENES=$(wc -l < /scratch/leduc.an/AAS_Evo/ANALYSIS/gene_list_for_msa.txt)
sbatch --array=1-${NUM_GENES}%10 scripts/mutation_analysis/submit_msa_generation.sh
```

MSA files named by UniProt accession (`P04637.a3m`). Pre-existing MSAs are auto-detected and skipped.

### 6. Coevolution Analysis

Predicts compensatory translation errors using EVcouplings/Direct Coupling Analysis (DCA). Given a destabilizing missense mutation at position i, uses the Potts model coupling tensor J(i,a;j,b) to identify covarying positions and predict which amino acid substitution could compensate.

**Method**: Mean-field DCA computes evolutionary couplings between all position pairs. The coupling parameters directly encode which amino acid combinations are evolutionarily preferred, enabling more accurate compensatory predictions than simple mutual information.

```bash
# After MSA generation completes (auto-finds gene list from ANALYSIS/):
sbatch scripts/mutation_analysis/submit_coevolution.sh
```

Output: `COEVOL/compensatory_predictions.tsv`

### 7. Generate Compensatory FASTAs

```bash
# After coevolution analysis completes:
sbatch scripts/fasta_gen/submit_compensatory_fastas.sh
```

Output: `FASTA/compensatory/all_compensatory.fasta` — tryptic peptides containing both the original destabilizing mutation and the predicted compensatory substitution.

### 8. Proteogenomics FASTA Generation

```bash
# Generate per-sample mutant FASTAs + per-plex search databases
# Automatically includes compensatory entries if FASTA/compensatory/ exists
sbatch scripts/fasta_gen/submit_proteogenomics.sh
```

Creates custom MS search databases per TMT plex: reference proteome + plex-specific mutant tryptic peptides + plex-specific compensatory peptides.

**Tryptic peptide approach**: Instead of adding full mutant proteins (~500 AA each), the pipeline extracts only tryptic peptides containing mutations (~15 AA avg). This minimizes database size and improves FDR statistics.

**Header format**:
```
# Observed mutations (from VEP):
>mut|P04637|TP53|R273H|genetic|C3L-00001|tumor
SVTCTYSPALNKMFCQLAK

# Predicted compensatory mutations (from coevolution):
>comp|P04637|TP53|R273H_G245S|predicted|C3L-00001|tumor
SVTCTYSPALNKMFCQLAK
```

Fields: `type|accession|gene|swap|source|patient|sample_type`

### 9. MS Database Search (FragPipe)

```bash
# Set up per-plex manifests and submit FragPipe searches
bash scripts/ms_search/run_ms_search.sh /path/to/template.workflow

# Or generate manifests first, then configure workflow separately:
bash scripts/ms_search/run_ms_search.sh
# -> Place workflow at MS_SEARCH/fragpipe.workflow, then:
NUM_PLEXES=$(wc -l < /scratch/leduc.an/AAS_Evo/MS_SEARCH/plex_list.txt)
sbatch --array=1-${NUM_PLEXES}%5 scripts/ms_search/submit_fragpipe.sh
```

Each plex is searched independently (different custom FASTA per plex). All fractions from a plex are grouped as one experiment. TMT channel annotations map channels to patient IDs.

## Reference Files

All reference files are downloaded and indexed by `scripts/setup/setup_seq_files.sh`. They live in `/scratch/leduc.an/AAS_Evo/SEQ_FILES/`:

| File | Source | Size |
|------|--------|------|
| `hg38.fa` + `.fai` | UCSC Genome Browser (chr-prefix, matching GDC BAMs) | ~3 GB |
| `cds.chr.bed` | GENCODE v46 CDS regions (merged, standard chromosomes) | ~2 MB |
| `uniprot_human_canonical.fasta` | UniProt reference proteome UP000005640 (reviewed, canonical) | ~25 MB |
| `AlphaMissense_hg38.tsv.gz` + `.tbi` | DeepMind AlphaMissense pathogenicity predictions | ~6 GB |
| `uniref50` (MMseqs2 db) | UniProt Reference Clusters at 50% identity | ~60 GB |
| VEP container + cache | Ensembl VEP Apptainer image | ~15 GB |

## Requirements

- Python 3.8+, numpy
- `gdc-client` (GDC Data Transfer Tool)
- `samtools`, `bcftools` (variant calling)
- Ensembl VEP via Apptainer (annotation + AlphaMissense plugin)
- `bedtools` (for generating CDS BED, if needed)
- `mmseqs2` (MSA generation)
- `htslib` / `tabix` (AlphaMissense indexing)
- FragPipe (MS database search)

## Data Sources

| Source | Data Type | Access | Files |
|--------|-----------|--------|-------|
| [GDC](https://portal.gdc.cancer.gov) | WXS BAM files | Controlled (dbGaP) | ~2,098 matched |
| [PDC](https://pdc.cancer.gov) | TMT RAW files | Open (signed URLs) | ~5,020 matched |

## License

MIT License
