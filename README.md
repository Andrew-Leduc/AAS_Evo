# AAS_Evo

Multi-omics data pipeline for CPTAC3 (Clinical Proteomic Tumor Analysis Consortium). Downloads and processes matched genomics (WXS BAM files from GDC) and proteomics (TMT-labeled RAW files from PDC) data for integrated analysis of amino acid substitutions.

## Repository Structure

```
AAS_Evo/
├── config/
│   └── paths.py                      # Environment-aware path config
├── scripts/
│   ├── download/
│   │   ├── gdc/
│   │   │   ├── submit_download.sh    # Split manifest & submit SLURM jobs
│   │   │   ├── download_chunk.sh     # Single-chunk SLURM download job
│   │   │   ├── download.py           # GDC download via gdc-client
│   │   │   ├── fetch_metadata.py     # Fetch sample metadata from GDC API
│   │   │   ├── filter_wxs_manifest.py
│   │   │   └── setup_chunks.sh       # Split manifest into 500-BAM chunks
│   │   ├── pdc/
│   │   │   ├── submit_download.sh    # SLURM job wrapper
│   │   │   ├── download.py           # PDC download (rate-limited)
│   │   │   └── consolidate_metadata.py
│   │   └── mapping_report.py         # GDC-PDC sample matching
│   ├── proc_bams/
│   │   ├── run_pipeline.sh           # Wrapper: auto-generates file lists + submits jobs
│   │   ├── submit_variant_call.sh    # SLURM array: variant calling
│   │   ├── submit_vep.sh             # SLURM array: VEP annotation
│   │   └── consolidate_missense.sh   # Merge missense mutations
│   └── proteogenomics/
│       ├── generate_mutant_fastas.py # Per-sample mutant FASTAs from VEP
│       ├── combine_plex_fastas.py    # Combine by TMT plex
│       └── submit_proteogenomics.sh  # SLURM wrapper
└── .claude/
    └── CLAUDE.md                     # Detailed project context
```

## Cluster Data Layout

```
/scratch/leduc.an/AAS_Evo/
├── BAMS/              # GDC BAM files (by UUID subdirectory)
├── RAW/               # PDC RAW files
├── VCF/               # Variant calls from BAM processing
├── VEP/               # VEP annotations + missense tables
├── SEQ_FILES/         # Reference files
│   ├── hg38.fa        # Human reference genome (GRCh38)
│   ├── hg38.fa.fai    # Reference index
│   ├── cds.chr.bed    # CDS regions for targeted variant calling
│   └── uniprot_human_canonical.fasta  # UniProt reviewed proteome
├── FASTA/             # Custom proteogenomics FASTAs
│   ├── per_sample/    # Per-sample mutant entries
│   └── per_plex/      # Reference + plex-specific mutants
└── logs/              # SLURM job logs
```

## Workflows

### 1. Data Download (Chunked)

BAMs are downloaded in chunks of ~500 to stay within storage limits:

```bash
# One-time: split manifest into chunks
bash scripts/download/gdc/setup_chunks.sh

# Per chunk (repeat for each chunk manifest):
bash scripts/download/gdc/submit_download.sh path/to/chunk_00.tsv
```

**PDC RAW files** (open access):
```bash
sbatch scripts/download/pdc/submit_download.sh
```

### 2. BAM Processing Pipeline

Processes all BAMs currently in `BAMS/`, outputs to `VCF/` and `VEP/`:

```bash
# Step 1: Variant calling (finds BAMs automatically)
bash scripts/proc_bams/run_pipeline.sh variant-call

# Step 2: VEP annotation (finds VCFs automatically)
bash scripts/proc_bams/run_pipeline.sh vep

# Step 3: Delete BAMs, download next chunk, repeat
rm -rf /scratch/leduc.an/AAS_Evo/BAMS/*

# After ALL chunks: consolidate missense mutations
bash scripts/proc_bams/consolidate_missense.sh
```

VCF/VEP outputs persist across chunks. Both scripts skip already-processed samples.

**Final output** (`all_missense_mutations.tsv`): sample_id, genomic position, gene symbol, protein change (HGVSp), amino acid swap, gnomAD population frequency, variant allele frequency.

### 3. Proteogenomics FASTA Generation

```bash
# One-time: download UniProt reference proteome
wget -O /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29"

# Generate per-sample mutant FASTAs + per-TMT-plex search databases
sbatch scripts/proteogenomics/submit_proteogenomics.sh
```

Creates custom MS search databases: reference proteome + sample-specific missense mutations, combined per TMT plex. Links VEP output → GDC UUID → TMT plex via sample metadata.

## Reference Files

Reference files stored in `/scratch/leduc.an/AAS_Evo/SEQ_FILES/`. To reproduce:

**Human reference genome (hg38)**:
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

**CDS regions (GENCODE)**:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
zcat gencode.v46.annotation.gtf.gz \
    | awk '$3 == "CDS" {print $1"\t"$4-1"\t"$5}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge > cds.chr.bed
```

## Requirements

- Python 3.8+
- `gdc-client` (GDC Data Transfer Tool)
- `samtools`, `bcftools` (variant calling)
- Ensembl VEP via Apptainer (annotation)
- `bedtools` (for generating CDS BED, if needed)

## Data Sources

| Source | Data Type | Access | Files |
|--------|-----------|--------|-------|
| [GDC](https://portal.gdc.cancer.gov) | WXS BAM files | Controlled (dbGaP) | ~2,098 matched |
| [PDC](https://pdc.cancer.gov) | TMT RAW files | Open (signed URLs) | ~5,020 matched |

## License

MIT License
