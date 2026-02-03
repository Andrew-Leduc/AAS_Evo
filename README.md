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
│   ├── fasta_gen/                    # Custom proteogenomics FASTAs
│   │   ├── generate_mutant_fastas.py # Per-sample mutant FASTAs from VEP
│   │   ├── combine_plex_fastas.py    # Combine by TMT plex
│   │   └── submit_proteogenomics.sh  # SLURM wrapper for FASTA generation
│   └── mutation_analysis/            # MSA generation & coevolution
│       ├── generate_msas.py          # MSA generation via MMseqs2
│       ├── submit_msa_generation.sh  # SLURM array: one gene per task
│       ├── coevolution_analysis.py   # MI+APC covariation & compensatory prediction
│       └── submit_coevolution.sh     # SLURM wrapper for coevolution analysis
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
├── MSA/               # Per-gene multiple sequence alignments (A3M)
├── COEVOL/            # Coevolution analysis output
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
sbatch scripts/fasta_gen/submit_proteogenomics.sh
```

Creates custom MS search databases: reference proteome + sample-specific missense mutations, combined per TMT plex. Links VEP output → GDC UUID → TMT plex via sample metadata.

### 4. MSA Generation & Coevolution Analysis

Predicts compensatory translation errors: given a destabilizing missense mutation, finds covarying positions via MI+APC and predicts which amino acid substitution could compensate.

```bash
# One-time: download and index UniRef90 (~28 GB compressed)
wget -O uniref90.fasta.gz \
    https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
mmseqs createdb uniref90.fasta.gz /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniref90

# Step 1: Generate gene list (genes with mutations that need MSAs)
python3 scripts/mutation_analysis/generate_msas.py --make-gene-list \
    --vep-tsv /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv \
    --msa-dir /scratch/leduc.an/AAS_Evo/MSA \
    --ref-fasta /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta \
    -o /scratch/leduc.an/AAS_Evo/gene_list.txt

# Step 2: Generate MSAs (SLURM array, one gene per task)
NUM_GENES=$(wc -l < /scratch/leduc.an/AAS_Evo/gene_list.txt)
sbatch --array=1-${NUM_GENES}%10 scripts/mutation_analysis/submit_msa_generation.sh

# Step 3: Run coevolution analysis (after MSA jobs complete)
sbatch scripts/mutation_analysis/submit_coevolution.sh
```

MSA files are named by UniProt accession (`P04637.a3m`). The gene list uses gene symbols from VEP; the scripts handle the mapping. Pre-existing MSAs (named by gene, accession, or entry name) are auto-detected and skipped.

## Reference Files

Reference files stored in `/scratch/leduc.an/AAS_Evo/SEQ_FILES/`. To reproduce:

**Human reference genome (hg38)**:
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

**CDS regions (GENCODE)** — GENCODE uses `chr` prefixes matching GDC BAMs:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz

# Extract CDS, filter to standard chromosomes, sort, merge overlapping intervals
zcat gencode.v46.annotation.gtf.gz \
    | awk '$3 == "CDS" {print $1"\t"$4-1"\t"$5}' \
    | grep -E '^chr([0-9]+|X|Y|M)\b' \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > cds.chr.bed
```
The grep removes alt/patch contigs not in GDC BAMs. The merge is required — without it, overlapping CDS intervals cause bcftools mpileup to emit unsorted positions.

## Requirements

- Python 3.8+, numpy
- `gdc-client` (GDC Data Transfer Tool)
- `samtools`, `bcftools` (variant calling)
- Ensembl VEP via Apptainer (annotation)
- `bedtools` (for generating CDS BED, if needed)
- `mmseqs2` (MSA generation)

## Data Sources

| Source | Data Type | Access | Files |
|--------|-----------|--------|-------|
| [GDC](https://portal.gdc.cancer.gov) | WXS BAM files | Controlled (dbGaP) | ~2,098 matched |
| [PDC](https://pdc.cancer.gov) | TMT RAW files | Open (signed URLs) | ~5,020 matched |

## License

MIT License
