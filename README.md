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
│   │   │   └── filter_wxs_manifest.py
│   │   ├── pdc/
│   │   │   ├── submit_download.sh    # SLURM job wrapper
│   │   │   ├── download.py           # PDC download (rate-limited)
│   │   │   └── consolidate_metadata.py
│   │   └── mapping_report.py         # GDC-PDC sample matching
│   └── proc_bams/
│       ├── submit_variant_call.sh    # SLURM array: variant calling
│       ├── submit_vep.sh             # SLURM array: VEP annotation
│       └── consolidate_missense.sh   # Merge missense mutations
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
│   └── cds.chr.bed    # CDS regions for targeted variant calling
└── logs/              # SLURM job logs
```

## Workflows

### 1. Data Download

**GDC BAM files** (controlled access, requires dbGaP approval):
```bash
# On login node - splits manifest into 20 chunks, submits each as a SLURM job
bash scripts/download/gdc/submit_download.sh
```

**PDC RAW files** (open access):
```bash
sbatch scripts/download/pdc/submit_download.sh
```

### 2. BAM Processing Pipeline

```bash
# Step 1: Variant calling (CDS regions only)
find /scratch/leduc.an/AAS_Evo/BAMS -name "*.bam" > /scratch/leduc.an/AAS_Evo/bam_list.txt
NUM_BAMS=$(wc -l < /scratch/leduc.an/AAS_Evo/bam_list.txt)
sbatch --array=1-${NUM_BAMS}%10 scripts/proc_bams/submit_variant_call.sh

# Step 2: VEP annotation (with gnomAD frequencies)
ls /scratch/leduc.an/AAS_Evo/VCF/*.vcf.gz > /scratch/leduc.an/AAS_Evo/vcf_list.txt
NUM_VCFS=$(wc -l < /scratch/leduc.an/AAS_Evo/vcf_list.txt)
sbatch --array=1-${NUM_VCFS}%10 scripts/proc_bams/submit_vep.sh

# Step 3: Consolidate all missense mutations
bash scripts/proc_bams/consolidate_missense.sh
```

**Final output** (`all_missense_mutations.tsv`): sample_id, genomic position, gene symbol, protein change (HGVSp), amino acid swap, gnomAD population frequency, variant allele frequency.

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
