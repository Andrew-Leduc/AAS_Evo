#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --job-name=setup_seq
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/setup_seq_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/setup_seq_%j.err

#
# Download and index all external reference files needed for the pipeline.
#
# This script is idempotent: it skips files that already exist.
# Can be run interactively or submitted as a SLURM job.
#
# What it downloads:
#   1. Human reference genome (hg38.fa) from UCSC           (~3 GB)
#   2. GENCODE CDS regions → cds.chr.bed                    (~2 MB)
#   3. UniProt human reference proteome (canonical)          (~25 MB)
#   4. AlphaMissense pathogenicity data + tabix index        (~6 GB)
#   5. UniRef30 ColabFold pre-built MMseqs2 database           (~68 GB + ~110 GB)
#   6. VEP Apptainer container + cache                       (~15 GB)
#
# Usage:
#   # Interactive (may be slow for large downloads):
#   bash scripts/setup/setup_seq_files.sh
#
#   # As SLURM job (recommended):
#   sbatch scripts/setup/setup_seq_files.sh
#
#   # Skip large downloads (UniRef90, VEP cache):
#   bash scripts/setup/setup_seq_files.sh --skip-large
#

set -euo pipefail

# --------- PATHS ----------
DATA_DIR="/scratch/leduc.an/AAS_Evo"
SEQ_DIR="${DATA_DIR}/SEQ_FILES"
TOOLS_DIR="/scratch/leduc.an/tools"
VEP_DIR="${TOOLS_DIR}/vep"
# ---------------------------

SKIP_LARGE=false
if [[ "${1:-}" == "--skip-large" ]]; then
    SKIP_LARGE=true
fi

mkdir -p "$SEQ_DIR"
mkdir -p "$TOOLS_DIR"
mkdir -p "${DATA_DIR}/logs"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Load modules if available (cluster)
module load samtools 2>/dev/null || true
module load mmseqs2 2>/dev/null || true

echo "========================================="
echo "Pipeline Reference File Setup"
echo "========================================="
echo "Target directory: $SEQ_DIR"
echo "Skip large files: $SKIP_LARGE"
echo ""

# ============================================================
# 1. Human Reference Genome (hg38.fa)
# ============================================================
# Source: UCSC Genome Browser. Uses chr-prefix matching GDC BAMs.
if [[ -f "${SEQ_DIR}/hg38.fa" && -f "${SEQ_DIR}/hg38.fa.fai" ]]; then
    log "SKIP hg38.fa (already exists with index)"
else
    log "Downloading hg38.fa from UCSC..."
    wget -q -O "${SEQ_DIR}/hg38.fa.gz" \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

    log "Decompressing hg38.fa.gz..."
    gunzip -f "${SEQ_DIR}/hg38.fa.gz"

    log "Indexing hg38.fa with samtools..."
    samtools faidx "${SEQ_DIR}/hg38.fa"

    log "DONE hg38.fa"
fi
echo ""

# ============================================================
# 2. CDS Regions BED (cds.chr.bed)
# ============================================================
# Source: GENCODE v46. CDS regions for targeted variant calling.
# Uses chr-prefix matching UCSC hg38 and GDC BAMs.
if [[ -f "${SEQ_DIR}/cds.chr.bed" ]]; then
    log "SKIP cds.chr.bed (already exists)"
else
    log "Downloading GENCODE v46 annotation..."
    GTF_GZ="${SEQ_DIR}/gencode.v46.annotation.gtf.gz"
    wget -q -O "$GTF_GZ" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"

    log "Extracting CDS regions..."
    # Extract CDS, filter to standard chromosomes, sort, merge overlapping
    zcat "$GTF_GZ" \
        | awk '$3 == "CDS" {print $1"\t"$4-1"\t"$5}' \
        | grep -E '^chr([0-9]+|X|Y|M)\b' \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > "${SEQ_DIR}/cds.chr.bed"

    rm -f "$GTF_GZ"
    REGIONS=$(wc -l < "${SEQ_DIR}/cds.chr.bed")
    log "DONE cds.chr.bed ($REGIONS regions)"
fi
echo ""

# ============================================================
# 3. UniProt Human Reference Proteome (canonical only)
# ============================================================
# Source: UniProt reference proteome UP000005640 (reviewed/Swiss-Prot).
# One canonical protein per gene. GN= field maps to VEP SYMBOL.
if [[ -f "${SEQ_DIR}/uniprot_human_canonical.fasta" ]]; then
    log "SKIP uniprot_human_canonical.fasta (already exists)"
else
    log "Downloading UniProt human reference proteome (reviewed, canonical)..."
    wget -q -O "${SEQ_DIR}/uniprot_human_canonical.fasta" \
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29"

    NPROTEINS=$(grep -c "^>" "${SEQ_DIR}/uniprot_human_canonical.fasta")
    log "DONE uniprot_human_canonical.fasta ($NPROTEINS proteins)"
fi
echo ""

# ============================================================
# 4. AlphaMissense Pathogenicity Data
# ============================================================
# Source: DeepMind AlphaMissense. Used by VEP plugin.
# Requires tabix (htslib) for indexing.
if [[ -f "${SEQ_DIR}/AlphaMissense_hg38.tsv.gz" && -f "${SEQ_DIR}/AlphaMissense_hg38.tsv.gz.tbi" ]]; then
    log "SKIP AlphaMissense (already exists with index)"
else
    log "Downloading AlphaMissense_hg38.tsv.gz (~6 GB)..."
    wget -q -O "${SEQ_DIR}/AlphaMissense_hg38.tsv.gz" \
        "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"

    log "Indexing with tabix..."
    # Try system tabix first, then local build
    TABIX_BIN="tabix"
    if ! command -v tabix &>/dev/null; then
        if [[ -f "/scratch/leduc.an/tools/htslib-1.21/tabix" ]]; then
            TABIX_BIN="/scratch/leduc.an/tools/htslib-1.21/tabix"
        else
            log "ERROR: tabix not found. Install htslib or build from source:"
            log "  wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2"
            log "  tar xjf htslib-1.21.tar.bz2 && cd htslib-1.21 && make"
            exit 1
        fi
    fi
    "$TABIX_BIN" -s 1 -b 2 -e 2 -S 1 "${SEQ_DIR}/AlphaMissense_hg38.tsv.gz"

    log "DONE AlphaMissense"
fi
echo ""

# ============================================================
# 5. UniRef30 Database (for MSA generation)
# ============================================================
# Source: ColabFold pre-built UniRef30 MMseqs2 database (Mirdita et al., 2022).
# Pre-indexed — no mmseqs createdb step needed.
# Download: ~68 GB compressed. Extracted: ~110 GB.
UNIREF30_DIR="${SEQ_DIR}/uniref30_2302"
if [[ "$SKIP_LARGE" == true ]]; then
    log "SKIP UniRef30 (--skip-large)"
elif [[ -f "${UNIREF30_DIR}/uniref30_2302.dbtype" || -f "${UNIREF30_DIR}/uniref30_2302" ]]; then
    log "SKIP UniRef30 MMseqs2 database (already exists)"
else
    UNIREF30_GZ="${SEQ_DIR}/uniref30_2302.tar.gz"

    if [[ ! -f "$UNIREF30_GZ" ]]; then
        log "Downloading uniref30_2302.tar.gz (~68 GB)..."
        wget --continue -O "$UNIREF30_GZ" \
            "https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2302.tar.gz"
    else
        log "uniref30_2302.tar.gz already downloaded, extracting..."
    fi

    log "Extracting UniRef30 database (this may take a while)..."
    tar -xzf "$UNIREF30_GZ" -C "$SEQ_DIR"

    log "DONE UniRef30 MMseqs2 database"
fi
echo ""

# ============================================================
# 6. VEP Apptainer Container + Cache
# ============================================================
# Source: Ensembl VEP Docker image, converted to Apptainer/Singularity.
if [[ "$SKIP_LARGE" == true ]]; then
    log "SKIP VEP container + cache (--skip-large)"
elif [[ -f "${VEP_DIR}/ensembl-vep.sif" ]]; then
    log "SKIP VEP container (already exists at ${VEP_DIR}/ensembl-vep.sif)"

    # Still check cache
    if [[ -d "${SEQ_DIR}/vep_cache/homo_sapiens" ]]; then
        log "SKIP VEP cache (already exists)"
    else
        log "Installing VEP cache (offline mode)..."
        mkdir -p "${SEQ_DIR}/vep_cache"
        apptainer exec "${VEP_DIR}/ensembl-vep.sif" \
            vep_install --AUTO c --ASSEMBLY GRCh38 --SPECIES homo_sapiens \
            --CACHEDIR "${SEQ_DIR}/vep_cache"
        log "DONE VEP cache"
    fi
else
    log "Pulling VEP Apptainer image..."
    mkdir -p "$VEP_DIR"
    apptainer pull "${VEP_DIR}/ensembl-vep.sif" \
        docker://ensemblorg/ensembl-vep

    log "Installing VEP cache (offline mode)..."
    mkdir -p "${SEQ_DIR}/vep_cache"
    apptainer exec "${VEP_DIR}/ensembl-vep.sif" \
        vep_install --AUTO c --ASSEMBLY GRCh38 --SPECIES homo_sapiens \
        --CACHEDIR "${SEQ_DIR}/vep_cache"

    log "DONE VEP container + cache"
fi
echo ""

# ============================================================
# Summary
# ============================================================
echo "========================================="
echo "Setup Summary"
echo "========================================="

check_file() {
    local path="$1"
    local label="$2"
    if [[ -f "$path" ]]; then
        local size
        size=$(du -sh "$path" 2>/dev/null | cut -f1)
        echo "  [OK]   $label ($size)"
    elif [[ -d "$path" ]]; then
        echo "  [OK]   $label (directory)"
    else
        echo "  [MISS] $label"
    fi
}

check_file "${SEQ_DIR}/hg38.fa"                       "hg38.fa"
check_file "${SEQ_DIR}/hg38.fa.fai"                   "hg38.fa.fai"
check_file "${SEQ_DIR}/cds.chr.bed"                    "cds.chr.bed"
check_file "${SEQ_DIR}/uniprot_human_canonical.fasta"  "uniprot_human_canonical.fasta"
check_file "${SEQ_DIR}/AlphaMissense_hg38.tsv.gz"      "AlphaMissense_hg38.tsv.gz"
check_file "${SEQ_DIR}/AlphaMissense_hg38.tsv.gz.tbi"  "AlphaMissense_hg38.tsv.gz.tbi"
check_file "${UNIREF30_DIR}/uniref30_2302.dbtype"       "UniRef30 MMseqs2 database"
check_file "${VEP_DIR}/ensembl-vep.sif"                "VEP Apptainer container"
check_file "${SEQ_DIR}/vep_cache/homo_sapiens"         "VEP cache"
echo ""
echo "All reference files are stored in: $SEQ_DIR"
