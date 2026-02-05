#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --job-name=protgen
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/proteogenomics_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/proteogenomics_%j.err

#
# Generate custom protein FASTA search databases for proteogenomics.
#
# Step 1: Per-sample mutant FASTAs from VEP missense output
# Step 2: Per-TMT-plex FASTAs (reference + plex-specific mutants)
#
# Prerequisites:
#   - VEP annotation completed (*.vep.tsv files in VEP/)
#   - Reference proteome downloaded (or run: bash scripts/setup/setup_seq_files.sh)
#
# Usage:
#   sbatch submit_proteogenomics.sh
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
# SLURM copies batch scripts to a spool directory, breaking dirname $0.
# Fall back to SLURM_SUBMIT_DIR when the expected Python script isn't found.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && ! -f "${SCRIPTS_DIR}/generate_mutant_fastas.py" ]]; then
    SCRIPTS_DIR="${SLURM_SUBMIT_DIR}/scripts/fasta_gen"
fi
META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo_meta"
DATA_DIR="/scratch/leduc.an/AAS_Evo"

REF_FASTA="${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta"
VEP_DIR="${DATA_DIR}/VEP"
FASTA_DIR="${DATA_DIR}/FASTA"

TMT_MAP="${META_DIR}/PDC_meta/pdc_file_tmt_map.tsv"
GDC_META="${META_DIR}/GDC_meta/gdc_meta_matched.tsv"
# --------------------------

module load python/3.8.1 2>/dev/null || true

mkdir -p "${FASTA_DIR}/per_sample"
mkdir -p "${FASTA_DIR}/per_plex"
mkdir -p "${DATA_DIR}/logs"

# Validate inputs
for f in "$REF_FASTA" "$TMT_MAP" "$GDC_META"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

if [[ ! -d "$VEP_DIR" ]]; then
    echo "ERROR: VEP directory not found: $VEP_DIR"
    exit 1
fi

echo "========================================="
echo "Proteogenomics FASTA Generation"
echo "========================================="
echo ""

# Step 1: Generate per-sample mutant FASTAs
echo "[$(date)] Step 1: Generating per-sample mutant FASTAs..."
python3 "${SCRIPTS_DIR}/generate_mutant_fastas.py" \
    --ref-fasta "$REF_FASTA" \
    --vep-dir "$VEP_DIR" \
    -o "${FASTA_DIR}/per_sample"

echo ""

# Step 2: Combine by TMT plex (with compensatory mutations if available)
COMP_FASTA="${FASTA_DIR}/compensatory/all_compensatory.fasta"
echo "[$(date)] Step 2: Combining per-plex FASTAs..."
COMP_FLAG=""
if [[ -f "$COMP_FASTA" ]]; then
    COMP_FLAG="--compensatory-fasta $COMP_FASTA"
    echo "  Including compensatory mutations from: $COMP_FASTA"
fi
python3 "${SCRIPTS_DIR}/combine_plex_fastas.py" \
    --ref-fasta "$REF_FASTA" \
    --sample-dir "${FASTA_DIR}/per_sample" \
    --tmt-map "$TMT_MAP" \
    --gdc-meta "$GDC_META" \
    $COMP_FLAG \
    -o "${FASTA_DIR}/per_plex"

echo ""
echo "[$(date)] Done."
echo "Per-sample FASTAs: ${FASTA_DIR}/per_sample/"
echo "Per-plex FASTAs:   ${FASTA_DIR}/per_plex/"
