#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --job-name=msa_gen
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/msa_gen_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/msa_gen_%A_%a.err

#
# Generate MSAs via MMseqs2 for coevolution analysis.
#
# One SLURM array task per gene. Each task searches the gene's protein
# sequence against UniRef30 and produces an A3M alignment file.
#
# Prerequisites:
#   - Gene list generated via filter_and_rank.py (writes ANALYSIS/gene_list_for_msa.txt)
#     or via generate_msas.py --make-gene-list (writes gene_list.txt)
#
#   - UniRef30 ColabFold MMseqs2 database in SEQ_FILES/uniref30_2302/
#       (download with: sbatch scripts/setup/download_uniref30.sh)
#
# Usage:
#   NUM_GENES=$(wc -l < /scratch/leduc.an/AAS_Evo/ANALYSIS/gene_list_for_msa.txt)
#   sbatch --array=1-${NUM_GENES}%10 submit_msa_generation.sh
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
# SLURM copies batch scripts to a spool directory, breaking dirname $0.
# Fall back to SLURM_SUBMIT_DIR when the expected Python script isn't found.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && ! -f "${SCRIPTS_DIR}/generate_msas.py" ]]; then
    SCRIPTS_DIR="${SLURM_SUBMIT_DIR}/scripts/mutation_analysis"
fi
DATA_DIR="/scratch/leduc.an/AAS_Evo"

# Gene list: prefer filter_and_rank.py output, fall back to legacy location
if [[ -f "${DATA_DIR}/ANALYSIS/gene_list_for_msa.txt" ]]; then
    GENE_LIST="${DATA_DIR}/ANALYSIS/gene_list_for_msa.txt"
else
    GENE_LIST="${DATA_DIR}/gene_list.txt"
fi
REF_FASTA="${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta"
TARGET_DB="${DATA_DIR}/SEQ_FILES/uniref30_2302/uniref30_2302"
MSA_DIR="${DATA_DIR}/MSA"
TMP_DIR="${DATA_DIR}/tmp/msa_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
# --------------------------

module load python/3.8.1 2>/dev/null || true
module load mmseqs2 2>/dev/null || true
export PATH="$HOME/bin/mmseqs/bin:$PATH"

mkdir -p "$MSA_DIR"
mkdir -p "${DATA_DIR}/logs"

# Validate inputs
if [[ ! -f "$GENE_LIST" ]]; then
    echo "ERROR: Gene list not found: $GENE_LIST"
    echo "Generate it first with: python3 generate_msas.py --make-gene-list ..."
    exit 1
fi

if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA"
    exit 1
fi

# Check target database exists (mmseqs dbs have multiple files)
if [[ ! -f "${TARGET_DB}.dbtype" ]] && [[ ! -f "${TARGET_DB}" ]]; then
    echo "ERROR: Target database not found: $TARGET_DB"
    echo "Download with: sbatch scripts/setup/download_uniref30.sh"
    exit 1
fi

echo "========================================="
echo "MSA Generation - Task ${SLURM_ARRAY_TASK_ID}"
echo "========================================="
echo ""

python3 "${SCRIPTS_DIR}/generate_msas.py" \
    --gene-list "$GENE_LIST" \
    --index "${SLURM_ARRAY_TASK_ID}" \
    --ref-fasta "$REF_FASTA" \
    --target-db "$TARGET_DB" \
    --msa-dir "$MSA_DIR" \
    --tmp-dir "$TMP_DIR" \
    --threads 4

echo ""
echo "[$(date)] Done."
