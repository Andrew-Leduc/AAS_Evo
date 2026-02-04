#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --job-name=coevol
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/coevolution_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/coevolution_%j.err

#
# Coevolutionary analysis: predict compensatory translation errors.
#
# For each gene with both missense mutations and an MSA, computes MI+APC
# coevolutionary couplings and predicts which amino acid substitution at a
# covarying position could compensate for a destabilizing mutation.
#
# Prerequisites:
#   - MSAs generated (run submit_msa_generation.sh first)
#   - Consolidated missense table: VEP/all_missense_mutations.tsv
#   - Reference proteome: SEQ_FILES/uniprot_human_canonical.fasta
#   - MSA directory with per-gene alignment files (A3M/FASTA/Stockholm)
#
# Usage:
#   sbatch submit_coevolution.sh
#
# To run on specific genes only:
#   sbatch submit_coevolution.sh TP53 BRCA1
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="/scratch/leduc.an/AAS_Evo"

REF_FASTA="${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta"
VEP_TSV="${DATA_DIR}/VEP/all_missense_mutations.tsv"
MSA_DIR="${DATA_DIR}/MSA"
# Gene list: prefer filter_and_rank.py output, fall back to legacy location
if [[ -f "${DATA_DIR}/ANALYSIS/gene_list_for_msa.txt" ]]; then
    GENE_LIST="${DATA_DIR}/ANALYSIS/gene_list_for_msa.txt"
else
    GENE_LIST="${DATA_DIR}/gene_list.txt"
fi
OUT_DIR="${DATA_DIR}/COEVOL"
# --------------------------

module load python/3.8.1

mkdir -p "$OUT_DIR"
mkdir -p "${DATA_DIR}/logs"

# Validate inputs
for f in "$REF_FASTA" "$VEP_TSV"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

if [[ ! -d "$MSA_DIR" ]]; then
    echo "ERROR: MSA directory not found: $MSA_DIR"
    exit 1
fi

echo "========================================="
echo "Coevolutionary Compensatory Analysis"
echo "========================================="
echo ""

# Build optional flags
EXTRA_FLAGS=""

# Use gene list if it exists (from MSA generation step)
if [[ -f "$GENE_LIST" ]]; then
    EXTRA_FLAGS="--gene-list $GENE_LIST"
    echo "Using gene list: $GENE_LIST"
fi

# Override with specific genes if provided as positional arguments
if [[ $# -gt 0 ]]; then
    EXTRA_FLAGS="--genes $*"
    echo "Limiting to genes: $*"
fi

echo ""
echo "[$(date)] Starting coevolution analysis..."
python3 "${SCRIPTS_DIR}/coevolution_analysis.py" \
    --msa-dir "$MSA_DIR" \
    --vep-tsv "$VEP_TSV" \
    --ref-fasta "$REF_FASTA" \
    -o "${OUT_DIR}/compensatory_predictions.tsv" \
    --min-neff 50 \
    --top-k-positions 10 \
    --top-k-substitutions 3 \
    $EXTRA_FLAGS

echo ""
echo "[$(date)] Done."
echo "Output: ${OUT_DIR}/"
