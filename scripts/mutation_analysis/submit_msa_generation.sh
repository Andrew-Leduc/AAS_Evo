#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --partition=short
#SBATCH --job-name=msa_gen
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/msa_gen_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/msa_gen_%A_%a.err

#
# Generate MSAs via MMseqs2 for coevolution analysis.
#
# One SLURM array task per gene. Each task searches the gene's protein
# sequence against UniRef90 and produces an A3M alignment file.
#
# Prerequisites:
#   - Gene list generated:
#       python3 scripts/proteogenomics/generate_msas.py --make-gene-list \
#           --vep-tsv /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv \
#           --msa-dir /scratch/leduc.an/AAS_Evo/MSA \
#           --ref-fasta /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta \
#           -o /scratch/leduc.an/AAS_Evo/gene_list.txt
#
#   - UniRef90 MMseqs2 database in SEQ_FILES/uniref90
#       (build with: mmseqs createdb uniref90.fasta SEQ_FILES/uniref90)
#
# Usage:
#   NUM_GENES=$(wc -l < /scratch/leduc.an/AAS_Evo/gene_list.txt)
#   sbatch --array=1-${NUM_GENES}%10 submit_msa_generation.sh
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="/scratch/leduc.an/AAS_Evo"

GENE_LIST="${DATA_DIR}/gene_list.txt"
REF_FASTA="${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta"
UNIREF90_DB="${DATA_DIR}/SEQ_FILES/uniref90"
MSA_DIR="${DATA_DIR}/MSA"
TMP_DIR="${DATA_DIR}/tmp/msa_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
# --------------------------

module load python/3.8.1
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

# Check UniRef90 database exists (mmseqs dbs have multiple files)
if [[ ! -f "${UNIREF90_DB}.dbtype" ]] && [[ ! -f "${UNIREF90_DB}" ]]; then
    echo "ERROR: UniRef90 database not found: $UNIREF90_DB"
    echo "Build with: mmseqs createdb uniref90.fasta ${UNIREF90_DB}"
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
    --uniref90-db "$UNIREF90_DB" \
    --msa-dir "$MSA_DIR" \
    --tmp-dir "$TMP_DIR" \
    --threads 4

echo ""
echo "[$(date)] Done."
