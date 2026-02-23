#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --job-name=spurs_ddg
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/spurs_ddg_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/spurs_ddg_%A_%a.err

set -euo pipefail

DATA_DIR="/scratch/leduc.an/AAS_Evo"
REF_FASTA="${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta"
SPURS_OUT="${DATA_DIR}/SPURS"

# This is the input you said the next step will use:
GENE_LIST="${DATA_DIR}/gene_list.txt"

mkdir -p "${DATA_DIR}/logs" "${SPURS_OUT}"

module load miniconda3 2>/dev/null || true
module load anaconda3 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${DATA_DIR}/conda_envs/spurs_env"

python3 "$(cd "$(dirname "$0")" && pwd)/spurs_predict_per_gene.py" \
  --gene-list "${GENE_LIST}" \
  --gene-index "${SLURM_ARRAY_TASK_ID}" \
  --ref-fasta "${REF_FASTA}" \
  --out-dir "${SPURS_OUT}"
