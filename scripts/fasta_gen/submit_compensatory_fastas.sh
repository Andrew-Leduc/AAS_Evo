#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --job-name=comp_fasta
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/comp_fasta_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/comp_fasta_%j.err

#
# Generate compensatory mutation FASTAs from coevolution predictions
#
# Usage:
#   sbatch submit_compensatory_fastas.sh
#
# Prerequisites:
#   1. Coevolution analysis completed (compensatory_predictions.tsv exists)
#   2. filter_and_rank.py completed (top_5000_mutations.tsv exists)
#

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
# SLURM copies batch scripts to a spool directory, breaking dirname $0.
# Fall back to SLURM_SUBMIT_DIR when the expected Python script isn't found.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && ! -f "${SCRIPTS_DIR}/generate_compensatory_fastas.py" ]]; then
    SCRIPTS_DIR="${SLURM_SUBMIT_DIR}/scripts/fasta_gen"
fi
DATA_DIR="/scratch/leduc.an/AAS_Evo"

mkdir -p /scratch/leduc.an/AAS_Evo/logs

python3 "${SCRIPTS_DIR}/generate_compensatory_fastas.py" \
    --predictions "${DATA_DIR}/COEVOL/compensatory_predictions.tsv" \
    --ref-fasta "${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta" \
    --ranked-mutations "${DATA_DIR}/ANALYSIS/top_5000_mutations.tsv" \
    -o "${DATA_DIR}/FASTA/compensatory"
