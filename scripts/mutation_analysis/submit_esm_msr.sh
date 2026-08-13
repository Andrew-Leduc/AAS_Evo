#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --job-name=esm_msr
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/esm_msr_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/esm_msr_%j.err

# ESM-MSR double-mutant / epistasis scoring (github.com/SKTeamLab/esm-msr).
# Usage:
#   export HF_TOKEN=hf_...            # HF read token for esm3-sm-open-v1
#   sbatch submit_esm_msr.sh                       # 1A0F smoke test
#   sbatch submit_esm_msr.sh <input.csv> <out.csv> # real run on our pairs
# The venv/repo live off the AAS repo (installed separately); paths below.

set -euo pipefail

ESM_DIR="${ESM_DIR:-/projects/slavov/andrew/esm/esm-msr}"
VENV="${VENV:-/scratch/leduc.an/msr_venv}"
CKPT="${CKPT:-LoRA_models/esm-msr-small/epoch=03-val_rho_combined_avg=0.816.ckpt}"
OUTDIR=/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap
mkdir -p /scratch/leduc.an/AAS_Evo/logs "$OUTDIR"

module load anaconda3/2024.06 cuda/13.2.0
source "$VENV/bin/activate"
export HF_HOME="${HF_HOME:-/projects/slavov/andrew/esm/hf_cache}"
export TMPDIR="${TMPDIR:-/scratch/leduc.an/tmp}"
: "${HF_TOKEN:?set HF_TOKEN (a HuggingFace read token) before sbatch}"

cd "$ESM_DIR"
python -c "import torch; print('cuda available:', torch.cuda.is_available())"

INPUT_CSV="${1:-}"
OUTPUT_CSV="${2:-$OUTDIR/esm_msr_scores_all.csv}"

if [[ -z "$INPUT_CSV" ]]; then
    echo "[smoke test] scoring 1A0F single mutants"
    python src/esm_msr/inference.py \
        --checkpoint_path "$CKPT" \
        --code 1A0F --chain A --mode singles --skip_reverse \
        --output_csv "$OUTDIR/esm_msr_smoke.csv"
    head "$OUTDIR/esm_msr_smoke.csv"
else
    echo "[run] scoring $INPUT_CSV -> $OUTPUT_CSV"
    python src/esm_msr/inference.py \
        --checkpoint_path "$CKPT" \
        --input_csv "$INPUT_CSV" --mode singles+doubles \
        --output_csv "$OUTPUT_CSV"
    echo "wrote $OUTPUT_CSV"
fi
