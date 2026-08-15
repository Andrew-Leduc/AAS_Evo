#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --job-name=esm_msr_arr
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/esm_msr_arr_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/esm_msr_arr_%A_%a.err

# ESM-MSR scoring, ONE structure per array task (inference.py requires it).
# Prereq: python3 esm_msr_shard.py split <input.csv> <SHARD_DIR>
# Submit:  NUM=$(wc -l < $SHARD_DIR/list.txt)
#          sbatch --array=1-${NUM}%20 submit_esm_msr_array.sh
# Then:    python3 esm_msr_shard.py concat <OUT_DIR> <combined.csv>

set -euo pipefail

ESM_DIR="${ESM_DIR:-/projects/slavov/andrew/esm/esm-msr}"
VENV="${VENV:-/scratch/leduc.an/msr_venv}"
CKPT="${CKPT:-LoRA_models/esm-msr-small/epoch=03-val_rho_combined_avg=0.816.ckpt}"
BASE=/scratch/leduc.an/AAS_Evo/ANALYSIS/contact_saap
SHARD_DIR="${SHARD_DIR:-$BASE/esm_shards}"
OUT_DIR="${OUT_DIR:-$BASE/esm_scores_shards}"
mkdir -p "$OUT_DIR" /scratch/leduc.an/AAS_Evo/logs

module load anaconda3/2024.06 cuda/13.2.0
source "$VENV/bin/activate"
export HF_HOME="${HF_HOME:-/projects/slavov/andrew/esm/hf_cache}"
export TMPDIR="${TMPDIR:-/scratch/leduc.an/tmp}"
PROXY="${HTTP_PROXY_URL:-http://10.99.0.130:3128}"
export https_proxy="$PROXY" http_proxy="$PROXY"
: "${HF_TOKEN:?set HF_TOKEN (a HuggingFace read token) before sbatch}"

CODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SHARD_DIR/list.txt")
[[ -z "$CODE" ]] && { echo "no code at line $SLURM_ARRAY_TASK_ID"; exit 1; }
OUT="$OUT_DIR/${CODE}.csv"
if [[ -f "$OUT" ]]; then echo "[done] $CODE"; exit 0; fi   # resumable

cd "$ESM_DIR"
python -c "import torch; print('cuda:', torch.cuda.is_available())"
python src/esm_msr/inference.py \
    --checkpoint_path "$CKPT" \
    --input_csv "$SHARD_DIR/${CODE}.csv" --mode singles+doubles \
    --output_csv "$OUT"
echo "wrote $OUT"
