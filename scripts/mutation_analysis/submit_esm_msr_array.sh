#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --job-name=esm_msr_arr
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/esm_msr_arr_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/esm_msr_arr_%A_%a.err

# ESM-MSR scoring. inference.py runs ONE structure per call, but to stay under the
# QOS submit-job limit each array task processes a CHUNK of proteins in a loop.
# Prereq: python3 esm_msr_shard.py split <input.csv> <SHARD_DIR>
# Submit:  NUM=$(wc -l < $SHARD_DIR/list.txt); CHUNK=15
#          NTASK=$(( (NUM + CHUNK - 1) / CHUNK ))
#          CHUNK=$CHUNK sbatch --array=1-${NTASK}%10 submit_esm_msr_array.sh
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

CHUNK="${CHUNK:-15}"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
end=$(( start + CHUNK - 1 ))
echo "task ${SLURM_ARRAY_TASK_ID}: shards ${start}..${end} (CHUNK=${CHUNK})"

cd "$ESM_DIR"
python -c "import torch; print('cuda:', torch.cuda.is_available())"

for ln in $(seq "$start" "$end"); do
    CODE=$(sed -n "${ln}p" "$SHARD_DIR/list.txt")
    [[ -z "$CODE" ]] && continue
    OUT="$OUT_DIR/${CODE}.csv"
    if [[ -f "$OUT" ]]; then echo "[done] $CODE"; continue; fi   # resumable
    echo "[$(date +%H:%M:%S)] scoring $CODE"
    python src/esm_msr/inference.py \
        --checkpoint_path "$CKPT" \
        --input_csv "$SHARD_DIR/${CODE}.csv" --mode "${MODE:-singles+doubles}" \
        --output_csv "$OUT" || echo "[FAIL] $CODE"
done
echo "task ${SLURM_ARRAY_TASK_ID} done"
