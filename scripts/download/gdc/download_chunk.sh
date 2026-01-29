#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --job-name=gdc_dl
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/gdc_dl_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/gdc_dl_%j.err

#
# Download a chunk of BAM files from GDC
#
# This script is submitted by submit_download.sh with a sub-manifest.
#
# Usage:
#   sbatch download_chunk.sh <sub_manifest> <token>
#

set -e

MANIFEST="$1"
TOKEN="$2"
OUTPUT_DIR="/scratch/leduc.an/AAS_Evo/BAMS"

if [[ -z "$MANIFEST" || -z "$TOKEN" ]]; then
    echo "Usage: sbatch download_chunk.sh <manifest> <token>"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

CHUNK_NAME=$(basename "$MANIFEST")
TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')

echo "========================================"
echo "GDC Download: $CHUNK_NAME"
echo "========================================"
echo "Started: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: $(hostname)"
echo "Files in chunk: $TOTAL"
echo ""

cd "$OUTPUT_DIR"
gdc-client download \
    -m "$MANIFEST" \
    -t "$TOKEN" \
    -n 4 \
    --retry-amount 5

echo ""
echo "========================================"
echo "Completed: $(date)"
echo "========================================"
