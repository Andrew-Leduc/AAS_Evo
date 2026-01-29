#!/usr/bin/env bash
#
# Split GDC manifest into chunks and submit each as a separate SLURM job.
# Run this directly on the login node (not via sbatch).
#
# Usage:
#   bash submit_download.sh [manifest] [token] [num_chunks]
#
# Defaults:
#   manifest:   $META_DIR/GDC_meta/manifests/manifest_wxs_bams.tsv
#   token:      $META_DIR/GDC_meta/gdc-user-token_AL.txt
#   num_chunks: 20
#

set -e

# Project paths
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo_meta"
DATA_DIR="/scratch/leduc.an/AAS_Evo"

# Arguments
DEFAULT_MANIFEST="${META_DIR}/GDC_meta/manifests/manifest_wxs_bams.tsv"
DEFAULT_TOKEN="${META_DIR}/GDC_meta/gdc-user-token_AL.txt"
DEFAULT_CHUNKS=20

MANIFEST="${1:-$DEFAULT_MANIFEST}"
TOKEN="${2:-$DEFAULT_TOKEN}"
NUM_CHUNKS="${3:-$DEFAULT_CHUNKS}"

CHUNK_DIR="${META_DIR}/GDC_meta/manifests/manifest_chunks"

# Validate inputs
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: $MANIFEST"
    exit 1
fi

if [[ ! -f "$TOKEN" ]]; then
    echo "ERROR: Token not found: $TOKEN"
    exit 1
fi

# Count total files (exclude header)
TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
LINES_PER_CHUNK=$(( (TOTAL + NUM_CHUNKS - 1) / NUM_CHUNKS ))

echo "========================================"
echo "GDC Batch Download"
echo "========================================"
echo "Manifest: $MANIFEST"
echo "Token: $TOKEN"
echo "Total files: $TOTAL"
echo "Chunks: $NUM_CHUNKS (~$LINES_PER_CHUNK files each)"
echo ""

# Save the header line
HEADER=$(head -1 "$MANIFEST")

# Create chunk directory
mkdir -p "$CHUNK_DIR"
mkdir -p "${DATA_DIR}/logs"

# Split manifest into chunks (skip header, then re-add it to each chunk)
tail -n +2 "$MANIFEST" | split -l "$LINES_PER_CHUNK" -d -a 2 - "${CHUNK_DIR}/chunk_"

# Submit each chunk as a separate job
SUBMITTED=0
for chunk in "${CHUNK_DIR}"/chunk_*; do
    # Prepend header to make it a valid manifest
    CHUNK_MANIFEST="${chunk}.tsv"
    { echo "$HEADER"; cat "$chunk"; } > "$CHUNK_MANIFEST"
    rm "$chunk"

    CHUNK_SIZE=$(tail -n +2 "$CHUNK_MANIFEST" | wc -l | tr -d ' ')
    CHUNK_NAME=$(basename "$CHUNK_MANIFEST")

    echo "Submitting $CHUNK_NAME ($CHUNK_SIZE files)..."
    sbatch "$SCRIPTS_DIR/download_chunk.sh" "$CHUNK_MANIFEST" "$TOKEN"

    SUBMITTED=$((SUBMITTED + 1))

    # Pause between submissions to avoid overwhelming scheduler
    if [[ $SUBMITTED -lt $NUM_CHUNKS ]]; then
        sleep 5
    fi
done

echo ""
echo "Submitted $SUBMITTED jobs."
echo "Monitor with: squeue -u \$USER -n gdc_dl"
echo "Check logs:   ls ${DATA_DIR}/logs/gdc_dl_*.out"
