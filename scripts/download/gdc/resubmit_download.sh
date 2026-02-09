#!/usr/bin/env bash
#
# Re-submit GDC downloads for BAMs that haven't been downloaded yet.
# Checks BAMS/ directory for existing UUID subdirectories with .bam files,
# builds a filtered manifest of only missing files, then splits and submits.
#
# Run this directly on the login node (not via sbatch).
#
# Usage:
#   bash resubmit_download.sh [manifest] [token] [num_chunks]
#

set -e

# Project paths
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo/metadata"
DATA_DIR="/scratch/leduc.an/AAS_Evo"
BAMS_DIR="${DATA_DIR}/BAMS"

# Arguments
DEFAULT_MANIFEST="${META_DIR}/GDC_meta/manifests/manifest_wxs_bams.tsv"
DEFAULT_TOKEN="${META_DIR}/GDC_meta/.gdc-user-token.txt"
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

if [[ ! -d "$BAMS_DIR" ]]; then
    echo "BAMS directory not found: $BAMS_DIR"
    echo "No downloads exist yet. Use submit_download.sh instead."
    exit 1
fi

# Header line from manifest
HEADER=$(head -1 "$MANIFEST")
TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')

# Check which UUIDs already have a complete .bam file downloaded
# GDC downloads each BAM into BAMS/{uuid}/{filename}.bam
# Manifest columns: id  filename  md5  size  state
echo "Checking existing downloads in $BAMS_DIR ..."

FILTERED_MANIFEST="${CHUNK_DIR}/missing_manifest.tsv"
mkdir -p "$CHUNK_DIR"

echo "$HEADER" > "$FILTERED_MANIFEST"

HAVE=0
MISSING=0
PARTIAL=0

while IFS=$'\t' read -r uuid filename md5sum expected_size state; do
    bam_path="${BAMS_DIR}/${uuid}/${filename}"

    if [[ -f "$bam_path" ]]; then
        actual_size=$(stat -c%s "$bam_path" 2>/dev/null || stat -f%z "$bam_path" 2>/dev/null || echo 0)
        if [[ "$actual_size" == "$expected_size" ]]; then
            HAVE=$((HAVE + 1))
        else
            # Partial download — remove the incomplete file so gdc-client retries
            rm -f "$bam_path"
            echo -e "${uuid}\t${filename}\t${md5sum}\t${expected_size}\t${state}" >> "$FILTERED_MANIFEST"
            PARTIAL=$((PARTIAL + 1))
            MISSING=$((MISSING + 1))
        fi
    else
        echo -e "${uuid}\t${filename}\t${md5sum}\t${expected_size}\t${state}" >> "$FILTERED_MANIFEST"
        MISSING=$((MISSING + 1))
    fi
done < <(tail -n +2 "$MANIFEST")

echo ""
echo "========================================"
echo "GDC Re-submit Download"
echo "========================================"
echo "Manifest:    $MANIFEST"
echo "Total:       $TOTAL"
echo "Complete:    $HAVE"
echo "Partial:     $PARTIAL (removed, will re-download)"
echo "Missing:     $MISSING"
echo ""

if [[ $MISSING -eq 0 ]]; then
    echo "All files already downloaded. Nothing to do."
    rm -f "$FILTERED_MANIFEST"
    exit 0
fi

# Adjust chunk count if fewer files than chunks
if [[ $MISSING -lt $NUM_CHUNKS ]]; then
    NUM_CHUNKS=$MISSING
fi

LINES_PER_CHUNK=$(( (MISSING + NUM_CHUNKS - 1) / NUM_CHUNKS ))

echo "Submitting $MISSING files in $NUM_CHUNKS chunks (~$LINES_PER_CHUNK each)"
echo ""

# Clean old chunks
rm -f "${CHUNK_DIR}"/chunk_*

# Split filtered manifest into chunks
tail -n +2 "$FILTERED_MANIFEST" | split -l "$LINES_PER_CHUNK" -d -a 2 - "${CHUNK_DIR}/chunk_"

mkdir -p "${DATA_DIR}/logs"

# Submit each chunk
SUBMITTED=0
for chunk in "${CHUNK_DIR}"/chunk_*; do
    CHUNK_MANIFEST="${chunk}.tsv"
    { echo "$HEADER"; cat "$chunk"; } > "$CHUNK_MANIFEST"
    rm "$chunk"

    CHUNK_SIZE=$(tail -n +2 "$CHUNK_MANIFEST" | wc -l | tr -d ' ')
    CHUNK_NAME=$(basename "$CHUNK_MANIFEST")

    echo "Submitting $CHUNK_NAME ($CHUNK_SIZE files)..."
    sbatch "$SCRIPTS_DIR/download_chunk.sh" "$CHUNK_MANIFEST" "$TOKEN"

    SUBMITTED=$((SUBMITTED + 1))

    if [[ $SUBMITTED -lt $NUM_CHUNKS ]]; then
        sleep 5
    fi
done

echo ""
echo "Submitted $SUBMITTED jobs."
echo "Monitor with: squeue -u \$USER -n gdc_dl"
echo "Check logs:   ls ${DATA_DIR}/logs/gdc_dl_*.out"
