#!/usr/bin/env bash
#
# Split the GDC manifest into persistent chunks for sequential processing.
#
# Each chunk is a valid GDC manifest (with header) that can be passed directly
# to download_chunk.sh. This avoids downloading all ~2K BAMs at once.
#
# Usage:
#   bash scripts/chunk_workflow/setup_chunks.sh            # 500 BAMs per chunk
#   bash scripts/chunk_workflow/setup_chunks.sh 300         # custom chunk size
#   bash scripts/chunk_workflow/setup_chunks.sh 500 /path/to/manifest.tsv
#

set -euo pipefail

# --------- PATHS ----------
META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo_meta"
DEFAULT_MANIFEST="${META_DIR}/GDC_meta/manifest_wxs_bams.tsv"
CHUNK_DIR="${META_DIR}/GDC_meta/manifests/chunks"
# ---------------------------

CHUNK_SIZE="${1:-500}"
MANIFEST="${2:-$DEFAULT_MANIFEST}"

# Validate
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: $MANIFEST"
    exit 1
fi

if ! [[ "$CHUNK_SIZE" =~ ^[0-9]+$ ]] || [[ "$CHUNK_SIZE" -lt 1 ]]; then
    echo "ERROR: Chunk size must be a positive integer, got: $CHUNK_SIZE"
    exit 1
fi

# Count total files (exclude header)
HEADER=$(head -1 "$MANIFEST")
TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')

if [[ "$TOTAL" -eq 0 ]]; then
    echo "ERROR: Manifest has no data rows"
    exit 1
fi

NUM_CHUNKS=$(( (TOTAL + CHUNK_SIZE - 1) / CHUNK_SIZE ))
LAST_CHUNK_SIZE=$(( TOTAL - (NUM_CHUNKS - 1) * CHUNK_SIZE ))

echo "========================================"
echo "GDC Manifest Chunk Setup"
echo "========================================"
echo "Manifest:       $MANIFEST"
echo "Total BAMs:     $TOTAL"
echo "Chunk size:     $CHUNK_SIZE"
echo "Chunks:         $NUM_CHUNKS"
echo "Last chunk:     $LAST_CHUNK_SIZE BAMs"
echo ""

# Clean previous chunks
mkdir -p "$CHUNK_DIR"
rm -f "${CHUNK_DIR}"/chunk_*.tsv

# Split into raw chunks (no header)
tail -n +2 "$MANIFEST" | split -l "$CHUNK_SIZE" -d -a 2 - "${CHUNK_DIR}/chunk_"

# Prepend header to each chunk and rename to .tsv
for chunk in "${CHUNK_DIR}"/chunk_*; do
    # Skip if already a .tsv (shouldn't happen after rm above)
    [[ "$chunk" == *.tsv ]] && continue

    CHUNK_TSV="${chunk}.tsv"
    { echo "$HEADER"; cat "$chunk"; } > "$CHUNK_TSV"
    rm "$chunk"

    CHUNK_NUM=$(basename "$CHUNK_TSV" .tsv | sed 's/chunk_//')
    CHUNK_COUNT=$(tail -n +2 "$CHUNK_TSV" | wc -l | tr -d ' ')
    echo "  chunk_${CHUNK_NUM}.tsv  ($CHUNK_COUNT BAMs)"
done

echo ""
echo "Chunks written to: $CHUNK_DIR"
echo ""
echo "To process chunk 1:"
echo "  bash scripts/chunk_workflow/process_chunk.sh 1"
echo ""
echo "Workflow per chunk:"
echo "  1. process_chunk.sh N download    # Download BAMs"
echo "  2. process_chunk.sh N prep-calls  # Write bam_list.txt"
echo "  3. sbatch variant calling         # (command printed by prep-calls)"
echo "  4. process_chunk.sh N prep-vep    # Write vcf_list.txt"
echo "  5. sbatch VEP annotation          # (command printed by prep-vep)"
echo "  6. process_chunk.sh N cleanup     # Delete BAMs, keep VCF/VEP"
