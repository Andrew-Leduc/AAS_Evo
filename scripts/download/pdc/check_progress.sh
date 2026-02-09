#!/usr/bin/env bash
#
# Check PDC download progress against manifest.
# Does NOT move or modify any files.
#
# Usage:
#   bash check_progress.sh [raw_dir] [manifest]
#

set -e

META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo/metadata"
RAW_DIR="${1:-/scratch/leduc.an/AAS_Evo/RAW}"
MANIFEST="${2:-${META_DIR}/PDC_meta/pdc_all_files.tsv}"

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: $MANIFEST"
    exit 1
fi

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: RAW directory not found: $RAW_DIR"
    exit 1
fi

# Get expected file names from manifest (column 2 = "File Name")
EXPECTED=$(tail -n +2 "$MANIFEST" | cut -f2 | sed 's/^ *//;s/ *$//' | sort -u)
EXPECTED_COUNT=$(echo "$EXPECTED" | wc -l | tr -d ' ')

# Get all downloaded .raw files (by name only)
FOUND=$(find "$RAW_DIR" -iname "*.raw" -type f -exec basename {} \; | sort -u)
FOUND_COUNT=$(echo "$FOUND" | grep -c . 2>/dev/null || echo 0)

# Cross-reference: which expected files exist on disk
HAVE_COUNT=0
MISSING_COUNT=0
MISSING_LIST=""

while IFS= read -r fname; do
    if echo "$FOUND" | grep -qxF "$fname"; then
        HAVE_COUNT=$((HAVE_COUNT + 1))
    else
        MISSING_COUNT=$((MISSING_COUNT + 1))
        MISSING_LIST="${MISSING_LIST}${fname}"$'\n'
    fi
done <<< "$EXPECTED"

PCT=$(awk "BEGIN {printf \"%.1f\", ($HAVE_COUNT / $EXPECTED_COUNT) * 100}")

echo "========================================"
echo "PDC Download Progress"
echo "========================================"
echo "Manifest:  $MANIFEST"
echo "RAW dir:   $RAW_DIR"
echo ""
echo "Expected:   $EXPECTED_COUNT files"
echo "Downloaded: $HAVE_COUNT files"
echo "Missing:    $MISSING_COUNT files"
echo "Completion: ${PCT}%"
echo "========================================"

if [[ $MISSING_COUNT -gt 0 && $MISSING_COUNT -le 50 ]]; then
    echo ""
    echo "Missing files:"
    echo "$MISSING_LIST" | head -50
elif [[ $MISSING_COUNT -gt 50 ]]; then
    echo ""
    echo "First 50 missing files:"
    echo "$MISSING_LIST" | head -50
    echo "  ... and $((MISSING_COUNT - 50)) more"
fi
