#!/bin/bash
#
# Download PDC raw files from pdc_all_files.tsv manifest
#
# Usage:
#   ./download_files.sh [manifest_file] [output_dir]
#
# Defaults:
#   manifest: $META_DIR/PDC_meta/pdc_all_files.tsv
#   output: $RAW_DIR
#
# NOTE: Download URLs expire! Re-export manifests from PDC if downloads fail.

set -e

# Default paths (adjust for your environment)
DEFAULT_MANIFEST="${HOME}/AAS_Evo_project/AAS_Evo/metadata/PDC_meta/pdc_all_files.tsv"
DEFAULT_OUTPUT="/scratch/leduc.an/AAS_Evo/RAW"

MANIFEST="${1:-$DEFAULT_MANIFEST}"
OUTPUT_DIR="${2:-$DEFAULT_OUTPUT}"

# Check manifest exists
if [[ ! -f "$MANIFEST" ]]; then
    echo "Error: Manifest not found: $MANIFEST"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count files
TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
echo "Downloading $TOTAL files to $OUTPUT_DIR"
echo ""

# Download each file
COUNT=0
FAILED=0

# Skip header, read TSV columns
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r file_id file_name run_id study_name pdc_study_id file_size md5sum tissue_folder download_url; do
    COUNT=$((COUNT + 1))

    # Skip if URL is empty
    if [[ -z "$download_url" ]]; then
        echo "[$COUNT/$TOTAL] Skipping $file_name - no download URL"
        continue
    fi

    # Skip if file already exists and has correct size
    OUTPUT_FILE="$OUTPUT_DIR/$file_name"
    if [[ -f "$OUTPUT_FILE" ]]; then
        EXISTING_SIZE=$(stat -f%z "$OUTPUT_FILE" 2>/dev/null || stat -c%s "$OUTPUT_FILE" 2>/dev/null)
        if [[ "$EXISTING_SIZE" == "$file_size" ]]; then
            echo "[$COUNT/$TOTAL] Skipping $file_name - already exists"
            continue
        fi
    fi

    echo "[$COUNT/$TOTAL] Downloading $file_name..."

    # Download with wget (curl as fallback)
    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "$OUTPUT_FILE" "$download_url" || {
            echo "  Failed: $file_name"
            FAILED=$((FAILED + 1))
            rm -f "$OUTPUT_FILE"
        }
    else
        curl -# -L -o "$OUTPUT_FILE" "$download_url" || {
            echo "  Failed: $file_name"
            FAILED=$((FAILED + 1))
            rm -f "$OUTPUT_FILE"
        }
    fi
done

echo ""
echo "Download complete. Failed: $FAILED"
