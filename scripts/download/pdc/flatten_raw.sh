#!/usr/bin/env bash
#
# Flatten PDC download directory structure and report progress.
#
# PDC downloads create nested folders like:
#   RAW/PDC000125/1/Raw Mass Spectra/.../Proprietary/file.raw
#
# This script moves all .raw files directly into RAW/ and removes
# the empty nested directories.
#
# Usage:
#   bash flatten_raw.sh [raw_dir] [manifest]
#
# Defaults:
#   raw_dir:  /scratch/leduc.an/AAS_Evo/RAW
#   manifest: $META_DIR/PDC_meta/pdc_all_files.tsv
#

set -e

META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo_meta"
RAW_DIR="${1:-/scratch/leduc.an/AAS_Evo/RAW}"
MANIFEST="${2:-${META_DIR}/PDC_meta/pdc_all_files.tsv}"

echo "========================================"
echo "PDC RAW File Report"
echo "========================================"
echo "RAW directory: $RAW_DIR"
echo ""

# Count expected files from manifest
if [[ -f "$MANIFEST" ]]; then
    # File Name is column 2 in the manifest
    EXPECTED=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
    echo "Expected files (from manifest): $EXPECTED"
else
    echo "WARNING: Manifest not found: $MANIFEST"
    echo "         Cannot report completion percentage."
    EXPECTED=0
fi

# Count all .raw files (case-insensitive) in nested structure
FOUND=$(find "$RAW_DIR" -iname "*.raw" -type f | wc -l | tr -d ' ')
echo "Downloaded .raw files found: $FOUND"

if [[ $EXPECTED -gt 0 ]]; then
    PCT=$(awk "BEGIN {printf \"%.1f\", ($FOUND / $EXPECTED) * 100}")
    echo "Completion: ${PCT}%"
fi
echo ""

# Check for files already flat in RAW/
FLAT=$(find "$RAW_DIR" -maxdepth 1 -iname "*.raw" -type f | wc -l | tr -d ' ')
NESTED=$((FOUND - FLAT))

if [[ $NESTED -eq 0 && $FLAT -eq $FOUND ]]; then
    echo "All files are already flat in $RAW_DIR. Nothing to do."
    exit 0
fi

echo "Files in nested subdirectories: $NESTED"
echo "Files already flat in RAW/: $FLAT"
echo ""

if [[ $NESTED -eq 0 ]]; then
    echo "No nested files to flatten."
    exit 0
fi

# Find missing files (if manifest available)
if [[ -f "$MANIFEST" && $EXPECTED -gt 0 ]]; then
    MISSING=$((EXPECTED - FOUND))
    if [[ $MISSING -gt 0 ]]; then
        echo "Missing files: $MISSING"
        echo "  (Re-run the PDC download to retry HTTP 202 / failed files)"
        echo ""
    fi
fi

# Flatten: move all nested .raw files to RAW/
echo "Flattening nested .raw files to $RAW_DIR ..."
MOVED=0
SKIPPED=0

find "$RAW_DIR" -mindepth 2 -iname "*.raw" -type f | while read -r filepath; do
    filename=$(basename "$filepath")
    dest="$RAW_DIR/$filename"

    if [[ -f "$dest" ]]; then
        # File already exists at top level - check if same size
        src_size=$(stat -c%s "$filepath" 2>/dev/null || stat -f%z "$filepath" 2>/dev/null)
        dst_size=$(stat -c%s "$dest" 2>/dev/null || stat -f%z "$dest" 2>/dev/null)
        if [[ "$src_size" == "$dst_size" ]]; then
            rm "$filepath"
            echo "  [DUP]  $filename (removed nested copy)"
        else
            echo "  [WARN] $filename exists with different size, keeping both"
        fi
    else
        mv "$filepath" "$dest"
        echo "  [MOVE] $filename"
    fi
done

# Remove empty directories
echo ""
echo "Removing empty directories..."
find "$RAW_DIR" -mindepth 1 -type d -empty -delete 2>/dev/null || true

# Final count
FINAL=$(find "$RAW_DIR" -maxdepth 1 -iname "*.raw" -type f | wc -l | tr -d ' ')
echo ""
echo "========================================"
echo "Done. $FINAL .raw files in $RAW_DIR"
if [[ $EXPECTED -gt 0 ]]; then
    PCT=$(awk "BEGIN {printf \"%.1f\", ($FINAL / $EXPECTED) * 100}")
    echo "Completion: ${PCT}%"
fi
echo "========================================"
