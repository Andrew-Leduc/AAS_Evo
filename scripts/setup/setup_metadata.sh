#!/usr/bin/env bash
#
# Programmatic metadata setup for AAS_Evo.
#
# Fetches all metadata from GDC and PDC APIs, then consolidates and
# generates matched manifests. After running this, you're ready to
# download data.
#
# Usage:
#   bash scripts/setup/setup_metadata.sh
#
#   # Skip PDC signed URLs (faster, refresh later with refresh_urls.py):
#   bash scripts/setup/setup_metadata.sh --skip-urls
#
# Prerequisites:
#   - Python 3.8+ with standard library (no extra packages needed)
#   - GDC token is only needed for downloads, not for this setup script
#

set -e

REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
META_DIR="${REPO_DIR}/metadata"
SKIP_URLS=""

# Parse args
for arg in "$@"; do
    case "$arg" in
        --skip-urls) SKIP_URLS="--skip-urls" ;;
        *) echo "Unknown arg: $arg"; exit 1 ;;
    esac
done

echo "============================================================"
echo "  AAS_Evo Metadata Setup"
echo "============================================================"
echo ""
echo "  Repo:     $REPO_DIR"
echo "  Meta dir: $META_DIR"
echo ""

# -------------------------------------------------------------------
# Step 0: Create directory structure
# -------------------------------------------------------------------
echo "[Step 0] Creating directory structure..."
mkdir -p "${META_DIR}/GDC_meta/manifests/chunks"
mkdir -p "${META_DIR}/PDC_meta"
echo "  Done."
echo ""

# -------------------------------------------------------------------
# Step 1: Fetch GDC WXS BAM manifest (via GDC REST API)
# -------------------------------------------------------------------
echo "[Step 1] Fetching GDC WXS BAM manifest..."
python3 "${REPO_DIR}/scripts/setup/fetch_gdc_manifest.py"
echo ""

# -------------------------------------------------------------------
# Step 2: Fetch GDC sample metadata (enrich manifest with case info)
# -------------------------------------------------------------------
echo "[Step 2] Fetching GDC sample metadata..."
python3 "${REPO_DIR}/scripts/download/gdc/fetch_metadata.py"
echo ""

# -------------------------------------------------------------------
# Step 3: Fetch PDC metadata (per-tissue CSVs via GraphQL API)
# -------------------------------------------------------------------
echo "[Step 3] Fetching PDC metadata from GraphQL API..."
python3 "${REPO_DIR}/scripts/setup/fetch_pdc_metadata.py" $SKIP_URLS
echo ""

# -------------------------------------------------------------------
# Step 4: Consolidate PDC metadata (merge per-tissue CSVs)
# -------------------------------------------------------------------
echo "[Step 4] Consolidating PDC metadata..."
python3 "${REPO_DIR}/scripts/download/pdc/consolidate_metadata.py"
echo ""

# -------------------------------------------------------------------
# Step 5: Generate GDC↔PDC matching report and pruned manifests
# -------------------------------------------------------------------
echo "[Step 5] Generating initial matching report..."
python3 "${REPO_DIR}/scripts/download/mapping_report.py"
echo ""

# -------------------------------------------------------------------
# Step 6: Recover unmatched PDC patients from GDC
# -------------------------------------------------------------------
echo "[Step 6] Searching GDC for BAMs matching unmatched PDC patients..."
GDC_META="${META_DIR}/GDC_meta/gdc_meta.tsv"
TMT_MAP="${META_DIR}/PDC_meta/pdc_file_tmt_map.tsv"
UNMATCHED_OUT="${META_DIR}/GDC_meta/gdc_meta_unmatched.tsv"

if [[ -f "$GDC_META" && -f "$TMT_MAP" ]]; then
    python3 "${REPO_DIR}/scripts/download/gdc/fetch_unmatched_bams.py" \
        --gdc-meta "$GDC_META" \
        --tmt-map "$TMT_MAP" \
        -o "$UNMATCHED_OUT"

    # If new BAMs were found, merge and re-run matching
    if [[ -f "$UNMATCHED_OUT" ]] && [[ $(wc -l < "$UNMATCHED_OUT") -gt 1 ]]; then
        echo ""
        echo "  Merging newly found BAMs into gdc_meta.tsv..."
        tail -n +2 "$UNMATCHED_OUT" >> "$GDC_META"
        echo "  Re-running matching report..."
        python3 "${REPO_DIR}/scripts/download/mapping_report.py"
    fi
else
    echo "  Skipping (gdc_meta.tsv or pdc_file_tmt_map.tsv not found)"
fi
echo ""

# -------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------
echo "============================================================"
echo "  Setup Complete"
echo "============================================================"
echo ""

# Count files in key outputs
if [[ -f "${META_DIR}/GDC_meta/manifests/manifest_wxs_bams.tsv" ]]; then
    GDC_COUNT=$(tail -n +2 "${META_DIR}/GDC_meta/manifests/manifest_wxs_bams.tsv" | wc -l | tr -d ' ')
    echo "  GDC WXS BAMs:     ${GDC_COUNT}"
fi
if [[ -f "${META_DIR}/GDC_meta/gdc_meta_matched.tsv" ]]; then
    GDC_MATCHED=$(tail -n +2 "${META_DIR}/GDC_meta/gdc_meta_matched.tsv" | wc -l | tr -d ' ')
    echo "  GDC matched:      ${GDC_MATCHED}"
fi
if [[ -f "${META_DIR}/PDC_meta/pdc_all_files.tsv" ]]; then
    PDC_COUNT=$(tail -n +2 "${META_DIR}/PDC_meta/pdc_all_files.tsv" | wc -l | tr -d ' ')
    echo "  PDC RAW files:    ${PDC_COUNT}"
fi
if [[ -f "${META_DIR}/PDC_meta/pdc_all_files_matched.tsv" ]]; then
    PDC_MATCHED=$(tail -n +2 "${META_DIR}/PDC_meta/pdc_all_files_matched.tsv" | wc -l | tr -d ' ')
    echo "  PDC matched:      ${PDC_MATCHED}"
fi

echo ""
echo "Next steps:"
echo "  1. Split GDC manifest into chunks:"
echo "       bash scripts/download/gdc/setup_chunks.sh"
echo "  2. Download data (GDC BAMs require a dbGaP token file):"
echo "       bash scripts/download/gdc/submit_download.sh <chunk> <token>"
echo "       sbatch scripts/download/pdc/submit_download.sh"
echo ""
