#!/usr/bin/env bash
#
# Programmatic metadata setup for AAS_Evo.
#
# PDC-first approach: fetches PDC metadata first (source of truth for all
# searched patients), then queries GDC only for those specific patients.
# Produces gdc_meta_matched.tsv directly — no separate matching step needed.
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
# Step 1: Fetch PDC metadata (per-tissue CSVs via GraphQL API)
# -------------------------------------------------------------------
echo "[Step 1] Fetching PDC metadata from GraphQL API..."
python3 "${REPO_DIR}/scripts/setup/fetch_pdc_metadata.py" $SKIP_URLS
echo ""

# -------------------------------------------------------------------
# Step 2: Consolidate PDC metadata (merge per-tissue CSVs)
# -------------------------------------------------------------------
echo "[Step 2] Consolidating PDC metadata..."
python3 "${REPO_DIR}/scripts/download/pdc/consolidate_metadata.py"
echo ""

# -------------------------------------------------------------------
# Step 3: Fetch GDC BAM metadata for all PDC patients (PDC-first)
# -------------------------------------------------------------------
echo "[Step 3] Fetching GDC WXS BAMs for PDC patients..."
python3 "${REPO_DIR}/scripts/setup/fetch_matched_bams.py"
echo ""

# -------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------
echo "============================================================"
echo "  Setup Complete"
echo "============================================================"
echo ""

if [[ -f "${META_DIR}/GDC_meta/gdc_meta_matched.tsv" ]]; then
    GDC_MATCHED=$(tail -n +2 "${META_DIR}/GDC_meta/gdc_meta_matched.tsv" | wc -l | tr -d ' ')
    echo "  GDC matched:      ${GDC_MATCHED}"
fi
if [[ -f "${META_DIR}/GDC_meta/manifests/manifest_wxs_bams_tissue.tsv" ]]; then
    TISSUE_COUNT=$(tail -n +2 "${META_DIR}/GDC_meta/manifests/manifest_wxs_bams_tissue.tsv" | wc -l | tr -d ' ')
    echo "  Tissue manifest:  ${TISSUE_COUNT} BAMs"
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
echo "  1. Split tissue manifest into download chunks:"
echo "       bash scripts/download/gdc/setup_chunks.sh 400 \\"
echo "           metadata/GDC_meta/manifests/manifest_wxs_bams_tissue.tsv"
echo "  2. Download data (GDC BAMs require a dbGaP token file):"
echo "       bash scripts/download/gdc/submit_download.sh <chunk> <token>"
echo "       sbatch scripts/download/pdc/submit_download.sh"
echo ""
