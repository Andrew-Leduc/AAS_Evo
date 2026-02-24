#!/usr/bin/env bash
#
# Set up and submit MaxQuant MS searches for all TMT plexes.
#
# Steps:
#   1. Unzip MaxQuant if not already done
#   2. Generate per-plex mqpar.xml files
#   3. Submit SLURM array job
#
# Prerequisites:
#   - MaxQuant zip downloaded to cluster (e.g. ~/bin/MaxQuant_2.x.x.x.zip)
#   - RAW files available in DATA_DIR/RAW/
#   - Per-plex FASTAs generated (sbatch scripts/fasta_gen/submit_proteogenomics.sh)
#   - .NET 8.0+ available (module load dotnet/8.0)
#
# Usage:
#   bash scripts/ms_search/maxquant/run_maxquant.sh [path/to/MaxQuant.zip]
#

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="/scratch/leduc.an/AAS_Evo"
REPO_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo"
META_DIR="${REPO_DIR}/metadata"

RAW_DIR="${DATA_DIR}/RAW"
FASTA_DIR="${DATA_DIR}/FASTA/per_plex"
TMT_MAP="${META_DIR}/PDC_meta/pdc_file_tmt_map.tsv"
SEARCH_DIR="${DATA_DIR}/MQ_SEARCH"

MAXQUANT_INSTALL="${MAXQUANT_DIR:-/home/leduc.an/bin/MaxQuant}"
DOTNET_INSTALL="${DOTNET_ROOT:-/home/leduc.an/bin/dotnet}"
MQ_ZIP="${1:-}"

echo "========================================="
echo "MaxQuant Search Pipeline Setup"
echo "========================================="
echo ""

# ── Step 1: Unzip MaxQuant if needed ─────────────────────────────────────────
if [[ ! -f "${MAXQUANT_INSTALL}/bin/MaxQuantCmd.exe" ]]; then
    if [[ -z "$MQ_ZIP" ]]; then
        # Try to find a zip automatically in ~/bin/
        MQ_ZIP=$(ls ~/bin/MaxQuant*.zip 2>/dev/null | head -1 || true)
    fi

    if [[ -z "$MQ_ZIP" || ! -f "$MQ_ZIP" ]]; then
        echo "ERROR: MaxQuantCmd.exe not found at ${MAXQUANT_INSTALL}"
        echo "  Either:"
        echo "    a) Provide the zip: bash $0 /path/to/MaxQuant_X.Y.Z.zip"
        echo "    b) Set MAXQUANT_DIR to where it's already unzipped"
        exit 1
    fi

    echo "[$(date)] Unzipping MaxQuant: $MQ_ZIP"
    mkdir -p "$MAXQUANT_INSTALL"
    unzip -q "$MQ_ZIP" -d "$MAXQUANT_INSTALL"

    # MaxQuant zips usually have a top-level folder inside — flatten if needed
    inner=$(ls "$MAXQUANT_INSTALL" | head -1)
    if [[ -f "${MAXQUANT_INSTALL}/${inner}/bin/MaxQuantCmd.exe" ]]; then
        mv "${MAXQUANT_INSTALL}/${inner}"/* "$MAXQUANT_INSTALL/"
        rmdir "${MAXQUANT_INSTALL}/${inner}" 2>/dev/null || true
    fi

    echo "  Installed to: $MAXQUANT_INSTALL"
else
    echo "[$(date)] MaxQuant already installed: $MAXQUANT_INSTALL"
fi

# ── Step 1b: Install .NET 8 locally if not present ───────────────────────────
export DOTNET_ROOT="$DOTNET_INSTALL"
export PATH="${DOTNET_ROOT}:${PATH}"

if ! command -v dotnet &>/dev/null || [[ "$(dotnet --version | cut -d. -f1)" -lt 8 ]]; then
    echo "[$(date)] Installing .NET 8 locally to ${DOTNET_INSTALL}..."
    mkdir -p "$DOTNET_INSTALL"
    # Official Microsoft install script — installs without root
    curl -sSL https://dot.net/v1/dotnet-install.sh -o /tmp/dotnet-install.sh
    bash /tmp/dotnet-install.sh --channel 8.0 --install-dir "$DOTNET_INSTALL"
    rm -f /tmp/dotnet-install.sh
    echo "  .NET installed: $(dotnet --version)"
else
    echo "[$(date)] .NET already available: $(dotnet --version)"
fi
echo ""

# ── Step 2: Generate mqpar.xml files ─────────────────────────────────────────
echo "[$(date)] Generating mqpar.xml files..."

mkdir -p "${SEARCH_DIR}"

python3 "${SCRIPTS_DIR}/generate_mqpar.py" \
    --tmt-map   "$TMT_MAP" \
    --raw-dir   "$RAW_DIR" \
    --fasta-dir "$FASTA_DIR" \
    --threads   "${SLURM_CPUS_PER_TASK:-16}" \
    -o          "$SEARCH_DIR"

echo ""

# ── Step 3: Submit SLURM array ────────────────────────────────────────────────
PLEX_LIST="${SEARCH_DIR}/plex_list.txt"
if [[ ! -f "$PLEX_LIST" ]]; then
    echo "ERROR: Plex list not generated: $PLEX_LIST"
    exit 1
fi

NUM_PLEXES=$(wc -l < "$PLEX_LIST" | tr -d ' ')
if [[ "$NUM_PLEXES" -eq 0 ]]; then
    echo "No plexes ready. Check RAW files in $RAW_DIR and FASTAs in $FASTA_DIR"
    exit 1
fi

mkdir -p "${DATA_DIR}/logs"

echo "[$(date)] Submitting MaxQuant jobs..."
echo "  Plexes to search: $NUM_PLEXES"
echo "  Concurrent jobs:  5"
echo ""

sbatch --array=1-${NUM_PLEXES}%5 "${SCRIPTS_DIR}/submit_maxquant.sh"

echo ""
echo "Monitor with: squeue -u \$USER -n maxquant"
echo ""
echo "Results will be in:"
echo "  ${SEARCH_DIR}/results/{plex_id}/combined/txt/"
echo ""
echo "Key output files per plex:"
echo "  msms.txt          — PSM-level results (analog of psm.txt)"
echo "  peptides.txt      — peptide-level"
echo "  proteinGroups.txt — protein-level (minPeptides=1, single-peptide ok)"
echo "  evidence.txt      — all peptide evidence with TMT reporter intensities"
