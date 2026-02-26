#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --job-name=maxquant
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/maxquant_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/maxquant_%A_%a.err

#
# Run MaxQuant on one TMT plex per SLURM array task.
#
# MaxQuant advantages over FragPipe for proteogenomics:
#   - minPeptidesProtein=1 in mqpar.xml: single-peptide proteins are reported
#     without ProteinProphet's multi-peptide probability penalty
#   - Decoys handled internally (reversed sequences) — no add_decoys step
#
# Prerequisites:
#   - run_maxquant.sh has been run (creates mqpar/ files and plex_list.txt)
#   - MaxQuant downloaded and unzipped at MAXQUANT_DIR
#   - .NET 8 available (module load or direct path)
#
# Usage:
#   NUM_PLEXES=$(wc -l < /scratch/leduc.an/AAS_Evo/MQ_SEARCH/plex_list.txt)
#   sbatch --array=1-${NUM_PLEXES}%5 submit_maxquant.sh
#

set -euo pipefail

# ── PATHS ─────────────────────────────────────────────────────────────────────
DATA_DIR="/scratch/leduc.an/AAS_Evo"
SEARCH_DIR="${DATA_DIR}/MQ_SEARCH"
PLEX_LIST="${SEARCH_DIR}/plex_list.txt"

# MaxQuant — adjust path to wherever you unzipped it
MAXQUANT_DIR="${MAXQUANT_DIR:-/home/leduc.an/bin/MaxQuant}"
MAXQUANT_EXE="${MAXQUANT_DIR}/bin/MaxQuantCmd.dll"

# .NET 8 — installed locally by run_maxquant.sh (no module available on Discovery)
DOTNET_ROOT="${DOTNET_ROOT:-/home/leduc.an/bin/dotnet}"
export DOTNET_ROOT
export PATH="${DOTNET_ROOT}:${PATH}"
# ──────────────────────────────────────────────────────────────────────────────

mkdir -p "${DATA_DIR}/logs"

# Get plex ID for this array task
PLEX_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PLEX_LIST")
if [[ -z "$PLEX_ID" ]]; then
    echo "ERROR: No plex at line ${SLURM_ARRAY_TASK_ID} of $PLEX_LIST"
    exit 1
fi

MQPAR="${SEARCH_DIR}/mqpar/${PLEX_ID}.xml"
OUTDIR="${SEARCH_DIR}/results/${PLEX_ID}"

# Validate
if [[ ! -f "$MQPAR" ]]; then
    echo "ERROR: mqpar.xml not found: $MQPAR"
    echo "  Run: bash scripts/ms_search/maxquant/run_maxquant.sh"
    exit 1
fi

if [[ ! -f "$MAXQUANT_EXE" ]]; then
    echo "ERROR: MaxQuantCmd.dll not found: $MAXQUANT_EXE"
    echo "  Set MAXQUANT_DIR or edit this script"
    exit 1
fi

# Skip if already completed (evidence.txt is always written by MaxQuant)
if [[ -f "${OUTDIR}/combined/txt/evidence.txt" ]]; then
    echo "[$(date)] Already completed: $PLEX_ID"
    exit 0
fi

mkdir -p "$OUTDIR"

# Verify dotnet is accessible
if ! command -v dotnet &>/dev/null; then
    echo "ERROR: dotnet not found at ${DOTNET_ROOT}"
    echo "  Run run_maxquant.sh first — it installs .NET 8 locally."
    exit 1
fi

DOTNET_VERSION=$(dotnet --version)
echo "[$(date)] Starting MaxQuant for plex: $PLEX_ID"
echo "  dotnet version: $DOTNET_VERSION"
echo "  mqpar:          $MQPAR"
echo "  output:         $OUTDIR/combined/"
echo "  threads:        ${SLURM_CPUS_PER_TASK}"
echo ""

# Required on RHEL8+ HPC clusters: .NET 6+ needs ICU or globalization invariant
# mode for culture-sensitive ops in the Thermo RAW reader. Without this,
# MaxQuant silently reads 0 MS scans from all RAW files.
export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1

# RAW files are symlinked into $OUTDIR by generate_mqpar.py, so MaxQuant
# creates combined/ at $OUTDIR/combined/ automatically.
dotnet "$MAXQUANT_EXE" "$MQPAR"

echo ""
echo "[$(date)] Done: $PLEX_ID"
echo "Results: ${OUTDIR}/combined/txt/"
