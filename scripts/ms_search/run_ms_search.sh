#!/usr/bin/env bash
#
# Set up and submit FragPipe MS searches for all TMT plexes.
#
# Steps:
#   1. Generate per-plex manifests from TMT mapping
#   2. Submit SLURM array job to run FragPipe per plex
#
# Prerequisites:
#   - RAW files downloaded and flattened (bash scripts/download/pdc/flatten_raw.sh)
#   - Per-plex FASTAs generated (sbatch scripts/fasta_gen/submit_proteogenomics.sh)
#   - FragPipe installed and FRAGPIPE_BIN set (or edit submit_fragpipe.sh)
#   - Workflow file exported from FragPipe GUI (TMT-11 closed search recommended)
#
# Usage:
#   bash scripts/ms_search/run_ms_search.sh [workflow_template]
#
# Arguments:
#   workflow_template  Optional path to a FragPipe .workflow file.
#                      If provided, per-plex workflows are generated with
#                      the correct FASTA path patched in.
#                      If omitted, place a default workflow at:
#                        MS_SEARCH/fragpipe.workflow
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo_meta"
DATA_DIR="/scratch/leduc.an/AAS_Evo"

RAW_DIR="${DATA_DIR}/RAW"
FASTA_DIR="${DATA_DIR}/FASTA/per_plex"
TMT_MAP="${META_DIR}/PDC_meta/pdc_file_tmt_map.tsv"
SEARCH_DIR="${DATA_DIR}/MS_SEARCH"
# ---------------------------

WORKFLOW_TEMPLATE="${1:-}"

module load python/3.8.1

echo "========================================="
echo "MS Search Pipeline Setup"
echo "========================================="
echo ""

# Validate inputs
for f in "$TMT_MAP"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

for d in "$RAW_DIR" "$FASTA_DIR"; do
    if [[ ! -d "$d" ]]; then
        echo "ERROR: Required directory not found: $d"
        exit 1
    fi
done

# Step 1: Generate manifests
echo "[$(date)] Generating FragPipe manifests..."

WORKFLOW_FLAG=""
if [[ -n "$WORKFLOW_TEMPLATE" ]]; then
    if [[ ! -f "$WORKFLOW_TEMPLATE" ]]; then
        echo "ERROR: Workflow template not found: $WORKFLOW_TEMPLATE"
        exit 1
    fi
    WORKFLOW_FLAG="--workflow-template $WORKFLOW_TEMPLATE"
fi

python3 "${SCRIPTS_DIR}/generate_manifests.py" \
    --tmt-map "$TMT_MAP" \
    --raw-dir "$RAW_DIR" \
    --fasta-dir "$FASTA_DIR" \
    $WORKFLOW_FLAG \
    -o "$SEARCH_DIR"

echo ""

# Check plex list
PLEX_LIST="${SEARCH_DIR}/plex_list.txt"
if [[ ! -f "$PLEX_LIST" ]]; then
    echo "ERROR: Plex list not generated"
    exit 1
fi

NUM_PLEXES=$(wc -l < "$PLEX_LIST" | tr -d ' ')
if [[ "$NUM_PLEXES" -eq 0 ]]; then
    echo "No plexes ready for search. Check that:"
    echo "  1. RAW files exist in $RAW_DIR"
    echo "  2. Per-plex FASTAs exist in $FASTA_DIR"
    exit 1
fi

# Check for workflow file
if [[ -z "$WORKFLOW_TEMPLATE" && ! -f "${SEARCH_DIR}/fragpipe.workflow" ]]; then
    echo "WARNING: No workflow file found."
    echo ""
    echo "Before submitting the search, provide a FragPipe workflow:"
    echo "  Option A: Export from FragPipe GUI -> ${SEARCH_DIR}/fragpipe.workflow"
    echo "  Option B: Re-run with template:"
    echo "    bash $0 /path/to/template.workflow"
    echo ""
    echo "Manifests have been generated. Plex list: $PLEX_LIST"
    echo "You can submit manually after providing the workflow:"
    echo "  sbatch --array=1-${NUM_PLEXES}%5 ${SCRIPTS_DIR}/submit_fragpipe.sh"
    exit 0
fi

# Step 2: Submit SLURM array job
echo "[$(date)] Submitting FragPipe jobs..."
echo "  Plexes to search: $NUM_PLEXES"
echo "  Concurrent jobs:  5"
echo ""

mkdir -p "${DATA_DIR}/logs"

sbatch --array=1-${NUM_PLEXES}%5 "${SCRIPTS_DIR}/submit_fragpipe.sh"

echo ""
echo "Monitor with: squeue -u \$USER -n fragpipe"
echo ""
echo "When complete, results will be in:"
echo "  ${SEARCH_DIR}/results/{plex_id}/"
