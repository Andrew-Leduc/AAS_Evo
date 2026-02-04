#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --job-name=fragpipe
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/fragpipe_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/fragpipe_%A_%a.err

#
# Run FragPipe on one TMT plex per SLURM array task.
#
# Each plex is searched against its own custom FASTA database containing
# the reference proteome + plex-specific mutations + compensatory entries.
#
# Prerequisites:
#   - generate_manifests.py has been run (creates manifests/, workflows/, etc.)
#   - FragPipe, MSFragger, Philosopher, and IonQuant are installed
#
# Usage:
#   NUM_PLEXES=$(wc -l < /scratch/leduc.an/AAS_Evo/MS_SEARCH/plex_list.txt)
#   sbatch --array=1-${NUM_PLEXES}%5 submit_fragpipe.sh
#

set -euo pipefail

# --------- PATHS ----------
DATA_DIR="/scratch/leduc.an/AAS_Evo"
SEARCH_DIR="${DATA_DIR}/MS_SEARCH"
PLEX_LIST="${SEARCH_DIR}/plex_list.txt"
FASTA_DIR="${DATA_DIR}/FASTA/per_plex"

# FragPipe installation paths — adjust these to your installation
FRAGPIPE_BIN="${FRAGPIPE_BIN:-/scratch/leduc.an/tools/fragpipe/bin/fragpipe}"
# --------------------------

mkdir -p "${DATA_DIR}/logs"

# Get plex ID for this array task
PLEX_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PLEX_LIST")

if [[ -z "$PLEX_ID" ]]; then
    echo "ERROR: No plex at line ${SLURM_ARRAY_TASK_ID} of $PLEX_LIST"
    exit 1
fi

MANIFEST="${SEARCH_DIR}/manifests/${PLEX_ID}.fp-manifest"
FASTA="${FASTA_DIR}/${PLEX_ID}.fasta"
OUTDIR="${SEARCH_DIR}/results/${PLEX_ID}"
ANNOTATION="${SEARCH_DIR}/annotations/${PLEX_ID}_annotation.tsv"

# Use per-plex workflow if generated, otherwise fall back to default
if [[ -f "${SEARCH_DIR}/workflows/${PLEX_ID}.workflow" ]]; then
    WORKFLOW="${SEARCH_DIR}/workflows/${PLEX_ID}.workflow"
elif [[ -f "${SEARCH_DIR}/fragpipe.workflow" ]]; then
    WORKFLOW="${SEARCH_DIR}/fragpipe.workflow"
else
    echo "ERROR: No workflow file found."
    echo "  Expected: ${SEARCH_DIR}/workflows/${PLEX_ID}.workflow"
    echo "  Or:       ${SEARCH_DIR}/fragpipe.workflow"
    echo ""
    echo "  To create one:"
    echo "    1. Open FragPipe GUI, configure TMT-11 closed search"
    echo "    2. Export workflow to ${SEARCH_DIR}/fragpipe.workflow"
    echo "    3. Or re-run generate_manifests.py with --workflow-template"
    exit 1
fi

# Validate inputs
for f in "$MANIFEST" "$FASTA" "$WORKFLOW"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

# Skip if already completed (check for combined_protein.tsv from Philosopher)
if [[ -f "${OUTDIR}/combined_protein.tsv" ]]; then
    echo "[$(date)] Already completed: $PLEX_ID"
    exit 0
fi

mkdir -p "$OUTDIR"

echo "[$(date)] Starting FragPipe search for plex: $PLEX_ID"
echo "  Manifest:  $MANIFEST"
echo "  FASTA:     $FASTA"
echo "  Workflow:  $WORKFLOW"
echo "  Output:    $OUTDIR"
echo "  Threads:   ${SLURM_CPUS_PER_TASK}"
echo ""

# Copy annotation to output directory (FragPipe TMT-Integrator looks for it there)
if [[ -f "$ANNOTATION" ]]; then
    cp "$ANNOTATION" "${OUTDIR}/experiment_annotation.tsv"
    echo "  Annotation copied to output directory"
fi

# Run FragPipe in headless mode
"$FRAGPIPE_BIN" \
    --headless \
    --workflow "$WORKFLOW" \
    --manifest "$MANIFEST" \
    --output "$OUTDIR" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --ram "$((SLURM_MEM_PER_NODE / 1024))" \
    2>&1

echo ""
echo "[$(date)] Done: $PLEX_ID"
echo "Results: $OUTDIR"
