#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=108:00:00
#SBATCH --partition=slavov
#SBATCH --job-name=pdc_dl
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/pdc_dl_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/pdc_dl_%j.err

#
# Submit PDC download job
#
# Usage:
#   sbatch submit_download.sh [manifest] [output_dir]
#
# Defaults:
#   manifest: $META_DIR/PDC_meta/pdc_all_files.tsv
#   output: /scratch/leduc.an/AAS_Evo/RAW
#

set -e

# Load Python module (required on compute nodes)
module load python/3.8.1 2>/dev/null || true

# Project paths
SCRIPTS_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo"
META_DIR="${SCRIPTS_DIR}/metadata"
DATA_DIR="/scratch/leduc.an/AAS_Evo"

# Default arguments
DEFAULT_MANIFEST="${META_DIR}/PDC_meta/pdc_all_files.tsv"
DEFAULT_OUTPUT="${DATA_DIR}/RAW"

MANIFEST="${1:-$DEFAULT_MANIFEST}"
OUTPUT_DIR="${2:-$DEFAULT_OUTPUT}"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "${DATA_DIR}/logs"

echo "========================================"
echo "PDC Download Job"
echo "========================================"
echo "Started: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: $(hostname)"
echo ""
echo "Manifest: $MANIFEST"
echo "Output: $OUTPUT_DIR"
echo ""

# Check manifest exists
if [[ ! -f "$MANIFEST" ]]; then
    echo "Error: Manifest not found: $MANIFEST"
    exit 1
fi

# Count files to download
TOTAL=$(tail -n +2 "$MANIFEST" | wc -l | tr -d ' ')
echo "Files in manifest: $TOTAL"
echo ""

# Ensure required Python packages are installed
python3 -m pip install --user --quiet requests urllib3

# Run download script
# Use -u for unbuffered output so logs appear in real-time
# Use python3 explicitly (cluster's 'python' may be Python 2)
cd "$SCRIPTS_DIR/scripts/download/pdc"
python3 -u download.py "$MANIFEST" -o "$OUTPUT_DIR"

echo ""
echo "========================================"
echo "Completed: $(date)"
echo "========================================"
