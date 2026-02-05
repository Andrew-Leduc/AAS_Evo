#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --job-name=dl_uniref30
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/dl_uniref30_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/dl_uniref30_%j.err

#
# Download ColabFold pre-built UniRef30 MMseqs2 database.
#
# This is a pre-indexed database — no mmseqs createdb step needed.
# Just download, extract, and point the MSA pipeline at it.
#
# Source: ColabFold project (Mirdita et al., 2022)
# URL: https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2302.tar.gz
#
# Size: ~68 GB compressed, ~110 GB extracted
#
# Usage:
#   sbatch scripts/setup/download_uniref30.sh
#

set -euo pipefail

SEQ_DIR="/scratch/leduc.an/AAS_Evo/SEQ_FILES"
DB_PREFIX="${SEQ_DIR}/uniref30_2302"
TARBALL="${SEQ_DIR}/uniref30_2302.tar.gz"
URL="https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2302.tar.gz"

mkdir -p "$SEQ_DIR"
mkdir -p /scratch/leduc.an/AAS_Evo/logs

# Check if already extracted (ColabFold extracts directly into SEQ_DIR, not a subdirectory)
if [[ -f "${DB_PREFIX}.tsv" || -f "${DB_PREFIX}.dbtype" ]]; then
    echo "[$(date)] UniRef30 database already exists at ${DB_PREFIX}"
    echo "  To re-download, remove the uniref30_2302* files and ${TARBALL}"
    exit 0
fi

# Download
if [[ -f "$TARBALL" ]]; then
    echo "[$(date)] Tarball already exists, skipping download: ${TARBALL}"
else
    echo "[$(date)] Downloading UniRef30 from ColabFold (~68 GB)..."
    echo "  URL: ${URL}"
    wget --continue -O "$TARBALL" "$URL"
    echo "[$(date)] Download complete."
fi

# Extract (ColabFold tarball extracts files directly, not into a subdirectory)
echo "[$(date)] Extracting (this will take a while)..."
tar -xzf "$TARBALL" -C "$SEQ_DIR"

# Verify
if [[ -f "${DB_PREFIX}.tsv" ]]; then
    echo "[$(date)] UniRef30 database extracted successfully."
    echo "  Database prefix: ${DB_PREFIX}"
    ls -lh "${SEQ_DIR}/" | grep uniref30
else
    echo "[$(date)] WARNING: Expected database files not found."
    echo "  Check extraction. Contents of ${SEQ_DIR}:"
    ls "${SEQ_DIR}/" | grep -i uniref30 || echo "  (no uniref30 files found)"
fi

# Optionally remove tarball to save space
echo ""
echo "[$(date)] Done."
echo "To save ~68 GB disk space, remove the tarball:"
echo "  rm ${TARBALL}"
