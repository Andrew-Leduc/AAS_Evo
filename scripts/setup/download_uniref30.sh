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
DB_DIR="${SEQ_DIR}/uniref30_2302"
TARBALL="${SEQ_DIR}/uniref30_2302.tar.gz"
URL="https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2302.tar.gz"

mkdir -p "$SEQ_DIR"
mkdir -p /scratch/leduc.an/AAS_Evo/logs

# Check if already extracted
if [[ -f "${DB_DIR}/uniref30_2302" || -f "${DB_DIR}/uniref30_2302.dbtype" ]]; then
    echo "[$(date)] UniRef30 database already exists at ${DB_DIR}"
    echo "  To re-download, remove: rm -rf ${DB_DIR} ${TARBALL}"
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

# Extract
echo "[$(date)] Extracting (this will take a while)..."
mkdir -p "$DB_DIR"
tar -xzf "$TARBALL" -C "$SEQ_DIR"

# Verify
if [[ -f "${DB_DIR}/uniref30_2302" || -f "${DB_DIR}/uniref30_2302.dbtype" ]]; then
    echo "[$(date)] UniRef30 database extracted successfully."
    echo "  Location: ${DB_DIR}"
    ls -lh "${DB_DIR}/"
else
    # Some tarballs extract into the current directory directly
    # Check if files ended up in SEQ_DIR instead
    if [[ -f "${SEQ_DIR}/uniref30_2302.dbtype" ]]; then
        echo "[$(date)] UniRef30 database extracted to ${SEQ_DIR}"
        echo "  Database prefix: ${SEQ_DIR}/uniref30_2302"
    else
        echo "[$(date)] WARNING: Expected database files not found."
        echo "  Check extraction path. Contents of ${SEQ_DIR}:"
        ls "${SEQ_DIR}/" | grep -i uniref30 || echo "  (no uniref30 files found)"
    fi
fi

# Optionally remove tarball to save space
echo ""
echo "[$(date)] Done."
echo "To save ~68 GB disk space, remove the tarball:"
echo "  rm ${TARBALL}"
