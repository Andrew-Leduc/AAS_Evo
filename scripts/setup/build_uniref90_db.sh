#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=short
#SBATCH --job-name=build_uniref90
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/build_uniref90_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/build_uniref90_%j.err

#
# Build MMseqs2 database from UniRef90 FASTA.
# Decompresses first to avoid "Generic" type issues with .gz input.
#
# Usage:
#   sbatch scripts/setup/build_uniref90_db.sh
#

set -euo pipefail

SEQ_DIR="/scratch/leduc.an/AAS_Evo/SEQ_FILES"
GZ_FILE="${SEQ_DIR}/uniref90.fasta.gz"
FASTA_FILE="${SEQ_DIR}/uniref90.fasta"
DB_PREFIX="${SEQ_DIR}/uniref90"

mkdir -p /scratch/leduc.an/AAS_Evo/logs

if [[ ! -f "$GZ_FILE" ]]; then
    echo "ERROR: uniref90.fasta.gz not found at $GZ_FILE"
    exit 1
fi

# Remove any broken database files from previous attempts
echo "[$(date)] Cleaning up previous database files..."
rm -f "${DB_PREFIX}" "${DB_PREFIX}".dbtype "${DB_PREFIX}".index \
      "${DB_PREFIX}".lookup "${DB_PREFIX}"_h "${DB_PREFIX}"_h.dbtype \
      "${DB_PREFIX}"_h.index

# Decompress (keep the .gz)
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "[$(date)] Decompressing uniref90.fasta.gz..."
    gunzip -k "$GZ_FILE"
else
    echo "[$(date)] uniref90.fasta already decompressed"
fi

echo "[$(date)] Building MMseqs2 database (this will take a while)..."
mmseqs createdb "$FASTA_FILE" "$DB_PREFIX" --threads 4

# Verify
if [[ -f "${DB_PREFIX}.dbtype" ]]; then
    DBTYPE=$(cat "${DB_PREFIX}.dbtype" | od -An -tu1 | tr -d ' ')
    echo "[$(date)] Database built successfully. dbtype=$DBTYPE"
    if [[ "$DBTYPE" == "0" ]]; then
        echo "  Type: Aminoacid (correct)"
    else
        echo "  WARNING: Expected type 0 (Aminoacid), got $DBTYPE"
    fi
else
    echo "[$(date)] ERROR: Database build failed - no .dbtype file"
    exit 1
fi

# Clean up decompressed FASTA to save space (~90 GB)
echo "[$(date)] Removing decompressed FASTA to save space..."
rm -f "$FASTA_FILE"

echo "[$(date)] Done."
echo "Database: ${DB_PREFIX}"
