#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=short
#SBATCH --job-name=build_uniref50
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/build_uniref50_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/build_uniref50_%j.err

#
# Download UniRef50 FASTA and build MMseqs2 database.
#
# UniRef50 is a good middle ground:
#   - ~70M sequences (vs 188M for UniRef90, vs 30M for UniRef30)
#   - 12 GB compressed download, ~24 GB uncompressed
#   - Fast MSA generation with good coverage
#
# Usage:
#   sbatch scripts/setup/build_uniref50_db.sh
#

set -euo pipefail

SEQ_DIR="/scratch/leduc.an/AAS_Evo/SEQ_FILES"
GZ_FILE="${SEQ_DIR}/uniref50.fasta.gz"
FASTA_FILE="${SEQ_DIR}/uniref50.fasta"
DB_PREFIX="${SEQ_DIR}/uniref50"

mkdir -p /scratch/leduc.an/AAS_Evo/logs

# Check if database already exists
if [[ -f "${DB_PREFIX}.dbtype" ]]; then
    echo "[$(date)] UniRef50 database already exists at ${DB_PREFIX}"
    echo "  To rebuild, remove: rm ${DB_PREFIX}* ${GZ_FILE}"
    exit 0
fi

# Download
if [[ -f "$GZ_FILE" ]]; then
    echo "[$(date)] uniref50.fasta.gz already downloaded"
else
    echo "[$(date)] Downloading uniref50.fasta.gz (~12 GB)..."
    wget --continue -O "$GZ_FILE" \
        "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
    echo "[$(date)] Download complete."
fi

# Remove any broken database files from previous attempts
echo "[$(date)] Cleaning up previous database files..."
rm -f "${DB_PREFIX}" "${DB_PREFIX}".dbtype "${DB_PREFIX}".index \
      "${DB_PREFIX}".lookup "${DB_PREFIX}"_h "${DB_PREFIX}"_h.dbtype \
      "${DB_PREFIX}"_h.index

# Decompress (keep the .gz)
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "[$(date)] Decompressing uniref50.fasta.gz..."
    gunzip -k "$GZ_FILE"
else
    echo "[$(date)] uniref50.fasta already decompressed"
fi

# Build MMseqs2 database
echo "[$(date)] Building MMseqs2 database (this will take a while)..."
module load mmseqs2 2>/dev/null || true
export PATH="$HOME/bin/mmseqs/bin:$PATH"

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

# Clean up decompressed FASTA to save space (~24 GB)
echo "[$(date)] Removing decompressed FASTA to save space..."
rm -f "$FASTA_FILE"

echo "[$(date)] Done."
echo "Database: ${DB_PREFIX}"
echo ""
echo "To use: update TARGET_DB in submit_msa_generation.sh to point to ${DB_PREFIX}"
