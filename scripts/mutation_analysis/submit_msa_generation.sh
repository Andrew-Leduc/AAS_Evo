#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --job-name=msa_gen
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/msa_gen_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/msa_gen_%A_%a.err

#
# Generate MSAs via MMseqs2 for coevolution analysis.
# Pure bash implementation - no Python overhead.
#
# Usage:
#   NUM_GENES=$(wc -l < /scratch/leduc.an/AAS_Evo/ANALYSIS/gene_list_for_msa.txt)
#   sbatch --array=1-${NUM_GENES}%10 submit_msa_generation.sh
#

set -euo pipefail

DATA_DIR="/scratch/leduc.an/AAS_Evo"

# Gene list: prefer filter_and_rank.py output, fall back to legacy location
if [[ -f "${DATA_DIR}/ANALYSIS/gene_list_for_msa.txt" ]]; then
    GENE_LIST="${DATA_DIR}/ANALYSIS/gene_list_for_msa.txt"
else
    GENE_LIST="${DATA_DIR}/gene_list.txt"
fi
REF_FASTA="${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta"
TARGET_DB="${DATA_DIR}/SEQ_FILES/uniref50"
MSA_DIR="${DATA_DIR}/MSA"
TMP_DIR="${DATA_DIR}/tmp/msa_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# MMseqs2 parameters
THREADS=4
NUM_ITERS=2
SENSITIVITY=5.7

module load mmseqs2 2>/dev/null || true
export PATH="$HOME/bin/mmseqs/bin:$PATH"

mkdir -p "$MSA_DIR"
mkdir -p "${DATA_DIR}/logs"

echo "========================================="
echo "MSA Generation - Task ${SLURM_ARRAY_TASK_ID}"
echo "========================================="

# Get gene symbol from list
GENE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENE_LIST")
if [[ -z "$GENE" ]]; then
    echo "ERROR: No gene at index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
echo "Gene: $GENE"

# Extract protein from reference FASTA using grep + awk
# UniProt headers: >sp|ACCESSION|ENTRY_HUMAN ... GN=GENE ...
HEADER=$(grep "GN=${GENE} " "$REF_FASTA" | head -1)
if [[ -z "$HEADER" ]]; then
    # Try without trailing space (gene at end of line)
    HEADER=$(grep "GN=${GENE}$" "$REF_FASTA" | head -1)
fi
if [[ -z "$HEADER" ]]; then
    echo "WARNING: Gene $GENE not found in reference FASTA. Skipping."
    exit 0
fi

# Parse accession from header: >sp|P04637|TP53_HUMAN ...
ACCESSION=$(echo "$HEADER" | sed 's/^>sp|//' | sed 's/^>tr|//' | cut -d'|' -f1)
echo "Accession: $ACCESSION"

# Check if MSA already exists
OUTPUT_A3M="${MSA_DIR}/${ACCESSION}.a3m"
if [[ -f "$OUTPUT_A3M" ]]; then
    NSEQS=$(grep -c "^>" "$OUTPUT_A3M" || echo 0)
    echo "Already exists: $OUTPUT_A3M ($NSEQS sequences)"
    exit 0
fi

# Extract sequence using awk (handles multi-line FASTA)
mkdir -p "$TMP_DIR"
QUERY_FASTA="${TMP_DIR}/query.fasta"
awk -v gene="GN=${GENE} " -v gene2="GN=${GENE}$" '
    /^>/ {
        if (match($0, gene) || match($0, gene2)) {
            printing=1; print; next
        } else {
            printing=0
        }
    }
    printing { print }
' "$REF_FASTA" > "$QUERY_FASTA"

if [[ ! -s "$QUERY_FASTA" ]]; then
    echo "ERROR: Failed to extract sequence for $GENE"
    rm -rf "$TMP_DIR"
    exit 1
fi

SEQ_LEN=$(grep -v "^>" "$QUERY_FASTA" | tr -d '\n' | wc -c)
echo "Sequence length: $SEQ_LEN aa"

# MMseqs2 pipeline
QUERY_DB="${TMP_DIR}/queryDB"
RESULT_DB="${TMP_DIR}/resultDB"
MSA_DB="${TMP_DIR}/msaDB"
MSA_OUT="${TMP_DIR}/msa_out"
MMSEQS_TMP="${TMP_DIR}/mmseqs_tmp"

mkdir -p "$MSA_OUT" "$MMSEQS_TMP"

echo ""
echo "[$(date)] Creating query database..."
mmseqs createdb "$QUERY_FASTA" "$QUERY_DB"

echo "[$(date)] Searching target DB (iters=$NUM_ITERS, sens=$SENSITIVITY)..."
mmseqs search "$QUERY_DB" "$TARGET_DB" "$RESULT_DB" "$MMSEQS_TMP" \
    --num-iterations "$NUM_ITERS" \
    -s "$SENSITIVITY" \
    --threads "$THREADS"

echo "[$(date)] Converting results to A3M..."
mmseqs result2msa "$QUERY_DB" "$TARGET_DB" "$RESULT_DB" "$MSA_DB" \
    --msa-format-mode 6

echo "[$(date)] Unpacking MSA database..."
mmseqs unpackdb "$MSA_DB" "$MSA_OUT" --unpack-name-mode 0

# Copy output
OUT_FILES=($(ls "$MSA_OUT" | grep -v '\.index$' || true))
if [[ ${#OUT_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No MSA output produced"
    rm -rf "$TMP_DIR"
    exit 1
fi

cp "${MSA_OUT}/${OUT_FILES[0]}" "$OUTPUT_A3M"
NSEQS=$(grep -c "^>" "$OUTPUT_A3M" || echo 0)
echo ""
echo "[$(date)] Done."
echo "Output: $OUTPUT_A3M"
echo "Sequences: $NSEQS"

# Cleanup
rm -rf "$TMP_DIR"
