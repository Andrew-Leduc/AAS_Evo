#!/usr/bin/env bash
# interrogate_bams.sh
#
# SLURM array job: one task per unique not-have UUID in targets.tsv.
# For each target position, runs samtools mpileup and counts REF vs ALT reads.
#
# Usage (after extract_plex_targets.py has run):
#
#   TARGETS=/scratch/leduc.an/AAS_Evo/MS_SEARCH/bam_interrogate/targets.tsv
#   N_UUIDS=$(awk 'NR>1{print $1}' "$TARGETS" | sort -u | wc -l)
#   sbatch --array=1-${N_UUIDS}%8 interrogate_bams.sh \
#       --targets "$TARGETS" \
#       --bam-dir /scratch/leduc.an/AAS_Evo/BAMS \
#       --ref     /scratch/leduc.an/AAS_Evo/SEQ_FILES/hg38.fa \
#       --out-dir /scratch/leduc.an/AAS_Evo/MS_SEARCH/bam_interrogate
#
# Output:
#   bam_support_{SLURM_ARRAY_TASK_ID}.tsv  (per-UUID, concatenated later)
#   bam_support.tsv                         (after all tasks complete via --merge)
#
# Run --merge after the array finishes:
#   bash interrogate_bams.sh --merge --out-dir /scratch/.../bam_interrogate
#
#SBATCH --job-name=interrogate_bams
#SBATCH --partition=slavov
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/interrogate_bams_%A_%a.out
#SBATCH --error=logs/interrogate_bams_%A_%a.err

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
TARGETS="/scratch/leduc.an/AAS_Evo/MS_SEARCH/bam_interrogate/targets.tsv"
BAM_DIR="/scratch/leduc.an/AAS_Evo/BAMS"
REF="/scratch/leduc.an/AAS_Evo/SEQ_FILES/hg38.fa"
OUT_DIR="/scratch/leduc.an/AAS_Evo/MS_SEARCH/bam_interrogate"
MERGE=0
MIN_MAPQ=20    # samtools -q
MIN_BQ=20      # samtools -Q

# ── Argument parsing ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --targets) TARGETS="$2"; shift 2 ;;
        --bam-dir) BAM_DIR="$2"; shift 2 ;;
        --ref)     REF="$2";     shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --merge)   MERGE=1;      shift   ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

mkdir -p "$OUT_DIR"

# ── Load samtools ─────────────────────────────────────────────────────────────
export PATH="/shared/EL9/explorer/samtools/1.21/bin:$PATH"
echo "samtools: $(which samtools)  version: $(samtools --version | head -1)"

# ── Merge mode ───────────────────────────────────────────────────────────────
if [[ "$MERGE" -eq 1 ]]; then
    echo "Merging per-UUID output files..."
    OUTFILE="$OUT_DIR/bam_support.tsv"
    # Write header from first chunk
    FIRST=$(ls "$OUT_DIR"/bam_support_*.tsv 2>/dev/null | head -1)
    if [[ -z "$FIRST" ]]; then
        echo "ERROR: No bam_support_*.tsv files found in $OUT_DIR"; exit 1
    fi
    head -1 "$FIRST" > "$OUTFILE"
    for f in "$OUT_DIR"/bam_support_*.tsv; do
        tail -n +2 "$f" >> "$OUTFILE"
    done
    N=$(tail -n +2 "$OUTFILE" | wc -l)
    echo "Merged: $OUTFILE  ($N rows)"
    exit 0
fi

# ── SLURM array task mode ────────────────────────────────────────────────────
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

echo "Task ${TASK_ID}: loading targets from $TARGETS"

# Collect unique not_have_uuids in sorted order (skip header)
mapfile -t ALL_UUIDS < <(awk -F'\t' 'NR>1{print $1}' "$TARGETS" | sort -u)
N_UUIDS=${#ALL_UUIDS[@]}
echo "  Total unique UUIDs: $N_UUIDS"

if [[ "$TASK_ID" -gt "$N_UUIDS" ]]; then
    echo "Task ID $TASK_ID exceeds number of UUIDs ($N_UUIDS). Exiting."
    exit 0
fi

# 1-based → 0-based index
UUID="${ALL_UUIDS[$((TASK_ID - 1))]}"
echo "  UUID: $UUID"

# ── Locate BAM ───────────────────────────────────────────────────────────────
BAM_SEARCH=$(find "$BAM_DIR/$UUID" -name "*.bam" 2>/dev/null | head -1 || true)
if [[ -z "$BAM_SEARCH" ]]; then
    echo "ERROR: No BAM found for UUID $UUID in $BAM_DIR/$UUID" >&2
    exit 1
fi
BAM="$BAM_SEARCH"
echo "  BAM: $BAM"

# ── Index BAM if needed ──────────────────────────────────────────────────────
if [[ ! -f "${BAM}.bai" ]] && [[ ! -f "${BAM%.bam}.bai" ]]; then
    echo "  Indexing BAM..."
    samtools index "$BAM"
fi

# ── Extract target positions for this UUID ───────────────────────────────────
# Columns: not_have_uuid(1) not_have_case(2) not_have_channel(3)
#          have_uuid(4) have_case(5) have_channel(6)
#          gene(7) accession(8) swap(9) peptide(10)
#          ratio(11) have_vaf(12) CHROM(13) POS(14) REF(15) ALT(16)
HEADER=$(head -1 "$TARGETS")
COL_UUID=$(echo "$HEADER"   | tr '\t' '\n' | grep -n "^not_have_uuid$"    | cut -d: -f1)
COL_CASE=$(echo "$HEADER"   | tr '\t' '\n' | grep -n "^not_have_case$"    | cut -d: -f1)
COL_CHAN=$(echo "$HEADER"   | tr '\t' '\n' | grep -n "^not_have_channel$" | cut -d: -f1)
COL_HUUID=$(echo "$HEADER"  | tr '\t' '\n' | grep -n "^have_uuid$"        | cut -d: -f1)
COL_HCASE=$(echo "$HEADER"  | tr '\t' '\n' | grep -n "^have_case$"        | cut -d: -f1)
COL_GENE=$(echo "$HEADER"   | tr '\t' '\n' | grep -n "^gene$"             | cut -d: -f1)
COL_ACC=$(echo "$HEADER"    | tr '\t' '\n' | grep -n "^accession$"        | cut -d: -f1)
COL_SWAP=$(echo "$HEADER"   | tr '\t' '\n' | grep -n "^swap$"             | cut -d: -f1)
COL_PEP=$(echo "$HEADER"    | tr '\t' '\n' | grep -n "^peptide$"          | cut -d: -f1)
COL_RATIO=$(echo "$HEADER"  | tr '\t' '\n' | grep -n "^ratio$"            | cut -d: -f1)
COL_HVAF=$(echo "$HEADER"   | tr '\t' '\n' | grep -n "^have_vaf$"         | cut -d: -f1)
COL_CHROM=$(echo "$HEADER"  | tr '\t' '\n' | grep -n "^CHROM$"            | cut -d: -f1)
COL_POS=$(echo "$HEADER"    | tr '\t' '\n' | grep -n "^POS$"              | cut -d: -f1)
COL_REF=$(echo "$HEADER"    | tr '\t' '\n' | grep -n "^REF$"              | cut -d: -f1)
COL_ALT=$(echo "$HEADER"    | tr '\t' '\n' | grep -n "^ALT$"              | cut -d: -f1)

# Extract rows for this UUID
TMPFILE=$(mktemp)
awk -F'\t' -v uuid="$UUID" -v c="$COL_UUID" 'NR>1 && $c==uuid' "$TARGETS" > "$TMPFILE"
N_TARGETS=$(wc -l < "$TMPFILE")
echo "  Target positions for this UUID: $N_TARGETS"

# ── Output file ──────────────────────────────────────────────────────────────
OUTFILE="$OUT_DIR/bam_support_${TASK_ID}.tsv"
printf '%s\n' "not_have_uuid\tnot_have_case\tnot_have_channel\thave_uuid\thave_case\tgene\taccession\tswap\tpeptide\tratio\thave_vaf\tCHROM\tPOS\tREF\tALT\ttotal_depth\talt_depth\talt_vaf" \
    > "$OUTFILE"

# ── Pileup and parse ─────────────────────────────────────────────────────────
parse_pileup() {
    # Args: pileup_string  alt_base
    local PILEUP_STR="$1"
    local ALT="$2"
    local total=0
    local alt_count=0

    # Count total bases (A/T/G/C/a/t/g/c/./,) — ignore indels, starts, ends
    # Strip indel notation: [+-][0-9]+[ACGTNacgtn]+
    local CLEAN
    CLEAN=$(echo "$PILEUP_STR" | sed 's/[+-][0-9]\+[ACGTNacgtn]\+//g; s/\^.//g; s/\$//g; s/[*<>]//g')

    total=$(echo "$CLEAN" | tr -cd 'ACGTacgt.,' | wc -c)
    alt_count=$(echo "$CLEAN" | tr -cd "${ALT}$(echo "$ALT" | tr '[:upper:]' '[:lower:]')" | wc -c)

    echo "$total $alt_count"
}

while IFS=$'\t' read -r -a COLS; do
    NOT_HAVE_UUID="${COLS[$((COL_UUID-1))]}"
    NOT_HAVE_CASE="${COLS[$((COL_CASE-1))]}"
    NOT_HAVE_CHAN="${COLS[$((COL_CHAN-1))]}"
    HAVE_UUID="${COLS[$((COL_HUUID-1))]}"
    HAVE_CASE="${COLS[$((COL_HCASE-1))]}"
    GENE="${COLS[$((COL_GENE-1))]}"
    ACC="${COLS[$((COL_ACC-1))]}"
    SWAP="${COLS[$((COL_SWAP-1))]}"
    PEP="${COLS[$((COL_PEP-1))]}"
    RATIO="${COLS[$((COL_RATIO-1))]}"
    HAVE_VAF="${COLS[$((COL_HVAF-1))]}"
    CHROM="${COLS[$((COL_CHROM-1))]}"
    POS="${COLS[$((COL_POS-1))]}"
    REF_BASE="${COLS[$((COL_REF-1))]}"   # reference allele (e.g. "A") — NOT the genome path
    ALT_BASE="${COLS[$((COL_ALT-1))]}"   # alt allele (e.g. "G")

    echo "  Pileup: ${CHROM}:${POS} REF=${REF_BASE} ALT=${ALT_BASE} (${GENE} ${SWAP})"

    # Run samtools mpileup at the exact position; $REF is the genome FASTA path
    PILEUP=$(samtools mpileup \
        -r "${CHROM}:${POS}-${POS}" \
        -q ${MIN_MAPQ} \
        -Q ${MIN_BQ} \
        -f "$REF" \
        "$BAM" 2>/dev/null || true)

    if [[ -z "$PILEUP" ]]; then
        # No coverage at this position
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0.0\n' \
            "$NOT_HAVE_UUID" "$NOT_HAVE_CASE" "$NOT_HAVE_CHAN" \
            "$HAVE_UUID" "$HAVE_CASE" \
            "$GENE" "$ACC" "$SWAP" "$PEP" \
            "$RATIO" "$HAVE_VAF" \
            "$CHROM" "$POS" "$REF_BASE" "$ALT_BASE" \
            >> "$OUTFILE"
        continue
    fi

    # Pileup columns: CHROM  POS  REF  DEPTH  BASES  QUALS
    PILEUP_DEPTH=$(echo "$PILEUP" | awk '{print $4}')
    PILEUP_BASES=$(echo "$PILEUP" | awk '{print $5}')

    # Count ALT occurrences in bases string
    ALT_UPPER="$ALT_BASE"
    ALT_LOWER=$(echo "$ALT_BASE" | tr '[:upper:]' '[:lower:]')
    # Remove indel notation before counting
    BASES_CLEAN=$(echo "$PILEUP_BASES" | sed 's/[+-][0-9]\+[ACGTNacgtn]\+//g; s/\^.//g; s/\$//g; s/[*<>]//g')
    ALT_DEPTH=$(echo "$BASES_CLEAN" | tr -cd "${ALT_UPPER}${ALT_LOWER}" | wc -c)
    TOTAL_DEPTH=$(echo "$BASES_CLEAN" | tr -cd 'ACGTacgt.,' | wc -c)

    if [[ "$TOTAL_DEPTH" -gt 0 ]]; then
        ALT_VAF=$(awk "BEGIN {printf \"%.4f\", ${ALT_DEPTH}/${TOTAL_DEPTH}}")
    else
        ALT_VAF="0.0"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$NOT_HAVE_UUID" "$NOT_HAVE_CASE" "$NOT_HAVE_CHAN" \
        "$HAVE_UUID" "$HAVE_CASE" \
        "$GENE" "$ACC" "$SWAP" "$PEP" \
        "$RATIO" "$HAVE_VAF" \
        "$CHROM" "$POS" "$REF_BASE" "$ALT_BASE" \
        "$TOTAL_DEPTH" "$ALT_DEPTH" "$ALT_VAF" \
        >> "$OUTFILE"

done < "$TMPFILE"

rm -f "$TMPFILE"

echo "Done. Output: $OUTFILE"
echo "  Positions processed: $N_TARGETS"
