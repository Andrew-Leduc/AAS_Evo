#!/usr/bin/env bash
#
# Run the BAM processing pipeline on all BAMs currently in BAMS/.
#
# Generates the file lists automatically and submits the SLURM array jobs.
# Run each step after the previous one completes.
#
# Usage:
#   bash scripts/proc_bams/run_pipeline.sh variant-call chunk_00
#   bash scripts/proc_bams/run_pipeline.sh vep chunk_00
#
# The chunk name (e.g. chunk_00) is required. VCF and VEP output will be
# stored in per-chunk subdirectories: VCF/chunk_00/, VEP/chunk_00/.
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="/scratch/leduc.an/AAS_Evo"
BAMS_DIR="${DATA_DIR}/BAMS"
BAM_LIST="${DATA_DIR}/bam_list.txt"
VCF_LIST="${DATA_DIR}/vcf_list.txt"
# ---------------------------

usage() {
    echo "Usage: bash run_pipeline.sh <variant-call|vep> <chunk_name>"
    echo ""
    echo "  variant-call chunk_00   Find all BAMs in BAMS/, submit variant calling jobs"
    echo "  vep          chunk_00   Find all VCFs in VCF/chunk_00/, submit VEP annotation jobs"
    echo ""
    echo "  chunk_name is required (e.g. chunk_00, chunk_01, ...)"
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
fi

ACTION="$1"
CHUNK_NAME="$2"

# Per-chunk output directories
VCF_DIR="${DATA_DIR}/VCF/${CHUNK_NAME}"
VEP_DIR="${DATA_DIR}/VEP/${CHUNK_NAME}"

do_variant_call() {
    echo "Scanning for BAM files in ${BAMS_DIR}/ ..."

    find "$BAMS_DIR" -name "*.bam" > "$BAM_LIST"
    NUM_BAMS=$(wc -l < "$BAM_LIST" | tr -d ' ')

    if [[ "$NUM_BAMS" -eq 0 ]]; then
        echo "No BAM files found in ${BAMS_DIR}/"
        rm -f "$BAM_LIST"
        exit 1
    fi

    # Count how many already have VCFs (will be skipped)
    ALREADY_DONE=0
    while read -r bam; do
        uuid=$(basename "$(dirname "$bam")")
        if [[ -f "${VCF_DIR}/${uuid}.vcf.gz" && -f "${VCF_DIR}/${uuid}.vcf.gz.tbi" ]]; then
            ALREADY_DONE=$((ALREADY_DONE + 1))
        fi
    done < "$BAM_LIST"

    TO_PROCESS=$((NUM_BAMS - ALREADY_DONE))

    echo ""
    echo "========================================"
    echo "Variant Calling"
    echo "========================================"
    echo "BAMs found:        $NUM_BAMS"
    echo "Already processed: $ALREADY_DONE (will be skipped)"
    echo "To process:        $TO_PROCESS"
    echo "BAM list:          $BAM_LIST"
    echo ""

    if [[ "$TO_PROCESS" -eq 0 ]]; then
        echo "All BAMs already have VCF output. Nothing to do."
        echo ""
        echo "Next step: VEP annotation"
        echo "  bash scripts/proc_bams/run_pipeline.sh vep"
        exit 0
    fi

    mkdir -p "$VCF_DIR"
    mkdir -p "${DATA_DIR}/logs"

    sbatch --export="ALL,CHUNK_NAME=${CHUNK_NAME}" \
           --array=1-${NUM_BAMS}%20 "${SCRIPTS_DIR}/submit_variant_call.sh"

    echo ""
    echo "Monitor with: squeue -u \$USER -n var_call"
    echo ""
    echo "When complete, run:"
    echo "  bash scripts/proc_bams/run_pipeline.sh vep ${CHUNK_NAME}"
}

do_vep() {
    echo "Scanning for VCF files in ${VCF_DIR}/ ..."

    find "${VCF_DIR}" -maxdepth 1 -name "*.vcf.gz" > "$VCF_LIST" 2>/dev/null || true
    NUM_VCFS=$(wc -l < "$VCF_LIST" | tr -d ' ')

    if [[ "$NUM_VCFS" -eq 0 ]]; then
        echo "No VCF files found in ${VCF_DIR}/"
        rm -f "$VCF_LIST"
        exit 1
    fi

    # Count how many already have VEP output (will be skipped)
    ALREADY_DONE=0
    while read -r vcf; do
        uuid=$(basename "$vcf" .vcf.gz)
        if [[ -f "${VEP_DIR}/${uuid}.vep.tsv" ]]; then
            ALREADY_DONE=$((ALREADY_DONE + 1))
        fi
    done < "$VCF_LIST"

    TO_PROCESS=$((NUM_VCFS - ALREADY_DONE))

    echo ""
    echo "========================================"
    echo "VEP Annotation"
    echo "========================================"
    echo "VCFs found:        $NUM_VCFS"
    echo "Already processed: $ALREADY_DONE (will be skipped)"
    echo "To process:        $TO_PROCESS"
    echo "VCF list:          $VCF_LIST"
    echo ""

    if [[ "$TO_PROCESS" -eq 0 ]]; then
        echo "All VCFs already have VEP output. Nothing to do."
        echo ""
        echo "Pipeline complete for this batch. You can now:"
        echo "  1. Delete BAMs:  rm -rf ${BAMS_DIR}/*"
        echo "  2. Download next chunk and repeat"
        exit 0
    fi

    mkdir -p "$VEP_DIR"
    mkdir -p "${DATA_DIR}/logs"

    sbatch --export="ALL,CHUNK_NAME=${CHUNK_NAME}" \
           --array=1-${NUM_VCFS}%20 "${SCRIPTS_DIR}/submit_vep.sh"

    echo ""
    echo "Monitor with: squeue -u \$USER -n vep"
    echo ""
    echo "When complete, this batch is done. You can:"
    echo "  1. Delete BAMs:  rm -rf ${BAMS_DIR}/*"
    echo "  2. Download next chunk and repeat"
    echo "  3. After all chunks: bash scripts/proc_bams/consolidate_missense.sh"
}

case "$ACTION" in
    variant-call) do_variant_call ;;
    vep)          do_vep ;;
    *)            usage ;;
esac
