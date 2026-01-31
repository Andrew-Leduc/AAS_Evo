#!/usr/bin/env bash
#
# Process a single chunk of BAMs through the genomics pipeline.
#
# Usage:
#   bash process_chunk.sh <chunk_number>                # Show status
#   bash process_chunk.sh <chunk_number> download       # Submit download job
#   bash process_chunk.sh <chunk_number> prep-calls     # Write bam_list.txt
#   bash process_chunk.sh <chunk_number> prep-vep       # Write vcf_list.txt
#   bash process_chunk.sh <chunk_number> cleanup        # Delete BAMs (keeps VCF/VEP)
#

set -euo pipefail

# --------- PATHS ----------
SCRIPTS_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
META_DIR="/home/leduc.an/AAS_Evo_project/AAS_Evo_meta"
DATA_DIR="/scratch/leduc.an/AAS_Evo"
CHUNK_DIR="${META_DIR}/GDC_meta/manifests/chunks"
TOKEN="${META_DIR}/GDC_meta/gdc-user-token_AL.txt"
BAMS_DIR="${DATA_DIR}/BAMS"
VCF_DIR="${DATA_DIR}/VCF"
VEP_DIR="${DATA_DIR}/VEP"
# ---------------------------

usage() {
    echo "Usage: bash process_chunk.sh <chunk_number> [download|prep-calls|prep-vep|cleanup]"
    echo ""
    echo "Subcommands:"
    echo "  (none)      Show chunk status"
    echo "  download    Submit SLURM download job for this chunk"
    echo "  prep-calls  Write bam_list.txt for variant calling"
    echo "  prep-vep    Write vcf_list.txt for VEP annotation"
    echo "  cleanup     Delete BAMs (only if VCF exists)"
    exit 1
}

if [[ $# -lt 1 ]]; then
    usage
fi

CHUNK_NUM="$1"
ACTION="${2:-status}"

# Zero-pad chunk number to match split output (e.g., 1 -> 00, 12 -> 12)
CHUNK_PAD=$(printf "%02d" "$((CHUNK_NUM - 1))")
CHUNK_MANIFEST="${CHUNK_DIR}/chunk_${CHUNK_PAD}.tsv"

if [[ ! -f "$CHUNK_MANIFEST" ]]; then
    echo "ERROR: Chunk manifest not found: $CHUNK_MANIFEST"
    echo ""
    echo "Available chunks:"
    ls "${CHUNK_DIR}"/chunk_*.tsv 2>/dev/null | while read -r f; do
        num=$(basename "$f" .tsv | sed 's/chunk_//')
        count=$(tail -n +2 "$f" | wc -l | tr -d ' ')
        echo "  Chunk $((10#$num + 1)): $count BAMs"
    done
    exit 1
fi

# Extract UUIDs from chunk manifest (column 1, skip header)
get_uuids() {
    tail -n +2 "$CHUNK_MANIFEST" | cut -f1
}

CHUNK_TOTAL=$(tail -n +2 "$CHUNK_MANIFEST" | wc -l | tr -d ' ')

# Count how many BAMs/VCFs/VEP files exist for this chunk
count_bams() {
    local count=0
    while read -r uuid; do
        if ls "${BAMS_DIR}/${uuid}/"*.bam &>/dev/null; then
            count=$((count + 1))
        fi
    done < <(get_uuids)
    echo "$count"
}

count_vcfs() {
    local count=0
    while read -r uuid; do
        if [[ -f "${VCF_DIR}/${uuid}.vcf.gz" ]]; then
            count=$((count + 1))
        fi
    done < <(get_uuids)
    echo "$count"
}

count_vep() {
    local count=0
    while read -r uuid; do
        if [[ -f "${VEP_DIR}/${uuid}.vep.tsv" ]]; then
            count=$((count + 1))
        fi
    done < <(get_uuids)
    echo "$count"
}

# ===================== SUBCOMMANDS =====================

do_status() {
    BAMS=$(count_bams)
    VCFS=$(count_vcfs)
    VEPS=$(count_vep)

    echo "========================================"
    echo "Chunk ${CHUNK_NUM}: ${CHUNK_TOTAL} BAMs"
    echo "========================================"
    echo "Manifest:    $CHUNK_MANIFEST"
    echo ""
    echo "  Downloaded:   ${BAMS}/${CHUNK_TOTAL}"
    echo "  VCFs:         ${VCFS}/${CHUNK_TOTAL}"
    echo "  VEP TSVs:     ${VEPS}/${CHUNK_TOTAL}"
    echo "  BAMs present: ${BAMS}"
    echo ""

    # Suggest next step
    if [[ "$BAMS" -eq 0 && "$VCFS" -eq 0 ]]; then
        echo "Next step: download"
        echo "  bash scripts/chunk_workflow/process_chunk.sh ${CHUNK_NUM} download"
    elif [[ "$BAMS" -gt 0 && "$VCFS" -lt "$CHUNK_TOTAL" ]]; then
        echo "Next step: prep-calls (variant calling)"
        echo "  bash scripts/chunk_workflow/process_chunk.sh ${CHUNK_NUM} prep-calls"
    elif [[ "$VCFS" -gt 0 && "$VEPS" -lt "$CHUNK_TOTAL" ]]; then
        echo "Next step: prep-vep"
        echo "  bash scripts/chunk_workflow/process_chunk.sh ${CHUNK_NUM} prep-vep"
    elif [[ "$VEPS" -gt 0 && "$BAMS" -gt 0 ]]; then
        echo "Next step: cleanup"
        echo "  bash scripts/chunk_workflow/process_chunk.sh ${CHUNK_NUM} cleanup"
    elif [[ "$VEPS" -gt 0 && "$BAMS" -eq 0 ]]; then
        echo "Status: COMPLETE (BAMs cleaned, VCF/VEP outputs retained)"
    fi
}

do_download() {
    local NUM_SUBJOBS="${3:-10}"

    if [[ ! -f "$TOKEN" ]]; then
        echo "ERROR: GDC token not found: $TOKEN"
        exit 1
    fi

    mkdir -p "${DATA_DIR}/logs"

    # Split this chunk's manifest into sub-chunks for parallel download
    local SUB_DIR="${CHUNK_DIR}/chunk_${CHUNK_PAD}_subs"
    mkdir -p "$SUB_DIR"
    rm -f "${SUB_DIR}"/sub_*

    # Adjust sub-job count if fewer BAMs than requested jobs
    if [[ "$CHUNK_TOTAL" -lt "$NUM_SUBJOBS" ]]; then
        NUM_SUBJOBS="$CHUNK_TOTAL"
    fi

    local LINES_PER_SUB=$(( (CHUNK_TOTAL + NUM_SUBJOBS - 1) / NUM_SUBJOBS ))
    local HEADER
    HEADER=$(head -1 "$CHUNK_MANIFEST")

    tail -n +2 "$CHUNK_MANIFEST" | split -l "$LINES_PER_SUB" -d -a 2 - "${SUB_DIR}/sub_"

    echo "========================================"
    echo "Download: Chunk ${CHUNK_NUM} (${CHUNK_TOTAL} BAMs)"
    echo "========================================"
    echo "Splitting into ${NUM_SUBJOBS} parallel jobs (~${LINES_PER_SUB} BAMs each)"
    echo ""

    local SUBMITTED=0
    for sub in "${SUB_DIR}"/sub_*; do
        [[ "$sub" == *.tsv ]] && continue

        local SUB_TSV="${sub}.tsv"
        { echo "$HEADER"; cat "$sub"; } > "$SUB_TSV"
        rm "$sub"

        local SUB_COUNT
        SUB_COUNT=$(tail -n +2 "$SUB_TSV" | wc -l | tr -d ' ')
        local SUB_NAME
        SUB_NAME=$(basename "$SUB_TSV")

        echo "  Submitting ${SUB_NAME} (${SUB_COUNT} BAMs)..."
        sbatch "${SCRIPTS_DIR}/scripts/download/gdc/download_chunk.sh" \
            "$SUB_TSV" "$TOKEN"

        SUBMITTED=$((SUBMITTED + 1))
        if [[ $SUBMITTED -lt $NUM_SUBJOBS ]]; then
            sleep 2
        fi
    done

    echo ""
    echo "Submitted ${SUBMITTED} download jobs."
    echo "Monitor with: squeue -u \$USER -n gdc_dl"
    echo "Check logs:   ls ${DATA_DIR}/logs/gdc_dl_*.out"
    echo ""
    echo "When complete, run:"
    echo "  bash scripts/chunk_workflow/process_chunk.sh ${CHUNK_NUM} prep-calls"
}

do_prep_calls() {
    BAM_LIST="${DATA_DIR}/bam_list.txt"

    # Back up existing bam_list.txt
    if [[ -f "$BAM_LIST" ]]; then
        cp "$BAM_LIST" "${BAM_LIST}.bak"
    fi

    echo "Scanning for downloaded BAMs in chunk ${CHUNK_NUM}..."

    local found=0
    local missing=0
    > "$BAM_LIST"

    while read -r uuid; do
        bam=$(ls "${BAMS_DIR}/${uuid}/"*.bam 2>/dev/null | head -1)
        if [[ -n "$bam" ]]; then
            echo "$bam" >> "$BAM_LIST"
            found=$((found + 1))
        else
            missing=$((missing + 1))
        fi
    done < <(get_uuids)

    echo ""
    echo "Written: ${BAM_LIST} (${found} BAMs)"
    if [[ "$missing" -gt 0 ]]; then
        echo "Warning: ${missing} BAMs not found (not yet downloaded?)"
    fi

    echo ""
    echo "Run:"
    echo "  NUM_BAMS=\$(wc -l < ${BAM_LIST})"
    echo "  sbatch --array=1-\${NUM_BAMS}%10 scripts/proc_bams/submit_variant_call.sh"
    echo ""
    echo "When complete, run:"
    echo "  bash scripts/chunk_workflow/process_chunk.sh ${CHUNK_NUM} prep-vep"
}

do_prep_vep() {
    VCF_LIST="${DATA_DIR}/vcf_list.txt"

    # Back up existing vcf_list.txt
    if [[ -f "$VCF_LIST" ]]; then
        cp "$VCF_LIST" "${VCF_LIST}.bak"
    fi

    echo "Scanning for VCFs in chunk ${CHUNK_NUM}..."

    local found=0
    local missing=0
    > "$VCF_LIST"

    while read -r uuid; do
        vcf="${VCF_DIR}/${uuid}.vcf.gz"
        if [[ -f "$vcf" ]]; then
            echo "$vcf" >> "$VCF_LIST"
            found=$((found + 1))
        else
            missing=$((missing + 1))
        fi
    done < <(get_uuids)

    echo ""
    echo "Written: ${VCF_LIST} (${found} VCFs)"
    if [[ "$missing" -gt 0 ]]; then
        echo "Warning: ${missing} VCFs not found (variant calling incomplete?)"
    fi

    echo ""
    echo "Run:"
    echo "  NUM_VCFS=\$(wc -l < ${VCF_LIST})"
    echo "  sbatch --array=1-\${NUM_VCFS}%10 scripts/proc_bams/submit_vep.sh"
    echo ""
    echo "When complete, run:"
    echo "  bash scripts/chunk_workflow/process_chunk.sh ${CHUNK_NUM} cleanup"
}

do_cleanup() {
    echo "Checking chunk ${CHUNK_NUM} for safe BAM cleanup..."
    echo ""

    local safe=0
    local incomplete=0
    local already_clean=0
    local total_bytes=0
    local incomplete_uuids=()

    while read -r uuid; do
        bam=$(ls "${BAMS_DIR}/${uuid}/"*.bam 2>/dev/null | head -1)

        if [[ -z "$bam" ]]; then
            already_clean=$((already_clean + 1))
            continue
        fi

        if [[ -f "${VCF_DIR}/${uuid}.vcf.gz" ]]; then
            # Safe: VCF exists, can delete BAM
            bam_size=$(du -sb "${BAMS_DIR}/${uuid}" 2>/dev/null | cut -f1)
            total_bytes=$((total_bytes + ${bam_size:-0}))
            safe=$((safe + 1))
        else
            incomplete=$((incomplete + 1))
            incomplete_uuids+=("$uuid")
        fi
    done < <(get_uuids)

    # Convert bytes to GB
    total_gb=$(echo "scale=1; ${total_bytes} / 1073741824" | bc 2>/dev/null || echo "?")

    echo "========================================"
    echo "Cleanup Summary: Chunk ${CHUNK_NUM}"
    echo "========================================"
    echo "  Safe to delete:     ${safe}"
    echo "  Incomplete (skip):  ${incomplete}"
    echo "  Already cleaned:    ${already_clean}"
    echo "  Space to recover:   ~${total_gb} GB"
    echo ""

    if [[ "$incomplete" -gt 0 ]]; then
        echo "WARNING: ${incomplete} BAMs have no VCF yet and will NOT be deleted:"
        for u in "${incomplete_uuids[@]}"; do
            echo "  $u"
        done
        echo ""
    fi

    if [[ "$safe" -eq 0 ]]; then
        echo "Nothing to clean up."
        exit 0
    fi

    # Confirm
    read -p "Delete ${safe} BAM directories (~${total_gb} GB)? [y/N] " confirm
    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
        echo "Aborted."
        exit 0
    fi

    # Delete
    local deleted=0
    while read -r uuid; do
        if [[ -d "${BAMS_DIR}/${uuid}" && -f "${VCF_DIR}/${uuid}.vcf.gz" ]]; then
            rm -rf "${BAMS_DIR}/${uuid}"
            deleted=$((deleted + 1))
        fi
    done < <(get_uuids)

    echo ""
    echo "Deleted ${deleted} BAM directories."
    echo ""
    echo "To process the next chunk:"
    echo "  bash scripts/chunk_workflow/process_chunk.sh $((CHUNK_NUM + 1))"
}

# ===================== DISPATCH =====================

case "$ACTION" in
    status)     do_status ;;
    download)   do_download ;;
    prep-calls) do_prep_calls ;;
    prep-vep)   do_prep_vep ;;
    cleanup)    do_cleanup ;;
    *)
        echo "Unknown action: $ACTION"
        usage
        ;;
esac
