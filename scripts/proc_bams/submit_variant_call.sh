#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --partition=short
#SBATCH --job-name=var_call
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/var_call_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/var_call_%A_%a.err
#SBATCH --array=1-100%10

#
# Batch variant calling for all BAM files
#
# Usage:
#   # First, generate the BAM list:
#   find /scratch/leduc.an/AAS_Evo/BAMS -name "*.bam" > /scratch/leduc.an/AAS_Evo/bam_list.txt
#
#   # Then submit with the correct array size:
#   NUM_BAMS=$(wc -l < /scratch/leduc.an/AAS_Evo/bam_list.txt)
#   sbatch --array=1-${NUM_BAMS}%10 submit_variant_call.sh
#
# The %10 limits to 10 concurrent jobs.
#

set -euo pipefail

# --------- PATHS ----------
BAM_LIST="/scratch/leduc.an/AAS_Evo/bam_list.txt"
REF="/scratch/leduc.an/AAS_Evo/SEQ_FILES/hg38.fa"
CDS_BED="/scratch/leduc.an/AAS_Evo/SEQ_FILES/cds.chr.bed"  # Coding regions only (speeds up ~10x)
OUTDIR="/scratch/leduc.an/AAS_Evo/VCF"
MIN_DP=10
# --------------------------

# Load modules
module load samtools || true
module load bcftools || true

# Create output directory
mkdir -p "$OUTDIR"
mkdir -p /scratch/leduc.an/AAS_Evo/logs

# Get the BAM file for this array task
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")

if [[ -z "$BAM" ]]; then
    echo "ERROR: No BAM file at line ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM"
    exit 1
fi

# Extract sample ID from BAM path (assumes GDC structure: .../uuid/filename.bam)
# Use the parent directory name (UUID) as sample ID
SAMPLE_ID=$(basename "$(dirname "$BAM")")
base=$(basename "$BAM" .bam)

LOG="$OUTDIR/${SAMPLE_ID}.log"
VCF_GZ="$OUTDIR/${SAMPLE_ID}.vcf.gz"
VCF_TMP="$OUTDIR/${SAMPLE_ID}.vcf.gz.tmp"
TSV="$OUTDIR/${SAMPLE_ID}.variants.tsv"

# Skip if already processed
if [[ -f "$VCF_GZ" && -f "${VCF_GZ}.tbi" ]]; then
    echo "[$(date)] Already processed: $SAMPLE_ID"
    exit 0
fi

# Clean up any partial output from a previous failed run
rm -f "$VCF_GZ" "${VCF_GZ}.tbi" "$VCF_TMP" "$TSV"

echo "[$(date)] Processing: $SAMPLE_ID" | tee "$LOG"
echo "[$(date)] BAM: $BAM" | tee -a "$LOG"
echo "[$(date)] Task ID: ${SLURM_ARRAY_TASK_ID}" | tee -a "$LOG"

# Ensure BAM is indexed
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "[$(date)] Indexing BAM..." | tee -a "$LOG"
    samtools index -@ 2 "$BAM" 2>>"$LOG"
fi

# Check reference
if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: Reference index not found: ${REF}.fai" | tee -a "$LOG"
    exit 1
fi

# Check CDS BED file
if [[ ! -f "$CDS_BED" ]]; then
    echo "ERROR: CDS BED file not found: $CDS_BED" | tee -a "$LOG"
    exit 1
fi

echo "[$(date)] Calling variants with bcftools (CDS regions only)..." | tee -a "$LOG"

# Call variants -> compressed VCF
# -R restricts to CDS regions only (much faster for WXS data)
# Write to temp file first; move to final path only on success
bcftools mpileup \
    -f "$REF" \
    -R "$CDS_BED" \
    -Ou \
    -q 20 -Q 20 \
    -a FORMAT/AD,FORMAT/DP \
    "$BAM" 2>>"$LOG" \
| bcftools call \
    -mv \
    -Ou 2>>"$LOG" \
| bcftools filter \
    -Ou \
    -i "QUAL>=20 && FORMAT/DP>=${MIN_DP}" 2>>"$LOG" \
| bcftools view \
    -Oz -o "$VCF_TMP" 2>>"$LOG"

bcftools index -t "$VCF_TMP" 2>>"$LOG"
mv "$VCF_TMP" "$VCF_GZ"
mv "${VCF_TMP}.tbi" "${VCF_GZ}.tbi"

echo "[$(date)] Writing TSV..." | tee -a "$LOG"

# Output TSV with VAF
{
    echo -e "sample_id\tCHROM\tPOS\tREF\tALT\tQUAL\tDP\tAD_ref\tAD_alt\tVAF"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%AD]\n' "$VCF_GZ" \
    | awk -v sid="$SAMPLE_ID" -F'\t' 'BEGIN{OFS="\t"}
        {
            split($7, a, ",");
            ad_ref = (a[1] ~ /^[0-9]+$/ ? a[1] : "NA");
            ad_alt = (a[2] ~ /^[0-9]+$/ ? a[2] : "NA");
            dp = ($6 ~ /^[0-9]+$/ ? $6 : "NA");
            vaf = "NA";
            if (dp != "NA" && ad_alt != "NA" && dp > 0) vaf = ad_alt / dp;
            print sid,$1,$2,$3,$4,$5,dp,ad_ref,ad_alt,vaf
        }'
} > "$TSV"

echo "[$(date)] Done: $SAMPLE_ID" | tee -a "$LOG"
echo "VCF: $VCF_GZ"
echo "TSV: $TSV"
