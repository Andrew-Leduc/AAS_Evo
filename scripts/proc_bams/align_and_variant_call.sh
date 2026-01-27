#!/usr/bin/env bash
set -euo pipefail

# --------- USER SETTINGS ----------
BAM="/scratch/leduc.an/singularity/BAMS/0043084a-5de1-41fe-bb22-af0b2223e197/d2f91b7f-3ba9-4f92-9999-e34aaf39e8b6_wxs_gdc_realn.bam"
REF="/scratch/leduc.an/singularity/hg38.fa"   # MUST be bg38/hg38 FASTA with .fai index alongside
OUTDIR="/scratch/leduc.an/singularity/variants_test"
MIN_DP=10
# ---------------------------------

mkdir -p "$OUTDIR"
base="$(basename "$BAM" .bam)"

LOG="$OUTDIR/${base}.log"
VCF_GZ="$OUTDIR/${base}.bcftools.vcf.gz"
TSV="$OUTDIR/${base}.variants.tsv"

echo "[$(date)] BAM=$BAM" | tee "$LOG"
echo "[$(date)] REF=$REF" | tee -a "$LOG"
echo "[$(date)] OUTDIR=$OUTDIR" | tee -a "$LOG"

# Load tools (adjust for your cluster)
module load samtools || true
module load bcftools || true

# Ensure BAM is indexed
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
  echo "[$(date)] Indexing BAM..." | tee -a "$LOG"
  samtools index -@ 4 "$BAM" 2>>"$LOG"
fi

# Ensure reference is indexed
if [[ ! -f "${REF}.fai" ]]; then
  echo "ERROR: Reference FASTA index not found: ${REF}.fai" | tee -a "$LOG"
  exit 1
fi

echo "[$(date)] Calling variants with bcftools..." | tee -a "$LOG"

# Call variants -> compressed VCF
# Notes:
# -Ou streams uncompressed BCF between steps (fast)
# -q/-Q filter mapping/base quality (reasonable defaults)
# -a FORMAT/AD,FORMAT/DP ensures allele depths & depth are present
bcftools mpileup \
  -f "$REF" \
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
  -Oz -o "$VCF_GZ" 2>>"$LOG"

bcftools index -t "$VCF_GZ" 2>>"$LOG"

echo "[$(date)] Writing TSV (CHROM POS REF ALT DP AD_ref AD_alt VAF)..." | tee -a "$LOG"

# Output columns:
# CHROM POS REF ALT QUAL DP AD_ref AD_alt VAF
# AD is "ref,alt" for biallelic calls; we split it safely.
{
  echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tAD_ref\tAD_alt\tVAF"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%AD]\n' "$VCF_GZ" \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      {
        # AD like "123,45" (ref,alt); handle missing/odd cases
        split($7, a, ",");
        ad_ref = (a[1] ~ /^[0-9]+$/ ? a[1] : "NA");
        ad_alt = (a[2] ~ /^[0-9]+$/ ? a[2] : "NA");
        dp = ($6 ~ /^[0-9]+$/ ? $6 : "NA");
        vaf = "NA";
        if (dp != "NA" && ad_alt != "NA" && dp > 0) vaf = ad_alt / dp;
        print $1,$2,$3,$4,$5,dp,ad_ref,ad_alt,vaf
      }'
} > "$TSV"

echo "[$(date)] Done." | tee -a "$LOG"
echo "VCF: $VCF_GZ"
echo "TSV: $TSV"
echo "LOG: $LOG"
