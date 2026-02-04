#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --partition=short
#SBATCH --job-name=vep
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/vep_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/vep_%A_%a.err
#SBATCH --array=1-100%10

#
# Batch VEP annotation for all VCF files
#
# Usage:
#   # First, generate the VCF list:
#   ls /scratch/leduc.an/AAS_Evo/VCF/*.vcf.gz > /scratch/leduc.an/AAS_Evo/vcf_list.txt
#
#   # Then submit with the correct array size:
#   NUM_VCFS=$(wc -l < /scratch/leduc.an/AAS_Evo/vcf_list.txt)
#   sbatch --array=1-${NUM_VCFS}%10 submit_vep.sh
#

set -euo pipefail

# --------- PATHS ----------
DATA_DIR="/scratch/leduc.an/AAS_Evo"
VCF_LIST="${DATA_DIR}/vcf_list.txt"
SIF="/scratch/leduc.an/tools/vep/ensembl-vep.sif"
VEP_CACHE="${DATA_DIR}/SEQ_FILES/vep_cache"
FASTA="${DATA_DIR}/SEQ_FILES/hg38.fa"
OUTDIR="${DATA_DIR}/VEP"
# --------------------------

mkdir -p "$OUTDIR"
mkdir -p /scratch/leduc.an/AAS_Evo/logs

# Get the VCF file for this array task
VCF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$VCF_LIST")

if [[ -z "$VCF" ]]; then
    echo "ERROR: No VCF at line ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$VCF" ]]; then
    echo "ERROR: VCF not found: $VCF"
    exit 1
fi

# Extract sample ID from VCF filename
SAMPLE_ID=$(basename "$VCF" .vcf.gz)

VEP_OUT="$OUTDIR/${SAMPLE_ID}.vep.vcf.gz"
VEP_TSV="$OUTDIR/${SAMPLE_ID}.vep.tsv"
LOG="$OUTDIR/${SAMPLE_ID}.vep.log"
WARN="$OUTDIR/${SAMPLE_ID}.vep.warnings.txt"

# Skip if already processed
if [[ -f "$VEP_TSV" ]]; then
    echo "[$(date)] Already processed: $SAMPLE_ID"
    exit 0
fi

echo "[$(date)] Running VEP on: $SAMPLE_ID" | tee "$LOG"

# Run VEP with gnomAD frequencies
# --af_gnomad adds gnomAD allele frequencies
# --pick selects one consequence per variant (canonical preferred)
apptainer exec \
    -B "${VEP_CACHE}":/cache \
    -B "$(dirname "$FASTA")":/ref \
    -B "$(dirname "$VCF")":/input \
    -B "$OUTDIR":/output \
    "${SIF}" \
    vep \
        -i "/input/$(basename "$VCF")" \
        -o "/output/${SAMPLE_ID}.vep.vcf.gz" \
        --vcf --compress_output bgzip \
        --cache --dir_cache /cache \
        --assembly GRCh38 --offline \
        --fasta "/ref/$(basename "$FASTA")" \
        --symbol --hgvs --protein --canonical \
        --af_gnomad \
        --pick \
        --fields "Consequence,SYMBOL,Gene,BIOTYPE,HGVSc,HGVSp,Protein_position,Amino_acids,Codons,gnomAD_AF" \
        --fork 8 \
        --warning_file "/output/${SAMPLE_ID}.vep.warnings.txt" \
        2>&1 | tee -a "$LOG"

echo "[$(date)] Extracting missense variants..." | tee -a "$LOG"

# Extract missense variants to TSV
# Parse VEP annotations from VCF INFO field
{
    echo -e "sample_id\tCHROM\tPOS\tREF\tALT\tConsequence\tSYMBOL\tGene\tHGVSp\tAmino_acids\tProtein_position\tgnomAD_AF\tDP\tAD_ref\tAD_alt\tVAF"

    zcat "$VEP_OUT" 2>/dev/null | grep -v "^#" | awk -v sid="$SAMPLE_ID" -F'\t' '
    BEGIN { OFS="\t" }
    {
        chrom = $1
        pos = $2
        ref = $4
        alt = $5
        info = $8

        # Extract CSQ field from INFO
        if (match(info, /CSQ=[^;]+/)) {
            csq = substr(info, RSTART+4, RLENGTH-4)

            # Split on comma for multiple annotations, take first (--pick should give one)
            split(csq, annotations, ",")
            ann = annotations[1]

            # Split annotation by pipe
            split(ann, fields, "|")

            consequence = fields[1]
            symbol = fields[2]
            gene = fields[3]
            biotype = fields[4]
            hgvsc = fields[5]
            hgvsp = fields[6]
            protein_pos = fields[7]
            amino_acids = fields[8]
            codons = fields[9]
            gnomad_af = fields[10]

            # Only output missense variants
            if (consequence ~ /missense/) {
                # Extract DP, AD from FORMAT/sample fields if present
                dp = "NA"
                ad_ref = "NA"
                ad_alt = "NA"
                vaf = "NA"
                if (NF >= 10) {
                    # Parse FORMAT and sample columns
                    split($9, fmt, ":")
                    split($10, vals, ":")
                    for (i=1; i<=length(fmt); i++) {
                        if (fmt[i] == "DP") {
                            dp = vals[i]
                        }
                        if (fmt[i] == "AD") {
                            split(vals[i], ad, ",")
                            ad_ref = ad[1]
                            ad_alt = ad[2]
                            if (ad[1] + ad[2] > 0) {
                                vaf = ad[2] / (ad[1] + ad[2])
                            }
                        }
                    }
                }

                print sid, chrom, pos, ref, alt, consequence, symbol, gene, hgvsp, amino_acids, protein_pos, gnomad_af, dp, ad_ref, ad_alt, vaf
            }
        }
    }'
} > "$VEP_TSV"

echo "[$(date)] Done: $SAMPLE_ID" | tee -a "$LOG"
echo "VEP VCF: $VEP_OUT"
echo "Missense TSV: $VEP_TSV"
