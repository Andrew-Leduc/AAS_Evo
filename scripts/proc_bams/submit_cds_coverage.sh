#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --job-name=cds_cov
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/cds_cov_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/cds_cov_%A_%a.err
#SBATCH --array=1-100%10

#
# Compute per-sample, per-gene CDS coverage from BAM files.
#
# Reports which genes have sufficient sequencing depth for variant
# calling, so you know where mutations could be missed due to
# poor coverage in the custom proteogenomics FASTAs.
#
# Prerequisite:
#   bash make_cds_gene_bed.sh   # One-time: creates cds_genes.bed
#
# Usage:
#   NUM_BAMS=$(wc -l < /scratch/leduc.an/AAS_Evo/bam_list.txt)
#   sbatch --array=1-${NUM_BAMS}%10 submit_cds_coverage.sh
#
# Output per sample:
#   /scratch/leduc.an/AAS_Evo/coverage/{sample_id}.gene_coverage.tsv
#
# Then consolidate with:
#   bash consolidate_coverage.sh
#

set -euo pipefail

# --------- PATHS ----------
BAM_LIST="/scratch/leduc.an/AAS_Evo/bam_list.txt"
CDS_GENE_BED="/scratch/leduc.an/AAS_Evo/SEQ_FILES/cds_genes.bed"
OUTDIR="/scratch/leduc.an/AAS_Evo/coverage"
MIN_COV_PCT=80   # Genes with <80% bases at >=10x are flagged LOW_COVERAGE
# --------------------------

module load samtools || true

mkdir -p "$OUTDIR"
mkdir -p /scratch/leduc.an/AAS_Evo/logs

# Check for gene-annotated BED
if [[ ! -f "$CDS_GENE_BED" ]]; then
    echo "ERROR: Gene-annotated CDS BED not found: $CDS_GENE_BED"
    echo "       Run: bash make_cds_gene_bed.sh"
    exit 1
fi

# Get BAM for this array task
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")

if [[ -z "$BAM" ]]; then
    echo "ERROR: No BAM at line ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

SAMPLE_ID=$(basename "$(dirname "$BAM")")
GENE_TSV="$OUTDIR/${SAMPLE_ID}.gene_coverage.tsv"

# Skip if already done
if [[ -f "$GENE_TSV" ]]; then
    echo "[$(date)] Already processed: $SAMPLE_ID"
    exit 0
fi

echo "[$(date)] Computing per-gene CDS coverage: $SAMPLE_ID"

# Ensure BAM is indexed
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "[$(date)] Indexing BAM..."
    samtools index "$BAM"
fi

# samtools bedcov: reports total read bases covering each BED interval
# Output: chr  start  end  gene  total_read_bases
# We use this to compute mean coverage per interval, then aggregate by gene.
#
# For per-base threshold analysis we use samtools depth on the gene BED.
# samtools depth -a includes zero-coverage positions.

# Step 1: Per-base depth across all CDS regions (gene-annotated BED)
# -a: report all positions including zero depth
# -q 20: min mapping quality  -Q 20: min base quality (matches variant calling)
# Output: chr  pos  depth
#
# We join this with the gene BED to get per-gene, per-base depth,
# then aggregate per gene.

# Create a temporary file with per-base depths
TMP_DEPTH=$(mktemp "${OUTDIR}/${SAMPLE_ID}.depth.XXXXXX")
trap "rm -f '$TMP_DEPTH'" EXIT

samtools depth -a -b "$CDS_GENE_BED" -q 20 -Q 20 "$BAM" > "$TMP_DEPTH"

# Use bedtools intersect (if available) or awk-based approach to map
# positions back to genes. Since we have the gene BED sorted, we can
# use a streaming approach.

# Approach: for each interval in CDS_GENE_BED, extract depths for
# positions in that interval from the depth file, and tally per gene.
# This is done with awk by reading both files.

echo "[$(date)] Aggregating per-gene coverage..."

{
    echo -e "sample_id\tgene\tcds_bases\tbases_0x\tbases_ge1x\tbases_ge10x\tbases_ge20x\tmean_depth\tpct_ge10x\tstatus"

    # Read BED intervals and depth file together
    # For each BED interval, count bases at each threshold
    awk -v sid="$SAMPLE_ID" -v min_pct="$MIN_COV_PCT" '
    BEGIN { OFS = "\t" }

    # First pass: read depth file into array (chr:pos -> depth)
    # This can be memory-intensive for large genomes but CDS is ~35MB of bases
    NR == FNR {
        key = $1 ":" $2
        depth[key] = $3
        next
    }

    # Second pass: read BED file, look up depths for each position
    {
        chr = $1; start = $2; end = $3; gene = $4
        for (pos = start + 1; pos <= end; pos++) {
            key = chr ":" pos
            d = (key in depth) ? depth[key] : 0
            gene_total[gene]++
            gene_depth_sum[gene] += d
            if (d == 0) gene_0x[gene]++
            if (d >= 1) gene_1x[gene]++
            if (d >= 10) gene_10x[gene]++
            if (d >= 20) gene_20x[gene]++
        }
    }

    END {
        for (gene in gene_total) {
            total = gene_total[gene]
            b0  = gene_0x[gene] + 0
            b1  = gene_1x[gene] + 0
            b10 = gene_10x[gene] + 0
            b20 = gene_20x[gene] + 0
            md  = gene_depth_sum[gene] / total
            pct = (b10 / total) * 100

            status = (pct >= min_pct) ? "OK" : "LOW_COVERAGE"

            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\t%s\n",
                sid, gene, total, b0, b1, b10, b20, md, pct, status
        }
    }' "$TMP_DEPTH" "$CDS_GENE_BED" | sort -k2,2

} > "$GENE_TSV"

# Print summary
TOTAL_GENES=$(tail -n +2 "$GENE_TSV" | wc -l)
LOW_COV=$(tail -n +2 "$GENE_TSV" | grep -c "LOW_COVERAGE" || true)

echo "[$(date)] Done: $SAMPLE_ID"
echo "  Total genes: $TOTAL_GENES"
echo "  Low coverage (<${MIN_COV_PCT}% at >=10x): $LOW_COV"
echo "  Output: $GENE_TSV"
