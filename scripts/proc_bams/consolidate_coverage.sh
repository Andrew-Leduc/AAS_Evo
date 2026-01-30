#!/usr/bin/env bash
#
# Consolidate per-sample gene coverage into summary reports.
#
# Usage:
#   bash consolidate_coverage.sh
#
# Output:
#   coverage/all_gene_coverage.tsv       - All per-sample per-gene stats
#   coverage/sample_summary.tsv          - Per-sample overall CDS coverage
#   coverage/poorly_covered_genes.tsv    - Genes with low coverage across samples
#

set -euo pipefail

COV_DIR="/scratch/leduc.an/AAS_Evo/coverage"
ALL_TSV="$COV_DIR/all_gene_coverage.tsv"
SAMPLE_SUMMARY="$COV_DIR/sample_summary.tsv"
POOR_GENES="$COV_DIR/poorly_covered_genes.tsv"

echo "[$(date)] Consolidating gene coverage reports..."

# 1. Merge all per-sample gene coverage TSVs
echo -e "sample_id\tgene\tcds_bases\tbases_0x\tbases_ge1x\tbases_ge10x\tbases_ge20x\tmean_depth\tpct_ge10x\tstatus" > "$ALL_TSV"

count=0
for tsv in "$COV_DIR"/*.gene_coverage.tsv; do
    if [[ -f "$tsv" ]]; then
        tail -n +2 "$tsv" >> "$ALL_TSV"
        ((count++))
    fi
done

echo "  Merged $count sample files -> $ALL_TSV"

# 2. Per-sample summary (overall CDS coverage)
{
    echo -e "sample_id\ttotal_genes\tlow_coverage_genes\ttotal_cds_bases\tbases_ge10x\tpct_ge10x_overall\tmean_depth_overall"

    awk -F'\t' 'NR > 1 {
        sid = $1
        sample_genes[sid]++
        if ($10 == "LOW_COVERAGE") sample_low[sid]++
        sample_cds[sid] += $3
        sample_10x[sid] += $6
        sample_depth_sum[sid] += $8 * $3  # weighted mean depth
    }
    END {
        for (sid in sample_genes) {
            total = sample_cds[sid]
            b10 = sample_10x[sid]
            pct = (total > 0) ? (b10 / total) * 100 : 0
            md = (total > 0) ? sample_depth_sum[sid] / total : 0
            low = sample_low[sid] + 0
            printf "%s\t%d\t%d\t%d\t%d\t%.1f\t%.1f\n",
                sid, sample_genes[sid], low, total, b10, pct, md
        }
    }' "$ALL_TSV" | sort -k6,6n

} > "$SAMPLE_SUMMARY"

# 3. Genes frequently poorly covered across samples
{
    echo -e "gene\tsamples_total\tsamples_low_coverage\tpct_samples_low\tmedian_pct_ge10x"

    awk -F'\t' 'NR > 1 {
        gene = $2
        gene_count[gene]++
        if ($10 == "LOW_COVERAGE") gene_low[gene]++
        # Collect pct_ge10x values for median calculation
        gene_pcts[gene] = gene_pcts[gene] "," $9
    }
    END {
        for (gene in gene_count) {
            total = gene_count[gene]
            low = gene_low[gene] + 0
            pct_low = (low / total) * 100

            # Only report genes that are low in >20% of samples
            if (pct_low < 20) continue

            # Calculate median pct_ge10x
            n = split(gene_pcts[gene], vals, ",")
            # Remove empty first element from leading comma
            j = 0
            for (i = 1; i <= n; i++) {
                if (vals[i] != "") {
                    j++
                    arr[j] = vals[i] + 0
                }
            }
            # Simple sort for median
            for (a = 1; a <= j; a++)
                for (b = a+1; b <= j; b++)
                    if (arr[a] > arr[b]) { tmp = arr[a]; arr[a] = arr[b]; arr[b] = tmp }
            median = (j % 2 == 1) ? arr[int(j/2)+1] : (arr[j/2] + arr[j/2+1]) / 2

            printf "%s\t%d\t%d\t%.1f\t%.1f\n", gene, total, low, pct_low, median
        }
    }' "$ALL_TSV" | sort -k4,4rn

} > "$POOR_GENES"

# Print summary
TOTAL_SAMPLES=$(tail -n +2 "$SAMPLE_SUMMARY" | wc -l)
TOTAL_POOR=$(tail -n +2 "$POOR_GENES" | wc -l)

echo ""
echo "========================================="
echo "Coverage Consolidation Summary"
echo "========================================="
echo "Samples:  $TOTAL_SAMPLES"
echo "Genes frequently poorly covered (>20% of samples): $TOTAL_POOR"
echo ""
echo "Output files:"
echo "  $ALL_TSV"
echo "  $SAMPLE_SUMMARY"
echo "  $POOR_GENES"
echo ""

if [[ $TOTAL_POOR -gt 0 ]]; then
    echo "Top 20 most frequently poorly covered genes:"
    head -21 "$POOR_GENES" | column -t -s$'\t'
fi
