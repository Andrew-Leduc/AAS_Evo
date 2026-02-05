#!/usr/bin/env bash
#
# Aggregate all per-MSA quality statistics into a single summary TSV.
#
# Usage:
#   bash aggregate_msa_stats.sh
#
# Output:
#   /scratch/leduc.an/AAS_Evo/MSA_STATS/msa_quality_summary.tsv
#

set -euo pipefail

DATA_DIR="/scratch/leduc.an/AAS_Evo"
STATS_DIR="${DATA_DIR}/MSA_STATS"
MSA_DIR="${DATA_DIR}/MSA"
OUTPUT="${STATS_DIR}/msa_quality_summary.tsv"

if [[ ! -d "$STATS_DIR" ]]; then
    echo "ERROR: Stats directory not found: $STATS_DIR"
    echo "Run MSA generation first."
    exit 1
fi

echo "Aggregating MSA quality statistics..."

# Header
echo -e "accession\tgene\tn_sequences\tneff\tquery_length\tavg_identity\tavg_coverage\tcolumn_occupancy_mean\tcolumn_occupancy_min\tquality_flags" > "$OUTPUT"

# Load reference for accession -> gene mapping
REF_FASTA="${DATA_DIR}/SEQ_FILES/uniprot_human_canonical.fasta"

# Process each stats file
for stats_file in "$STATS_DIR"/*.stats.json; do
    [[ -f "$stats_file" ]] || continue

    accession=$(basename "$stats_file" .stats.json)

    # Get gene name from reference FASTA
    gene=$(grep "|${accession}|" "$REF_FASTA" 2>/dev/null | head -1 | sed 's/.*GN=\([^ ]*\).*/\1/' || echo "NA")

    # Extract stats from JSON
    python3 -c "
import json
import sys

try:
    with open('$stats_file') as f:
        d = json.load(f)

    cols = [
        '$accession',
        '$gene',
        d.get('n_sequences', 'NA'),
        d.get('neff', 'NA'),
        d.get('query_length', 'NA'),
        d.get('avg_identity', 'NA'),
        d.get('avg_coverage', 'NA'),
        d['column_occupancy'].get('mean', 'NA') if 'column_occupancy' in d else 'NA',
        d['column_occupancy'].get('min', 'NA') if 'column_occupancy' in d else 'NA',
        ','.join(d.get('quality_flags', [])) or 'OK'
    ]
    print('\t'.join(str(c) for c in cols))
except Exception as e:
    print(f'$accession\t$gene\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tERROR:{e}', file=sys.stderr)
" >> "$OUTPUT" 2>/dev/null || echo -e "${accession}\t${gene}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tPARSE_ERROR" >> "$OUTPUT"

done

# Count statistics
TOTAL=$(tail -n +2 "$OUTPUT" | wc -l)
LOW_NEFF=$(tail -n +2 "$OUTPUT" | grep -c "LOW_NEFF" || echo 0)
MARGINAL_NEFF=$(tail -n +2 "$OUTPUT" | grep -c "MARGINAL_NEFF" || echo 0)
OK=$(tail -n +2 "$OUTPUT" | grep -c $'\tOK$' || echo 0)

echo ""
echo "========================================="
echo "MSA Quality Summary"
echo "========================================="
echo "Total MSAs:              $TOTAL"
echo "Good quality (OK):       $OK"
echo "Marginal Neff (<100):    $MARGINAL_NEFF"
echo "Low Neff (<30):          $LOW_NEFF"
echo ""
echo "Output: $OUTPUT"
echo ""
echo "To view problematic MSAs:"
echo "  grep -E 'LOW_NEFF|LOW_COVERAGE' $OUTPUT"
