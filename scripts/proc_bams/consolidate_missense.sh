#!/usr/bin/env bash
#
# Consolidate all missense mutation TSVs into a single table
#
# Usage:
#   bash consolidate_missense.sh
#
# Output:
#   /scratch/leduc.an/AAS_Evo/VEP/all_missense_mutations.tsv
#

set -euo pipefail

VEP_DIR="/scratch/leduc.an/AAS_Evo/VEP"
OUTPUT="$VEP_DIR/all_missense_mutations.tsv"

echo "[$(date)] Consolidating missense mutations..."

# Write header once
echo -e "sample_id\tCHROM\tPOS\tREF\tALT\tConsequence\tSYMBOL\tGene\tHGVSp\tAmino_acids\tProtein_position\tgnomADe_AF\tam_pathogenicity\tam_class\tDP\tAD_ref\tAD_alt\tVAF" > "$OUTPUT"

# Append all TSV files (skipping headers)
count=0
for tsv in "$VEP_DIR"/*.vep.tsv; do
    if [[ -f "$tsv" ]]; then
        # Skip header line, append data
        tail -n +2 "$tsv" >> "$OUTPUT"
        ((++count))
    fi
done

# Count total mutations
total=$(tail -n +2 "$OUTPUT" | wc -l)

echo "[$(date)] Done."
echo "Processed: $count samples"
echo "Total missense mutations: $total"
echo "Output: $OUTPUT"

# Summary statistics
echo ""
echo "=== Summary ==="
echo "Top 20 most frequently mutated genes:"
tail -n +2 "$OUTPUT" | cut -f7 | sort | uniq -c | sort -rn | head -20

echo ""
echo "Mutations by consequence type:"
tail -n +2 "$OUTPUT" | cut -f6 | sort | uniq -c | sort -rn
