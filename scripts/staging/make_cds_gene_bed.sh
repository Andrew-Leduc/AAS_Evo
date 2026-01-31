#!/usr/bin/env bash
#
# Create a gene-annotated CDS BED file from GENCODE GTF.
# This is a one-time setup step. The output is used by the
# coverage script to report per-gene sequencing coverage.
#
# Usage:
#   bash make_cds_gene_bed.sh
#
# Output:
#   /scratch/leduc.an/AAS_Evo/SEQ_FILES/cds_genes.bed
#   Format: chr  start  end  gene_name
#
# The GTF is downloaded if not already present.
#

set -euo pipefail

SEQ_DIR="/scratch/leduc.an/AAS_Evo/SEQ_FILES"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
GTF_GZ="${SEQ_DIR}/gencode.v46.annotation.gtf.gz"
OUT_BED="${SEQ_DIR}/cds_genes.bed"

if [[ -f "$OUT_BED" ]]; then
    echo "Gene-annotated CDS BED already exists: $OUT_BED"
    echo "  $(wc -l < "$OUT_BED") intervals"
    exit 0
fi

mkdir -p "$SEQ_DIR"

# Download GTF if needed
if [[ ! -f "$GTF_GZ" ]]; then
    echo "Downloading GENCODE GTF..."
    wget -q -O "$GTF_GZ" "$GTF_URL"
fi

echo "Extracting CDS intervals with gene names..."

# Extract CDS entries with gene_name from GTF
# Output: chr  start(0-based)  end  gene_name
zcat "$GTF_GZ" | awk -F'\t' '
$3 == "CDS" {
    # Parse gene_name from attributes field
    n = split($9, attrs, ";")
    gene = ""
    for (i = 1; i <= n; i++) {
        gsub(/^ +/, "", attrs[i])
        if (attrs[i] ~ /^gene_name/) {
            split(attrs[i], kv, "\"")
            gene = kv[2]
        }
    }
    if (gene != "") print $1, $4-1, $5, gene
}' OFS="\t" | sort -k1,1 -k2,2n -k4,4 > "$OUT_BED"

INTERVALS=$(wc -l < "$OUT_BED")
GENES=$(cut -f4 "$OUT_BED" | sort -u | wc -l)

echo "Done: $OUT_BED"
echo "  $INTERVALS CDS intervals"
echo "  $GENES unique genes"
