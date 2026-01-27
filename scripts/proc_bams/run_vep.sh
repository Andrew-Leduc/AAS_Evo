#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --job-name=test_align
#SBATCH --mem=16G
#SBATCH --partition=short

SIF=/scratch/leduc.an/tools/vep/ensembl-vep.sif
VEP_CACHE=/scratch/leduc.an/tools/vep/cache
FASTA=/scratch/leduc.an/singularity/hg38.fa   # <-- change to your real path

LOG=d2f91b7f-3ba9-4f92-9999-e34aaf39e8b6.vep.runlog.txt
WARN=d2f91b7f-3ba9-4f92-9999-e34aaf39e8b6.vep.warnings.txt

apptainer exec -B "${VEP_CACHE}":/cache -B "$(dirname "$FASTA")":/ref "${SIF}" \
  vep -i d2f91b7f-3ba9-4f92-9999-e34aaf39e8b6_wxs_gdc_realn.bcftools.vcf.gz \
      -o d2f91b7f-3ba9-4f92-9999-e34aaf39e8b6_wxs_gdc_realn.vep.vcf.gz \
      --vcf --cache --dir_cache /cache --assembly GRCh38 --offline \
      --fasta /ref/$(basename "$FASTA") \
      --symbol --hgvs --protein --canonical --fork 8 \
      --stats_text --warning_file "$WARN" --verbose \
      2>&1 | tee "$LOG"