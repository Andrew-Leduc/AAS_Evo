#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --partition=short
#SBATCH --job-name=contact_maps
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/contact_maps_%A_%a.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/contact_maps_%A_%a.err

# Regenerate AlphaFold contact maps (min heavy-atom residue-residue distance).
# Reconstructed after the original contact_make.py was lost with scratch.
# Submit as a sharded array, e.g. 20 shards:
#   mkdir -p /scratch/leduc.an/AAS_Evo/logs
#   NSHARDS=20 sbatch --array=1-20%20 scripts/mutation_analysis/submit_contact_maps.sh
# Each task processes accs[(task_id-1)::NSHARDS]. Resumable (skips existing .npy).

set -euo pipefail
module load python/3.13.5

REPO_DIR="${REPO_DIR:-/home/leduc.an/AAS_Evo_project/AAS_Evo}"
cd "$REPO_DIR"

NSHARDS="${NSHARDS:-1}"

python3 scripts/mutation_analysis/make_contact_maps.py \
    --missense /projects/slavov/andrew/AAS_EVO/all_missense_mutations.tsv \
    --ref-fasta /scratch/leduc.an/AAS_Evo/SEQ_FILES/uniprot_human_canonical.fasta \
    --out-dir /scratch/leduc.an/AAS_Evo/SPURS/contact_maps \
    --pdb-cache /scratch/leduc.an/AAS_Evo/SPURS/pdb_cache \
    --nshards "$NSHARDS"
