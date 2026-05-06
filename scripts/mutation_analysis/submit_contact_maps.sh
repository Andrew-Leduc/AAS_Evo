#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --job-name=contact_maps
#SBATCH --output=/scratch/leduc.an/AAS_Evo/logs/contact_maps_%j.out
#SBATCH --error=/scratch/leduc.an/AAS_Evo/logs/contact_maps_%j.err

set -euo pipefail

source /projects/marubi/software/miniconda3/bin/activate
conda activate evcouplings2

cd /scratch/leduc.an/AAS_Evo/SPURS/contact_maps

python contact_make.py
