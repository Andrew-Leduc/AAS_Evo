#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="/scratch/leduc.an/AAS_Evo"
ENV_DIR="${DATA_DIR}/conda_envs/spurs_env"
PKG_DIR="${DATA_DIR}/conda_pkgs"
HF_CACHE="${DATA_DIR}/SPURS/hf_cache"

mkdir -p "${DATA_DIR}/conda_envs" "${PKG_DIR}" "${HF_CACHE}"

module load miniconda3 2>/dev/null || true
module load anaconda3 2>/dev/null || true

command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found"; exit 1; }

export CONDA_PKGS_DIRS="${PKG_DIR}"

if [[ ! -x "${ENV_DIR}/bin/python" ]]; then
  conda create -y -p "${ENV_DIR}" python=3.9 pip
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_DIR}"

# pip/setuptools pin: needed for old pytorch-lightning metadata used by SPURS
python -m pip install --no-cache-dir "pip==24.0" "setuptools<70" wheel

# CPU-only torch (recommended for this pipeline)
python -m pip install --no-cache-dir --index-url https://download.pytorch.org/whl/cpu \
  torch==1.12.0+cpu torchvision==0.13.0+cpu torchaudio==0.12.0+cpu

# SPURS deps + SPURS itself
python -m pip install --no-cache-dir omegaconf "hydra-core==1.2.0"
python -m pip install --no-cache-dir "pytorch-lightning==1.7.3" "torchmetrics==0.11.4"
python -m pip install --no-cache-dir "git+https://github.com/facebookresearch/esm.git"
python -m pip install --no-cache-dir spurs

python -c "import torch, spurs; print('OK torch', torch.__version__, '| spurs', spurs.__file__)"
echo "DONE: ${ENV_DIR}"
