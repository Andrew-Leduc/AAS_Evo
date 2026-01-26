# AAS_Evo

Genomics and proteomics data analysis pipeline for downloading and processing data from GDC (Genomic Data Commons) and PDC (Proteomics Data Commons) portals.

## Overview

This repository contains scripts for:
- Downloading genomics data (BAM files) from GDC
- Downloading proteomics data (RAW files) from PDC
- Analysis pipelines for multi-omics integration

## Directory Structure

### Local Development
```
AAS_Evo/
├── config/
│   └── paths.py          # Path configuration (local vs cluster)
├── scripts/
│   ├── download/         # Data download scripts
│   │   ├── gdc/          # GDC download utilities
│   │   └── pdc/          # PDC download utilities
│   └── analysis/         # Analysis scripts
├── utils/                # Shared utilities
└── README.md
```

### Cluster Layout
```
/home/leduc.an/AAS_Evo/           # Scripts (this repo)
/scratch/leduc.an/AAS_Evo/
├── BAMS/                          # Genomics data (BAM files)
└── RAW/                           # Proteomics data (RAW files)
```

## Setup

### Local Development
```bash
git clone https://github.com/Andrew-Leduc/AAS_Evo.git
cd AAS_Evo
pip install -r requirements.txt
```

### Cluster Deployment
```bash
cd /home/leduc.an
git clone https://github.com/Andrew-Leduc/AAS_Evo.git
cd AAS_Evo
pip install -r requirements.txt

# Create data directories
mkdir -p /scratch/leduc.an/AAS_Evo/BAMS
mkdir -p /scratch/leduc.an/AAS_Evo/RAW
```

## Configuration

The project automatically detects whether it's running locally or on the cluster. To override, set the environment variable:

```bash
export AAS_ENV=cluster  # or 'local'
```

## Usage

```bash
# Download GDC data
python scripts/download/gdc/download.py --manifest manifest.txt

# Download PDC data
python scripts/download/pdc/download.py --manifest manifest.txt
```

## Requirements

- Python 3.8+
- gdc-client (for GDC downloads)
- PDC API access

## License

MIT License
