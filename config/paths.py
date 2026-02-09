"""
Path configuration for AAS_Evo project.
Automatically detects local vs cluster environment.
"""

import os
from pathlib import Path

def get_environment():
    """Detect whether running on cluster or locally."""
    # Check for explicit override
    env = os.environ.get('AAS_ENV')
    if env:
        return env.lower()

    # Auto-detect based on path existence
    if Path('/scratch/leduc.an').exists():
        return 'cluster'
    return 'local'

ENV = get_environment()

# Base paths
REPO_DIR = Path(__file__).resolve().parent.parent

if ENV == 'cluster':
    # Cluster paths
    SCRIPTS_DIR = Path('/home/leduc.an/AAS_Evo_project/AAS_Evo')
    DATA_DIR = Path('/scratch/leduc.an/AAS_Evo')
    BAMS_DIR = DATA_DIR / 'BAMS'
    RAW_DIR = DATA_DIR / 'RAW'
else:
    # Local development paths
    SCRIPTS_DIR = REPO_DIR
    DATA_DIR = SCRIPTS_DIR / 'data'
    BAMS_DIR = DATA_DIR / 'BAMS'
    RAW_DIR = DATA_DIR / 'RAW'

# Metadata lives inside the repo (same path on all environments)
META_DIR = REPO_DIR / 'metadata'

# Derived paths
DOWNLOAD_DIR = SCRIPTS_DIR / 'scripts' / 'download'
ANALYSIS_DIR = SCRIPTS_DIR / 'scripts' / 'analysis'
UTILS_DIR = SCRIPTS_DIR / 'utils'

# GDC/PDC specific
GDC_MANIFEST_DIR = SCRIPTS_DIR / 'manifests' / 'gdc'
PDC_MANIFEST_DIR = SCRIPTS_DIR / 'manifests' / 'pdc'

def ensure_dirs():
    """Create necessary directories if they don't exist."""
    dirs = [BAMS_DIR, RAW_DIR, GDC_MANIFEST_DIR, PDC_MANIFEST_DIR]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)

def print_config():
    """Print current configuration for debugging."""
    print(f"Environment: {ENV}")
    print(f"Scripts directory: {SCRIPTS_DIR}")
    print(f"Data directory: {DATA_DIR}")
    print(f"Meta directory: {META_DIR}")
    print(f"BAMS directory: {BAMS_DIR}")
    print(f"RAW directory: {RAW_DIR}")

if __name__ == '__main__':
    print_config()
