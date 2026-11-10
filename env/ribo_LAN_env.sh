#!/bin/bash
# ribo_LAN_env.sh: environment activator for ribo‑seq & RNA‑seq pipeline
# Usage: source ribo_LAN_env.sh

CONDA_ENV_NAME="ribo_simple"

# Auto‑detect conda init script
if ! command -v conda &> /dev/null; then
    echo "[ERROR] conda command not found in current shell."
    return 1
fi

CONDA_BASE=$(conda info --base)
CONDA_INIT_PATH="${CONDA_BASE}/etc/profile.d/conda.sh"

if [ ! -f "${CONDA_INIT_PATH}" ]; then
    echo "[ERROR] Cannot find conda init script at ${CONDA_INIT_PATH}"
    echo "Check your conda installation."
    return 1
fi

source "$CONDA_INIT_PATH"

if conda info --envs | grep -q "^$CONDA_ENV_NAME\s"; then
    conda activate "$CONDA_ENV_NAME"
    echo "[OK] Conda environment '$CONDA_ENV_NAME' activated."
else
    echo "[ERROR] Conda environment '$CONDA_ENV_NAME' not found."
    echo "Please run: conda env create -f environment.yml at repository root."
    return 1
fi

# Check global command symlink
if command -v ribo_LAN >/dev/null 2>&1 && command -v rna_LAN >/dev/null 2>&1; then
    echo "[OK] Global commands 'rna_LAN' 'ribo_LAN' available."
else
    echo "[WARNING] global commands missing, please re‑run bash install.sh under repository root."
fi

echo -e "\n🎉 simple‑LAN environment ready."
echo "Usage examples:"
echo "  rna_LAN --step1 sample.csv"
echo "  ribo_LAN --step1 ribo_sample.csv"
