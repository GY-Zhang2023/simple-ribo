#!/bin/bash
# rna_LAN_env.sh: environment activator for bulk RNA‑seq only
# Usage: source rna_LAN_env.sh

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
    echo "Please run: conda env create -f environment.yml"
    return 1
fi

if command -v rna_LAN >/dev/null 2>&1; then
    echo "[OK] Global command 'rna_LAN' available."
else
    echo "[WARNING] global command missing, re‑run bash install.sh"
fi

echo -e "\n🎉 simple‑LAN RNA‑seq environment ready."
echo "Usage example:"
echo "  rna_LAN --step1 sample.csv"
