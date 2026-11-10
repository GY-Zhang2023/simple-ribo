#!/bin/bash
# install.sh (proposed replacement for simple-ribo)
# - resolves the repo from its own location (no hardcoded paths)
# - creates the two tool dirs the launchers expect, via SYMLINKS into the repo
#   (single source of truth: edit config in repo, no copy drift)
# - installs global commands in ~/.local/bin
# usage: bash install.sh [--uninstall]

set -euo pipefail

REPO=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)
TOOLS="${REPO}/TOOLS"
LOCAL_BIN="${HOME}/.local/bin"
RIBO_HOME="${HOME}/ribo_simple_LAN"
RNA_HOME="${HOME}/rna_simple_LAN"

[[ -d "${TOOLS}" ]] || { echo "[ERROR] TOOLS missing under ${REPO}"; exit 1; }
[[ "$(uname -s)" == Linux ]] || { echo "[ERROR] Linux only"; exit 1; }

if [[ "${1:-}" == "--uninstall" ]]; then
    rm -f "${LOCAL_BIN}/ribo_LAN" "${LOCAL_BIN}/rna_LAN"
    rm -rf "${RIBO_HOME}" "${RNA_HOME}"
    echo "[OK] removed tool dirs and global commands (repo untouched)"
    exit 0
fi

mkdir -p "${LOCAL_BIN}" "${RIBO_HOME}/TOOLS" "${RNA_HOME}/TOOLS"

echo "[1/4] link ribo scripts -> ${RIBO_HOME}/TOOLS"
for f in ribo_LAN step1_filter_clean.sh step2_align.sh step3_riboWaltz.R \
         step4_A.P.E_around_motif.R step5_personalize.R step6_ribo_group_diff.R; do
    ln -sf "${TOOLS}/${f}" "${RIBO_HOME}/TOOLS/${f}"
done
ln -sf "${REPO}/env/ribo_LAN_env.sh" "${RIBO_HOME}/ribo_LAN_env.sh"

echo "[2/4] link rna scripts -> ${RNA_HOME}/TOOLS"
for f in rna_LAN RNAstep1_filter_clean.sh RNAstep2_align.sh RNAstep3_countPCA.R \
         RNAstep4_TE.R RNAstep5_TE_diff_analysis.R; do
    ln -sf "${TOOLS}/${f}" "${RNA_HOME}/TOOLS/${f}"
done
ln -sf "${REPO}/env/rna_LAN_env.sh" "${RNA_HOME}/rna_LAN_env.sh"

echo "[3/4] permissions"
chmod +x "${TOOLS}"/*.sh "${TOOLS}"/ribo_LAN "${TOOLS}"/rna_LAN "${REPO}/env"/*.sh

echo "[4/4] global commands"
ln -sf "${TOOLS}/ribo_LAN" "${LOCAL_BIN}/ribo_LAN"
ln -sf "${TOOLS}/rna_LAN"  "${LOCAL_BIN}/rna_LAN"

if ! grep -qs 'HOME/.local/bin' "${HOME}/.bashrc"; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "${HOME}/.bashrc"
    echo "[INFO] appended ~/.local/bin to ~/.bashrc (re-login or source it)"
fi

cat <<EOF

Done.
Next:
  1) create the env:      mamba env create -f ${REPO}/environment.yml
  2) activate:            source ${RIBO_HOME}/ribo_LAN_env.sh
  3) preflight check:     bash ${REPO}/tools/check_env.sh
  4) self test:           ribo_LAN -h ; rna_LAN -h
EOF
