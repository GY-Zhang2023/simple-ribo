#!/bin/bash
# check_env.sh - preflight check for simple-ribo
# verifies: conda env, CLI tools, R packages, reference files/indexes, hisat2 path
# usage: bash tools/check_env.sh

set -uo pipefail

REPO=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)
TOOLS="${REPO}/TOOLS"
FAIL=0
WARN=0

ok()   { printf "  [OK]   %s\n" "$1"; }
bad()  { printf "  [FAIL] %s\n" "$1"; FAIL=$((FAIL+1)); }
warn() { printf "  [WARN] %s\n" "$1"; WARN=$((WARN+1)); }

echo "== 1. conda environment =="
if [[ "${CONDA_DEFAULT_ENV:-}" == "ribo_simple" ]]; then
    ok "conda env 'ribo_simple' active (${CONDA_PREFIX})"
else
    bad "conda env 'ribo_simple' NOT active (launchers will exit). run: source ~/ribo_simple_LAN/ribo_LAN_env.sh"
fi

echo "== 2. command line tools =="
for t in cutadapt bowtie2 STAR samtools hisat2 stringtie qualimap bedtools fastqc trim_galore Rscript dos2unix; do
    if command -v "$t" >/dev/null 2>&1; then ok "$t -> $(command -v "$t")"; else bad "$t missing"; fi
done

echo "== 3. R packages =="
if command -v Rscript >/dev/null 2>&1; then
    Rscript -e '
pkgs <- c("data.table","ggplot2","pheatmap","ggrepel","future","future.apply","furrr",
          "dplyr","tidyr","purrr","scales","cowplot","gridExtra","RColorBrewer","R.utils",
          "factoextra","FactoMineR","ggseqlogo","BiocParallel","S4Vectors","IRanges",
          "GenomicRanges","GenomicFeatures","GenomicAlignments","Biostrings","txdbmaker",
          "DESeq2","clusterProfiler","riboWaltz","org.Hs.eg.db","org.Mm.eg.db")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) { cat("  [FAIL] missing: ", paste(miss, collapse=", "), "\n", sep=""); quit(status=1) }
cat("  [OK]   all required R packages present\n")' || FAIL=$((FAIL+1))
else
    bad "Rscript missing, skip package check"
fi

echo "== 4. reference files (parsed from the config blocks of each script) =="

# variables used inside the config blocks of RNAstep2_align.sh
HOME_REF="$HOME/reference"
HOME_TOOLS="$HOME/rna_simple_LAN/TOOLS"
HOME_BIN="$HOME/tools"

check_value() {
    local key="$1" val="$2" tag="$3" force="${4:-}"
    [[ -z "$val" ]] && { bad "${tag}.${key}: empty path"; return; }

    if [[ "$force" == "bt2" ]]; then
        if compgen -G "${val}.1.bt2*" >/dev/null; then
            ok "${tag}.${key} -> ${val}"
        else
            bad "${tag}.${key} bowtie2 index missing: ${val}.1.bt2"
        fi
        return
    fi

    case "${key}" in
        *hisat2_index*)
            warn "${tag}.${key}: skip HISAT2 index check, configured path=${val}"
            ;;
        *star_index*|*genome_index*)
            warn "${tag}.${key}: skip STAR index check, configured path=${val}"
            ;;
        *)
            if [[ -f "$val" ]]; then
                ok "${tag}.${key} -> ${val}"
            else
                bad "${tag}.${key} missing: ${val}"
            fi ;;
    esac
}


# print "key|value" for every  ASSOC["key"]="value"  assignment found in a file
parse_array() {
    local file="$1" var="$2"
    [[ -f "$file" ]] || return 1
    local line k v
    while IFS= read -r line; do
        if [[ "$line" =~ ${var}\[\"([^\"]*)\"\]=\"(.*)\"$ ]]; then
            k="${BASH_REMATCH[1]}"; v="${BASH_REMATCH[2]}"
            set +u                       # config blocks may reference not-yet-defined vars
            printf '%s|%s\n' "$k" "$(eval echo "\"$v\"")"
            set -u
        fi
    done < <(grep -E "^[[:space:]]*${var}\[\"" "$file" 2>/dev/null)
}

load_plain_vars() {
    local file="$1" line name val
    [[ -f "$file" ]] || return 0
    while IFS= read -r line; do
        if [[ "$line" =~ ^([A-Z_][A-Z0-9_]*)=\"(.*)\"$ ]]; then
            name="${BASH_REMATCH[1]}"
            val="$(eval echo "\"${BASH_REMATCH[2]}\"")"
            export "${name}=${val}"
        fi
    done < <(grep -E '^[A-Z_][A-Z0-9_]+="' "$file" 2>/dev/null)
}

check_block() {
    local file="$1" var="$2" tag="$3" force="${4:-}"
    if [[ ! -f "$file" ]]; then bad "${tag}: script not found ($(basename "$file"))"; return; fi
    local n=0
    while IFS='|' read -r k v; do
        [[ -z "${k:-}" ]] && continue
        n=$((n+1)); check_value "$k" "$v" "$tag" "$force"
    done < <(parse_array "$file" "$var")
    [[ $n -eq 0 ]] && bad "${tag}: no ${var} entries parsed in $(basename "$file")"
}

check_block "${TOOLS}/step2_align.sh"          SPECIES_CONFIG          ribo_config
check_block "${TOOLS}/step1_filter_clean.sh"   SPECIES_REF_MAP         ribo_rrna bt2
load_plain_vars "${TOOLS}/RNAstep1_filter_clean.sh"
check_block "${TOOLS}/RNAstep2_align.sh"       SPECIES_CONFIG          rna_config
check_block "${TOOLS}/RNAstep1_filter_clean.sh" SPECIES_MT_INDEX_MAP   rna_mtDNA bt2
check_block "${TOOLS}/RNAstep1_filter_clean.sh" SPECIES_RRRNA_TRNA_INDEX_MAP rna_rrna bt2

echo "== 5. rna aligner path (hardcoded in RNAstep2_align.sh) =="
if [[ -x "${HOME}/tools/hisat2/hisat2" ]]; then
    ok "~/tools/hisat2/hisat2 present"
else
    warn "~/tools/hisat2/hisat2 missing -> ln -sf \"\$(which hisat2)\" ~/tools/hisat2/hisat2  (or edit RNAstep2_align.sh HOME_BIN)"
fi

echo "== 6. resources =="
printf "  [INFO] cpu=%s  mem=%sG  free_disk_on_home=%sG\n" \
    "$(nproc)" "$(free -g | awk '/Mem:/{print $2}')" \
    "$(df -h "$HOME" 2>/dev/null | awk 'NR==2{print $4}')"

echo
if [[ $FAIL -gt 0 ]]; then
    echo "RESULT: ${FAIL} blocking issue(s), ${WARN} warning(s)."
    exit 1
fi
echo "RESULT: all checks passed (${WARN} warning(s))."
