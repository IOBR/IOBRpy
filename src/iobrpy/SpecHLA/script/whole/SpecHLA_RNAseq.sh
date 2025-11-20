#!/usr/bin/env bash
# SpecHLA_RNAseq.sh
# Thin wrapper around SpecHLA.sh for RNA-seq / exon typing.
# Keeps the original SpecHLA.sh CLI (e.g. -n, -1, -2, -o, -j),
# and only appends "-u 1" to force RNA/exon mode.
# It also tries to build SpecHap and ExtractHAIRs if they are missing.

set -euo pipefail

# Directory of this script: .../SpecHLA/script/whole
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# SpecHLA project root: .../SpecHLA
SPEC_HLA_ROOT="$( cd "${SCRIPT_DIR}/../.." && pwd )"
SPEC_HLA_MAIN="${SCRIPT_DIR}/SpecHLA.sh"

ensure_spec_tools() {
    # Ensure SpecHap and ExtractHAIRs are built.
    local root="$1"
    local spechap_bin="${root}/bin/SpecHap/build/SpecHap"
    local extract_bin="${root}/bin/extractHairs/build/ExtractHAIRs"

    if [[ -x "${spechap_bin}" && -x "${extract_bin}" ]]; then
        echo "[SpecHLA_RNAseq] Found SpecHap and ExtractHAIRs in ${root}/bin"
        return 0
    fi

    echo "[SpecHLA_RNAseq] SpecHap or ExtractHAIRs not found, trying to build them..."

    # Prefer a dedicated installer if it exists (file exists is enough)
    if [[ -f "${root}/install_spechap.sh" ]]; then
        (
            cd "${root}"
            bash ./install_spechap.sh
        )
    # Fallback: use the main index script if available
    elif [[ -f "${root}/index.sh" ]]; then
        (
            cd "${root}"
            bash ./index.sh
        )
    else
        echo "[SpecHLA_RNAseq] ERROR: install_spechap.sh or index.sh not found in ${root}." >&2
        echo "[SpecHLA_RNAseq] Please build SpecHap and ExtractHAIRs manually." >&2
        return 1
    fi

    # Re-check after build attempt
    if [[ ! -x "${spechap_bin}" || ! -x "${extract_bin}" ]]; then
        echo "[SpecHLA_RNAseq] ERROR: Failed to build SpecHap or ExtractHAIRs." >&2
        echo "[SpecHLA_RNAseq] You may need to run install_spechap.sh or index.sh manually and check the build logs." >&2
        return 1
    fi

    echo "[SpecHLA_RNAseq] Successfully built SpecHap and ExtractHAIRs."
}

# If user only wants help, do not try to build anything, just delegate.
if [[ $# -eq 1 && ( "$1" = "-h" || "$1" = "--help" ) ]]; then
    if [[ ! -f "${SPEC_HLA_MAIN}" ]]; then
        echo "[SpecHLA_RNAseq] ERROR: Cannot find SpecHLA.sh at ${SPEC_HLA_MAIN}" >&2
        exit 1
    fi
    exec bash "${SPEC_HLA_MAIN}" "$1"
fi

# Check that the main SpecHLA.sh exists
if [[ ! -f "${SPEC_HLA_MAIN}" ]]; then
    echo "[SpecHLA_RNAseq] ERROR: Cannot find SpecHLA.sh at ${SPEC_HLA_MAIN}" >&2
    exit 1
fi

# Try to ensure SpecHap / ExtractHAIRs are available.
# This is safe even if you do not use PacBio/Nanopore/HiC/10x,
# it only makes sure the binaries exist for scripts like link_fragment.py.
ensure_spec_tools "${SPEC_HLA_ROOT}"

echo "[SpecHLA_RNAseq] Delegating to SpecHLA.sh in RNA/exon mode (-u 1)"
echo "[SpecHLA_RNAseq] SpecHLA.sh: ${SPEC_HLA_MAIN}"
echo "[SpecHLA_RNAseq] SpecHLA root: ${SPEC_HLA_ROOT}"

# Keep the original arguments exactly as they are,
# and append "-u 1" to force the RNA/exon pipeline type.
exec bash "${SPEC_HLA_MAIN}" "$@" -u 1