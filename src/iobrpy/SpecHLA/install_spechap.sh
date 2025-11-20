#!/usr/bin/env bash
set -euo pipefail

# Helper script to build SpecHap and ExtractHAIRs for SpecHLA.
#  - Use the directory of this script as the SpecHLA root (relative path, no hard-coded absolute paths).
#  - Prefer the current conda environment prefix for CMAKE_PREFIX_PATH.
#  - Automatically create build directories.
#  - Use multi-core parallel compilation when possible.

echo "[install_spechap] Starting SpecHap / ExtractHAIRs build..."

# -----------------------------
# Resolve SpecHLA root directory (directory of this script)
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPEC_HLA_ROOT="${SCRIPT_DIR}"
echo "[install_spechap] SPEC_HLA_ROOT: ${SPEC_HLA_ROOT}"

# -----------------------------
# Determine CMAKE_PREFIX_PATH (from current conda environment)
# -----------------------------
CMAKE_PREFIX=""
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    CMAKE_PREFIX="${CONDA_PREFIX}"
    echo "[install_spechap] CONDA_PREFIX: ${CMAKE_PREFIX}"
else
    echo "[install_spechap] WARNING: CONDA_PREFIX is not set; CMake will use system paths only."
fi

# -----------------------------
# Determine make parallel jobs
# -----------------------------
if [[ -n "${MAKE_JOBS:-}" ]]; then
    JOBS="${MAKE_JOBS}"
else
    if command -v nproc >/dev/null 2>&1; then
        JOBS="$(nproc)"
    else
        JOBS=4
    fi
fi
echo "[install_spechap] Using make -j ${JOBS}"

# -----------------------------
# Helper: wrap cmake invocation
# -----------------------------
run_cmake() {
    if [[ -n "${CMAKE_PREFIX}" ]]; then
        cmake .. -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX}"
    else
        cmake ..
    fi
}

# -----------------------------
# Helper: ensure ARPACK is available
# -----------------------------
ensure_arpack() {
    echo "[install_spechap] Checking for ARPACK library..."

    # ARPACK is expected to be installed into the current conda env:
    #   ${CONDA_PREFIX}/lib/libarpack.*
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
        echo "[install_spechap] WARNING: CONDA_PREFIX is not set; cannot auto-install ARPACK."
        echo "[install_spechap]          If SpecHap fails to link with ssaupd_/dseupd_ symbols,"
        echo "[install_spechap]          please manually install ARPACK, e.g.:"
        echo "              conda install -c conda-forge arpack"
        return
    fi

    local lib_dir="${CONDA_PREFIX}/lib"
    local has_arpack="false"

    # Use a glob test instead of plain 'ls' to avoid exiting under 'set -e' when no files match.
    if compgen -G "${lib_dir}/libarpack.*" >/dev/null 2>&1; then
        has_arpack="true"
    fi

    if [[ "${has_arpack}" == "true" ]]; then
        echo "[install_spechap] Found ARPACK library under ${lib_dir}; skip installation."
        return
    fi

    echo "[install_spechap] ARPACK library not detected under ${lib_dir}."

    # Locate conda executable (prefer CONDA_EXE, fall back to 'conda' in PATH)
    local conda_exe=""
    if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
        conda_exe="${CONDA_EXE}"
    elif command -v conda >/dev/null 2>&1; then
        conda_exe="$(command -v conda)"
    fi

    if [[ -z "${conda_exe}" ]]; then
        echo "[install_spechap] WARNING: 'conda' command not found; cannot auto-install ARPACK."
        echo "[install_spechap]          Please install it manually in this environment, e.g.:"
        echo "              conda install -c conda-forge arpack"
        return
    fi

    echo "[install_spechap] Installing ARPACK via conda (arpack from conda-forge)..."
    "${conda_exe}" install -y -c conda-forge arpack

    # Re-check after installation
    if compgen -G "${lib_dir}/libarpack.*" >/dev/null 2>&1; then
        echo "[install_spechap] ARPACK successfully installed under ${lib_dir}."
    else
        echo "[install_spechap] WARNING: ARPACK still not found after conda install."
        echo "[install_spechap]          SpecHap link may still fail with undefined ARPACK symbols."
    fi
}

# -----------------------------
# Ensure ARPACK before building SpecHap / ExtractHAIRs
# -----------------------------
ensure_arpack

# -----------------------------
# Build SpecHap
# -----------------------------
echo "[install_spechap] Building SpecHap..."

SPECHAP_SRC="${SPEC_HLA_ROOT}/bin/SpecHap"
SPECHAP_BUILD="${SPECHAP_SRC}/build"

if [[ ! -d "${SPECHAP_SRC}" ]]; then
    echo "[install_spechap] ERROR: SpecHap source directory not found at ${SPECHAP_SRC}"
    exit 1
fi

mkdir -p "${SPECHAP_BUILD}"
cd "${SPECHAP_BUILD}"

run_cmake
make -j "${JOBS}"

echo "[install_spechap] SpecHap build finished."

# -----------------------------
# Build ExtractHAIRs
# -----------------------------
echo "[install_spechap] Building ExtractHAIRs..."

EXTRACT_SRC="${SPEC_HLA_ROOT}/bin/extractHairs"
EXTRACT_BUILD="${EXTRACT_SRC}/build"

if [[ ! -d "${EXTRACT_SRC}" ]]; then
    echo "[install_spechap] ERROR: ExtractHAIRs source directory not found at ${EXTRACT_SRC}"
    exit 1
fi

mkdir -p "${EXTRACT_BUILD}"
cd "${EXTRACT_BUILD}"

run_cmake
make -j "${JOBS}"

# Some versions may produce a binary with a slightly different name.
# Create a symlink named 'ExtractHAIRs' to make downstream scripts robust.
if [[ ! -x "${EXTRACT_BUILD}/ExtractHAIRs" ]]; then
    for alt in extractHairs extractHAIRs; do
        if [[ -x "${EXTRACT_BUILD}/${alt}" ]]; then
            ln -sf "${EXTRACT_BUILD}/${alt}" "${EXTRACT_BUILD}/ExtractHAIRs"
            echo "[install_spechap] Linked ${alt} -> ExtractHAIRs"
            break
        fi
    done
fi

echo "[install_spechap] ExtractHAIRs build finished."

# -----------------------------
# Ensure SpecHLA_RNAseq.sh is executable
# -----------------------------
RNA_SCRIPT="${SPEC_HLA_ROOT}/script/whole/SpecHLA_RNAseq.sh"
if [[ -f "${RNA_SCRIPT}" ]]; then
    chmod +x "${RNA_SCRIPT}"
    echo "[install_spechap] Marked ${RNA_SCRIPT} as executable."
else
    echo "[install_spechap] WARNING: SpecHLA_RNAseq.sh not found at ${RNA_SCRIPT}"
fi

echo "[install_spechap] Done."