#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Wrapper script for running SpecHLA RNAseq mode from iobrpy.

Features:
- Detect SpecHLA root automatically.
- Make sure SpecHap and ExtractHAIRs are built (run install_spechap.sh if needed).
- Build Bowtie2 index for the DRB reference only when *.bt2 files are missing.
- Automatically check and (when possible) install Python and external dependencies.
- Delegate the actual workflow to SpecHLA_RNAseq.sh.

This file is intended to be wired to the CLI entry point `iobrpy spechla`.
"""

import os
import sys
import argparse
import subprocess
from shutil import which
from iobrpy.utils.print_colorful_message import print_colorful_message

# ------------------------------------------------------------
# Generic helpers
# ------------------------------------------------------------
def run_cmd(cmd, cwd=None):
    """
    Run a shell command and raise an error if it fails.

    Parameters
    ----------
    cmd : list[str]
        Command with arguments, e.g. ["bash", "script.sh", "arg1"].
    cwd : str or None
        Working directory for the command.
    """
    print(f"[SpecHLA] Running command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[SpecHLA] ERROR: command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)


def prepend_to_path(path):
    """
    Prepend a directory to PATH if it is not already there.

    Parameters
    ----------
    path : str
        Directory to prepend.
    """
    if not path or not os.path.isdir(path):
        return
    current = os.environ.get("PATH", "")
    paths = current.split(os.pathsep) if current else []
    if path in paths:
        return
    os.environ["PATH"] = path + os.pathsep + current if current else path
    print(f"[SpecHLA] Prepending '{path}' to PATH")


# ------------------------------------------------------------
# Detect SpecHLA root
# ------------------------------------------------------------
def detect_spec_hla_root():
    """
    Detect the root directory of SpecHLA.

    Strategy
    --------
    - Assume this script lives under .../iobrpy/SpecHLA/ or similar.
    - Take the directory of this file as SPEC_HLA_ROOT.

    Returns
    -------
    str
        Absolute path to SpecHLA root directory.
    """
    this_file = os.path.abspath(__file__)
    spec_hla_root = os.path.dirname(this_file)
    print(f"[SpecHLA] Using SpecHLA root: {spec_hla_root}")
    return spec_hla_root


# ------------------------------------------------------------
# Python dependencies (pysam, biopython)
# ------------------------------------------------------------
def ensure_python_module(mod_name, pip_name=None, version=None):
    """
    Ensure a Python module can be imported, optionally enforcing a specific version.

    Parameters
    ----------
    mod_name : str
        Name used in `import mod_name`.
    pip_name : str or None
        Package name for pip install. If None, use mod_name.
    version : str or None
        Required module version (e.g. "0.19.0"). If provided, a mismatch will
        trigger installation of this exact version.
    """
    install_name = pip_name or mod_name
    needs_install = False

    try:
        m = __import__(mod_name)
        if version is not None:
            installed_version = getattr(m, "__version__", None)
            if installed_version != version:
                print(
                    f"[SpecHLA] Python module '{mod_name}' version mismatch: "
                    f"found {installed_version}, require {version}."
                )
                needs_install = True
            else:
                print(
                    f"[SpecHLA] Python module '{mod_name}' is available "
                    f"with required version {version}."
                )
                return
        else:
            print(f"[SpecHLA] Python module '{mod_name}' is available.")
            return
    except ImportError:
        print(f"[SpecHLA] Python module '{mod_name}' not found.")
        needs_install = True

    if not needs_install:
        return

    # Build install spec with pinned version if requested
    if version is not None and "==" not in install_name:
        install_spec = f"{install_name}=={version}"
    else:
        install_spec = install_name

    print(
        f"[SpecHLA] Trying to install '{install_spec}' via pip in current environment..."
    )
    cmd = [sys.executable, "-m", "pip", "install", install_spec]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(
            f"[SpecHLA] ERROR: failed to install '{install_spec}' (exit {e.returncode}).\n"
            f"          Please install it manually in the current conda environment, e.g.\n"
            f"          conda install -c bioconda {install_name}={version if version else ''}\n"
            f"          or: pip install {install_spec}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Re-try import after installation
    try:
        m = __import__(mod_name)
        if version is not None:
            installed_version = getattr(m, "__version__", None)
            if installed_version != version:
                print(
                    f"[SpecHLA] ERROR: module '{mod_name}' still not at required "
                    f"version {version} (found {installed_version}).",
                    file=sys.stderr,
                )
                sys.exit(1)
        print(f"[SpecHLA] Python module '{mod_name}' successfully installed and imported.")
    except ImportError:
        print(
            f"[SpecHLA] ERROR: module '{mod_name}' still cannot be imported after installation.\n"
            f"          Please check that you are using the intended python / conda environment.",
            file=sys.stderr,
        )
        sys.exit(1)


def ensure_python_deps():
    """
    Ensure core Python dependencies for SpecHLA scripts are available.

    - pysam  : used by assign_reads_to_genes.py / phase_variants.py
    - Bio    : from biopython, used by g_group_annotation.py
    """
    # Pin pysam to 0.19.0 to match SpecHLA's expected API / htslib behavior
    ensure_python_module("pysam", "pysam")
    # biopython does not require a specific version here
    ensure_python_module("Bio", "biopython")


# ------------------------------------------------------------
# Conda helpers for external tools and libraries
# ------------------------------------------------------------
def detect_conda_exe():
    """
    Try to detect the conda executable.

    Returns
    -------
    str or None
        Path to conda executable, or None if not found.
    """
    # Best source: CONDA_EXE when running inside a conda environment
    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe and os.path.isfile(conda_exe):
        return conda_exe

    # Fallback: search in PATH
    path = which("conda")
    if path:
        return path

    return None

def is_conda_package_installed(conda_exe, package_name):
    """
    Check whether a given package is already installed in the current
    conda environment.

    This is a best-effort heuristic based on `conda list package_name`.
    It returns True if the package appears in the output, otherwise False.
    If anything goes wrong, it falls back to False.
    """
    try:
        result = subprocess.run(
            [conda_exe, "list", package_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except Exception as e:
        print(
            f"[SpecHLA] WARNING: failed to query conda for package '{package_name}': {e}",
        )
        return False

    if result.returncode != 0:
        return False

    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        # Typical `conda list` line: name  version  build  channel
        cols = line.split()
        if cols and cols[0] == package_name:
            return True

    return False

def conda_install_packages(conda_exe, packages):
    """
    Install a list of packages into the current conda environment.

    Parameters
    ----------
    conda_exe : str
        Path to the conda executable.
    packages : list[str]
        Package names to install.
    """
    if not packages:
        return

    unique = sorted(set(packages))
    print(
        "[SpecHLA] Missing components detected; trying to install via conda:\n"
        f"          {' '.join(unique)}"
    )

    # Use bioconda + conda-forge to cover most bioinformatics tools.
    cmd = [conda_exe, "install", "-y", "-c", "bioconda", "-c", "conda-forge"] + unique
    run_cmd(cmd)

def get_conda_package_version(conda_exe, package_name):
    """
    Return installed version of a conda package in the current environment, or None.

    This uses `conda list <pkg>` and parses the plain-text table:
    name  version  build  channel
    """
    try:
        result = subprocess.run(
            [conda_exe, "list", package_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except Exception as e:
        print(
            f"[SpecHLA] WARNING: failed to query conda for package '{package_name}': {e}",
            file=sys.stderr,
        )
        return None

    if result.returncode != 0:
        return None

    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        # Typical: name  version  build  channel
        if cols and cols[0] == package_name and len(cols) >= 2:
            return cols[1]

    return None

def get_cmake_version(cmake_path):
    """
    Query the version of a cmake executable.

    Parameters
    ----------
    cmake_path : str
        Path to the cmake executable (e.g. from shutil.which("cmake")).

    Returns
    -------
    tuple[int, int, int] or None
        Version as (major, minor, patch) if it can be parsed, otherwise None.
    """
    try:
        result = subprocess.run(
            [cmake_path, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except Exception as e:
        print(
            f"[SpecHLA] WARNING: failed to query cmake version from '{cmake_path}': {e}",
            file=sys.stderr,
        )
        return None

    if result.returncode != 0 or not result.stdout:
        return None

    # First line usually looks like: "cmake version 3.29.2"
    first_line = result.stdout.splitlines()[0].strip()
    tokens = first_line.split()

    version_str = None
    for tok in tokens:
        if tok[0].isdigit():
            version_str = tok
            break

    if version_str is None:
        return None

    parts = []
    for p in version_str.split("."):
        # Strip non-digit suffixes, e.g. "3.29.2-rc1"
        digits = "".join(ch for ch in p if ch.isdigit())
        if digits:
            parts.append(int(digits))

    # Normalize to (major, minor, patch)
    while len(parts) < 3:
        parts.append(0)

    return tuple(parts[:3])

def ensure_vcflib_after_freebayes():
    """
    Ensure that vcflib has the required version *after* resolving freebayes.

    Logic
    -----
    - Use `conda list vcflib` to read the installed version.
    - Require vcflib==1.0.10.
    - If missing or mismatched, run:
        conda install -c bioconda -c conda-forge vcflib=1.0.10
    """
    required_version = "1.0.10"
    conda_exe = detect_conda_exe()

    if conda_exe is None:
        print(
            "[SpecHLA] ERROR: 'conda' executable not found; cannot check vcflib version via 'conda list'.\n"
            f"[SpecHLA]        Please make sure vcflib={required_version} is installed in the current environment.",
            file=sys.stderr,
        )
        sys.exit(1)

    installed_version = get_conda_package_version(conda_exe, "vcflib")
    if installed_version == required_version:
        print(f"[SpecHLA] vcflib is available with required version {required_version}.")
        return

    # Not installed or wrong version -> try to fix via conda
    spec = f"vcflib={required_version}"
    print(
        f"[SpecHLA] vcflib is missing or has incompatible version "
        f"(found {installed_version or 'none'}, require {required_version}).\n"
        f"[SpecHLA] Trying to install '{spec}' via conda...",
        file=sys.stderr,
    )
    cmd = [
        conda_exe,
        "install",
        "-y",
        "-c",
        "bioconda",
        "-c",
        "conda-forge",
        spec,
    ]
    run_cmd(cmd)

    # Re-check after installation
    installed_version = get_conda_package_version(conda_exe, "vcflib")
    if installed_version != required_version:
        print(
            f"[SpecHLA] ERROR: vcflib is still not at required version {required_version} "
            f"after conda installation (found {installed_version or 'none'}).",
            file=sys.stderr,
        )
        sys.exit(1)

# ------------------------------------------------------------
# External tools (bwa, samtools, freebayes, bgzip, tabix, bowtie2, bcftools, blastn)
# ------------------------------------------------------------
def ensure_external_tools(spec_hla_root):
    """
    Ensure external command-line tools are available.

    This function:
    - Prepends SpecHLA/bin and SpecHLA/bin/fermikit/fermi.kit to PATH,
      so bundled tools (bcftools, bwa, etc.) can be picked up.
    - Checks for required tools and tries to auto-install missing ones
      into the current conda environment using `conda install`.
    - For bcftools and blastn, always prefers the bundled binaries under
      SpecHLA/bin instead of installing them via conda.
    """

    # 1) Make sure SpecHLA/bin and fermikit kit are on PATH.
    bin_dir = os.path.join(spec_hla_root, "bin")
    fermikit_dir = os.path.join(bin_dir, "fermikit", "fermi.kit")

    prepend_to_path(bin_dir)
    prepend_to_path(fermikit_dir)

    # --- Make sure libwfa2.so.0 in SpecHLA/bin is visible to freebayes ---
    wfa_lib = os.path.join(bin_dir, "libwfa2.so.0")
    if os.path.exists(wfa_lib):
        ld_library_path = os.environ.get("LD_LIBRARY_PATH", "")
        ld_paths = ld_library_path.split(os.pathsep) if ld_library_path else []
        if bin_dir not in ld_paths:
            # Prepend SpecHLA/bin to LD_LIBRARY_PATH so the dynamic linker
            # can find libwfa2.so.0 when launching freebayes.
            new_ld = bin_dir if not ld_library_path else bin_dir + os.pathsep + ld_library_path
            os.environ["LD_LIBRARY_PATH"] = new_ld
            print(f"[SpecHLA] Prepending '{bin_dir}' to LD_LIBRARY_PATH for libwfa2.so.0")
    else:
        print(
            "[SpecHLA] WARNING: libwfa2.so.0 not found in SpecHLA/bin; "
            "freebayes may still fail to start.",
            file=sys.stderr,
        )

    # Required tools and the corresponding conda package providing them.
    # NOTE:
    # - blastn is provided as a bundled binary in SpecHLA/bin and is NOT
    #   auto-installed via conda here.
    tool_to_conda_pkg = {
        "samtools": "samtools=1.21",
        "bwa": "bwa",
        "bowtie2": "bowtie2",
        "freebayes": "freebayes=1.3.8",
        # htslib provides bgzip/tabix with the versions SpecHLA expects
        "bgzip": "htslib=1.21",
        "tabix": "htslib=1.21",
    }

    missing_tools = []

    for tool, pkg in tool_to_conda_pkg.items():
        path = which(tool)
        if path is None:
            print(f"[SpecHLA] External program '{tool}' not found in PATH.")
            missing_tools.append(tool)
        else:
            print(f"[SpecHLA] Found {tool}: {path}")

    # ---------------- bundled bcftools ----------------
    bcftools_bundled = os.path.join(bin_dir, "bcftools")

    if os.path.exists(bcftools_bundled):
        if os.access(bcftools_bundled, os.X_OK):
            print(f"[SpecHLA] Using bundled bcftools at {bcftools_bundled}.")
            os.environ["BCFTOOLS"] = bcftools_bundled
        else:
            print(
                f"[SpecHLA] Found bundled bcftools at {bcftools_bundled} but it is not executable.\n"
                f"[SpecHLA] Trying to add execute permission (chmod +x)..."
            )
            try:
                st = os.stat(bcftools_bundled)
                os.chmod(bcftools_bundled, st.st_mode | 0o111)
            except Exception as e:
                print(
                    f"[SpecHLA] ERROR: failed to add execute permission to {bcftools_bundled}: {e}",
                    file=sys.stderr,
                )
                print(
                    "[SpecHLA]        Please run 'chmod +x' on this file manually and try again.",
                    file=sys.stderr,
                )
                sys.exit(1)

            if os.access(bcftools_bundled, os.X_OK):
                print(f"[SpecHLA] Successfully made bcftools executable at {bcftools_bundled}.")
                os.environ["BCFTOOLS"] = bcftools_bundled
            else:
                print(
                    f"[SpecHLA] ERROR: bcftools at {bcftools_bundled} is still not executable after chmod.",
                    file=sys.stderr,
                )
                print(
                    "[SpecHLA]        Please check file permissions and filesystem mount options.",
                    file=sys.stderr,
                )
                sys.exit(1)
    else:
        print(
            f"[SpecHLA] ERROR: bundled bcftools not found at {bcftools_bundled}.",
            file=sys.stderr,
        )
        print(
            "[SpecHLA]        Please make sure SpecHLA/bin/bcftools exists, "
            "for example by copying it from the original SpecHLA package.",
            file=sys.stderr,
        )
        sys.exit(1)

    # ---------------- bundled blastn ----------------
    # Always prefer the blastn binary shipped with SpecHLA under SpecHLA/bin.
    blastn_bundled = os.path.join(bin_dir, "blastn")

    if os.path.exists(blastn_bundled):
        if os.access(blastn_bundled, os.X_OK):
            print(f"[SpecHLA] Using bundled blastn at {blastn_bundled}.")
        else:
            print(
                f"[SpecHLA] Found bundled blastn at {blastn_bundled} but it is not executable.\n"
                f"[SpecHLA] Trying to add execute permission (chmod +x)..."
            )
            try:
                st = os.stat(blastn_bundled)
                os.chmod(blastn_bundled, st.st_mode | 0o111)
            except Exception as e:
                print(
                    f"[SpecHLA] ERROR: failed to add execute permission to {blastn_bundled}: {e}",
                    file=sys.stderr,
                )
                print(
                    "[SpecHLA]        Please run 'chmod +x' on this file manually and try again.",
                    file=sys.stderr,
                )
                sys.exit(1)

            if os.access(blastn_bundled, os.X_OK):
                print(f"[SpecHLA] Successfully made blastn executable at {blastn_bundled}.")
            else:
                print(
                    f"[SpecHLA] ERROR: blastn at {blastn_bundled} is still not executable after chmod.",
                    file=sys.stderr,
                )
                print(
                    "[SpecHLA]        Please check file permissions and filesystem mount options.",
                    file=sys.stderr,
                )
                sys.exit(1)
    else:
        print(
            f"[SpecHLA] ERROR: bundled blastn not found at {blastn_bundled}.",
            file=sys.stderr,
        )
        print(
            "[SpecHLA]        Please make sure SpecHLA/bin/blastn exists, "
            "for example by copying it from the original SpecHLA package.",
            file=sys.stderr,
        )
        sys.exit(1)

    # If all conda-managed tools are available, we are done here.
    if not missing_tools:
        ensure_vcflib_after_freebayes()
        return

    # At least one conda-managed tool is missing: try to install them via conda.
    conda_exe = detect_conda_exe()
    if conda_exe is None:
        print(
            "[SpecHLA] ERROR: Some external tools are missing and 'conda' could not be found.\n"
            "          Please install the following tools manually in your environment:\n"
            f"          {', '.join(sorted(missing_tools))}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Map tools -> packages and install the unique set.
    pkgs = [tool_to_conda_pkg[t] for t in missing_tools]
    conda_install_packages(conda_exe, pkgs)

    # Re-check after installation.
    still_missing = []
    for tool in missing_tools:
        path = which(tool)
        if path is None:
            still_missing.append(tool)
        else:
            print(f"[SpecHLA] After conda install, found {tool}: {path}")

    if still_missing:
        print(
            "[SpecHLA] ERROR: The following tools are still missing even after conda installation:\n"
            f"          {', '.join(sorted(still_missing))}\n"
            "          Please check your conda environment and install them manually.",
            file=sys.stderr,
        )
        sys.exit(1)

    ensure_vcflib_after_freebayes()

# ------------------------------------------------------------
# SpecHap / ExtractHAIRs & Bowtie2 index
# ------------------------------------------------------------
def ensure_spechap_built(spec_hla_root, threads):
    """
    Make sure SpecHap and ExtractHAIRs are built under SpecHLA/bin.

    Logic
    -----
    - If both build directories exist:
        <SPEC_HLA_ROOT>/bin/SpecHap/build
        <SPEC_HLA_ROOT>/bin/extractHairs/build
      then assume they are already built and do nothing.
    - Otherwise:
      - Only if ARPACK is not detected in the current conda environment,
        try to install ARPACK via conda.
      - Ensure CMake is available and has version <3.30
        (install cmake<3.30 via conda if needed).
      - Then run install_spechap.sh to build SpecHap and ExtractHAIRs.

    Parameters
    ----------
    spec_hla_root : str
        Path to SpecHLA root.
    threads : int
        Number of threads to pass to install_spechap.sh (if it uses it).
    """
    spec_hap_build_dir = os.path.join(spec_hla_root, "bin", "SpecHap", "build")
    extract_hairs_build_dir = os.path.join(spec_hla_root, "bin", "extractHairs", "build")

    spec_hap_ok = os.path.isdir(spec_hap_build_dir)
    extract_ok = os.path.isdir(extract_hairs_build_dir)

    # If both build directories exist, we assume SpecHap / ExtractHAIRs are
    # already compiled and ready to use -> no need to install ARPACK or CMake.
    if spec_hap_ok and extract_ok:
        print("[SpecHLA] Found SpecHap and ExtractHAIRs build directories:")
        print(f"          {spec_hap_build_dir}")
        print(f"          {extract_hairs_build_dir}")
        return

    print("[SpecHLA] SpecHap or ExtractHAIRs build directory not found;")
    print("[SpecHLA] will ensure dependencies (ARPACK / CMake<3.30) and then build them with install_spechap.sh...")

    # Step 1: detect conda and prepare a list of packages to install
    conda_exe = detect_conda_exe()
    pkgs_to_install = []

    # --- ARPACK: only auto-install when build dirs are missing AND ARPACK
    #             is not already present in the current conda environment. ---
    if conda_exe is not None:
        if is_conda_package_installed(conda_exe, "arpack"):
            print("[SpecHLA] Detected existing ARPACK installation; skip installing it.")
        else:
            print("[SpecHLA] ARPACK not detected; will install 'arpack' via conda.")
            pkgs_to_install.append("arpack")
    else:
        print(
            "[SpecHLA] WARNING: Could not find 'conda'; skipping automatic ARPACK installation.\n"
            "          If the build fails due to missing ARPACK, please install it manually "
            "in the current environment (e.g. via conda-forge).",
            file=sys.stderr,
        )

    # --- CMake: require a version strictly less than 3.30 ---
    required_cmake_max = (3, 30, 0)
    cmake_path = which("cmake")

    def _require_cmake_via_conda():
        if conda_exe is not None:
            pkgs_to_install.append("cmake<3.30")
        else:
            print(
                "[SpecHLA] ERROR: 'cmake' (version <3.30) is required to build SpecHap/ExtractHAIRs,\n"
                "          but neither a suitable cmake nor conda could be found.\n"
                "          Please install cmake<3.30 manually in your environment (e.g. via conda-forge)\n"
                "          and re-run this command.",
                file=sys.stderr,
            )
            sys.exit(1)

    if cmake_path is None:
        print("[SpecHLA] 'cmake' not found in PATH.")
        _require_cmake_via_conda()
    else:
        cmake_ver = get_cmake_version(cmake_path)
        if cmake_ver is None:
            print(
                f"[SpecHLA] WARNING: could not determine cmake version from '{cmake_path}'.\n"
                "[SpecHLA]          Will try to install cmake<3.30 via conda to satisfy SpecHLA's requirement.",
                file=sys.stderr,
            )
            _require_cmake_via_conda()
        else:
            if cmake_ver < required_cmake_max:
                print(
                    f"[SpecHLA] Found cmake {cmake_ver[0]}.{cmake_ver[1]}.{cmake_ver[2]} "
                    f"(OK: version < 3.30)."
                )
            else:
                print(
                    f"[SpecHLA] Detected cmake {cmake_ver[0]}.{cmake_ver[1]}.{cmake_ver[2]}, "
                    f"but SpecHLA requires cmake<3.30."
                )
                _require_cmake_via_conda()

    # Step 2: install missing dependencies via conda (if available)
    if conda_exe is not None and pkgs_to_install:
        # conda will skip packages that are already installed with a compatible version
        conda_install_packages(conda_exe, pkgs_to_install)

        # Re-check cmake after installation to ensure it is now available
        cmake_path = which("cmake")
        if cmake_path is None:
            print(
                "[SpecHLA] ERROR: 'cmake' is still not available after conda installation.\n"
                "          Please check your conda environment and install cmake<3.30 manually.",
                file=sys.stderr,
            )
            sys.exit(1)

        cmake_ver = get_cmake_version(cmake_path)
        if cmake_ver is None or not (cmake_ver < required_cmake_max):
            print(
                "[SpecHLA] ERROR: After conda installation, cmake does not satisfy the version requirement (<3.30).\n"
                f"          Detected version: {cmake_ver if cmake_ver is not None else 'unknown'}.\n"
                "          Please ensure that a cmake version <3.30 is installed and visible in PATH.",
                file=sys.stderr,
            )
            sys.exit(1)
        else:
            print(
                f"[SpecHLA] Using cmake {cmake_ver[0]}.{cmake_ver[1]}.{cmake_ver[2]} "
                f"after conda installation (OK: <3.30)."
            )

    # Step 3: run install_spechap.sh to build SpecHap / ExtractHAIRs
    install_script = os.path.join(spec_hla_root, "install_spechap.sh")
    if not os.path.exists(install_script):
        print(f"[SpecHLA] ERROR: install_spechap.sh not found at {install_script}", file=sys.stderr)
        sys.exit(1)

    run_cmd(["bash", install_script, str(threads)])
    print("[SpecHLA] Finished running install_spechap.sh.")

def ensure_bowtie2_index(spec_hla_root, bowtie2_build_path, drb_ref_relpath):
    """
    Ensure Bowtie2 index for DRB reference exists; build it only if missing.

    Parameters
    ----------
    spec_hla_root : str
        Path to SpecHLA root.
    bowtie2_build_path : str
        Path to the bowtie2-build executable.
    drb_ref_relpath : str
        Relative path from SpecHLA root to DRB reference FASTA
        (e.g. 'db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta').
    """
    # Absolute path to the DRB reference FASTA
    ref_fasta = os.path.join(spec_hla_root, drb_ref_relpath)
    if not os.path.exists(ref_fasta):
        print(f"[SpecHLA] ERROR: DRB reference fasta not found at {ref_fasta}", file=sys.stderr)
        sys.exit(1)

    # Bowtie2 index prefix is usually the same as the FASTA path
    index_prefix = ref_fasta

    # Expected Bowtie2 index files
    bt2_files = [
        f"{index_prefix}.1.bt2",
        f"{index_prefix}.2.bt2",
        f"{index_prefix}.3.bt2",
        f"{index_prefix}.4.bt2",
        f"{index_prefix}.rev.1.bt2",
        f"{index_prefix}.rev.2.bt2",
    ]

    if all(os.path.exists(f) for f in bt2_files):
        # Index already exists -> reuse it
        print("[SpecHLA] Detected existing Bowtie2 index for DRB reference, skip building.")
        return

    # Index missing -> build it
    print(f"[SpecHLA] Building Bowtie2 index for {os.path.basename(ref_fasta)}...")
    if bowtie2_build_path is None:
        print("[SpecHLA] ERROR: bowtie2-build not found in PATH or environment.", file=sys.stderr)
        sys.exit(1)

    run_cmd([bowtie2_build_path, ref_fasta, index_prefix])
    print("[SpecHLA] Successfully created Bowtie2 index for DRB reference.")


def detect_bowtie2_build():
    """
    Try to detect bowtie2-build executable.

    Returns
    -------
    str or None
        Path to bowtie2-build, or None if not found.
    """
    # 1) Environment variable (allows explicit override)
    env_path = os.environ.get("BOWTIE2_BUILD")
    if env_path and os.path.exists(env_path):
        return env_path

    # 2) Search in PATH
    return which("bowtie2-build")


# ------------------------------------------------------------
# Delegate to SpecHLA_RNAseq.sh
# ------------------------------------------------------------
def run_spechla_rnaseq(spec_hla_root, sample_name, read1, read2, outdir, threads):
    """
    Delegate to SpecHLA_RNAseq.sh with the given arguments.

    Parameters
    ----------
    spec_hla_root : str
        Path to SpecHLA root.
    sample_name : str
        Sample name (-n).
    read1 : str
        Path to R1 FASTQ.gz (-1).
    read2 : str
        Path to R2 FASTQ.gz (-2).
    outdir : str
        Output directory (-o).
    threads : int
        Number of threads (-j).
    """
    rnaseq_script = os.path.join(spec_hla_root, "script", "whole", "SpecHLA_RNAseq.sh")
    if not os.path.exists(rnaseq_script):
        print(f"[SpecHLA] ERROR: SpecHLA_RNAseq.sh not found at {rnaseq_script}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "bash",
        rnaseq_script,
        "-n",
        sample_name,
        "-1",
        read1,
        "-2",
        read2,
        "-o",
        outdir,
        "-j",
        str(threads),
    ]
    run_cmd(cmd)

def merge_hla_results(outdir, merged_filename="hla_result_merged.txt"):
    """
    Merge per-sample SpecHLA `hla.result.txt` files into a single table.

    This function scans all immediate subdirectories of `outdir`. For each
    subdirectory, it looks for a file named `hla.result.txt`. The header
    line is taken from the first file and reused for all samples. Data
    lines from all files are concatenated and written to
    `hla_result_merged.txt` in `outdir`.

    Parameters
    ----------
    outdir : str
        Root output directory used for SpecHLA (the same as `-o`).
    merged_filename : str, optional
        Name of the merged output file. Default is `hla_result_merged.txt`.
    """
    # We assume `os` and `sys` are already imported at the top of this file.
    if not os.path.isdir(outdir):
        # If the output directory does not exist, there is nothing to merge.
        print(f"[SpecHLA] WARNING: Output directory does not exist: {outdir}", file=sys.stderr)
        return

    # Collect all immediate subdirectories under the root output directory.
    # In the standard SpecHLA layout, each folder corresponds to a sample ID.
    sample_dirs = [
        d for d in sorted(os.listdir(outdir))
        if os.path.isdir(os.path.join(outdir, d))
    ]

    header = None
    data_lines = []

    for sample_dir in sample_dirs:
        # Each sample folder is expected to contain `hla.result.txt`.
        result_path = os.path.join(outdir, sample_dir, "hla.result.txt")
        if not os.path.isfile(result_path):
            # Skip this sample if no result file is found.
            continue

        with open(result_path, "r") as f:
            for i, line in enumerate(f):
                # Remove the trailing newline for consistent handling.
                line = line.rstrip("\n")

                # Skip completely empty lines (if any).
                if not line:
                    continue

                if i == 0:
                    # The first line is the header. We keep it only once.
                    if header is None:
                        header = line
                    elif header != line:
                        # If headers are different between files, print a warning
                        # but still continue merging the data rows.
                        print(
                            f"[SpecHLA] WARNING: Header mismatch in {result_path}",
                            file=sys.stderr,
                        )
                    # Do not append header lines to `data_lines`.
                    continue

                # All non-header, non-empty lines are treated as data rows.
                data_lines.append(line)

    if header is None:
        # If we never saw any header, it means no `hla.result.txt` files were found.
        print(
            f"[SpecHLA] WARNING: No 'hla.result.txt' files found under {outdir}; "
            f"skip generating merged result.",
            file=sys.stderr,
        )
        return

    merged_path = os.path.join(outdir, merged_filename)
    with open(merged_path, "w") as out_f:
        # Write the header once
        out_f.write(header + "\n")
        # Then write all collected data rows
        for line in data_lines:
            out_f.write(line + "\n")

    print(
        f"[SpecHLA] Wrote merged HLA typing table with {len(data_lines)} row(s) to: {merged_path}"
    )

# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
def parse_args(argv=None):
    """
    Parse command-line arguments for the spechla wrapper.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Run SpecHLA (RNAseq mode) through iobrpy wrapper."
    )
    parser.add_argument(
        "-n",
        "--name",
        required=True,
        help="Sample name.",
    )
    parser.add_argument(
        "-1",
        "--read1",
        required=True,
        help="Path to read1 FASTQ.gz.",
    )
    parser.add_argument(
        "-2",
        "--read2",
        required=True,
        help="Path to read2 FASTQ.gz.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use (default: 8).",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """
    Main entry point for the SpecHLA wrapper.
    """
    args = parse_args(argv)

    # 1. Detect SpecHLA root
    spec_hla_root = detect_spec_hla_root()

    # 2. Ensure Python dependencies (pysam / biopython)
    ensure_python_deps()

    # 3. Ensure external tools (bwa / samtools / freebayes / bgzip / tabix / bowtie2 / bcftools / blastn)
    ensure_external_tools(spec_hla_root)

    # 4. Make sure SpecHap / ExtractHAIRs are built
    ensure_spechap_built(spec_hla_root, args.threads)

    # 5. Make sure Bowtie2 index for DRB reference exists
    bowtie2_build_path = detect_bowtie2_build()
    drb_ref_relpath = os.path.join("db", "ref", "hla_gen.format.filter.extend.DRB.no26789.fasta")
    ensure_bowtie2_index(spec_hla_root, bowtie2_build_path, drb_ref_relpath)

    # 6. Run SpecHLA_RNAseq.sh
    run_spechla_rnaseq(
        spec_hla_root=spec_hla_root,
        sample_name=args.name,
        read1=args.read1,
        read2=args.read2,
        outdir=args.outdir,
        threads=args.threads,
    )

    try:
        merge_hla_results(args.outdir)
    except Exception as e:
        # If merging fails for any reason, we just print a warning
        # and do not stop the whole pipeline.
        print(
            f"[SpecHLA] WARNING: Failed to merge HLA result tables: {e}",
            file=sys.stderr,
        )

    # ---------------- IOBRpy banner ----------------
    print("   ")
    print_colorful_message("#########################################################", "blue")
    print_colorful_message(" IOBRpy: Immuno-Oncology Biological Research using Python ", "cyan")
    print_colorful_message(" If you encounter any issues, please report them at ", "cyan")
    print_colorful_message(" https://github.com/IOBR/IOBRpy/issues ", "cyan")
    print_colorful_message("#########################################################", "blue")
    print(" Author: Haonan Huang, Dongqiang Zeng")
    print(" Email: interlaken@smu.edu.cn ")
    print_colorful_message("#########################################################", "blue")
    print("   ")


if __name__ == "__main__":
    main()