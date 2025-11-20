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
def ensure_python_module(mod_name, pip_name=None):
    """
    Ensure a Python module can be imported. If not, try to install it via pip.

    Parameters
    ----------
    mod_name : str
        Name used in `import mod_name`.
    pip_name : str or None
        Package name for pip install. If None, use mod_name.
    """
    try:
        __import__(mod_name)
        print(f"[SpecHLA] Python module '{mod_name}' is available.")
        return
    except ImportError:
        install_name = pip_name or mod_name
        print(f"[SpecHLA] Python module '{mod_name}' not found.")
        print(f"[SpecHLA] Trying to install '{install_name}' via pip in current environment...")

        cmd = [sys.executable, "-m", "pip", "install", install_name]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"[SpecHLA] ERROR: failed to install '{install_name}' (exit {e.returncode}).\n"
                f"          Please install it manually in the current conda environment, e.g.\n"
                f"          conda install -c bioconda {install_name}\n"
                f"          or: pip install {install_name}",
                file=sys.stderr,
            )
            sys.exit(1)

    # Re-try import after installation
    try:
        __import__(mod_name)
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
    ensure_python_module("pysam", "pysam")
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


# ------------------------------------------------------------
# External tools (bwa, samtools, freebayes, bgzip, tabix, bowtie2, bcftools)
# ------------------------------------------------------------
def ensure_external_tools(spec_hla_root):
    """
    Ensure external command-line tools are available.

    This function:
    - Prepends SpecHLA/bin and SpecHLA/bin/fermikit/fermi.kit to PATH,
      so bundled tools (bcftools, bwa, etc.) can be picked up.
    - Checks for required tools and tries to auto-install missing ones
      into the current conda environment using `conda install`.
    """

    # 1) Make sure SpecHLA/bin and fermikit kit are on PATH.
    bin_dir = os.path.join(spec_hla_root, "bin")
    fermikit_dir = os.path.join(bin_dir, "fermikit", "fermi.kit")

    prepend_to_path(bin_dir)
    prepend_to_path(fermikit_dir)

    # Required tools and the corresponding conda package providing them.
    tool_to_conda_pkg = {
        "samtools": "samtools",
        "bwa": "bwa",
        "bowtie2": "bowtie2",
        "freebayes": "freebayes",
        "bgzip": "htslib",   # bgzip is provided by htslib
        "tabix": "htslib",   # tabix is also from htslib
    }

    missing_tools = []

    for tool, pkg in tool_to_conda_pkg.items():
        path = which(tool)
        if path is None:
            print(f"[SpecHLA] External program '{tool}' not found in PATH.")
            missing_tools.append(tool)
        else:
            print(f"[SpecHLA] Found {tool}: {path}")

    # Check bcftools separately: SpecHLA prefers the bundled one under bin/.
    bcftools_bundled = os.path.join(bin_dir, "bcftools")
    if os.path.exists(bcftools_bundled) and os.access(bcftools_bundled, os.X_OK):
        print(f"[SpecHLA] Found bundled bcftools at {bcftools_bundled} (SpecHLA/bin/bcftools).")
    else:
        bcftools_path = which("bcftools")
        if bcftools_path:
            print(f"[SpecHLA] Found bcftools in PATH: {bcftools_path}")
        else:
            print(
                "[SpecHLA] WARNING: bcftools was not detected.\n"
                "          SpecHLA will use SpecHLA/bin/bcftools if present.\n"
                "          If you prefer a system-wide bcftools, install it via:\n"
                "          conda install -c bioconda bcftools",
                file=sys.stderr,
            )

    if not missing_tools:
        # Everything is already available.
        return

    # At least one tool is missing: try to install them via conda.
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


# ------------------------------------------------------------
# SpecHap / ExtractHAIRs & Bowtie2 index
# ------------------------------------------------------------
def ensure_spechap_built(spec_hla_root, threads):
    """
    Make sure SpecHap and ExtractHAIRs are built under SpecHLA/bin.

    If either directory is missing, run install_spechap.sh to build them.

    Parameters
    ----------
    spec_hla_root : str
        Path to SpecHLA root.
    threads : int
        Number of threads to pass to install_spechap.sh (if it uses it).
    """
    # Expected layout:
    #   <SPEC_HLA_ROOT>/bin/SpecHap
    #   <SPEC_HLA_ROOT>/bin/extractHairs
    spec_hap_dir = os.path.join(spec_hla_root, "bin", "SpecHap")
    extract_hairs_dir = os.path.join(spec_hla_root, "bin", "extractHairs")

    spec_hap_ok = os.path.isdir(spec_hap_dir)
    extract_ok = os.path.isdir(extract_hairs_dir)

    if spec_hap_ok and extract_ok:
        print(f"[SpecHLA] Found SpecHap and ExtractHAIRs in {os.path.join(spec_hla_root, 'bin')}")
        return

    print("[SpecHLA] SpecHap or ExtractHAIRs not found under bin/, trying to build them...")
    install_script = os.path.join(spec_hla_root, "install_spechap.sh")
    if not os.path.exists(install_script):
        print(f"[SpecHLA] ERROR: install_spechap.sh not found at {install_script}", file=sys.stderr)
        sys.exit(1)

    # Note: ARPACK auto-install is now handled inside install_spechap.sh.
    # Note: even if 'threads' is ignored by the script, passing it is harmless.
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

    # 3. Ensure external tools (bwa / samtools / freebayes / bgzip / tabix / bowtie2 / bcftools)
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