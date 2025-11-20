#!/usr/bin/env python
"""
Python wrapper for SpecHLA's ExtractHLAread.sh.

This script:
  1. Parses the same CLI arguments as ExtractHLAread.sh (-s, -b, -r, -o).
  2. Checks whether required system tools (samtools, bam from bamUtil) are
     available in the current conda environment.
  3. Tries to auto-install missing tools via conda/mamba (bioconda + conda-forge).
  4. Calls ExtractHLAread.sh via bash.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List
from iobrpy.utils.print_colorful_message import print_colorful_message

# Mapping from required binary name -> conda package name.
# ExtractHLAread.sh uses "samtools" and "bam" (from bamUtil). :contentReference[oaicite:0]{index=0}
REQUIRED_TOOLS: Dict[str, str] = {
    "samtools": "samtools",
    "bam": "bamutil",
}


def _running_in_conda() -> bool:
    """
    Return True if the current Python process looks like it is running
    inside a conda environment.
    """
    return bool(os.environ.get("CONDA_PREFIX") or os.environ.get("CONDA_DEFAULT_ENV"))


def _detect_conda_executable() -> str | None:
    """
    Try to find a conda-like executable in PATH (mamba, conda, micromamba).
    Returns the path to the first one that is found, or None if none exist.
    """
    for exe in ("mamba", "conda", "micromamba"):
        path = shutil.which(exe)
        if path is not None:
            return path
    return None


def ensure_dependencies(auto_install: bool = True) -> None:
    """
    Check whether required command-line tools are available.

    If some are missing and auto_install is True, try to install the matching
    conda packages into the current environment using bioconda + conda-forge.

    Raises RuntimeError if tools cannot be installed or are still missing.
    """
    missing_binaries: List[str] = [
        bin_name for bin_name in REQUIRED_TOOLS if shutil.which(bin_name) is None
    ]
    if not missing_binaries:
        return

    if not auto_install:
        missing_str = ", ".join(missing_binaries)
        raise RuntimeError(
            f"Missing required tools: {missing_str}. "
            "Please install them (for example: "
            "'conda install -c bioconda -c conda-forge samtools bamutil')."
        )

    if not _running_in_conda():
        missing_str = ", ".join(missing_binaries)
        raise RuntimeError(
            "Missing required tools and current process does not appear to run "
            "inside a conda environment. Auto-installation is disabled.\n"
            f"Please install these tools manually: {missing_str}.\n"
            "Example:\n"
            "  conda install -c bioconda -c conda-forge samtools bamutil"
        )

    conda_exe = _detect_conda_executable()
    if conda_exe is None:
        missing_str = ", ".join(missing_binaries)
        raise RuntimeError(
            "Missing required tools but no conda/mamba executable was found in PATH.\n"
            f"Please install these tools manually: {missing_str}."
        )

    # Map binaries to conda packages
    packages = sorted({REQUIRED_TOOLS[bin_name] for bin_name in missing_binaries})

    cmd = [
        conda_exe,
        "install",
        "-y",
        "-c",
        "bioconda",
        "-c",
        "conda-forge",
    ] + packages

    print(
        f"[extract_hla_read] Installing missing tools via: {' '.join(cmd)}",
        file=sys.stderr,
    )
    subprocess.run(cmd, check=True)

    # Re-check after installation
    still_missing = [bin_name for bin_name in missing_binaries if shutil.which(bin_name) is None]
    if still_missing:
        raise RuntimeError(
            "Failed to install some required tools automatically. "
            f"Still missing: {', '.join(still_missing)}"
        )


def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the local argparse parser that mirrors ExtractHLAread.sh.

    Original shell usage: :contentReference[oaicite:1]{index=1}
      <PATH-TO>/ExtractHLAread.sh -s <sample_id> -b <bamfile> -r <refGenome> -o <outdir>
    """
    parser = argparse.ArgumentParser(
        prog="extract_hla_read",
        description=(
            "Python wrapper for SpecHLA's ExtractHLAread.sh to extract HLA-related reads "
            "from a BAM/CRAM file and convert them to FASTQ."
        ),
    )
    parser.add_argument(
        "-s",
        "--sample",
        dest="sample_id",
        required=True,
        help="Desired sample name (for example: NA12878).",
    )
    parser.add_argument(
        "-b",
        "--bam",
        dest="bam_path",
        required=True,
        help="Sorted and indexed BAM or CRAM file.",
    )
    parser.add_argument(
        "-r",
        "--ref",
        dest="ref",
        required=True,
        choices=["hg38", "hg19"],
        help="Reference genome: hg38 or hg19.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help="Output folder where extracted reads will be written.",
    )
    parser.add_argument(
        "--no-auto-install",
        dest="no_auto_install",
        action="store_true",
        help=(
            "Disable automatic installation of missing tools (samtools, bamutil) "
            "even if they are not found in the current conda environment."
        ),
    )
    return parser


def run_extraction(sample_id: str, bam_path: Path, ref: str, outdir: Path) -> None:
    """
    Locate ExtractHLAread.sh under iobrpy/SpecHLA/script and run it via bash
    with the provided arguments.
    """
    # This file lives in iobrpy/SpecHLA; ExtractHLAread.sh is in SpecHLA/script. :contentReference[oaicite:2]{index=2}
    script_dir = Path(__file__).resolve().parent
    sh_path = script_dir / "script" / "ExtractHLAread.sh"

    if not sh_path.is_file():
        raise FileNotFoundError(
            f"Cannot find ExtractHLAread.sh at {sh_path}. "
            "Please make sure it is included in the installed iobrpy package."
        )

    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bash",
        str(sh_path),
        "-s",
        sample_id,
        "-b",
        str(bam_path),
        "-r",
        ref,
        "-o",
        str(outdir),
    ]

    # Print the command for debugging
    print("[extract_hla_read] Running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> None:
    """
    Entry point for the wrapper.

    When called from iobrpy's main CLI, `argv` will be the list of arguments
    after the subcommand name (e.g. ['-s', 'SAMPLE', '-b', 'sample.bam', ...]).
    When executed as a standalone script, it uses sys.argv[1:].
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        ensure_dependencies(auto_install=not args.no_auto_install)
    except RuntimeError as exc:
        # Use parser.error so we get a clean error message and exit code 2
        parser.error(str(exc))

    bam_path = Path(os.path.expanduser(args.bam_path)).resolve()
    outdir = Path(os.path.expanduser(args.outdir)).resolve()

    run_extraction(
        sample_id=args.sample_id,
        bam_path=bam_path,
        ref=args.ref,
        outdir=outdir,
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