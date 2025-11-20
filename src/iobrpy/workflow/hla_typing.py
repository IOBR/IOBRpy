#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch HLA typing workflow.

This script wires the two SpecHLA helpers shipped with iobrpy:

1) Run ExtractHLAread.sh on every BAM file in a directory.
2) Run SpecHLA_RNAseq.sh on every resulting FASTQ pair.
3) Merge per-sample HLA result tables into one combined file.

Typical usage (through the main iobrpy CLI):

    iobrpy hla_typing -b /path/to/bam_dir -r hg38 -o /path/to/output -j 8

Arguments
---------
-b / --bam-dir
    Directory that contains one or more *.bam files.
    * If a filename ends with "_Aligned.sortedByCoord.out.bam",
      the sample id is the prefix before that suffix.
    * Otherwise the sample id is the prefix before ".bam".
    The inferred sample id will be used for both -s in ExtractHLAread
    and -n in SpecHLA.
-r / --ref
    Reference genome for ExtractHLAread, "hg19" or "hg38".
-o / --outdir
    Root working directory. Two subdirectories will be created:
        <outdir>/ExtractHLAread/<sample_id>  (ExtractHLAread output)
        <outdir>/SpecHLA/<sample_id>        (SpecHLA output, assumed layout)
-j / --threads
    Number of threads for SpecHLA (-j). Default: 8.

Layout for SpecHLA results
--------------------------
For each sample <id>, SpecHLA is assumed to write an HLA result table to

    <outdir>/SpecHLA/<id>/hla.result.txt   (preferred)
or  <outdir>/SpecHLA/<id>/hla.results.txt  (fallback)

The format follows the example provided by the user and starts with
a version comment line, then a header line, then one or more data lines.

This script will append all per-sample result tables into a single
merged file:

    <outdir>/hla_result_merged.txt
"""

from __future__ import annotations

import argparse
import os
import sys
from tqdm import tqdm
from pathlib import Path
from typing import List, Tuple
from iobrpy.utils.print_colorful_message import print_colorful_message

from iobrpy.SpecHLA.extract_hla_read import (
    ensure_dependencies as ensure_extract_deps,
    run_extraction,
)
from iobrpy.SpecHLA.SpecHLA import (
    detect_spec_hla_root,
    ensure_python_deps,
    ensure_external_tools,
    ensure_spechap_built,
    ensure_bowtie2_index,
    detect_bowtie2_build,
    run_spechla_rnaseq,
)

# ---------------------------------------------------------------------------
# Sample collection & name inference
# ---------------------------------------------------------------------------

STAR_SUFFIX = "_Aligned.sortedByCoord.out.bam"


def infer_sample_id(bam_path: Path) -> str:
    """
    Infer sample id from a BAM file name.

    Rules
    -----
    1) If the name ends with "_Aligned.sortedByCoord.out.bam",
       strip that suffix.
    2) Otherwise, strip the ".bam" suffix.

    The resulting id will be used as:
        -s for ExtractHLAread.sh
        -n for SpecHLA_RNAseq.sh
    """
    name = bam_path.name
    if name.endswith(STAR_SUFFIX):
        return name[: -len(STAR_SUFFIX)]
    if name.endswith(".bam"):
        return name[: -len(".bam")]

    # In practice we only look at *.bam, so this should not happen.
    raise ValueError(f"Invalid BAM file name (not ending with .bam): {bam_path}")


def collect_samples(bam_dir: Path) -> List[Tuple[str, Path]]:
    """
    Scan the BAM directory and return a list of (sample_id, bam_path).
    """
    bam_files = sorted(bam_dir.glob("*.bam"))
    samples: List[Tuple[str, Path]] = []
    for bam in bam_files:
        sample_id = infer_sample_id(bam)
        samples.append((sample_id, bam.resolve()))
    return samples

# ---------------------------------------------------------------------------
# ExtractHLAread phase
# ---------------------------------------------------------------------------

def run_extract_phase(
    samples: List[Tuple[str, Path]],
    ref: str,
    extract_root: Path,
) -> None:
    """
    Run ExtractHLAread.sh for all samples.

    Parameters
    ----------
    samples :
        List of (sample_id, bam_path).
    ref :
        Reference genome for ExtractHLAread (hg19 or hg38).
    extract_root :
        Output root for ExtractHLAread; each sample gets a subdirectory.
    """
    if not samples:
        print("[hla_typing] No BAM files found; nothing to do.")
        return

    print(f"[hla_typing] Found {len(samples)} samples for ExtractHLAread.")
    extract_root.mkdir(parents=True, exist_ok=True)

    # Dependency check and optional auto-install. Mirrors the CLI wrapper.
    try:
        ensure_extract_deps(auto_install=True)
    except RuntimeError as exc:
        print(
            "[hla_typing] ERROR while checking/installing dependencies for "
            f"ExtractHLAread:\n{exc}",
            file=sys.stderr,
        )
        sys.exit(1)

    total = len(samples)
    print("[hla_typing] Starting ExtractHLAread for all samples...")

    for idx, (sample_id, bam_path) in enumerate(
        tqdm(samples, desc="[ExtractHLAread]", unit="sample", total=total),
        start=1,
    ):
        sample_outdir = extract_root / sample_id
        sample_outdir.mkdir(parents=True, exist_ok=True)
 
        tqdm.write(f"[hla_typing] [ExtractHLAread] Sample {sample_id} ({idx}/{total})")
        run_extraction(
            sample_id=sample_id,
            bam_path=bam_path,
            ref=ref,
            outdir=sample_outdir,
        )

# ---------------------------------------------------------------------------
# SpecHLA phase: locate FASTQs and call run_spechla_rnaseq
# ---------------------------------------------------------------------------

def find_fastqs_for_sample(sample_dir: Path, sample_id: str) -> Tuple[Path, Path]:
    """
    Find R1 and R2 FASTQ.gz files for a sample under ExtractHLAread/<sample_id>.

    The function searches for files named "*_1.fq.gz" and "*_2.fq.gz"
    and prefers matches whose basename starts with the sample id. If no
    such preference matches exist, it falls back to the first candidate.
    """
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"Sample directory not found: {sample_dir}")

    r1_candidates = sorted(sample_dir.glob("*_1.fq.gz"))
    r2_candidates = sorted(sample_dir.glob("*_2.fq.gz"))

    if not r1_candidates or not r2_candidates:
        raise FileNotFoundError(
            f"Could not find *_1.fq.gz / *_2.fq.gz under {sample_dir} "
            f"(sample {sample_id})."
        )

    def choose(cands):
        # Prefer files whose name starts with the sample id.
        for p in cands:
            if p.name.startswith(sample_id):
                return p
        return cands[0]

    r1 = choose(r1_candidates)
    r2 = choose(r2_candidates)
    return r1, r2


def run_spechla_phase(
    samples: List[Tuple[str, Path]],
    threads: int,
    extract_root: Path,
    spechla_outdir: Path,
) -> None:
    """
    Run SpecHLA in RNA-seq mode for all samples.

    Parameters
    ----------
    samples :
        List of (sample_id, bam_path).
    threads :
        Number of threads for SpecHLA (-j).
    extract_root :
        Root directory that holds ExtractHLAread/<sample_id>/... with FASTQs.
    spechla_outdir :
        Root directory for SpecHLA results.
        The underlying SpecHLA workflow is assumed to create one folder
        per sample named after -n (sample id).
    """
    if not samples:
        print("[hla_typing] No samples to process in SpecHLA phase.")
        return

    spechla_outdir.mkdir(parents=True, exist_ok=True)

    # Set up SpecHLA environment once.
    spec_hla_root = detect_spec_hla_root()
    ensure_python_deps()
    ensure_external_tools(spec_hla_root)
    ensure_spechap_built(spec_hla_root, threads)

    bowtie2_build_path = detect_bowtie2_build()
    drb_ref_relpath = os.path.join(
        "db", "ref", "hla_gen.format.filter.extend.DRB.no26789.fasta"
    )
    ensure_bowtie2_index(spec_hla_root, bowtie2_build_path, drb_ref_relpath)

    total = len(samples)
    print(f"[hla_typing] Starting SpecHLA for {total} samples...")

    for idx, (sample_id, _bam_path) in enumerate(
        tqdm(samples, desc="[SpecHLA]", unit="sample", total=total),
        start=1,
    ):
        sample_dir = extract_root / sample_id
        r1, r2 = find_fastqs_for_sample(sample_dir, sample_id)

        tqdm.write(f"[hla_typing] [SpecHLA] Sample {sample_id} ({idx}/{total})")
        # SpecHLA_RNAseq.sh is expected to put results into
        #   <spechla_outdir>/<sample_id>/
        # (this is the default layout of the upstream workflow).
        run_spechla_rnaseq(
            spec_hla_root=spec_hla_root,
            sample_name=sample_id,        # -n
            read1=str(r1),                # -1
            read2=str(r2),                # -2
            outdir=str(spechla_outdir),   # -o
            threads=threads,              # -j
        )

# ---------------------------------------------------------------------------
# Merge per-sample HLA result tables
# ---------------------------------------------------------------------------

def merge_hla_results(
    samples: List[Tuple[str, Path]],
    spechla_outdir: Path,
    outdir: Path,
) -> None:
    """
    Merge per-sample HLA result tables into a single file.

    For each sample_id, we look for:
        <spechla_outdir>/<sample_id>/hla.result.txt
    and, if missing,
        <spechla_outdir>/<sample_id>/hla.results.txt

    The file format is assumed to follow the example provided:

        # version: IPD-IMGT/HLA 3.38.0
        Sample  HLA_A_1 HLA_A_2 ...

    The merged file will contain:
        - The version line and header line from the first sample that has
          a non-empty result file.
        - All data lines from all samples (non-empty, non-comment lines
          after the header).

    The output is written as:
        <outdir>/hla_result_merged.txt
    """
    sample_result_files: List[Tuple[str, Path]] = []

    for sample_id, _bam_path in samples:
        sample_dir = spechla_outdir / sample_id
        candidates = [
            sample_dir / "hla.result.txt",
            sample_dir / "hla.results.txt",
        ]
        found = None
        for path in candidates:
            if path.is_file():
                found = path
                break
        if found is None:
            print(
                f"[hla_typing] WARNING: no hla.result(s).txt found for sample "
                f"{sample_id} under {sample_dir}; skipping this sample.",
                file=sys.stderr,
            )
            continue
        sample_result_files.append((sample_id, found))

    if not sample_result_files:
        print(
            "[hla_typing] WARNING: no HLA result files were found; "
            "merged table will not be created.",
            file=sys.stderr,
        )
        return

    merged_path = outdir / "hla_result_merged.txt"
    print(f"[hla_typing] Merging HLA result tables into: {merged_path}")

    with merged_path.open("w", encoding="utf-8") as out_f:
        wrote_header = False

        for sample_id, hla_file in sample_result_files:
            with hla_file.open("r", encoding="utf-8") as f:
                # Strip trailing newlines but keep empty-line information via strip check
                lines = [line.rstrip("\n") for line in f]

            # Remove leading empty lines
            lines = [ln for ln in lines if ln.strip() != ""]
            if not lines:
                print(
                    f"[hla_typing] WARNING: HLA result file for sample "
                    f"{sample_id} is empty: {hla_file}",
                    file=sys.stderr,
                )
                continue

            idx = 0
            local_version = None
            # Optional version/comment lines at the top (usually exactly one)
            while idx < len(lines) and lines[idx].startswith("#"):
                if local_version is None:
                    local_version = lines[idx]
                idx += 1

            if idx >= len(lines):
                print(
                    f"[hla_typing] WARNING: HLA result file for sample "
                    f"{sample_id} does not contain a header line: {hla_file}",
                    file=sys.stderr,
                )
                continue

            local_header = lines[idx]
            idx += 1

            if not wrote_header:
                # Write version line (if any) and header once.
                if local_version is not None:
                    out_f.write(local_version + "\n")
                out_f.write(local_header + "\n")
                wrote_header = True

            # Remaining lines are data lines.
            for data_line in lines[idx:]:
                if data_line.strip() == "":
                    continue
                # Avoid accidentally duplicating header lines
                if data_line.startswith("Sample\t") or data_line.startswith("Sample "):
                    continue
                out_f.write(data_line + "\n")

    print(f"[hla_typing] HLA result tables merged successfully.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the CLI argument parser for the hla_typing subcommand.
    """
    parser = argparse.ArgumentParser(
        prog="hla_typing",
        description=(
            "Batch HLA typing: run ExtractHLAread + SpecHLA on all BAM files "
            "in a directory, then merge per-sample results."
        ),
    )
    parser.add_argument(
        "-b",
        "--bam-dir",
        dest="bam_dir",
        required=True,
        help="Directory containing BAM files.",
    )
    parser.add_argument(
        "-r",
        "--ref",
        dest="ref",
        required=True,
        choices=["hg19", "hg38"],
        help="Reference genome passed to ExtractHLAread (hg19 or hg38).",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        help=(
            "Root output directory. "
            "Will create subfolders 'ExtractHLAread' and 'SpecHLA'."
        ),
    )
    parser.add_argument(
        "-j",
        "--threads",
        dest="threads",
        type=int,
        default=8,
        help="Number of threads for SpecHLA (-j). Default: 8.",
    )
    return parser


def main(argv: List[str] | None = None) -> None:
    """
    Entry point for the hla_typing subcommand.

    Parameters
    ----------
    argv :
        Optional list of arguments; if None, sys.argv[1:] is used.
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = build_arg_parser()
    args = parser.parse_args(argv)

    bam_dir = Path(os.path.expanduser(args.bam_dir)).resolve()
    if not bam_dir.is_dir():
        parser.error(f"BAM directory does not exist: {bam_dir}")

    outdir = Path(os.path.expanduser(args.outdir)).resolve()
    extract_root = outdir / "ExtractHLAread"
    spechla_outdir = outdir / "SpecHLA"

    samples = collect_samples(bam_dir)
    if not samples:
        parser.error(f"No .bam files found under directory: {bam_dir}")

    print(f"[hla_typing] Using BAM dir : {bam_dir}")
    print(f"[hla_typing] Output root  : {outdir}")
    print(f"[hla_typing] Reference    : {args.ref}")
    print(f"[hla_typing] Threads      : {args.threads}")
    print(f"[hla_typing] Detected {len(samples)} sample(s).")

    # 1) Run ExtractHLAread for all samples.
    run_extract_phase(samples, args.ref, extract_root)

    # 2) Run SpecHLA for all samples.
    run_spechla_phase(samples, args.threads, extract_root, spechla_outdir)

    # 3) Merge per-sample HLA result tables.
    merge_hla_results(samples, spechla_outdir, outdir)

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