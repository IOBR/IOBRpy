# src/iobrpy/workflow/trust4.py
# -*- coding: utf-8 -*-
"""
Wrapper for TRUST4 so that you can run:
    iobrpy trust4 [TRUST4 options...]

Features:
- Reuses -b for both single and batch:
    * -b <file.bam>          -> single-run (TRUST4 native semantics)
    * -b <directory_of_bams> -> batch over all *.bam in that directory (non-recursive)
- --fqdir: batch over paired FASTQ (prefix + _1/_2.fastq.gz).
- In batch modes (-b <dir> or --fqdir), CLI -o is treated as an OUTPUT ROOT directory.
  Each sample runs with:
      --od <OUT_ROOT>/<sample>
      -o  TRUST_<sample-prefix>
  so all files of one sample end up in <OUT_ROOT>/<sample>/ with prefix TRUST_<sample-prefix>.*.

After all TRUST4 runs finish (including single-run), this wrapper:
- Locates all *_report.tsv files under the chosen output directory
  (including subfolders, e.g. per-sample folders for bamdir/fqdir mode).
- Calls iobrpy.utils.process_immune_data_batch.process_immune_data_batch(...)
  to generate two CSV files in the output directory:
      trust4_immdata.csv
      trust4_immune_indices.csv
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import atexit
from typing import Dict, List, Optional, Tuple

from importlib.resources import files

# Optional immune post-processing (TRUST4 report -> immune indices)
try:  # pragma: no cover - optional dependency
    from iobrpy.utils.process_immune_data_batch import process_immune_data_batch
except Exception:  # pragma: no cover
    process_immune_data_batch = None  # type: ignore

# tqdm progress bar (with graceful fallback if not installed)
try:
    from tqdm import tqdm  # type: ignore
except Exception:  # pragma: no cover
    class _SimpleTqdm:
        """Minimal fallback when tqdm is unavailable."""

        def __init__(self, total=None, desc=None, unit=None, dynamic_ncols=None):
            self.total = total or 0
            self.desc = desc or "progress"
            self.n = 0

        def __enter__(self):
            print(f"[{self.desc}] total={self.total}")
            return self

        def update(self, n: int = 1):
            self.n += n
            print(f"[{self.desc}] {self.n}/{self.total}", flush=True)

        def __exit__(self, exc_type, exc, tb):
            pass

    def tqdm(*args, **kwargs):  # type: ignore
        return _SimpleTqdm(*args, **kwargs)


DEFAULT_F = "hg38_bcrtcr.fa"
DEFAULT_REF = "human_IMGT+C.fa"
RES_PKG = "iobrpy.resources"

_TRUST4_TMP_DIRS: List[str] = []

def _cleanup_tmp_dirs():
    for d in list(_TRUST4_TMP_DIRS):
        try:
            if d and os.path.isdir(d) and os.path.basename(d).startswith("iobrpy_trust4_"):
                shutil.rmtree(d, ignore_errors=True)
        except Exception:
            pass
    _TRUST4_TMP_DIRS.clear()

atexit.register(_cleanup_tmp_dirs)

def _extract_resource_to_tmp(name: str) -> Optional[str]:
    """Read a resource file from iobrpy.resources and write it to a temp file."""
    try:
        data = files(RES_PKG).joinpath(name).read_bytes()
    except Exception:
        return None
    tmpdir = tempfile.mkdtemp(prefix="iobrpy_trust4_")
    _TRUST4_TMP_DIRS.append(tmpdir)
    outp = os.path.join(tmpdir, name)
    with open(outp, "wb") as f:
        f.write(data)
    return outp

def _append_opt(cmd: List[str], flag: str, value, is_bool: bool = False):
    """Append a CLI option to command if value is present."""
    if is_bool:
        if value:
            cmd.append(flag)
    else:
        if value is not None:
            cmd.extend([flag, str(value)])


def _list_bams(bamdir: str) -> List[str]:
    """List *.bam files (non-recursive) in bamdir."""
    if not os.path.isdir(bamdir):
        raise FileNotFoundError(f"-b path is a directory but not found: {bamdir}")
    out: List[str] = []
    for fn in sorted(os.listdir(bamdir)):
        if fn.endswith(".bam"):
            p = os.path.join(bamdir, fn)
            if os.path.isfile(p):
                out.append(p)
    return out


def _infer_sample_from_fastq_name(filename: str) -> Optional[Tuple[str, str]]:
    """
    Return (sample_prefix, end_type) where end_type is 'R1' or 'R2',
    based on suffixes '_1.fastq.gz' / '_2.fastq.gz'.
    If not matched, return None.
    """
    if filename.endswith("_1.fastq.gz"):
        return (filename[: -len("_1.fastq.gz")], "R1")
    if filename.endswith("_2.fastq.gz"):
        return (filename[: -len("_2.fastq.gz")], "R2")
    return None


def _pair_fastqs(fqdir: str) -> Dict[str, Tuple[str, str]]:
    """
    Find paired FASTQs in fqdir and return mapping:
        sample_prefix -> (R1_path, R2_path)

    Rules:
    - Exactly one R1 and one R2 per sample prefix.
    - Raises ValueError if a sample has multiple R1 or multiple R2 (ambiguous),
      or if unmatched pairs remain.
    """
    if not os.path.isdir(fqdir):
        raise FileNotFoundError(f"--fqdir not found: {fqdir}")

    r1_map: Dict[str, List[str]] = {}
    r2_map: Dict[str, List[str]] = {}

    for fn in sorted(os.listdir(fqdir)):
        if not fn.endswith(".fastq.gz"):
            continue
        parsed = _infer_sample_from_fastq_name(fn)
        if parsed is None:
            continue
        prefix, which = parsed
        full = os.path.join(fqdir, fn)
        if which == "R1":
            r1_map.setdefault(prefix, []).append(full)
        else:
            r2_map.setdefault(prefix, []).append(full)

    pairs: Dict[str, Tuple[str, str]] = {}
    samples = sorted(set(r1_map.keys()) | set(r2_map.keys()))
    missing: List[Tuple[str, int, int]] = []
    multi: List[Tuple[str, int, int]] = []
    for s in samples:
        r1s = r1_map.get(s, [])
        r2s = r2_map.get(s, [])
        if len(r1s) != 1 or len(r2s) != 1:
            if len(r1s) == 0 or len(r2s) == 0:
                missing.append((s, len(r1s), len(r2s)))
            else:
                multi.append((s, len(r1s), len(r2s)))
            continue
        pairs[s] = (r1s[0], r2s[0])

    errs: List[str] = []
    if missing:
        errs.append(
            "Unmatched sample pairs (need exactly one R1 and one R2): "
            + ", ".join([f"{s}[R1={n1},R2={n2}]" for s, n1, n2 in missing])
        )
    if multi:
        errs.append(
            "Ambiguous samples with multiple R1 or R2 (please merge lanes or reduce to one pair per sample): "
            + ", ".join([f"{s}[R1={n1},R2={n2}]" for s, n1, n2 in multi])
        )
    if errs:
        raise ValueError("; ".join(errs))

    return pairs


def _build_common_options(args, f_path: str, ref_path: Optional[str]) -> List[str]:
    """Build the subset of TRUST4 options common to all runs (threads, refs, barcode etc.)."""
    cmd: List[str] = []
    _append_opt(cmd, "-f", f_path)
    if ref_path:
        _append_opt(cmd, "--ref", ref_path)
    _append_opt(cmd, "-t", args.t)
    _append_opt(cmd, "-k", args.k)

    _append_opt(cmd, "--barcode", args.barcode)
    _append_opt(cmd, "--barcodeLevel", args.barcodeLevel)
    _append_opt(cmd, "--barcodeWhitelist", args.barcodeWhitelist)
    _append_opt(cmd, "--barcodeTranslate", args.barcodeTranslate)
    _append_opt(cmd, "--UMI", args.UMI)
    _append_opt(cmd, "--readFormat", args.readFormat)

    _append_opt(cmd, "--repseq", args.repseq, is_bool=True)
    _append_opt(cmd, "--contigMinCov", args.contigMinCov)
    _append_opt(cmd, "--minHitLen", args.minHitLen)
    _append_opt(cmd, "--mateIdSuffixLen", args.mateIdSuffixLen)
    _append_opt(cmd, "--skipMateExtension", args.skipMateExtension, is_bool=True)
    _append_opt(cmd, "--abnormalUnmapFlag", args.abnormalUnmapFlag, is_bool=True)
    _append_opt(cmd, "--assembleWithRef", args.assembleWithRef, is_bool=True)
    _append_opt(cmd, "--noExtraction", args.noExtraction, is_bool=True)
    _append_opt(cmd, "--outputReadAssignment", args.outputReadAssignment, is_bool=True)

    _append_opt(cmd, "--stage", args.stage)
    _append_opt(cmd, "--clean", args.clean)
    return cmd


def _infer_single_output_root(args) -> str:
    """
    Infer the directory where TRUST4 wrote its outputs for single-run mode.

    Preference:
    1. args.od (explicit output directory)
    2. args.o if it is an existing directory
    3. directory component of args.o if that directory exists
    4. current working directory
    """
    if getattr(args, "od", None):
        return args.od
    o_val = getattr(args, "o", None)
    if o_val and os.path.isdir(o_val):
        return o_val
    if o_val:
        o_dir = os.path.dirname(o_val)
        if o_dir and os.path.isdir(o_dir):
            return o_dir
    return os.getcwd()


def _collect_report_files(root: str) -> List[str]:
    """
    Recursively collect all *_report.tsv files under root (including subdirectories).
    """
    report_files: List[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.endswith("_report.tsv"):
                report_files.append(os.path.join(dirpath, fn))
    report_files.sort()
    return report_files


def _derive_sample_name_from_path(src: str, output_root: str) -> str:
    """
    Derive a clean sample name from a TRUST4 *_report.tsv path.

    Rules:
    - If the report is in a subdirectory of output_root, use the top-level
      subdirectory name as the sample name, e.g.
          /out/SRR3474759_T01/TRUST_..._report.tsv -> SRR3474759_T01
    - If the report is directly under output_root, use the file name stem.
    - Strip a leading 'TRUST_' prefix and trailing alignment-related suffixes
      such as '_Aligned' or '.Aligned.sortedByCoord.out'.
    """
    rel = os.path.relpath(src, output_root)
    parts = rel.split(os.sep)
    if len(parts) > 1:
        # Sample lives in a per-sample subdirectory
        sample = parts[0]
    else:
        base = os.path.basename(src)
        if base.endswith("_report.tsv"):
            base = base[: -len("_report.tsv")]
        sample = base

    # Remove TRUST_ prefix if present
    if sample.startswith("TRUST_"):
        sample = sample[len("TRUST_") :]

    # Remove suffixes related to STAR-style alignment names
    sample = re.sub(
        r"([._]Aligned)\.sortedByCoord\.out$",
        "",
        sample,
        flags=re.IGNORECASE,
    )
    sample = re.sub(
        r"\.sortedByCoord\.out$",
        "",
        sample,
        flags=re.IGNORECASE,
    )

    return sample


def _run_immune_postprocessing(output_root: Optional[str]) -> None:
    """
    Run process_immune_data_batch on all *_report.tsv under output_root.

    For bamdir/fqdir batch modes, this includes reports in per-sample subfolders.
    To satisfy process_immune_data_batch(path_to_result, ...), which expects a
    single directory containing *_report.tsv files, we:
      1) find all *_report.tsv under output_root (recursively),
      2) stage them into a temporary directory as <sample>_report.tsv, where
         <sample> is derived from the path using _derive_sample_name_from_path,
      3) run process_immune_data_batch on that temporary directory.

    The resulting CSV files are written into output_root.
    """
    if not output_root:
        return

    if process_immune_data_batch is None:
        # process_immune_data_batch is optional; skip if unavailable
        print(
            "[IOBRpy|trust4] process_immune_data_batch not available; skip immune post-processing.",
            file=sys.stderr,
        )
        return

    output_root = os.path.abspath(output_root)
    if not os.path.isdir(output_root):
        print(
            f"[IOBRpy|trust4] Output root for immune post-processing is not a directory: {output_root}",
            file=sys.stderr,
        )
        return

    report_files = _collect_report_files(output_root)
    if not report_files:
        print(
            f"[IOBRpy|trust4] No *_report.tsv files found under {output_root!r}; skip immune post-processing.",
            file=sys.stderr,
        )
        return

    # Create a temporary directory containing symlinks/copies to all report TSVs.
    # File names in this staging directory are <sample>_report.tsv so that
    # downstream tools use the clean sample names in their outputs (both
    # trust4_immdata.csv and trust4_immune_indices.csv).
    staging_dir = tempfile.mkdtemp(prefix="iobrpy_trust4_reports_")
    name_counts: Dict[str, int] = {}
    for src in report_files:
        sample = _derive_sample_name_from_path(src, output_root)
        base_name = f"{sample}_report.tsv"

        # Avoid name collisions by appending a numeric suffix if needed
        count = name_counts.get(base_name, 0)
        name_counts[base_name] = count + 1
        if count == 0:
            dst_name = base_name
        else:
            dst_name = f"{sample}_{count}_report.tsv"

        dst = os.path.join(staging_dir, dst_name)

        try:
            os.symlink(src, dst)
        except Exception:
            try:
                shutil.copy2(src, dst)
            except Exception as e:
                print(
                    f"[IOBRpy|trust4] WARNING: could not stage {src} -> {dst}: {e}",
                    file=sys.stderr,
                )

    immdata_out_csv = os.path.join(output_root, "trust4_immdata.csv")
    indices_out_csv = os.path.join(output_root, "trust4_immune_indices.csv")

    try:
        print(
            f"[IOBRpy|trust4] Running immune post-processing on TRUST4 reports under: {output_root}",
            flush=True,
        )
        process_immune_data_batch(staging_dir, immdata_out_csv, indices_out_csv)
        print(
            f"[IOBRpy|trust4] Immune clone-level data saved to: {immdata_out_csv}",
            flush=True,
        )
        print(
            f"[IOBRpy|trust4] Immune diversity indices saved to: {indices_out_csv}",
            flush=True,
        )
    except FileNotFoundError as e:
        print(
            f"[IOBRpy|trust4] No *_report.tsv files found for immune post-processing: {e}",
            file=sys.stderr,
        )
    except Exception as e:
        print(
            f"[IOBRpy|trust4] Failed to run immune post-processing: {e}",
            file=sys.stderr,
        )
    finally:
        # Clean up staging directory
        try:
            shutil.rmtree(staging_dir)
        except Exception:
            pass


def main(argv: Optional[List[str]] = None):
    runner = os.environ.get("TRUST4_BIN", "run-trust4")
    if shutil.which(runner) is None:
        print(
            f"[IOBRpy|trust4] Cannot find TRUST4 executable '{runner}'. "
            f"Install TRUST4 (e.g. conda install -c bioconda trust4) or set TRUST4_BIN.",
            file=sys.stderr,
        )
        sys.exit(127)

    p = argparse.ArgumentParser(
        prog="iobrpy trust4",
        description="Run TRUST4 (TCR/BCR reconstruction) with IOBRpy defaults for -f/--ref if unspecified.",
        add_help=True,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Inputs
    p.add_argument(
        "-b",
        dest="bam",
        help="Path to BAM file OR a directory of *.bam (batch mode)",
    )
    p.add_argument("-1", dest="r1", help="Paired-end read 1 FASTQ")
    p.add_argument("-2", dest="r2", help="Paired-end read 2 FASTQ")
    p.add_argument("-u", dest="ru", help="Single-end read FASTQ")

    # Batch (FASTQ) mode
    p.add_argument(
        "--fqdir",
        dest="fqdir",
        help="Directory containing paired FASTQs (*.fastq.gz) with suffixes _1/_2",
    )

    # Reference files
    p.add_argument(
        "-f",
        dest="f",
        help="FASTA with coordinates/sequences of V/D/J/C genes",
    )
    p.add_argument(
        "--ref",
        dest="ref",
        help="Detailed V/D/J/C reference FASTA (e.g., IMGT)",
    )

    # General
    p.add_argument(
        "-o",
        dest="o",
        help="Single-run: output prefix; Batch: OUTPUT ROOT directory",
    )
    p.add_argument(
        "--od",
        dest="od",
        help="Output directory (TRUST4 native; ignored in batch modes)",
    )
    p.add_argument("-t", dest="t", type=int, help="Threads")
    p.add_argument("-k", dest="k", type=int, help="Starting k-mer size")

    # Single-cell / barcodes / UMI
    p.add_argument("--barcode", dest="barcode", help="Barcode field/file")
    p.add_argument(
        "--barcodeLevel",
        dest="barcodeLevel",
        choices=["cell", "molecule"],
        help="Barcode level",
    )
    p.add_argument(
        "--barcodeWhitelist",
        dest="barcodeWhitelist",
        help="Barcode whitelist",
    )
    p.add_argument(
        "--barcodeTranslate",
        dest="barcodeTranslate",
        help="Barcode translate file",
    )
    p.add_argument("--UMI", dest="UMI", help="UMI field/file")
    p.add_argument(
        "--readFormat",
        dest="readFormat",
        help="Format spec for reads/barcodes/UMI",
    )

    # Modes & fine controls
    p.add_argument(
        "--repseq",
        action="store_true",
        help="Bulk non-UMI-based TCR/BCR-seq",
    )
    p.add_argument(
        "--contigMinCov",
        type=int,
        help="Min bases covered by reads for contigs",
    )
    p.add_argument(
        "--minHitLen",
        type=int,
        help="Minimal hit length for valid overlap",
    )
    p.add_argument(
        "--mateIdSuffixLen",
        type=int,
        help="Mate suffix length in read id",
    )
    p.add_argument(
        "--skipMateExtension",
        action="store_true",
        help="Do not extend assemblies with mate info",
    )
    p.add_argument(
        "--abnormalUnmapFlag",
        action="store_true",
        help="Unmapped read-pair flag is nonconcordant",
    )
    p.add_argument(
        "--assembleWithRef",
        action="store_true",
        help="Assemble with --ref file",
    )
    p.add_argument(
        "--noExtraction",
        action="store_true",
        help="Directly use provided FASTQ files",
    )
    p.add_argument(
        "--outputReadAssignment",
        action="store_true",
        help="Output read assignment to *_assign.out",
    )

    # Pipeline stage / cleanup
    p.add_argument(
        "--stage",
        type=int,
        choices=[0, 1, 2, 3],
        help="0:begin; 1:assembly; 2:annotation; 3:report",
    )
    p.add_argument(
        "--clean",
        type=int,
        choices=[0, 1, 2],
        help="0:keep all; 1:clean intermediates; 2:only keep AIRR files",
    )

    args, unknown = p.parse_known_args(argv)

    # Determine mode
    batch_fq = args.fqdir is not None
    batch_bamdir = False
    single_bam = False

    if args.bam is not None:
        if os.path.isdir(args.bam):
            batch_bamdir = True
        elif os.path.isfile(args.bam):
            single_bam = True
        else:
            p.error(f"-b path not found: {args.bam}")

    has_pe = args.r1 is not None and args.r2 is not None
    has_se = args.ru is not None

    # Conflicts
    if batch_fq and (single_bam or batch_bamdir or has_pe or has_se):
        p.error("Conflicting inputs: --fqdir cannot be combined with -b/-1/-2/-u.")
    if batch_bamdir and (has_pe or has_se):
        p.error(
            "Conflicting inputs: -b <dir> (batch mode) cannot be combined with -1/-2/-u."
        )

    # Validate at least one input mode
    if not (batch_fq or batch_bamdir or single_bam or has_pe or has_se):
        p.error(
            "Provide one of: -b <BAM> | -b <DIR_OF_BAM> | -1 R1 -2 R2 | -u READS | --fqdir DIR"
        )

    # Resolve -f and --ref defaults if not provided
    f_path = args.f
    ref_path = args.ref
    if f_path is None:
        f_path = _extract_resource_to_tmp(DEFAULT_F)
        if f_path is None:
            p.error(
                f"Default -f resource '{DEFAULT_F}' not found in {RES_PKG}. "
                f"Place the file under iobrpy/resources/ or pass -f explicitly."
            )
    if ref_path is None:
        ref_path = _extract_resource_to_tmp(DEFAULT_REF)
        if ref_path is None:
            print(
                f"[IOBRpy|trust4] Warning: default --ref '{DEFAULT_REF}' not found in {RES_PKG}. "
                f"Continuing without --ref (less recommended).",
                file=sys.stderr,
            )

    # Common options (shared for all runs) -- DO NOT include --od here
    common_opts = _build_common_options(args, f_path, ref_path)
    common_opts.extend(unknown)  # pass-through any extra options

    # -------------------------
    # Batch over directory of BAMs via -b <dir>
    # -------------------------
    if batch_bamdir:
        if not args.o:
            p.error("When using -b <DIR>, please provide -o as an OUTPUT ROOT directory.")
        bams = _list_bams(args.bam)
        if not bams:
            p.error(f"No *.bam found under directory: {args.bam}")

        print(f"[IOBRpy|trust4] Found {len(bams)} BAM samples in: {args.bam}")
        overall_rc = 0
        with tqdm(
            total=len(bams),
            desc="TRUST4 (BAM batch)",
            unit="sample",
            dynamic_ncols=True,
        ) as pbar:
            for bam in bams:
                sample = os.path.splitext(os.path.basename(bam))[0]
                # Folder name: drop trailing _Aligned.sortedByCoord.out (or .Aligned...); case-insensitive
                folder = re.sub(
                    r"([._]Aligned)\.sortedByCoord\.out$",
                    "",
                    sample,
                    flags=re.IGNORECASE,
                )
                sample_dir = os.path.join(args.o, folder)
                os.makedirs(sample_dir, exist_ok=True)

                # Per-sample done flag: if this file exists, skip re-running TRUST4
                done_flag = os.path.join(sample_dir, f"{folder}.TRUST4.done")
                if os.path.exists(done_flag):
                    print(
                        f"[IOBRpy|trust4] Skip sample {folder}: found done flag {done_flag}",
                        flush=True,
                    )
                    # Still advance the progress bar so that the batch summary is correct
                    pbar.update(1)
                    continue

                # File prefix: keep _Aligned, only drop .sortedByCoord.out
                prefix_base = re.sub(
                    r"\.sortedByCoord\.out$",
                    "",
                    sample,
                    flags=re.IGNORECASE,
                )
                prefix = f"TRUST_{prefix_base}"

                cmd: List[str] = [runner]
                _append_opt(cmd, "-b", bam)
                cmd.extend(common_opts)
                _append_opt(cmd, "--od", sample_dir)
                _append_opt(cmd, "-o", prefix)

                print(
                    f"[IOBRpy|trust4] Running ({folder}):",
                    " ".join(map(str, cmd)),
                    flush=True,
                )
                try:
                    rc = subprocess.run(cmd, check=False).returncode
                except KeyboardInterrupt:
                    sys.exit(130)

                # If TRUST4 finishes successfully for this sample, create the done flag
                if rc == 0:
                    try:
                        with open(done_flag, "w") as f:
                            f.write(
                                f"SUCCESS: TRUST4 finished for sample {folder}\n"
                            )
                    except Exception as e:
                        # Do not fail the whole batch just because we could not write the flag
                        print(
                            f"[IOBRpy|trust4] WARNING: could not write done flag for {folder}: {e}",
                            file=sys.stderr,
                        )
                else:
                    overall_rc = rc

                if rc != 0:
                    overall_rc = rc
                pbar.update(1)

        # After all TRUST4 BAM batch runs, run immune post-processing on -o root
        _run_immune_postprocessing(args.o)
        sys.exit(overall_rc)

    # -------------------------
    # Batch over FASTQ directory
    # -------------------------
    if batch_fq:
        if not args.o:
            p.error("When using --fqdir, please provide -o as an OUTPUT ROOT directory.")
        try:
            pairs = _pair_fastqs(args.fqdir)
        except (FileNotFoundError, ValueError) as e:
            p.error(str(e))
        if not pairs:
            p.error(
                f"No paired FASTQ (*.fastq.gz with _1/_2) found in --fqdir: {args.fqdir}"
            )

        samples = sorted(pairs.keys())
        print(f"[IOBRpy|trust4] Found {len(samples)} FASTQ samples in: {args.fqdir}")
        overall_rc = 0
        with tqdm(
            total=len(samples),
            desc="TRUST4 (FASTQ batch)",
            unit="sample",
            dynamic_ncols=True,
        ) as pbar:
            for sample in samples:
                r1p, r2p = pairs[sample]
                # Folder name: drop trailing _Aligned.sortedByCoord.out if present
                folder = re.sub(
                    r"([._]Aligned)\.sortedByCoord\.out$",
                    "",
                    sample,
                    flags=re.IGNORECASE,
                )
                sample_dir = os.path.join(args.o, folder)
                os.makedirs(sample_dir, exist_ok=True)

                # Per-sample done flag: if this file exists, skip re-running TRUST4
                done_flag = os.path.join(sample_dir, f"{folder}.TRUST4.done")
                if os.path.exists(done_flag):
                    print(
                        f"[IOBRpy|trust4] Skip sample {folder}: found done flag {done_flag}",
                        flush=True,
                    )
                    pbar.update(1)
                    continue

                # File prefix: keep _Aligned, only drop .sortedByCoord.out
                prefix_base = re.sub(
                    r"\.sortedByCoord\.out$",
                    "",
                    sample,
                    flags=re.IGNORECASE,
                )
                prefix = f"TRUST_{prefix_base}"

                cmd = [runner]
                _append_opt(cmd, "-1", r1p)
                _append_opt(cmd, "-2", r2p)
                cmd.extend(common_opts)
                _append_opt(cmd, "--od", sample_dir)
                _append_opt(cmd, "-o", prefix)

                print(
                    f"[IOBRpy|trust4] Running ({folder}):",
                    " ".join(map(str, cmd)),
                    flush=True,
                )
                try:
                    rc = subprocess.run(cmd, check=False).returncode
                except KeyboardInterrupt:
                    sys.exit(130)

                # If TRUST4 finishes successfully for this sample, create the done flag
                if rc == 0:
                    try:
                        with open(done_flag, "w") as f:
                            f.write(
                                f"SUCCESS: TRUST4 finished for sample {folder}\n"
                            )
                    except Exception as e:
                        print(
                            f"[IOBRpy|trust4] WARNING: could not write done flag for {folder}: {e}",
                            file=sys.stderr,
                        )
                else:
                    overall_rc = rc
                
                pbar.update(1)

        # After all TRUST4 FASTQ batch runs, run immune post-processing on -o root
        _run_immune_postprocessing(args.o)
        sys.exit(overall_rc)

    # -------------------------
    # Single-run mode (BAM or FASTQ directly)
    # -------------------------
    cmd: List[str] = [runner]
    if single_bam:
        _append_opt(cmd, "-b", args.bam)
    if has_pe:
        _append_opt(cmd, "-1", args.r1)
        _append_opt(cmd, "-2", args.r2)
    if has_se:
        _append_opt(cmd, "-u", args.ru)

    cmd.extend(common_opts)
    _append_opt(cmd, "-o", args.o)
    _append_opt(cmd, "--od", args.od)  # single-run only; optional

    print("[IOBRpy|trust4] Running:", " ".join(map(str, cmd)), flush=True)
    try:
        proc = subprocess.run(cmd, check=False)
        rc = proc.returncode

        out_root = _infer_single_output_root(args)
        _run_immune_postprocessing(out_root)

        sys.exit(rc)
    except KeyboardInterrupt:
        sys.exit(130)


if __name__ == "__main__":
    # Support `python -m iobrpy.workflow.trust4 ...` or running the file directly.
    main()