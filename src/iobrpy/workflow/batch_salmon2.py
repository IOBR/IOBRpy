#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch RNA-seq quantification with Salmon.
- Scans a FASTQ directory for paired-end R1 files (suffix1), infers R2 by swapping "1"->"2" in the suffix.
- Runs `salmon quant` per sample with optional GTF.
- Skips a sample if both `task.complete` exists and `quant.sf` is non-empty.
- Shows a global progress bar across samples (one step per finished sample).
- Prints the final `quant.sf` path for each processed/skipped sample.
"""

import os
import random
import subprocess
from multiprocessing import Pool

# Progress bar (optional dependency)
try:
    from tqdm.auto import tqdm
except ImportError:  # graceful no-op if tqdm is unavailable
    def tqdm(x, **kwargs):
        return x

def _print_iobrpy_banner():
    """Pretty trailer banner at the very end of the run."""
    print("   ")
    try:
        from iobrpy.utils.print_colorful_message import print_colorful_message
        print_colorful_message("#########################################################", "blue")
        print_colorful_message(" IOBRpy: Immuno-Oncology Biological Research using Python ", "cyan")
        print_colorful_message(" If you encounter any issues, please report them at ", "cyan")
        print_colorful_message(" https://github.com/IOBR/IOBRpy/issues ", "cyan")
        print_colorful_message("#########################################################", "blue")
    except Exception:
        print("#########################################################")
        print(" IOBRpy: Immuno-Oncology Biological Research using Python ")
        print(" If you encounter any issues, please report them at ")
        print(" https://github.com/IOBR/IOBRpy/issues ")
        print("#########################################################")
    print(" Author: Haonan Huang, Dongqiang Zeng")
    print(" Email: interlaken@smu.edu.cn ")
    try:
        from iobrpy.utils.print_colorful_message import print_colorful_message as _pcm
        _pcm("#########################################################", "blue")
    except Exception:
        print("#########################################################")
    print("   ")

def process_sample(f1, path_out, index, suffix1, num_threads, gtf=None):
    """
    Process a single sample using Salmon for quantification.
    Always ends by printing: 'Saved to: <quant.sf>' for this sample.
    Returns a small status tuple for the parent progress bar.
    """
    # Infer R2 filename by swapping "1"->"2" in the provided suffix (e.g., _1.fastq.gz -> _2.fastq.gz)
    suffix2 = suffix1.replace("1", "2")
    f2 = f1.replace(suffix1, suffix2)

    # Derive sample ID from file name
    sample_id = os.path.basename(f1).replace(suffix1, "")
    sample_dir = os.path.join(path_out, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    quant_path = os.path.join(sample_dir, "quant.sf")
    complete_flag = os.path.join(sample_dir, "task.complete")

    # Skip logic
    if os.path.exists(complete_flag) and os.path.exists(quant_path) and os.path.getsize(quant_path) > 0:
        print(f"Skipping {sample_id} - task.complete and quant.sf exist and are not empty.")
        print(f"Saved to: {quant_path}", flush=True)
        return (sample_id, "skipped")

    # Build salmon command
    cmd = [
        "salmon", "quant",
        "-i", str(index),
        "-l", "ISF",
        "--gcBias",
        "-1", str(f1),
        "-2", str(f2),
        "-p", str(num_threads),
        "-o", str(sample_dir),
        "--validateMappings",
    ]
    if gtf:
        cmd.extend(["-g", str(gtf)])

    # Execute
    print(f"Processing {sample_id}...")
    subprocess.run(" ".join(cmd), shell=True, check=True)

    # Mark completion and print final path
    open(complete_flag, "w").close()
    print(f"Saved to: {quant_path}", flush=True)
    return (sample_id, "done")

def _pool_wrapper(args):
    """Helper so Pool can pass a single tuple argument."""
    return process_sample(*args)

def main():
    """Parse arguments, dispatch Salmon jobs with a progress bar, print trailer banner at the end."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Batch-run Salmon quantification over paired-end FASTQ files."
    )
    parser.add_argument(
        "--index", type=str, required=True,
        help="Path to Salmon index directory (for `salmon quant -i`)."
    )
    parser.add_argument(
        "--path_fq", type=str, required=True,
        help="Directory containing FASTQ files to scan."
    )
    parser.add_argument(
        "--path_out", type=str, required=True,
        help="Output directory; each sample will have its own subdirectory."
    )
    parser.add_argument(
        "--suffix1", type=str, default="_1.fastq.gz",
        help="Suffix of R1 FASTQ files (R2 is inferred by replacing '1' with '2')."
    )
    parser.add_argument(
        "--batch_size", type=int, default=1,
        help="Number of samples to process concurrently (process count)."
    )
    parser.add_argument(
        "--num_threads", type=int, default=8,
        help="Threads per Salmon process (`salmon quant -p`)."
    )
    parser.add_argument(
        "--gtf", type=str, default=None,
        help="Optional GTF file path passed to Salmon with `-g`."
    )
    args = parser.parse_args()

    # Prepare I/O
    os.makedirs(args.path_out, exist_ok=True)
    files_1 = sorted(
        os.path.join(args.path_fq, f)
        for f in os.listdir(args.path_fq)
        if f.endswith(args.suffix1)
    )
    random.shuffle(files_1)

    if not files_1:
        print("No FASTQ files found matching the given suffix.")
        _print_iobrpy_banner()
        return

    # Run with a sample-level progress bar (one update per finished sample)
    with Pool(processes=args.batch_size) as pool:
        tup_iter = ((f1, args.path_out, args.index, args.suffix1, args.num_threads, args.gtf) for f1 in files_1)
        for _ in tqdm(pool.imap_unordered(_pool_wrapper, tup_iter), total=len(files_1), desc="Samples", unit="sample"):
            pass

    _print_iobrpy_banner()

if __name__ == "__main__":
    main()
