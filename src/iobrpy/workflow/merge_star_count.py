#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge STAR gene count files into a single matrix.
- Expects multiple files ending with `_ReadsPerGene.out.tab` in a given folder.
- Reads column 0 (gene ID) and column 1 (counts), names the count column by sample (file stem).
- Concatenates columns across samples to form a wide matrix.
- Saves as <project>.STAR.count.tsv.gz in the same folder.
- Shows a progress bar while reading files concurrently.
"""

import os
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip

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

def _read_star_file(file_path: str) -> pd.DataFrame:
    """
    Read a STAR ReadsPerGene.out.tab file and return a Series/DataFrame with ID as index and counts as one column.
    The column name is the sample name derived from the file name (strip `_ReadsPerGene.out.tab`).
    """
    file_name = os.path.basename(file_path)
    sample_name = file_name.replace("_ReadsPerGene.out.tab", "")
    df = pd.read_csv(file_path, sep="\t", usecols=[0, 1], header=None, names=["ID", sample_name])
    return df.set_index("ID")

def _list_star_files(path: str):
    """List all files in `path` that end with `_ReadsPerGene.out.tab`."""
    return [
        os.path.join(path, f)
        for f in os.listdir(path)
        if f.endswith("_ReadsPerGene.out.tab")
    ]

def main():
    """Parse arguments, read/merge STAR outputs with a progress bar, save gzipped TSV, print trailer banner at the end."""
    parser = argparse.ArgumentParser(
        description="Merge STAR *_ReadsPerGene.out.tab files into a single gzipped TSV matrix."
    )
    parser.add_argument(
        "--path", metavar="PATH", type=str, required=True,
        help="Directory containing the STAR outputs to merge."
    )
    parser.add_argument(
        "--project", metavar="PROJECT", type=str, required=True,
        help="Output name prefix (creates <PROJECT>.STAR.count.tsv.gz under PATH)."
    )
    args = parser.parse_args()

    files = _list_star_files(args.path)

    if not files:
        print("No files ending with '_ReadsPerGene.out.tab' were found in the given path.")
        _print_iobrpy_banner()
        return

    # Concurrent reading with progress bar (I/O bound -> threads are effective)
    results = []
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(_read_star_file, fp): fp for fp in files}
        for fut in tqdm(as_completed(futures), total=len(files), desc="Reading STAR count files", unit="file"):
            results.append(fut.result())

    # Merge by index
    merged_df = pd.concat(results, axis=1)

    # Save gzipped TSV
    output_file = os.path.join(args.path, f"{args.project}.STAR.count.tsv.gz")
    print("Saving merged matrix...")
    with gzip.open(output_file, "wt", compresslevel=5) as f:
        merged_df.to_csv(f, sep="\t")

    # Brief info and final path
    print("Head of merged file:")
    print(merged_df.head())
    print("Number of rows:", merged_df.shape[0])
    print("Number of columns:", merged_df.shape[1])
    print(f"Saved to: {output_file}")

    _print_iobrpy_banner()

if __name__ == "__main__":
    main()
