#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge Salmon `quant.sf` files into two matrices (TPM and NumReads).

- Recursively finds all `quant.sf` under --path_salmon.
- Reads columns `Name`, `TPM`, `NumReads` from each `quant.sf`.
- Column names are sample names taken from the parent directory of each `quant.sf`.
- Merges with an outer join across all samples.
- Saves to `<project>_salmon_tpm.tsv.gz` and `<project>_salmon_count.tsv.gz` in --path_salmon.
- Shows a progress bar while loading files concurrently.
- Prints output file locations *BEFORE* the IOBRpy banner (per request).
"""

import os
import argparse
import gzip
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

# Optional progress bar
try:
    from tqdm.auto import tqdm
except ImportError:  # graceful no-op if tqdm is unavailable
    def tqdm(x, **kwargs):
        return x

def _print_iobrpy_banner():
    """Print the IOBRpy banner (after output paths)."""
    print(" ")
    try:
        from iobrpy.utils.print_colorful_message import print_colorful_message
        print_colorful_message("#########################################################", "blue")
        print_colorful_message(" IOBRpy: Immuno-Oncology Biological Research using Python ", "cyan")
        print_colorful_message(" If you encounter any issues, please report them at ", "cyan")
        print_colorful_message(" https://github.com/IOBR/IOBRpy/issues ", "cyan")
        print_colorful_message("#########################################################", "blue")
        print(" Author: Haonan Huang, Dongqiang Zeng")
        print(" Email: interlaken@smu.edu.cn ")
        print_colorful_message("#########################################################", "blue")
    except Exception:
        # Fallback without colors
        print("#########################################################")
        print(" IOBRpy: Immuno-Oncology Biological Research using Python ")
        print(" If you encounter any issues, please report them at ")
        print(" https://github.com/IOBR/IOBRpy/issues ")
        print("#########################################################")
        print(" Author: Haonan Huang, Dongqiang Zeng")
        print(" Email: interlaken@smu.edu.cn ")
        print("#########################################################")
    print(" ")

def _find_quant_files(root: str):
    """Recursively find all `quant.sf` under `root`."""
    hits = []
    for dirpath, dirnames, filenames in os.walk(root):
        if "quant.sf" in filenames:
            hits.append(os.path.join(dirpath, "quant.sf"))
    return sorted(hits)

def _load_one_quant(quant_path: str):
    """
    Load a single quant.sf and return two Series (TPM, Count) with index=Name and name=sample.
    - Salmon quant.sf header columns: Name Length EffectiveLength TPM NumReads
    """
    sample = os.path.basename(os.path.dirname(quant_path))
    df = pd.read_csv(quant_path, sep="\t", header=0, usecols=["Name", "TPM", "NumReads"])
    tpm = df.set_index("Name")["TPM"].rename(sample)
    cnt = df.set_index("Name")["NumReads"].rename(sample)
    return tpm, cnt

def main():
    parser = argparse.ArgumentParser(
        description="Merge Salmon quant.sf files into TPM and NumReads matrices."
    )
    parser.add_argument(
        "--path_salmon", type=str, required=True,
        help="Root directory containing per-sample Salmon outputs (searched recursively for quant.sf)."
    )
    parser.add_argument(
        "--project", type=str, required=True,
        help="Output file prefix. Will produce <project>_salmon_tpm.tsv.gz and <project>_salmon_count.tsv.gz."
    )
    parser.add_argument(
        "--num_processes", type=int, default=None,
        help="Number of threads used to load quant.sf files (I/O bound). Default: number of CPUs."
    )
    args = parser.parse_args()

    # Discover files
    quants = _find_quant_files(args.path_salmon)
    if not quants:
        print("No quant.sf files were found under the given --path_salmon.")
        _print_iobrpy_banner()
        return

    # Concurrent loading with progress bar
    tpm_list = []
    cnt_list = []
    max_workers = args.num_processes or os.cpu_count() or 4
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [ex.submit(_load_one_quant, q) for q in quants]
        for fut in tqdm(as_completed(futs), total=len(futs), desc="Loading quant.sf", unit="file"):
            tpm, cnt = fut.result()
            tpm_list.append(tpm)
            cnt_list.append(cnt)

    # Merge to wide DataFrames (outer join on index)
    tpm_df = pd.concat(tpm_list, axis=1)
    cnt_df = pd.concat(cnt_list, axis=1)

    # Save files (gzipped TSV) in the same root folder
    out_tpm = os.path.join(args.path_salmon, f"{args.project}_salmon_tpm.tsv.gz")
    out_cnt = os.path.join(args.path_salmon, f"{args.project}_salmon_count.tsv.gz")

    # Optional preview in console
    print("TPM head:")
    print(tpm_df.head())
    print("Count head:")
    print(cnt_df.head())

    # Persist results
    print("Saving TPM matrix...")
    with gzip.open(out_tpm, "wt", compresslevel=5) as f:
        tpm_df.to_csv(f, sep="\t")
    print("Saving Count matrix...")
    with gzip.open(out_cnt, "wt", compresslevel=5) as f:
        cnt_df.to_csv(f, sep="\t")

    # IMPORTANT: Output paths FIRST...
    print(f"Saved to (TPM): {out_tpm}", flush=True)
    print(f"Saved to (Count): {out_cnt}", flush=True)
    # ...THEN the banner
    _print_iobrpy_banner()

if __name__ == "__main__":
    main()
