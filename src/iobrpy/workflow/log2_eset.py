#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
log2_eset.py
Apply log2(x+1) to a gene expression matrix (rows=genes, cols=samples).
- Auto-detects the input delimiter (supports ',', '\\t', ';', '|'; also *.gz).
- Chooses output delimiter by output extension; otherwise mirrors input delimiter.

Usage:
  python log2_eset.py -i input.tsv -o output.tsv
  python log2_eset.py -i input.csv -o output.csv
  python log2_eset.py -i input.tsv.gz -o output.txt.gz
"""

import argparse
import sys
import os
import csv
import gzip
import numpy as np
import pandas as pd


CANDIDATE_SEPS = [",", "\t", ";", "|"]


def _open_text(path: str):
    """Open plain or gzipped text file in UTF-8 with replacement for bad bytes."""
    if path.lower().endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace", newline="")
    return open(path, mode="r", encoding="utf-8", errors="replace", newline="")


def _ext_lower(path: str) -> str:
    """Return lowercase extension without .gz suffix, e.g., .csv for *.csv.gz."""
    p = path.lower()
    if p.endswith(".gz"):
        p = p[:-3]
    _, ext = os.path.splitext(p)
    return ext  # includes leading dot, e.g., ".csv", ".tsv", ".txt", or ""


def _detect_input_sep(path: str) -> str:
    """
    Detect delimiter using csv.Sniffer first; fall back to simple frequency heuristic
    over the first N non-empty lines.
    """
    sample = ""
    try:
        with _open_text(path) as f:
            # Read up to ~64 KiB for sniffing
            sample = f.read(65536)
    except Exception as e:
        print(f"❌ Failed to open '{path}' for delimiter detection: {e}", file=sys.stderr)
        sys.exit(1)

    # Try csv.Sniffer
    if sample:
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="".join(CANDIDATE_SEPS))
            if getattr(dialect, "delimiter", None) in CANDIDATE_SEPS:
                return dialect.delimiter
        except Exception:
            pass  # fall through to heuristic

    # Heuristic: count occurrences of candidates across first ~10 non-empty lines
    counts = {sep: 0 for sep in CANDIDATE_SEPS}
    lines = [ln for ln in sample.splitlines() if ln.strip()]
    for ln in lines[:10]:
        for sep in CANDIDATE_SEPS:
            counts[sep] += ln.count(sep)

    # Choose the sep with the highest count; break ties by a preference order
    best_sep = max(CANDIDATE_SEPS, key=lambda s: (counts[s], -CANDIDATE_SEPS.index(s)))
    if counts[best_sep] > 0:
        return best_sep

    # As a last resort, use extension hint
    ext = _ext_lower(path)
    if ext == ".csv":
        return ","
    if ext == ".tsv":
        return "\t"

    # Default to tab (TSV) if totally ambiguous
    return "\t"


def _choose_output_sep(out_path: str, in_sep: str) -> str:
    """
    Decide output delimiter:
      - *.csv or *.csv.gz -> ','
      - *.tsv or *.tsv.gz -> '\t'
      - otherwise -> mirror input delimiter
    """
    ext = _ext_lower(out_path)
    if ext == ".csv":
        return ","
    if ext == ".tsv":
        return "\t"
    # Mirror input if extension is ambiguous (e.g., .txt, no extension, etc.)
    return in_sep or "\t"


def main():
    parser = argparse.ArgumentParser(
        description="Apply log2(x+1) to a gene expression matrix (rows=genes, cols=samples)."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input matrix (csv/tsv/txt; supports .gz). First column must be gene IDs (index)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output path. '.csv/.csv.gz' -> CSV, '.tsv/.tsv.gz' -> TSV, else mirror input delimiter."
    )
    args = parser.parse_args()

    in_path = args.input
    out_path = args.output

    # Auto-detect input delimiter
    in_sep = _detect_input_sep(in_path)
    sep_name = {"\t": "TAB", ",": "COMMA", ";": "SEMICOLON", "|": "PIPE"}.get(in_sep, repr(in_sep))
    print(f"ℹ️ Detected input delimiter: {sep_name}", file=sys.stderr)

    # Read matrix with detected delimiter; keep first column as gene index
    try:
        df = pd.read_csv(in_path, sep=in_sep, index_col=0)
    except Exception as e:
        print(f"❌ Failed to read input file '{in_path}' with sep={repr(in_sep)}: {e}", file=sys.stderr)
        sys.exit(1)

    if df.empty:
        print("❌ Input matrix is empty after reading. Check delimiter and index column.", file=sys.stderr)
        sys.exit(1)

    # Convert all columns to numeric (coerce errors to NaN) while preserving shape
    df_before_na = df.isna().sum().sum()
    df = df.apply(pd.to_numeric, errors="coerce")
    coerced_cells = int(df.isna().sum().sum() - df_before_na)
    if coerced_cells > 0:
        print(
            f"⚠️ Detected non-numeric entries; coerced to NaN in {coerced_cells} cells.",
            file=sys.stderr
        )

    # Domain check for log2(x+1)
    min_val = df.min().min(skipna=True)
    if pd.notna(min_val) and min_val < -1:
        print(f"❌ Found values < -1 (min={min_val}). log2(x+1) is undefined for x < -1.", file=sys.stderr)
        sys.exit(1)
    if pd.notna(min_val) and min_val < 0:
        print(f"⚠️ Found negative values (min={min_val}). Proceeding with log2(x+1); verify this is expected.", file=sys.stderr)

    # Transform
    try:
        out_df = np.log2(df + 1.0)
    except Exception as e:
        print(f"❌ Failed during log2(x+1) transformation: {e}", file=sys.stderr)
        sys.exit(1)

    # Choose output delimiter
    out_sep = _choose_output_sep(out_path, in_sep)
    out_sep_name = {"\t": "TAB", ",": "COMMA", ";": "SEMICOLON", "|": "PIPE"}.get(out_sep, repr(out_sep))
    print(f"ℹ️ Output delimiter: {out_sep_name}", file=sys.stderr)

    # Write output
    try:
        out_df.to_csv(out_path, sep=out_sep, index=True)
    except Exception as e:
        print(f"❌ Failed to write output file '{out_path}': {e}", file=sys.stderr)
        sys.exit(1)

    print(f"✅ Done. Wrote log2(x+1) matrix to: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()