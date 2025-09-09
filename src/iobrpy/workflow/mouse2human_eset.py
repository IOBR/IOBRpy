#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Converting mouse gene symbols to human gene symbols for an expression set.

- Uses local resources (iobrpy.resources/mus_human.pkl) for mapping.
- When --is_matrix is provided: rows are mouse gene symbols; columns are samples.
- When --is_matrix is omitted: input must contain a gene-symbol column specified by --column_of_symbol.
  The table will be de-duplicated first and then treated as a matrix with symbols as index.
- Uses anno_eset() from iobrpy.workflow.anno_eset to perform final mapping
  and duplicate handling consistent with the R workflow (method='mean').

CLI:
    python -m iobrpy.workflow.mouse2human_eset \
        -i input.csv -o output.csv \
        --is_matrix \
        --column_of_symbol SYMBOL

Input/Output separator rules:
- .tsv, .tsv.gz, .tab, .tab.gz, .txt, .txt.gz -> '\t'
- otherwise -> ','
- You can override by --sep (input) and --out_sep (output).
"""

import argparse
import gzip
import pickle
from pathlib import Path
from typing import Optional

import pandas as pd

# Optional progress bar
try:
    from tqdm.auto import tqdm
except Exception:  # pragma: no cover
    def tqdm(iterable=None, **kwargs):
        return iterable if iterable is not None else range(0)

# Built-in resources loader
from importlib.resources import files as ir_files

# Use your existing anno_eset from workflow
from iobrpy.workflow.anno_eset import anno_eset

# -------- Colorful printing utility (drop-in, no hard deps) --------
import os
import sys

# ANSI color codes covering common and bright variants
_ANSI_COLOR_CODES = {
    "black": "30", "red": "31", "green": "32", "yellow": "33",
    "blue": "34", "magenta": "35", "cyan": "36", "white": "37",
    "bright_black": "90", "bright_red": "91", "bright_green": "92",
    "bright_yellow": "93", "bright_blue": "94", "bright_magenta": "95",
    "bright_cyan": "96", "bright_white": "97",
}

def _enable_windows_ansi(stream) -> bool:
    """
    Try to enable ANSI escape sequences on Windows consoles.
    Returns True if ANSI is (now or already) supported; False otherwise.
    Strategy:
      1) Try enabling VT100 processing on the console handle via Win32 API.
      2) Fallback: if 'colorama' is installed, initialize it.
    """
    if os.name != "nt":
        return True

    # (1) Attempt to turn on VT processing (Windows 10+)
    try:
        import ctypes
        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        # Prefer the stream's handle if possible; otherwise use STD_OUTPUT_HANDLE
        handle = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
        mode = ctypes.c_uint32()
        if kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
            kernel32.SetConsoleMode(handle, mode.value | ENABLE_VIRTUAL_TERMINAL_PROCESSING)
            return True
    except Exception:
        pass

    # (2) Optional fallback to colorama if available
    try:
        import colorama  # type: ignore
        colorama.just_fix_windows_console()
        return True
    except Exception:
        return False

def _supports_color(stream) -> bool:
    """
    Return True if the given stream is a TTY that likely supports colors.
    - On non-Windows: any TTY is assumed to support ANSI.
    - On Windows: require successful ANSI enablement (see above).
    - If the stream is not a TTY (e.g., piped to file), return False.
    """
    if not hasattr(stream, "isatty") or not stream.isatty():
        return False
    if os.name == "nt":
        return _enable_windows_ansi(stream)
    return True

def print_colorful_message(text: str,
                           color: str = "cyan",
                           bold: bool = False,
                           stream=None) -> None:
    """
    Print colored text to the terminal, gracefully degrading to plain text
    when color is not supported (non-TTY, old terminals, or redirected output).

    Args:
        text:  The message to print.
        color: Color name (e.g., 'blue', 'cyan', 'red', 'bright_blue', ...).
        bold:  If True, apply bold style.
        stream: Target stream (defaults to sys.stdout).

    Behavior:
        - Uses ANSI escape sequences when supported.
        - Unknown color names fall back to 'cyan'.
        - Automatically appends a newline if 'text' does not end with one.
    """
    stream = stream or sys.stdout
    use_color = _supports_color(stream)
    newline = "" if text.endswith("\n") else "\n"

    if use_color:
        code = _ANSI_COLOR_CODES.get(color.lower(), _ANSI_COLOR_CODES["cyan"])
        prefix = f"\033[{'1;' if bold else ''}{code}m"
        suffix = "\033[0m"
        stream.write(f"{prefix}{text}{suffix}{newline}")
    else:
        stream.write(text + newline)
# -------- End of colorful printing utility --------

# ---------------------- Utilities ----------------------
def _infer_sep(path: str, explicit: Optional[str] = None) -> str:
    """Infer separator from extension unless explicitly provided.

    Rules:
    - .tsv, .tsv.gz, .tab, .tab.gz, .txt, .txt.gz  -> '\t'
    - otherwise                                     -> ','
    """
    if explicit is not None:
        return explicit
    p = str(path).lower()
    tab_exts = (".tsv", ".tsv.gz", ".tab", ".tab.gz", ".txt", ".txt.gz")
    if p.endswith(tab_exts):
        return "\t"
    return ","


def _read_matrix(path: str, sep: Optional[str] = None) -> pd.DataFrame:
    """Read expression matrix with row index as the first column."""
    sep = _infer_sep(path, sep)
    df = pd.read_csv(path, sep=sep, index_col=0, compression="infer")
    return df


def _read_table_with_symbol(path: str, sep: Optional[str] = None) -> pd.DataFrame:
    """Read non-matrix table (gene symbol is a normal column)."""
    sep = _infer_sep(path, sep)
    df = pd.read_csv(path, sep=sep, compression="infer")
    return df


def _remove_duplicate_genes(
    df: pd.DataFrame,
    column_of_symbol: str,
    method: str = "mean"
) -> pd.DataFrame:
    """
    Remove duplicate gene symbols by selecting one row per symbol.
    The selection follows the R logic:
      - compute a row 'score' using the aggregation method across numeric columns
      - sort by score (descending), keep the first occurrence per symbol

    Parameters
    ----------
    df : DataFrame
        Input table containing a symbol column and expression columns.
    column_of_symbol : str
        Column name that holds gene symbols.
    method : {'mean','sd','sum'}
        Aggregation for the scoring.

    Returns
    -------
    DataFrame
        DataFrame indexed by gene symbol with expression values only.
    """
    df = df.copy()

    if column_of_symbol not in df.columns:
        raise KeyError(
            f"Column '{column_of_symbol}' not found in the input table. "
            f"Available columns: {list(df.columns)}"
        )

    sym_col = column_of_symbol
    value_cols = [c for c in df.columns if c != sym_col]

    # Keep only numeric columns for scoring (non-numeric are ignored)
    numeric_cols = [c for c in value_cols if pd.api.types.is_numeric_dtype(df[c])]
    if not numeric_cols:
        raise ValueError("No numeric columns found for duplicate resolution.")

    # Compute score per row
    if method == "mean":
        score = df[numeric_cols].mean(axis=1, skipna=True)
    elif method == "sd":
        score = df[numeric_cols].std(axis=1, skipna=True)
    elif method == "sum":
        score = df[numeric_cols].sum(axis=1, skipna=True)
    else:
        score = df[numeric_cols].mean(axis=1, skipna=True)

    df["_score"] = score
    # Sort by score descending so the first duplicate kept is the highest score
    df.sort_values("_score", ascending=False, inplace=True)
    df = df.drop(columns=["_score"])

    # Drop duplicates keeping the first (highest score)
    df = df.drop_duplicates(subset=[sym_col], keep="first")

    # Set index to symbols and keep only expression columns
    df = df.set_index(sym_col)
    df = df[value_cols]

    return df


def _load_mus_human_df() -> pd.DataFrame:
    """
    Load the mouse->human mapping from built-in resources (mus_human.pkl).

    The pickle is expected to contain either:
        - a pandas DataFrame, or
        - a dict with a single DataFrame (or a key that contains a DataFrame).
    The DataFrame must include columns:
        'gene_symbol_mus' and 'gene_symbol_human'
    """
    resource_pkg = "iobrpy.resources"
    resource_path = ir_files(resource_pkg).joinpath("mus_human.pkl")

    with resource_path.open("rb") as f:
        obj = pickle.load(f)

    if isinstance(obj, pd.DataFrame):
        df = obj.copy()
    elif isinstance(obj, dict):
        # choose the only DataFrame or raise if ambiguous
        dfs = [v for v in obj.values() if isinstance(v, pd.DataFrame)]
        if len(dfs) == 1:
            df = dfs[0].copy()
        elif len(dfs) > 1:
            print("Warning: multiple DataFrames found in mus_human.pkl; using the first one.")
            df = dfs[0].copy()
        else:
            df = pd.DataFrame(obj)
    else:
        df = pd.DataFrame(obj)

    needed = {"gene_symbol_mus", "gene_symbol_human"}
    if not needed.issubset(set(df.columns)):
        raise KeyError(
            f"mus_human.pkl does not contain required columns {needed}. "
            f"Available columns: {list(df.columns)}"
        )
    return df


def _open_text_writer(path: Path):
    """Open a text writer; use gzip if path ends with .gz."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return path.open("w", newline="")


def _write_with_progress(df: pd.DataFrame, out_path: Path, sep: str, show_progress: bool = False):
    """
    Write a DataFrame to CSV/TSV/TXT with an optional progress bar.

    - Writes a blank top-left header cell (common in expression matrices).
    - Replaces NaN with empty strings during writing.
    - Uses gzip when the filename ends with .gz.
    """
    cols = list(map(str, df.columns.tolist()))
    total = df.shape[0]

    with _open_text_writer(out_path) as f:
        # header: blank cell + column names
        f.write(sep.join([""] + cols) + "\n")

        iterator = df.itertuples(index=True, name=None)  # (index, v1, v2, ...)
        if show_progress:
            iterator = tqdm(iterator, total=total, unit="row", desc="Saving")

        for row in iterator:
            idx = "" if row[0] is None else str(row[0])
            # Replace NaN with empty string and convert to str
            values = [("" if pd.isna(v) else v) for v in row[1:]]
            f.write(idx + sep + sep.join(map(str, values)) + "\n")


# ---------------------- Core ----------------------
def mouse2human_eset(
    eset: pd.DataFrame,
    is_matrix: bool = False,
    column_of_symbol: Optional[str] = None,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Convert mouse gene symbols to human gene symbols in an expression dataset.

    Parameters
    ----------
    eset : DataFrame
        Expression matrix (rows=genes, cols=samples) or a table with a symbol column.
    is_matrix : bool
        If True, 'eset' row index is mouse gene symbols. If False, 'column_of_symbol'
        must be set and duplicates will be resolved before conversion.
    column_of_symbol : Optional[str]
        Column name containing gene symbols when is_matrix=False.
    verbose : bool
        Print additional information.

    Returns
    -------
    DataFrame
        Expression matrix indexed by human gene symbols.
    """
    if not is_matrix:
        if not column_of_symbol:
            raise ValueError("When is_matrix=False, --column_of_symbol must be specified.")
        # Remove duplicates first (default scoring by mean to match R)
        eset = _remove_duplicate_genes(eset, column_of_symbol=column_of_symbol, method="mean")

    # Load mapping from local resources
    probe_data = _load_mus_human_df()
    if verbose:
        print("Loaded mus_human mapping:", probe_data.shape)

    # Call anno_eset (expects 'probe' to match row index of eset)
    # Here: probe -> mouse symbol, symbol -> human symbol
    result = anno_eset(
        eset_df=eset,
        annotation=probe_data,
        symbol="gene_symbol_human",
        probe="gene_symbol_mus",
        method="mean"
    )
    return result


# ---------------------- CLI ----------------------
def main():
    parser = argparse.ArgumentParser(
        description="Convert mouse gene symbols to human symbols for an expression set (local mapping only)."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Path to input file (CSV/TSV/TXT, optionally .gz).")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to save converted matrix (CSV/TSV/TXT, optionally .gz).")

    # Parity with R parameters (renamed with underscores)
    parser.add_argument("--is_matrix", action="store_true",
                        help="Treat input as a matrix (rows=genes, cols=samples). If omitted, expects a symbol column.")
    parser.add_argument("--column_of_symbol", default=None,
                        help="Column name containing gene symbols when not using --is_matrix.")
    parser.add_argument("--verbose", action="store_true",
                        help="Verbose output.")

    # Optional: allow custom input/output separator overrides
    parser.add_argument("--sep", default=None,
                        help="Input separator (',' or '\\t'). If omitted, infer by input extension.")
    parser.add_argument("--out_sep", default=None,
                        help="Output separator (',' or '\\t'). If omitted, infer by output extension.")

    # Progress bar toggle for saving
    parser.add_argument("--progress", action="store_true",default=True,
                        help="Show a progress bar during saving.")

    args = parser.parse_args()

    # Read input
    if args.is_matrix:
        df = _read_matrix(args.input, sep=args.sep)
    else:
        df = _read_table_with_symbol(args.input, sep=args.sep)

    if args.verbose:
        print(f"Loaded input shape: {df.shape}")
        print(f"is_matrix={args.is_matrix}, column_of_symbol={args.column_of_symbol}")

    # Run conversion
    out_df = mouse2human_eset(
        eset=df,
        is_matrix=args.is_matrix,
        column_of_symbol=args.column_of_symbol,
        verbose=args.verbose
    )

    # Save with inferred/overridden separator and compression (manual writer + optional progress)
    out_sep = _infer_sep(args.output, args.out_sep)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    _write_with_progress(out_df, out_path, sep=out_sep, show_progress=args.progress)

    print(f"[iobrpy] Converted matrix saved to: {out_path.resolve()}")
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
