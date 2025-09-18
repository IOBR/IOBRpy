#!/usr/bin/env python3
import argparse
import pickle
import pandas as pd
import numpy as np
from importlib.resources import files
from pathlib import Path
from typing import Optional
try:
    from tqdm.auto import tqdm
except Exception:
    def tqdm(iterable, **kwargs):
        return iterable
from iobrpy.utils.print_colorful_message import print_colorful_message
from iobrpy.utils.remove_version import strip_versions_in_dataframe, deduplicate_after_stripping

def _auto_sep_for_path(path: str, explicit_sep: Optional[str] = None):
    """
    Decide separator & engine for pandas.read_csv.
    - If explicit_sep is given (e.g., via --sep), use it with engine='c'.
    - If file ends with .tsv or .tsv.gz -> '\t', engine='c'
    - If file ends with .csv or .csv.gz -> ',', engine='c'
    - Otherwise (including unusual .gz like .expr.gz) -> sep=None, engine='python' (sniff)
    """
    if explicit_sep is not None:
        return explicit_sep, 'c'
    low = path.lower()
    if low.endswith(('.tsv', '.tsv.gz')):
        return '\t', 'c'
    if low.endswith(('.csv', '.csv.gz')):
        return ',', 'c'
    # Unusual extensions (including .gz): sniff
    return None, 'python'

def _load_external_annotation(path: Path, key: str = None):
    # 1) Pickle: support DataFrame or dict-of-DataFrame (select via --annotation-key)
    ext = path.suffix.lower()
    if ext in ('.pkl', '.pickle'):
        with path.open('rb') as f:
            obj = pickle.load(f)
        if isinstance(obj, dict):
            if key:
                if key not in obj:
                    raise KeyError(f"Key '{key}' not found in pkl. Available keys: {list(obj.keys())}")
                df = obj[key]
            else:
                keys = list(obj.keys())
                if len(keys) == 1:
                    df = obj[keys[0]]
                else:
                    raise KeyError(f"Pickle contains multiple keys: {keys}. Provide --annotation-key to choose one.")
    # 2) Excel: xls/xlsx as-is
    elif ext in ('.xls', '.xlsx'):
        df = pd.read_excel(path, index_col=None)
    else:
        # 3) Delimited text (csv/tsv/txt and ANY .gz suffix)
        #    Use auto sep sniffing for unusual .gz (e.g., .expr.gz)
        sep_infer, eng = _auto_sep_for_path(str(path), explicit_sep=None)
        try:
            read_kwargs = dict(
                filepath_or_buffer=path,
                sep=sep_infer,
                engine=eng,
                index_col=None,
                compression='infer',
                chunksize=200_000
            )
            if eng != 'python':
                read_kwargs['low_memory'] = False

            chunks = []
            for chunk in tqdm(pd.read_csv(**read_kwargs), desc="Loading annotation", unit="rows"):
                chunks.append(chunk)
            df = pd.concat(chunks, axis=0)
        except Exception as e:
            # Fallback: try pickle if user handed a pickled object with odd suffix
            try:
                with path.open('rb') as f:
                    obj = pickle.load(f)
                if isinstance(obj, dict):
                    if key:
                        if key not in obj:
                            raise KeyError(f"Key '{key}' not found in fallback pkl. Keys: {list(obj.keys())}")
                        df = obj[key]
                    else:
                        keys = list(obj.keys())
                        if len(keys) == 1:
                            df = obj[keys[0]]
                        else:
                            raise KeyError(f"Fallback pickle has multiple keys: {keys}. Provide --annotation-key.")
                elif isinstance(obj, pd.DataFrame):
                    df = obj
                else:
                    raise TypeError("Fallback pickle is not a DataFrame or dict of DataFrames.")
            except Exception as e2:
                raise TypeError(f"Unsupported/failed to parse annotation file '{path}'. "
                                f"Tried CSV/TSV (auto sep) and pickle fallback. Last errors: {e!r} / {e2!r}")

    if not isinstance(df, pd.DataFrame):
        raise TypeError("Loaded annotation is not a pandas DataFrame.")
    return df

def anno_eset(eset_df: pd.DataFrame,
              annotation,               # either a str key (for built-in) or a pd.DataFrame (external)
              symbol: str = "symbol",
              probe: str = "id",
              method: str = "mean") -> pd.DataFrame:
    """
    annotation: either a string key (lookup in built-in pkl) or a pandas.DataFrame provided by user.
    symbol/probe: column names in the annotation dataframe indicating gene symbol and probe id.
    """
    # === Prepare annotation (robust handling for different pkl structures) ===
    if isinstance(annotation, pd.DataFrame):
        annotation_df = annotation.copy()
    else:
        # load built-in resource
        resource_pkg = 'iobrpy.resources'
        resource_path = files(resource_pkg).joinpath('anno_eset.pkl')
        with resource_path.open('rb') as f:
            anno_dict = pickle.load(f)

        # if user passed a list-like annotation key by mistake, give clear error
        if isinstance(annotation, (list, tuple)):
            raise TypeError(f"Unexpected type for annotation argument: {type(annotation)}. "
                            f"Expected a string key for built-in resources or a DataFrame for external annotation.")

        if annotation not in anno_dict:
            raise KeyError(f"Annotation '{annotation}' not found in built-in resources. Available keys: {list(anno_dict.keys())}")

        raw = anno_dict[annotation]

        # normalize common shapes: DataFrame, list-of-DataFrame, dict-like, or other
        if isinstance(raw, pd.DataFrame):
            annotation_df = raw.copy()
        elif isinstance(raw, dict):
            # if dict, try to find a DataFrame value, or convert dict->DataFrame
            df_candidates = [v for v in raw.values() if isinstance(v, pd.DataFrame)]
            if len(df_candidates) == 1:
                annotation_df = df_candidates[0].copy()
            elif len(df_candidates) > 1:
                # ambiguous: pick first but warn
                print("Warning: built-in annotation entry is a dict with multiple DataFrames; taking the first one.")
                annotation_df = df_candidates[0].copy()
            else:
                # convert dict to DataFrame (fallback)
                annotation_df = pd.DataFrame(raw)
        elif isinstance(raw, (list, tuple)):
            # try to pick first DataFrame in list, or convert list->DataFrame
            df_candidates = [x for x in raw if isinstance(x, pd.DataFrame)]
            if len(df_candidates) >= 1:
                annotation_df = df_candidates[0].copy()
                if len(df_candidates) > 1:
                    print("Warning: annotation entry is a list with multiple DataFrames; using the first DataFrame found.")
            else:
                # try to coerce to DataFrame (e.g., list of tuples or rows)
                try:
                    annotation_df = pd.DataFrame(raw)
                    print("Note: converted list-like annotation entry to DataFrame.")
                except Exception as e:
                    raise TypeError(f"Unable to coerce annotation entry (type list) to DataFrame: {e}")
        else:
            # last resort: try to convert to DataFrame
            try:
                annotation_df = pd.DataFrame(raw)
                print(f"Note: coerced annotation object of type {type(raw)} into DataFrame.")
            except Exception:
                raise TypeError(f"Unsupported annotation object type: {type(raw)}. Expect DataFrame, dict, or list containing a DataFrame.")

    # quick sanity check
    if not isinstance(annotation_df, pd.DataFrame):
        raise TypeError(f"annotation_df is not a DataFrame after processing (type={type(annotation_df)}).")
    print("Annotation DataFrame shape:", annotation_df.shape)

    # rename provided symbol/probe to standardized names
    if symbol not in annotation_df.columns or probe not in annotation_df.columns:
        raise KeyError(f"Annotation does not contain specified columns. Expected symbol column '{symbol}' and probe column '{probe}'. Available columns: {list(annotation_df.columns)}")

    annotation_df = annotation_df.rename(columns={symbol: "symbol", probe: "probe_id"})
    annotation_df = annotation_df[["probe_id", "symbol"]]

    # filter out bad symbols
    annotation_df = annotation_df[annotation_df["symbol"] != "NA_NA"]
    annotation_df = annotation_df[annotation_df["symbol"].notna()]

    # Logging original count
    print(f"Row number of original eset: {eset_df.shape[0]}")

    # Annotate count
    probes_in = eset_df.index.isin(annotation_df["probe_id"]).sum()
    total_probes = eset_df.shape[0]
    print(f"Probes matched annotation: {probes_in} / {total_probes}")
    print(f"{100 * (probes_in / total_probes if total_probes else 0):.2f}% of probes were annotated")

    # Filter to annotated probes (preserve order of annotation_df)
    annotation_filtered = annotation_df[annotation_df["probe_id"].isin(eset_df.index)].copy()
    # reorder eset to match annotation_filtered probe_id order
    eset_filtered = eset_df.reindex(annotation_filtered["probe_id"]).copy()

    # Merge annotation (probe_id becomes a column)
    eset_reset = eset_filtered.reset_index().rename(columns={eset_filtered.index.name or 'index': 'probe_id'})
    merged = pd.merge(annotation_filtered, eset_reset, on="probe_id", how="inner")
    # drop probe_id column (we use symbol as index)
    merged.drop(columns=["probe_id"], inplace=True)

    # Handle duplicates: collapse by symbol using chosen method
    total_rows = merged.shape[0]
    unique_symbols = merged['symbol'].nunique()
    dups = total_rows - unique_symbols

    if dups > 0:
        data_cols = merged.columns.difference(['symbol'])
        if method == 'mean':
            merged['_score'] = merged[data_cols].mean(axis=1, skipna=True)
        elif method == 'sd':
            merged['_score'] = merged[data_cols].std(axis=1, skipna=True)
        elif method == 'sum':
            merged['_score'] = merged[data_cols].sum(axis=1, skipna=True)
        else:
            merged['_score'] = merged[data_cols].mean(axis=1, skipna=True)

        # keep highest scoring row per symbol
        merged.sort_values('_score', ascending=False, inplace=True)
        merged.drop(columns=['_score'], inplace=True)
        merged.drop_duplicates(subset=['symbol'], keep='first', inplace=True)

    result = merged.set_index('symbol')

    # Filter out rows all zero or all NA or NA in first column
    result = result.loc[~(result == 0).all(axis=1)]
    result = result.loc[~result.isna().all(axis=1)]
    if result.shape[1] > 0:
        first_col = result.columns[0]
        result = result.loc[result[first_col].notna()]

    print(f"Row number after filtering: {result.shape[0]}")
    return result


def main():
    parser = argparse.ArgumentParser(description="Annotate gene expression matrix and remove duplicated genes.")
    parser.add_argument('--input', '-i', dest='eset', required=True,
                        help='path to input matrix')
    parser.add_argument('--output', required=True, help='Path to save annotated matrix')
    parser.add_argument('--annotation', required=True,
                        help="Annotation key to use from built-in resources (one of: anno_hug133plus2, anno_rnaseq, anno_illumina, anno_grch38) OR ignored if --annotation-file is provided")
    parser.add_argument('--annotation-file', default=None,
                        help="Path to external annotation file (pkl/csv/tsv/xlsx). If provided, it overrides built-in annotation.")
    parser.add_argument('--annotation-key', default=None,
                        help="If external pkl contains multiple dataframes (a dict), select which key to use.")
    parser.add_argument('--symbol', default='symbol', help='Annotation symbol column name (in the annotation dataframe)')
    parser.add_argument('--probe', default='id', help='Annotation probe column name (in the annotation dataframe)')
    parser.add_argument('--method', default='mean', choices=['mean','sd','sum'], help='Dup handling method')
    parser.add_argument('--sep', default=None, help='Input sep; if omitted, infer from extension')
    parser.add_argument('--remove_version', action='store_true',
                        help='Strip Ensembl version suffix from row index before annotation')
    args = parser.parse_args()

    # Determine input separator
    sep = args.sep if args.sep is not None else ('\t' if args.eset.endswith(('.tsv','.tsv.gz')) else ',')

    # Read input
    sep_infer, eng = _auto_sep_for_path(args.eset, getattr(args, 'sep', None))
    read_kwargs = dict(
        filepath_or_buffer=args.eset,
        sep=sep_infer,
        engine=eng,
        index_col=0,
        compression='infer',
        chunksize=100_000
    )
    if eng != 'python':
        read_kwargs['low_memory'] = False

    chunks = []
    for chunk in tqdm(pd.read_csv(**read_kwargs), desc="Loading input", unit="rows"):
        chunks.append(chunk)
    df = pd.concat(chunks, axis=0)

    print(f"Loaded df shape: {df.shape}")

    if args.remove_version:
        df, n_stripped = strip_versions_in_dataframe(df, on_index=True, also_refseq=False)
        print(f"[anno_eset] stripped {n_stripped} versioned IDs from index")


    # Prepare annotation object (either DataFrame or built-in key)
    annotation_obj = None
    if args.annotation_file:
        p = Path(args.annotation_file)
        if not p.exists():
            raise FileNotFoundError(f"Annotation file not found: {p}")
        annotation_obj = _load_external_annotation(p, key=args.annotation_key)
    else:
        # use built-in key (string)
        annotation_obj = args.annotation

    # Annotate
    result = anno_eset(df, annotation_obj, symbol=args.symbol, probe=args.probe, method=args.method)

    # Save output: manual write to ensure first cell blank and headers aligned (CSV)
    cols = result.columns.tolist()
    with open(args.output, 'w', newline='') as out_f:
        # write header line: blank, then column names
        out_f.write(',' + ','.join(cols) + '\n')
        # write each row: index (symbol) and values
        for idx, row in tqdm(result.iterrows(), total=result.shape[0], desc="Writing rows", unit="row"):
            out_f.write(str(idx) + ',' + ','.join(map(str, row.values)) + '\n')
    print(f"Annotated matrix saved to {args.output}")

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

if __name__=='__main__':
    main()