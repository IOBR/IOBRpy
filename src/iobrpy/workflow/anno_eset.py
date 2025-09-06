#!/usr/bin/env python3
import argparse
import pickle
import pandas as pd
import numpy as np
from importlib.resources import files
from pathlib import Path
try:
    from tqdm.auto import tqdm
except Exception:
    def tqdm(iterable, **kwargs):
        return iterable

def _load_external_annotation(path: Path, key: str = None):
    ext = path.suffix.lower()
    if ext in ('.pkl', '.pickle'):
        with path.open('rb') as f:
            obj = pickle.load(f)
        # If it's a dict-like pkl, choose key
        if isinstance(obj, dict):
            if key:
                if key not in obj:
                    raise KeyError(f"Key '{key}' not found in pkl. Available keys: {list(obj.keys())}")
                df = obj[key]
            else:
                # if only one df inside, pick it; otherwise error to avoid ambiguity
                keys = list(obj.keys())
                if len(keys) == 1:
                    df = obj[keys[0]]
                else:
                    raise KeyError(f"Pickle contains multiple keys: {keys}. Provide --annotation-key to choose one.")
        elif isinstance(obj, pd.DataFrame):
            df = obj
        else:
            raise TypeError("Pickle does not contain a pandas DataFrame or dict of DataFrames.")
    elif ext in ('.csv', '.txt'):
        df = pd.read_csv(path, sep=',', index_col=None)
    elif ext in ('.tsv', '.gz'):
        # let pandas infer compression; try tab first
        try:
            df = pd.read_csv(path, sep='\t', index_col=None, compression='infer')
        except Exception:
            df = pd.read_csv(path, sep=',', index_col=None, compression='infer')
    elif ext in ('.xls', '.xlsx'):
        df = pd.read_excel(path, index_col=None)
    else:
        # try pickle as fallback
        with path.open('rb') as f:
            obj = pickle.load(f)
        if isinstance(obj, dict):
            if key:
                df = obj[key]
            else:
                keys = list(obj.keys())
                if len(keys) == 1:
                    df = obj[keys[0]]
                else:
                    raise KeyError(f"File ambiguous; provide --annotation-key. Keys: {keys}")
        elif isinstance(obj, pd.DataFrame):
            df = obj
        else:
            raise TypeError(f"Unsupported annotation file format: {ext}")
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
    args = parser.parse_args()

    # Determine input separator
    sep = args.sep if args.sep is not None else ('\t' if args.eset.endswith(('.tsv','.tsv.gz')) else ',')
    print(f"Input sep: '{sep}'")

    # Read input
    df = pd.read_csv(args.eset, sep=sep, index_col=0, compression='infer')
    print(f"Loaded df shape: {df.shape}")

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

if __name__=='__main__':
    main()
