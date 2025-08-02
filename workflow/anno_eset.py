import argparse
import pickle
import pandas as pd
import numpy as np
import importlib_resources as pkg_res

def anno_eset(eset_df: pd.DataFrame,
              annotation: str,
              symbol: str = "symbol",
              probe: str = "probe_id",
              method: str = "mean") -> pd.DataFrame:
    # Prepare annotation
    resource_pkg = 'iobrpy.resources'
    resource_path = pkg_res.files(resource_pkg).joinpath('anno_eset.pkl')
    with resource_path.open('rb') as f:
        anno_dict = pickle.load(f)
    if annotation not in anno_dict:
        raise KeyError(f"Annotation '{annotation}' not found in resources")
    annotation_df = anno_dict[annotation]
    print("[DEBUG] Annotation DataFrame shape:", annotation_df.shape)
    annotation = annotation_df.copy()
    annotation = annotation.rename(columns={symbol: "symbol", probe: "probe_id"})
    annotation = annotation[["probe_id", "symbol"]]
    annotation = annotation[annotation["symbol"] != "NA_NA"]
    annotation = annotation[annotation["symbol"].notna()]

    # Logging original count
    print(f"Row number of original eset: {eset_df.shape[0]}")

    # Annotate count
    probes_in = eset_df.index.isin(annotation["probe_id"]).sum()
    total_probes = eset_df.shape[0]
    print(f"[DEBUG] Probes matched annotation: {probes_in} / {total_probes}")
    print(f"{100 * (probes_in / total_probes if total_probes else 0):.2f}% of probes were annotated")

    # Filter to annotated probes
    annotation_filtered = annotation[annotation["probe_id"].isin(eset_df.index)]
    eset_filtered = eset_df.loc[annotation_filtered["probe_id"]].copy()

    # Merge annotation
    eset_reset = eset_filtered.reset_index().rename(columns={eset_filtered.index.name or 'index': 'probe_id'})
    merged = pd.merge(annotation_filtered, eset_reset, on="probe_id", how="inner")
    merged.drop(columns=["probe_id"], inplace=True)

    # Handle duplicates
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
    parser.add_argument('--annotation', required=True, choices=['anno_hug133plus2','anno_rnaseq','anno_illumina'],help="Annotation key to use (one of: anno_hug133plus2, anno_rnaseq, anno_illumina)")
    parser.add_argument('--symbol', default='symbol', help='Annotation symbol column')
    parser.add_argument('--probe', default='probe_id', help='Annotation probe column')
    parser.add_argument('--method', default='mean', choices=['mean','sd','sum'], help='Dup handling method')
    parser.add_argument('--sep', default=None, help='Input sep; if omitted, infer from extension')
    args = parser.parse_args()

    # Determine input separator
    sep = args.sep if args.sep is not None else ('\t' if args.eset.endswith(('.tsv','.tsv.gz')) else ',')
    print(f"[DEBUG] Input sep: '{sep}'")

    # Read input
    df = pd.read_csv(args.eset, sep=sep, index_col=0, compression='infer')
    print(f"[DEBUG] Loaded df shape: {df.shape}")

    # Annotate
    result = anno_eset(df, args.annotation, symbol=args.symbol, probe=args.probe, method=args.method)

    # Save output: manual write to ensure first cell blank and headers aligned
    cols = result.columns.tolist()
    with open(args.output, 'w', newline='') as out_f:
        # write header line: blank, then column names
        out_f.write(',' + ','.join(cols) + '\n')
        # write each row: index (symbol) and values
        for idx, row in result.iterrows():
            out_f.write(idx + ',' + ','.join(map(str, row.values)) + '\n')
    print(f"Annotated matrix saved to {args.output}")

if __name__=='__main__':
    main()
