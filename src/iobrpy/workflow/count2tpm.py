import argparse
import os
import pandas as pd
import numpy as np
from importlib.resources import files
from iobrpy.utils.print_colorful_message import print_colorful_message
from iobrpy.utils.remove_version import strip_versions_in_dataframe, deduplicate_after_stripping

try:
    from tqdm.auto import tqdm
except Exception:
    # fallback: if tqdm is not installed, tqdm(...) just returns the iterable unchanged
    def tqdm(iterable, **kwargs):
        return iterable

def feature_manipulation(data: pd.DataFrame,
                         feature: list = None,
                         is_matrix: bool = False,
                         print_result: bool = False) -> list:
    """
    Filter features (rows or columns) based on NA, non-numeric, infinite, and zero variance.
    Returns the list of valid feature names.
    """
    df = data.copy()
    if is_matrix:
        # Treat rows as features
        df = df.T
        feature = list(df.columns)
    else:
        if feature is None:
            feature = list(df.columns)
    # Remove features with any NA
    na_mask = df[feature].isna().any(axis=0)
    feature = [f for f in feature if not na_mask[f]]
    # Remove non-numeric
    feature = [f for f in feature if pd.api.types.is_numeric_dtype(df[f])]
    # Remove infinite values
    inf_mask = df[feature].apply(lambda col: np.isinf(col).any(), axis=0)
    feature = [f for f, bad in inf_mask.items() if not bad]
    # Remove zero-variance
    zero_var = df[feature].std(axis=0) == 0
    feature = [f for f in feature if not zero_var.get(f, False)]
    return feature

def remove_duplicate_genes(df: pd.DataFrame,
                           column_of_symbol: str = "symbol",
                           method: str = "mean",
                           show_progress: bool = True) -> pd.DataFrame:
    """
    Aggregate rows with duplicate gene symbols.
    method: "mean", "sd" (std), or "sum".
    Returns DataFrame indexed by unique gene symbol.
    """
    if column_of_symbol not in df.columns:
        raise ValueError(f"Column '{column_of_symbol}' not found in DataFrame.")
    symbols = df[column_of_symbol]
    data = df.drop(columns=[column_of_symbol])

    # No duplicates: set symbol as index
    if symbols.duplicated().sum() == 0:
        return df.set_index(column_of_symbol)

    # Prepare temp
    temp = data.copy()
    temp[column_of_symbol] = symbols

    # Choose agg function
    if method == "mean": agg_func = 'mean'
    elif method == "sd": agg_func = 'std'
    elif method == "sum": agg_func = 'sum'
    else: raise ValueError("method must be 'mean', 'sd', or 'sum'")

    # Order by agg descending (optional)
    numeric_cols = temp.select_dtypes(include=[np.number]).columns.tolist()
    if numeric_cols:
        order_series = getattr(temp[numeric_cols], agg_func)(axis=1)
        temp['_order'] = order_series
        temp = temp.sort_values('_order', ascending=False).drop(columns=['_order'])

    # Group and aggregate with progress bar
    groups = temp.groupby(column_of_symbol, as_index=True)
    iterator = tqdm(groups, total=getattr(groups, 'ngroups', None),
                    desc="Deduplicating", unit="gene") if show_progress else groups

    rows, idx = [], []
    for key, grp in iterator:
        rows.append(getattr(grp[numeric_cols], agg_func)() if numeric_cols else pd.Series(dtype=float))
        idx.append(key)

    aggregated = pd.DataFrame(rows, index=idx)
    return aggregated


def count2tpm(count_mat: pd.DataFrame,
              anno_grch38: pd.DataFrame,
              anno_gc_vm32: pd.DataFrame,
              idType: str = "Ensembl",
              org: str = "hsa",
              source: str = "local",
              effLength_df: pd.DataFrame = None,
              id_col: str = "id",
              gene_symbol_col: str = "symbol",
              length_col: str = "eff_length",
              check_data: bool = False,
              remove_version: bool = False) -> pd.DataFrame:
    # Debug: entry
    print(f"Enter count2tpm with idType={idType}, org={org}, source={source}")

    if org not in ["hsa", "mmus"]:
        raise ValueError("`org` must be 'hsa' or 'mmus'.")

    countMat = count_mat.copy().astype(float)
    print(f"initial countMat.shape: {countMat.shape}")
    # Optionally strip Ensembl version suffixes in index (ENSG / ENSMUSG)
    if remove_version:
        idx = countMat.index.astype(str)
        looks_ensembl = (idx.str.startswith('ENSG') | idx.str.startswith('ENSMUSG')).any()
        if looks_ensembl:
            countMat, n_stripped = strip_versions_in_dataframe(countMat, on_index=True, also_refseq=False)
            print(f"[count2tpm] stripped {n_stripped} Ensembl-versioned IDs (ENSG/ENSMUSG)")
    if countMat.isna().sum().sum() > 0 or check_data:
        n_na = countMat.isna().sum().sum()
        print(f"There are {n_na} missing or bad values; filtering genes via feature_manipulation.")
        features = feature_manipulation(countMat, is_matrix=True)
        countMat = countMat.loc[features]
        print(f"countMat.shape after feature_manipulation: {countMat.shape}")

    def _as_df(obj, name):
        if isinstance(obj, pd.DataFrame):
            return obj
        try:
            return pd.DataFrame(obj)
        except Exception as e:
            raise TypeError(f"{name} must be a DataFrame; got {type(obj).__name__}") from e

    if not (isinstance(anno_grch38, pd.DataFrame) and isinstance(anno_gc_vm32, pd.DataFrame)):
        resource_pkg = 'iobrpy.resources'
        resource_path = files(resource_pkg).joinpath('count2tpm_data.pkl')
        with resource_path.open('rb') as f:
            anno_data = pd.read_pickle(f, compression=None)
        anno_grch38 = _as_df(anno_data['anno_grch38'], 'anno_grch38')
        anno_gc_vm32 = _as_df(anno_data['anno_gc_vm32'], 'anno_gc_vm32')
    else:
        anno_grch38 = _as_df(anno_grch38, 'anno_grch38')
        anno_gc_vm32 = _as_df(anno_gc_vm32, 'anno_gc_vm32')

    len_series = None
    symbol_map = None

    if effLength_df is not None:
        print("Using manual effLength_df provided.")
        eff = effLength_df.rename(columns={id_col: 'id', length_col: 'eff_length', gene_symbol_col: 'gene_symbol'})
        eff = eff.drop_duplicates(subset='id').set_index('id')
        common = countMat.index.intersection(eff.index)
        countMat = countMat.loc[common]
        eff = eff.loc[common]
        len_series = eff['eff_length']
        symbol_map = eff['gene_symbol']
    else:
        print("Loading lengths from annotation tables.")
        if source != 'local':
            raise NotImplementedError("source='biomart' not implemented.")
        if org == 'hsa' and idType.lower() == 'ensembl':
            countMat.index = countMat.index.str[:15]
            df_ref = anno_grch38[['id', 'eff_length', 'symbol']].copy().set_index('id')
        elif org == 'hsa' and idType.lower() == 'entrez':
            print("Fuzzy calculation for Entrez IDs.")
            df_ref = anno_grch38[['entrez', 'eff_length', 'symbol']].copy()
            df_ref.columns = ['id', 'eff_length', 'symbol']
            df_ref = df_ref.drop_duplicates(subset='id').set_index('id')
        elif org == 'hsa' and idType.lower() == 'symbol':
            print("Fuzzy calculation for gene symbols.")
            df_ref = anno_grch38[['symbol', 'eff_length', 'gc']].copy()
            df_ref.columns = ['id', 'eff_length', 'gc']
            df_ref = df_ref.drop_duplicates(subset='id').set_index('id')
        elif org == 'mmus' and idType.lower() == 'ensembl':
            print("Using mmus annotation for Ensembl.")
            df_ref = anno_gc_vm32[['id', 'eff_length', 'symbol']].copy().set_index('id')
        elif org == 'mmus' and idType.lower() == 'mgi':
            print("Fuzzy calculation for MGI IDs.")
            df_ref = anno_gc_vm32[['mgi_id', 'eff_length', 'symbol']].copy()
            df_ref.columns = ['id', 'eff_length', 'symbol']
            df_ref = df_ref.drop_duplicates(subset='id').set_index('id')
        elif org == 'mmus' and idType.lower() == 'symbol':
            print("Fuzzy calculation for gene symbols mmus.")
            df_ref = anno_gc_vm32[['symbol', 'eff_length', 'gc']].copy()
            df_ref.columns = ['id', 'eff_length', 'gc']
            df_ref = df_ref.drop_duplicates(subset='id').set_index('id')
        else:
            raise ValueError(f"Unsupported idType for {org}: {idType}")
        common = countMat.index.intersection(df_ref.index)
        df_ref = df_ref.loc[common]
        countMat = countMat.loc[common]
        len_series = df_ref['eff_length']
        symbol_map = df_ref['symbol']

    # Drop duplicates in len_series/symbol_map
    if len_series.index.duplicated().any():
        dup = len_series.index.duplicated().sum()
        print(f"Dropping {dup} duplicates in len_series/symbol_map.")
        mask = ~len_series.index.duplicated(keep='first')
        len_series = len_series[mask]
        symbol_map = symbol_map[mask]

    # Align to countMat
    len_series = len_series.reindex(countMat.index)
    symbol_map = symbol_map.reindex(countMat.index)

    # Filter NA lengths
    valid = ~len_series.isna()
    countMat = countMat.loc[valid]
    len_series = len_series[valid]
    symbol_map = symbol_map[valid]

    # Aggregate duplicates in countMat
    if countMat.index.duplicated().any():
        countMat = countMat.groupby(countMat.index).sum()
        len_series = len_series.loc[countMat.index]
        symbol_map = symbol_map.loc[countMat.index]

    # Compute TPM with progress bars (per-sample operations)
    print("[INFO] Computing TPM (showing progress per sample)...")
    # prepare containers
    rpk = pd.DataFrame(index=countMat.index, columns=countMat.columns, dtype=float)
    divisor = (len_series / 1000)  # kb per gene

    # 1) compute RPK (reads per kilobase) per sample
    for col in tqdm(countMat.columns, desc="RPK per sample", unit="sample"):
        # elementwise division (aligned by index)
        rpk[col] = countMat[col] / divisor

    # 2) sum RPK per sample
    denom = pd.Series(index=rpk.columns, dtype=float)
    for col in tqdm(rpk.columns, desc="Sum RPK", unit="sample"):
        denom[col] = rpk[col].sum()

    # 3) normalize to TPM per sample
    tpm = pd.DataFrame(index=rpk.index, columns=rpk.columns, dtype=float)
    for col in tqdm(rpk.columns, desc="TPM normalize", unit="sample"):
        if denom[col] > 0:
            tpm[col] = rpk[col] / denom[col] * 1e6
        else:
            # avoid division by zero
            tpm[col] = 0.0


    # Insert symbol as first column
    tpm.insert(0, 'symbol', symbol_map.loc[tpm.index].values)

    # Deduplicate by symbol with progress
    print("[INFO] Removing duplicate genes...")
    tpm = remove_duplicate_genes(tpm, 'symbol')  # 不要再包 tqdm 了

    return tpm

def main():
    parser = argparse.ArgumentParser(description="Convert count matrix to TPM.")
    parser.add_argument('--input', type=str, required=True,
                        help="Path to input count matrix (CSV/TSV/GZ) with gene IDs as index.")
    parser.add_argument('--output', type=str, required=True,
                        help="Path (including filename) to save TPM output, e.g. /path/to/TPM.csv")
    parser.add_argument('--idtype', type=str, default="ensembl",
                        choices=["ensembl","entrez","symbol","mgi"],
                        help="Gene ID type.")
    parser.add_argument('--org', type=str, default="hsa", choices=["hsa","mmus"],
                        help="Organism: hsa or mmus.")
    parser.add_argument('--source', type=str, default="local", choices=["local","biomart"],
                        help="Source of gene lengths.")
    parser.add_argument('--effLength_csv', type=str,
                        help="Optional CSV with id, eff_length, and gene_symbol columns.")
    parser.add_argument('--id', dest='id_col', type=str, default="id",
                        help="Column name for gene ID in effLength CSV.")
    parser.add_argument('--length', dest='length_col', type=str, default="eff_length",
                        help="Column name for gene length in effLength CSV.")
    parser.add_argument('--gene_symbol', dest='gene_symbol_col', type=str, default="symbol",
                        help="Column name for gene symbol in effLength CSV.")
    parser.add_argument('--check_data', action='store_true',
                        help="Check and remove missing values in count matrix.")
    parser.add_argument('--remove_version', action='store_true',
                        help='Strip gene version suffix from row index before processing.')
    args = parser.parse_args()

    # Load count matrix
    if args.input.endswith('.gz'):
        count_mat = pd.read_csv(args.input, sep='\t', index_col=0, compression='gzip')
    else:
        sep = '\t' if args.input.endswith('.tsv') else ','
        count_mat = pd.read_csv(args.input, sep=sep, index_col=0)

    # Load manual effLength if provided
    eff_df = pd.read_csv(args.effLength_csv) if args.effLength_csv else None

    # Run conversion
    tpm_df = count2tpm(count_mat,None,None,
                        idType=args.idtype, org=args.org,
                        source=args.source, effLength_df=eff_df,
                        id_col=args.id_col, gene_symbol_col=args.gene_symbol_col,
                        length_col=args.length_col, check_data=args.check_data,remove_version=args.remove_version)

    cols = list(tpm_df.columns)
    cols[0] = ""
    tpm_df.columns = cols

    print(f"Output matrix shape: {tpm_df.shape[0]} rows × {tpm_df.shape[1]} columns")
    print(tpm_df.iloc[:5, :5])
    
    out_file = args.output
    tpm_df.to_csv(out_file, index=True)
    print(f"Saved TPM matrix to {out_file}")

    # ---- Friendly banner (English comments) ----
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
    # ---- End banner ----

if __name__ == '__main__':
    main()