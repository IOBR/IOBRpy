import argparse
import pickle
import pandas as pd
import numpy as np
import os
from importlib.resources import files, as_file
from tqdm import tqdm
from .count2tpm import count2tpm


def detect_sep(filename: str) -> str:
    ext = os.path.splitext(filename)[1].lower()
    if ext in ['.tsv', '.txt']:
        return '\t'
    return ','


def feature_manipulation(data: pd.DataFrame,
                         feature: list = None,
                         is_matrix: bool = False,
                         print_result: bool = False) -> list:
    df = data.copy()
    if is_matrix:
        df = df.T
        features = list(df.columns)
    else:
        features = feature if feature is not None else list(df.columns)
        features = [f for f in features if f in df.columns]
    # Remove NA
    complete = df[features].dropna(axis=1, how='any').columns.tolist()
    if print_result:
        print(f"After NA filter: {len(complete)}/{len(features)}")
    # Non-numeric
    numeric = [f for f in complete if pd.api.types.is_numeric_dtype(df[f])]
    if print_result:
        print(f"After non-numeric filter: {len(numeric)}/{len(complete)}")
    # Infinite
    finite = [f for f in numeric if np.isfinite(df[f]).all()]
    if print_result:
        print(f"After infinite filter: {len(finite)}/{len(numeric)}")
    # Zero variance
    if df.shape[0] > 1:
        valid = [f for f in finite if df[f].std() != 0]
    else:
        valid = finite
    if print_result:
        print(f"After zero-variance filter: {len(valid)}/{len(finite)}")
    return valid


def compute_LR_pairs(RNA_tpm: pd.DataFrame,
                     cancer_type: str,
                     intercell_networks: dict,
                     group_lrpairs: list,
                     verbose: bool = False) -> pd.DataFrame:
    """
    Compute LR interaction scores mirroring R logic exactly, using explicit position-based grouping.
    """
    # --- Apply feature_manipulation filter as in R wrapper ---
    valid_genes = feature_manipulation(RNA_tpm, is_matrix=True, print_result=verbose)
    RNA_tpm = RNA_tpm.loc[valid_genes]
    if verbose:
        print(f"[LOG] After feature_manipulation, genes: {len(valid_genes)}")

    # 1. Log2 transform TPM
    gene_expr = np.log2(RNA_tpm + 1)
    samples = gene_expr.columns.tolist()
    genes = gene_expr.index.tolist()
    if verbose:
        print(f"[LOG] gene_expr shape: {gene_expr.shape}")

    # 2. Build raw pairs list in order
    net = intercell_networks[cancer_type]
    lr_raw = []
    for i in range(len(net['ligands'])):
        a = str(net['ligands'][i])
        b = str(net['receptors'][i])
        lr_raw.append(f"{a}_{b}")
    # unique preserving order
    seen = set(); pairs = []
    for p in lr_raw:
        if p not in seen:
            seen.add(p); pairs.append(p)
    if verbose:
        print(f"[LOG] unique pairs: {len(pairs)}")

    # 3. Compute min expression for each pair -> build values matrix
    data = []
    for p in tqdm(pairs, desc="Computing LR pairs", unit="pair"):
        a, b = p.split('_')
        if a in genes and b in genes:
            vals = gene_expr.loc[[a, b]].min(axis=0).values
        else:
            vals = np.array([np.nan] * len(samples))
        data.append(vals)
    # data shape: [num_pairs x num_samples] -> transpose
    mat = np.vstack(data).T  # samples x pairs

    # 4. Apply grouping by positions
    cols = pairs.copy()
    for grp in group_lrpairs:
        main = grp['main']
        remove = grp.get('involved_pairs', [])
        combo = grp.get('combo_name', main)
        # rename main entries
        for idx, c in enumerate(cols):
            if c == main:
                cols[idx] = combo
        # drop remove entries by first match
        for ip in remove:
            if ip in cols:
                j = cols.index(ip)
                cols.pop(j)
                mat = np.delete(mat, j, axis=1)
    # remove duplicate columns created by renaming
    unique_cols = []
    unique_idx = []
    for idx, c in enumerate(cols):
        if c not in unique_cols:
            unique_cols.append(c)
            unique_idx.append(idx)
    mat = mat[:, unique_idx]
    cols = unique_cols
    if verbose:
        print(f"[LOG] after grouping columns: {len(cols)}")

    # 5. Drop columns where any sample is NaN (sum NA behavior)
    mask = ~np.isnan(mat).any(axis=0)
    final_cols = [c for c, m in zip(cols, mask) if m]
    mat = mat[:, mask]
    if verbose:
        print(f"[LOG] after drop any NA: {len(final_cols)}")

    # return DataFrame
    df_clean = pd.DataFrame(mat, index=samples, columns=final_cols)
    return df_clean


def LR_cal(input_file: str,
           output_file: str,
           data_type: str = 'count',
           id_type: str = 'ensembl',
           cancer_type: str = 'pancan',
           verbose: bool = False):
    sep_in = detect_sep(input_file)
    df = pd.read_csv(input_file, sep=sep_in, index_col=0)
    if verbose:
        print(f"Loaded {df.shape}")
    if data_type == 'count':
        df = count2tpm(df, idType=id_type, org='hsa', source='local')
        if verbose:
            print(f"TPM {df.shape}")
    resource_path = files('iobrpy').joinpath('resources', 'lr_data.pkl')
    with as_file(resource_path) as pkl_path:
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
    intercell_networks = data['intercell_networks']
    group_lrpairs = data['group_lrpairs']
    result = compute_LR_pairs(df, cancer_type, intercell_networks, group_lrpairs, verbose)
    sep_out = detect_sep(output_file)
    result.insert(0, 'ID', result.index)
    result.to_csv(output_file, sep=sep_out, index=False)
    if verbose:
        print(f"Saved {result.shape}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, help='Path to input expression matrix (genes x samples)')
    parser.add_argument('-o','--output', required=True, help='Path to save LR scores')
    parser.add_argument('--data_type', choices=['count','tpm'], default='tpm',help='Type of input data: count or tpm')
    parser.add_argument('--id_type', default='ensembl', help='Gene ID type')
    parser.add_argument('--cancer_type', default='pancan', help='Cancer type network')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    args = parser.parse_args()
    LR_cal(args.input, args.output, args.data_type,args.id_type, args.cancer_type,args.verbose)

if __name__ == '__main__':
    main()
