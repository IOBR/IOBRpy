import argparse
import pickle
import numpy as np
import pandas as pd
from scipy.optimize import nnls  # fast non-negative least squares
import statsmodels.api as sm
from tqdm import tqdm
import os
from importlib.resources import files
from iobrpy.utils.print_colorful_message import print_colorful_message

# -----------------------------
# Fast helpers (vectorized)
# -----------------------------

def infer_sep(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == '.csv':
        return ','
    elif ext in ('.tsv', '.txt'):
        return '\t'
    else:
        return None

def make_qn(df: pd.DataFrame) -> pd.DataFrame:
    """
    Quantile normalize columns of a DataFrame (fast, vectorized).

    Algorithm:
    1) For each column, sort its values (ascending).
    2) Take the row-wise mean across all sorted columns (mean distribution).
    3) For each column, assign the mean distribution back in the positions
       of that column's sorted order.

    This avoids Python-level nested loops.
    """
    if df.shape[1] == 0 or df.shape[0] == 0:
        return df.copy()

    vals = df.to_numpy(dtype=float, copy=False)
    order = np.argsort(vals, axis=0)                 # indices that would sort each column
    sorted_vals = np.take_along_axis(vals, order, axis=0)
    mean_sorted = sorted_vals.mean(axis=1)           # average at each rank across columns

    # Place the averaged ranks back to original positions
    out = np.empty_like(vals)
    # For each column j, positions order[:, j] should receive mean_sorted
    for j in range(vals.shape[1]):
        out[order[:, j], j] = mean_sorted

    return pd.DataFrame(out, index=df.index, columns=df.columns)

def _build_hgnc_alias_map(hgnc: pd.DataFrame):
    """
    Build a dict mapping any alias (ApprovedSymbol, PreviousSymbols, Synonyms) to ApprovedSymbol.
    Uses simple string splitting; robust to missing columns.
    """
    alias_to_approved = {}

    if hgnc is None or hgnc.empty:
        return alias_to_approved

    cols = set(hgnc.columns)
    has_prev = 'PreviousSymbols' in cols
    has_syn  = 'Synonyms' in cols
    has_approved = 'ApprovedSymbol' in cols

    if not has_approved:
        return alias_to_approved

    # Iterate rows once; O(n) and typically fast
    for _, row in hgnc[['ApprovedSymbol'] + ([ 'PreviousSymbols'] if has_prev else []) + ([ 'Synonyms'] if has_syn else [])].fillna('').iterrows():
        approved = str(row['ApprovedSymbol']).strip()
        if approved:
            # map the approved symbol to itself
            if approved not in alias_to_approved:
                alias_to_approved[approved] = approved

        if has_prev:
            for token in str(row.get('PreviousSymbols', '')).split(','):
                alias = token.strip()
                if alias and alias not in alias_to_approved:
                    alias_to_approved[alias] = approved

        if has_syn:
            for token in str(row.get('Synonyms', '')).split(','):
                alias = token.strip()
                if alias and alias not in alias_to_approved:
                    alias_to_approved[alias] = approved

    return alias_to_approved

def map_genes(df: pd.DataFrame, hgnc: pd.DataFrame) -> pd.DataFrame:
    """
    Map gene symbols to Approved HGNC symbols and aggregate duplicates by median.
    This is a faster, vectorized rewrite of the previous multi-pass logic.
    """
    if df is None or df.empty:
        return df.copy()

    alias_map = _build_hgnc_alias_map(hgnc)
    if not alias_map:
        # No mapping available; return as-is
        return df.copy()

    # Map index via a vectorized Series.map with a fallback to "keep original if already approved"
    idx = df.index.to_series()
    mapped = idx.map(lambda g: alias_map.get(g, alias_map.get(str(g).strip(), None)))
    # If still missing, keep original ONLY if it is an approved symbol
    approved_set = {v for v in alias_map.values()}
    mapped = mapped.where(mapped.notna(), idx.where(idx.isin(approved_set), np.nan))

    # Drop genes we cannot map; aggregate duplicates by median to match R behavior
    keep = mapped.notna()
    if not keep.any():
        # In pathological cases, return empty frame with original columns
        return pd.DataFrame(index=[], columns=df.columns, dtype=float)

    df2 = df.loc[keep].copy()
    df2.index = mapped[keep].astype(str).values
    return df2.groupby(df2.index).median()

def fix_mixture(mix: pd.DataFrame, hgnc: pd.DataFrame, arrays: bool = False) -> pd.DataFrame:
    """
    1) HGNC mapping (fast).
    2) If values look log-scaled (max < 50), convert back with 2**x.
    3) Optional quantile normalization for array data (fast, vectorized).
    4) Column-wise TPM-like scaling to 1e6 for numerical stability.
    """
    mix2 = map_genes(mix, hgnc)
    if mix2.size == 0:
        return mix2

    if np.nanmax(mix2.values) < 50:  # heuristic: log-scale check
        mix2 = pd.DataFrame(np.power(2.0, mix2.values), index=mix2.index, columns=mix2.columns)

    if arrays:
        mix2 = make_qn(mix2)

    col_sums = mix2.sum(axis=0).to_numpy(dtype=float, copy=False)
    # Avoid division by zero
    col_sums[col_sums == 0] = 1.0
    mix2 = mix2.div(col_sums, axis=1) * 1e6
    return mix2

# -----------------------------
# Core solvers
# -----------------------------

def DClsei(b, A, scaling):
    """
    Fast constrained least squares using NNLS (non-negative least squares).

    We solve:  min ||A x - b||_2  subject to x >= 0
    Then apply the quanTIseq mRNA scaling correction in the same way as before,
    and finally compute the 'Other' fraction as max(0, 1 - sum(x)).

    This replaces the per-sample SLSQP minimize() call with nnls(), which is
    typically orders of magnitude faster while preserving the non-negativity
    constraint and the downstream normalization behavior.
    """
    # R-like scaling to stabilize
    sc = np.linalg.norm(A, 2)
    if not np.isfinite(sc) or sc == 0:
        sc = 1.0
    A2 = A / sc
    b2 = b / sc

    # Lawson–Hanson NNLS
    est, _ = nnls(A2, b2)

    # Keep the original sum (tot), then apply mRNA scaling and re-normalize to tot
    tot = est.sum()
    if scaling is not None and len(scaling) == est.shape[0]:
        safe = np.where(np.asarray(scaling) > 0, np.asarray(scaling), 1.0)
        est = est / safe
        s = est.sum()
        if s > 0:
            est = est / s * tot

    other = max(0.0, 1.0 - est.sum())
    return np.concatenate([est, [other]])

def DCrr(b, A, method, scaling):
    """
    Robust regression (unchanged), but with small vectorization micro-optimizations.
    """
    exog = sm.add_constant(A, prepend=True, has_constant='add')
    m_norm = {
        'hampel': sm.robust.norms.Hampel(1.5, 3.5, 8),
        'huber': sm.robust.norms.HuberT(),
        'bisquare': sm.robust.norms.TukeyBiweight()
    }
    M = m_norm.get(method)
    if M is None:
        raise ValueError(f"Unknown method {method}")
    rlm = sm.RLM(b, exog, M=M)
    res = rlm.fit(maxiter=1000)
    est = np.clip(res.params[1:], 0, None)  # drop intercept, enforce non-negativity

    tot0 = est.sum()
    if tot0 > 0:
        est = est / tot0
        if scaling is not None and len(scaling) == est.shape[0]:
            safe = np.where(np.asarray(scaling) > 0, np.asarray(scaling), 1.0)
            est = est / safe
            s = est.sum()
            if s > 0:
                est = est / s * tot0
    return est

def quanTIseq(sig: pd.DataFrame, mix: pd.DataFrame, scaling, method: str) -> pd.DataFrame:
    """
    Deconvolute per sample.
    Keeps the original CLI semantics; only speeds up internals.
    """
    genes = sig.index.intersection(mix.index)
    if genes.empty:
        raise ValueError("No overlapping genes between signature and mixture matrix.")
    A = sig.loc[genes].to_numpy(dtype=float, copy=False)
    B = mix.loc[genes].to_numpy(dtype=float, copy=False)

    n_samples = B.shape[1]
    # Pre-allocate output
    if method == 'lsei':
        out = np.empty((n_samples, A.shape[1] + 1), dtype=float)  # +1 for 'Other'
    else:
        out = np.empty((n_samples, A.shape[1]), dtype=float)

    for i in tqdm(range(n_samples), desc="Deconvoluting samples"):
        b = B[:, i]
        if method == 'lsei':
            out[i, :] = DClsei(b, A, scaling)
        else:
            out[i, :] = DCrr(b, A, method, scaling)

    cols = list(sig.columns) + (['Other'] if method == 'lsei' else [])
    return pd.DataFrame(out, index=mix.columns, columns=cols)

# -----------------------------
# High-level entry
# -----------------------------

def deconvolute_quantiseq_default(mix, data, arrays=False, signame='TIL10', tumor=False, mRNAscale=True, method='lsei', rmgenes='unassigned'):
    print("Running quanTIseq deconvolution module")
    if rmgenes == 'unassigned':
        rmgenes = 'none' if arrays else 'default'
    if signame != 'TIL10':
        raise ValueError("Only TIL10 supported currently")

    sig = data['TIL10_signature'].copy()
    mRNA_df = data['TIL10_mRNA_scaling']
    lrm = data['TIL10_rmgenes']

    # extract mRNA scaling factors
    if mRNAscale:
        if isinstance(mRNA_df, pd.DataFrame):
            if 'celltype' in mRNA_df.columns and 'scaling' in mRNA_df.columns:
                series = mRNA_df.set_index('celltype')['scaling']
            elif mRNA_df.shape[1] == 2:
                series = mRNA_df.set_index(mRNA_df.columns[0])[mRNA_df.columns[1]]
            else:
                series = pd.Series(mRNA_df.iloc[:, 1].values, index=mRNA_df.iloc[:, 0].values)
        elif isinstance(mRNA_df, dict):
            series = pd.Series(mRNA_df)
        else:
            raise ValueError("Unrecognized mRNA scaling structure")
        mRNA = series.reindex(sig.columns).fillna(1).to_numpy()
    else:
        mRNA = np.ones(len(sig.columns), dtype=float)

    print(f"Gene expression normalization and re-annotation (arrays: {arrays})")
    mix2 = fix_mixture(mix, data.get('HGNC_genenames_20170418', pd.DataFrame()), arrays)

    if rmgenes != 'none':
        n1 = sig.shape[0]
        sig = sig.drop(index=lrm, errors='ignore')
        n2 = sig.shape[0]
        print(f"Removing {n1-n2} noisy genes")

    if tumor:
        ab = data.get('TIL10_TCGA_aberrant_immune_genes', [])
        n1 = sig.shape[0]
        sig = sig.drop(index=ab, errors='ignore')
        n2 = sig.shape[0]
        print(f"Removing {n1-n2} genes with high expression in tumors")

    ns = sig.shape[0]
    us = len(sig.index.intersection(mix2.index))
    print(f"Signature genes found in data set: {us}/{ns} ({(us*100.0/ns):.2f}%)")
    print(f"Mixture deconvolution (method: {method})")
    res1 = quanTIseq(sig, mix2, mRNA, method)

    # correct low Tregs cases (unchanged logic)
    if method == 'lsei' and {'Tregs', 'T.cells.CD4'}.issubset(set(sig.columns)):
        minT = 0.02
        i_cd4 = sig.columns.get_loc('T.cells.CD4')
        sig2 = sig.drop(columns=['T.cells.CD4'])
        m2 = np.delete(mRNA, i_cd4)
        r2 = quanTIseq(sig2, mix2, m2, method)
        mask = res1['Tregs'] < minT
        avgT = (r2.loc[mask, 'Tregs'] + res1.loc[mask, 'Tregs']) / 2
        res1.loc[mask, 'Tregs'] = avgT
        res1.loc[mask, 'T.cells.CD4'] = np.maximum(0, res1.loc[mask, 'T.cells.CD4'] - avgT)

    # Final normalization and formatting
    res = res1.div(res1.sum(axis=1), axis=0).reset_index()
    res.insert(0, 'Sample', res.pop('index'))
    print("Deconvolution successful!")
    return res

def main():
    parser = argparse.ArgumentParser(
        description="Deconvolute cell-type fractions from bulk RNA-seq using quanTIseq algorithm"
    )
    parser.add_argument(
        '-i', '--input',
        dest='mix',
        required=True,
        help="Path to the input mixture matrix TSV file (genes x samples)"
    )
    parser.add_argument(
        '-o', '--output',
        dest='out',
        required=True,
        help="Path where the deconvolution results TSV will be written"
    )
    parser.add_argument(
        '--arrays',
        action='store_true',
        help="Perform quantile normalization on array data before deconvolution"
    )
    parser.add_argument(
        '--signame',
        default='TIL10',
        help="Name of the signature set to use (default: TIL10)"
    )
    parser.add_argument(
        '--tumor',
        action='store_true',
        help="Remove genes with high expression in tumor samples"
    )
    parser.add_argument(
        '--scale_mrna',
        dest='mRNAscale',
        action='store_true',
        help="Enable mRNA scaling; use raw signature proportions otherwise"
    )
    parser.add_argument(
        '--method',
        choices=['lsei','hampel','huber','bisquare'],
        default='lsei',
        help="Robust regression method to use: lsei (least squares), or robust norms (hampel, huber, bisquare)"
    )
    parser.add_argument(
        '--rmgenes',
        default='unassigned',
        help="List of genes to remove (e.g., 'default', 'none', or comma-separated identifiers)"
    )
    args = parser.parse_args()

    data_file = files('iobrpy.resources').joinpath('quantiseq_data.pkl')
    with data_file.open('rb') as fh:
        data = pickle.load(fh)

    in_sep = infer_sep(args.mix)
    mix = pd.read_csv(args.mix, sep=in_sep, index_col=0)

    res = deconvolute_quantiseq_default(
        mix, data,
        arrays=args.arrays,
        signame=args.signame,
        tumor=args.tumor,
        mRNAscale=args.mRNAscale,
        method=args.method,
        rmgenes=args.rmgenes
    )

    out_sep = infer_sep(args.out) or '\t'
    res.columns = [
        'ID' if col == 'Sample'
        else f"{col.replace('.', '_')}_quantiseq"
        for col in res.columns
    ]
    eps = 1e-8
    num_cols = res.columns.drop('ID')
    res[num_cols] = res[num_cols].mask(res[num_cols].abs() < eps, 0)
    res.to_csv(args.out, sep=out_sep, index=False)

    print(f"Results saved to {args.out}")
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

if __name__ == '__main__':
    main()
