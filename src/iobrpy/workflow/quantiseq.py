import argparse
import pickle
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import statsmodels.api as sm
from tqdm import tqdm
import os
from importlib.resources import files

def infer_sep(path):
        ext = os.path.splitext(path)[1].lower()
        if ext == '.csv':
            return ','
        elif ext in ('.tsv', '.txt'):
            return '\t'
        else:
            return None

def make_qn(df: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalize columns of a DataFrame."""
    mat = df.values.astype(float)
    sorted_mat = np.sort(mat, axis=0)
    mean_vals = np.mean(sorted_mat, axis=1)
    ranks = np.argsort(np.argsort(mat, axis=0), axis=0)
    qn = np.zeros_like(mat)
    for i in tqdm(range(mat.shape[0]),desc="Quantile normalizing rows",total=mat.shape[0],unit="row"):
        for j in range(mat.shape[1]):
            qn[i, j] = mean_vals[ranks[i, j]]
    return pd.DataFrame(qn, index=df.index, columns=df.columns)


def map_genes(df: pd.DataFrame, hgnc: pd.DataFrame) -> pd.DataFrame:
    """Map gene symbols to approved HGNC symbols, aggregate duplicates by median."""
    curgenes = list(df.index)
    # Use R-like mapping if detailed HGNC info is available
    required_cols = {'ApprovedSymbol', 'Status', 'ApprovedName', 'PreviousSymbols', 'Synonyms'}
    if required_cols.issubset(hgnc.columns):
        # Initialize mappings
        newgenes = [None] * len(curgenes)
        newgenes2 = [None] * len(curgenes)
        hgnc_index = hgnc.set_index('ApprovedSymbol')
        # First pass: exact matches and withdrawn logic
        for idx, gene in tqdm(enumerate(curgenes),desc="Resolving primary symbols",total=len(curgenes),unit="gene"):
            if gene in hgnc_index.index:
                row = hgnc_index.loc[gene]
                status = row['Status']
                if status == 'Approved':
                    newgenes[idx] = gene
                elif status == 'EntryWithdrawn':
                    continue
                else:
                    newname = row['ApprovedName'].replace('symbolwithdrawn,see', '')
                    newgenes2[idx] = newname
        # Second pass: previous symbols and synonyms
        ps_map = hgnc_index['PreviousSymbols'].dropna().to_dict()
        sy_map = hgnc_index['Synonyms'].dropna().to_dict()
        for idx, gene in tqdm(enumerate(curgenes),desc="Checking previous/synonyms",total=len(curgenes),unit="gene"):
            if newgenes[idx] is None:
                # PreviousSymbols
                for symb, ps in ps_map.items():
                    if isinstance(ps, str) and gene in [g.strip() for g in ps.split(',')]:
                        newgenes2[idx] = symb
                        break
                # Synonyms
                if newgenes2[idx] is None:
                    for symb, sy in sy_map.items():
                        if isinstance(sy, str) and gene in [g.strip() for g in sy.split(',')]:
                            newgenes2[idx] = symb
                            break
        # Remove duplicates
        for i, val in enumerate(newgenes2):
            if val is not None and val in newgenes:
                newgenes2[i] = None
        # Fill primary names
        for i in range(len(curgenes)):
            if newgenes[i] is None and newgenes2[i] is not None:
                newgenes[i] = newgenes2[i]
        # Filter and aggregate
        keep = [i for i, v in enumerate(newgenes) if v is not None]
        df2 = df.iloc[keep].copy()
        df2.index = [newgenes[i] for i in keep]
        return df2.groupby(df2.index).median()
    else:
        # Fallback: simple mapping using available columns
        curgenes = list(df.index)
        newgenes = [None] * len(curgenes)
        if 'ApprovedSymbol' in hgnc.columns:
            hgnc_indexed = hgnc.set_index('ApprovedSymbol')
            prev_map = hgnc_indexed['PreviousSymbols'].dropna().to_dict() if 'PreviousSymbols' in hgnc_indexed.columns else {}
            syn_map = hgnc_indexed['Synonyms'].dropna().to_dict() if 'Synonyms' in hgnc_indexed.columns else {}
        else:
            hgnc_indexed, prev_map, syn_map = None, {}, {}
        # direct matches
        for idx, gene in tqdm(enumerate(curgenes),desc="Fallback: direct mapping",total=len(curgenes),unit="gene"):
            if hgnc_indexed is not None and gene in hgnc_indexed.index:
                newgenes[idx] = gene
        # previous and synonyms
        for idx, gene in tqdm(enumerate(curgenes),desc="Fallback: prev/synonyms",total=len(curgenes),unit="gene"):
            if newgenes[idx] is None:
                for symb, ps in prev_map.items():
                    if gene in [g.strip() for g in str(ps).split(',')]:
                        newgenes[idx] = symb
                        break
                if newgenes[idx] is None:
                    for symb, sy in syn_map.items():
                        if gene in [g.strip() for g in str(sy).split(',')]:
                            newgenes[idx] = symb
                            break
        keep = [i for i, v in enumerate(newgenes) if v is not None]
        df2 = df.iloc[keep].copy()
        df2.index = [newgenes[i] for i in keep]
        return df2.groupby(df2.index).median()


def fix_mixture(mix: pd.DataFrame, hgnc: pd.DataFrame, arrays: bool = False) -> pd.DataFrame:
    mix2 = map_genes(mix, hgnc)
    if mix2.values.max() < 50:
        mix2 = 2 ** mix2
    if arrays:
        mix2 = make_qn(mix2)
    mix2 = mix2.div(mix2.sum(axis=0), axis=1) * 1e6
    return mix2


def DClsei(b, A, scaling):
    # Constrained least squares with R-like scaling
    sc = np.linalg.norm(A, 2)
    A2 = A / sc
    b2 = b / sc
    def obj(x):
        return np.sum((A2.dot(x) - b2) ** 2)
    n = A2.shape[1]
    cons = [
        {'type': 'ineq', 'fun': lambda x: x},
        {'type': 'ineq', 'fun': lambda x: 1 - np.sum(x)}
    ]
    x0 = np.zeros(n)
    res = minimize(obj, x0, bounds=[(0, 1)] * n, constraints=cons)
    est = res.x
    tot = est.sum()
    est = np.where(scaling > 0, est / scaling, est)
    if est.sum() > 0:
        est = est / est.sum() * tot
    other = max(0.0, 1.0 - est.sum())
    return np.concatenate([est, [other]])


def DCrr(b, A, method, scaling):
    exog = sm.add_constant(A)
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
    est = np.clip(res.params[1:], 0, None)
    tot0 = est.sum()
    if tot0 > 0:
        est = est / tot0
        est = np.where(scaling > 0, est / scaling, est)
        if est.sum() > 0:
            est = est / est.sum() * tot0
    return est


def quanTIseq(sig, mix, scaling, method):
    genes = sig.index.intersection(mix.index)
    A = sig.loc[genes].values
    mix_arr = mix.loc[genes].values
    res_list = []
    for i in tqdm(range(mix_arr.shape[1]), desc="Deconvoluting samples"):
        b = mix_arr[:, i]
        if method == 'lsei':
            res = DClsei(b, A, scaling)
        else:
            res = DCrr(b, A, method, scaling)
        res_list.append(res)
    mat = np.vstack(res_list)
    cols = list(sig.columns) + (['Other'] if method == 'lsei' else [])
    return pd.DataFrame(mat, index=mix.columns, columns=cols)


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
        mRNA = series.reindex(sig.columns).fillna(1).values
    else:
        mRNA = np.ones(len(sig.columns))

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
    print(f"Signature genes found in data set: {us}/{ns} ({us*100/ns:.2f}%)")
    print(f"Mixture deconvolution (method: {method})")
    res1 = quanTIseq(sig, mix2, mRNA, method)

    # correct low Tregs cases
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

if __name__ == '__main__':
    main()