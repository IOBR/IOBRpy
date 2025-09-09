#!/usr/bin/env python3
"""
Python port of the R EPIC() function (Racle et al. 2017, eLife).
This implementation mirrors the R code step-by-step to ensure identical behavior, with added debug logging and solver/jitter options.
"""
import argparse
import os
import pickle
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm
from importlib.resources import files
from iobrpy.utils.print_colorful_message import print_colorful_message

# --------------- DEFAULT mRNA PER CELL -----------------
mRNA_cell_default = {
    'Bcells': 0.4016,
    'Macrophages': 1.4196,
    'Monocytes': 1.4196,
    'Neutrophils': 0.1300,
    'NKcells': 0.4396,
    'Tcells': 0.3952,
    'CD4_Tcells': 0.3952,
    'CD8_Tcells': 0.3952,
    'Thelper': 0.3952,
    'Treg': 0.3952,
    'otherCells': 0.4000,
    'default': 0.4000
}

# --------------- INTERNAL UTILITIES -------------------
def infer_sep(filepath: str) -> str:
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.csv':
        return ','
    if ext in ('.tsv', '.txt'):
        return '\t'
    return None

def _showwarning(message, category, filename, lineno, file=None, line=None):
    tqdm.write(warnings.formatwarning(message, category, filename, lineno))

warnings.showwarning = _showwarning

def merge_duplicates(mat: pd.DataFrame, in_type=None):
    dup = mat.index.duplicated()
    if dup.any():
        genes = mat.index[dup].unique()
        warnings.warn(
            f"There are {len(genes)} duplicated gene names"
            + (f" in the {in_type}" if in_type else "")
            + "; using median."
        )
        uni = mat[~dup]
        med = pd.DataFrame({g: mat.loc[g].median(axis=0) for g in genes}).T
        med.index.name = mat.index.name
        mat = pd.concat([uni, med], axis=0)
    return mat


def scale_counts(counts: pd.DataFrame,
                 sig_genes=None,
                 renorm_genes=None,
                 norm_fact=None):
    print(f"Scale_counts ENTRY: counts.shape={counts.shape},"
          f" sig_genes={len(sig_genes) if sig_genes is not None else None},"
          f" renorm_genes={len(renorm_genes) if renorm_genes is not None else None}")
    counts = counts.loc[:, ~counts.columns.duplicated()]
    if sig_genes is None:
        sig_genes = counts.index.tolist()
    if norm_fact is None:
        if renorm_genes is None:
            renorm_genes = counts.index.tolist()
        norm_fact = counts.loc[renorm_genes].sum(axis=0)
    sub = counts.loc[sig_genes]
    print(f"Scale_counts sub.shape={sub.shape}")
    scaled = sub.div(norm_fact, axis=1) * 1e6
    print(f"Scale_counts scaled.shape={scaled.shape}")
    return scaled, norm_fact


def EPIC(bulk: pd.DataFrame,
         reference: dict,
         mRNA_cell=None,
         mRNA_cell_sub=None,
         sig_genes=None,
         scale_exprs=True,
         with_other_cells=True,
         constrained_sum=True,
         range_based_optim=False,
         solver='SLSQP',
         init_jitter=0.0):
    # 1) remove NA genes
    na_all = bulk.isna().all(axis=1)
    if na_all.any():
        warnings.warn(f"{na_all.sum()} genes are NA in all bulk samples; removing.")
        bulk = bulk.loc[~na_all]
    # 2) merge duplicates
    bulk = merge_duplicates(bulk, "bulk samples")
    refP = merge_duplicates(reference['refProfiles'], "reference profiles").loc[:, ~reference['refProfiles'].columns.duplicated()]
    refV = (merge_duplicates(reference['refProfiles.var'], "reference profiles var").loc[:, refP.columns]
            if reference.get('var_present', False) else None)
    # 3) signature genes
    common = bulk.index.intersection(refP.index)
    if sig_genes is None:
        sig = [g for g in reference['sigGenes'] if g in common]
    else:
        sig = [g for g in sig_genes if g in common]
    sig = list(dict.fromkeys(sig))
    if len(sig) < refP.shape[1]:
        raise ValueError(f"Only {len(sig)} signature genes < {refP.shape[1]} cell types.")
    # 4) scaling
    if scale_exprs:
        if len(common) < 2000:
            warnings.warn(f"Few common genes ({len(common)}) may affect scaling.")
        bulk_s, norm_bulk = scale_counts(bulk, sig, common)
        ref_s, norm_ref = scale_counts(refP, sig, common)
        refV_s, _ = (scale_counts(refV, sig, common, norm_ref)
                     if refV is not None else (None, None))
    else:
        bulk_s = bulk.loc[sig]
        ref_s = refP.loc[sig]
        refV_s = refV.loc[sig] if refV is not None else None
    bulk_s = bulk_s.loc[sig]
    ref_s = ref_s.loc[sig]
    if refV_s is not None:
        refV_s = refV_s.loc[sig]
    # 5) mRNA_cell defaults
    if mRNA_cell is None:
        mRNA_cell = reference.get('mRNA_cell', mRNA_cell_default)
    if mRNA_cell_sub:
        mRNA_cell.update(mRNA_cell_sub)
    # 6) compute weights
    if refV_s is not None and not range_based_optim:
        w = (ref_s.div(refV_s + 1e-12)).sum(axis=1).values
        med_w = np.median(w[w > 0])
        w = np.minimum(w, 100 * med_w)
    else:
        w = np.ones(len(sig))
    # 7) constraints
    nC = ref_s.shape[1]
    cMin = 0 if with_other_cells else 0.99
    def make_constraints():
        cons = [{'type': 'ineq', 'fun': lambda x: x}]
        if constrained_sum:
            cons.append({'type': 'ineq', 'fun': lambda x: np.sum(x) - cMin})
            cons.append({'type': 'ineq', 'fun': lambda x: 1 - np.sum(x)})
        return cons
    # 8) optimization per sample
    mprops, gof_list = [], []
    with tqdm(total=len(bulk_s.columns),desc="EPIC deconvolution",unit="sample",leave=False) as pbar:
        for sample in bulk_s.columns:
            pbar.update(1)
            b = bulk_s[sample].values
            A = ref_s.values
            # objective
            if not range_based_optim:
                fun = lambda x: np.nansum(w * (A.dot(x) - b)**2)
            else:
                def fun(x):
                    vmax = (A + refV_s.values).dot(x) - b
                    vmin = (A - refV_s.values).dot(x) - b
                    err = np.zeros_like(b)
                    mask = np.sign(vmax) * np.sign(vmin) == 1
                    err[mask] = np.minimum(np.abs(vmax[mask]), np.abs(vmin[mask]))
                    return np.sum(err)
            # initial guess + jitter
            base = (1 - 1e-5) / nC
            if init_jitter > 0:
                jitter = init_jitter * (np.random.rand(nC) - 0.5)
                x0 = np.clip(base * (1 + jitter), 0, None)
            else:
                x0 = np.full(nC, base)
            # solver choice
            if solver == 'trust-constr':
                res = minimize(fun, x0, method='trust-constr', constraints=make_constraints(),
                                options={'gtol':1e-6, 'barrier_tol':1e-6})
            else:
                 res = minimize(fun, x0, method='SLSQP', constraints=make_constraints())
            x = res.x
            if not with_other_cells and constrained_sum:
                x = x / x.sum()
            mprops.append(x)
            # GOF metrics
            b_est = A.dot(x)
            sp = spearmanr(b, b_est)
            pe = pearsonr(b, b_est)
            try:
                a, b0 = np.polyfit(b, b_est, 1)
            except np.linalg.LinAlgError:
                a, b0 = np.nan, np.nan
            a0 = np.sum(b * b_est) / np.sum(b * b) if np.sum(b * b) else np.nan
            rmse = np.sqrt(fun(x) / len(sig))
            rmse0 = np.sqrt(fun(np.zeros_like(x)) / len(sig))
            gof_list.append({
                'convergeCode': res.status,
                'convergeMessage': res.message,
                'RMSE_weighted': rmse,
                'Root_mean_squared_geneExpr_weighted': rmse0,
                'spearmanR': sp.correlation, 'spearmanP': sp.pvalue,
                'pearsonR': pe.statistic, 'pearsonP': pe.pvalue,
                'regline_a_x': a,
                'regline_b': b0,
                'regline_a_x_through0': a0,
                'sum_mRNAProportions': np.sum(x)
            })
    cell_types = list(ref_s.columns)
    mRNA_df = pd.DataFrame(mprops, index=bulk_s.columns, columns=cell_types)
    if with_other_cells:
        mRNA_df['otherCells'] = 1 - mRNA_df.sum(axis=1)
    denom = [mRNA_cell.get(c, mRNA_cell.get('default', 1)) for c in mRNA_df.columns]
    cf = mRNA_df.div(denom, axis=1).div(mRNA_df.div(denom, axis=1).sum(axis=1), axis=0)
    gof_df = pd.DataFrame(gof_list, index=bulk_s.columns)
    return {'mRNAProportions': mRNA_df, 'cellFractions': cf, 'fit_gof': gof_df}

# ---------------- MAIN & I/O ----------------
def parse_keyval(arg: str) -> dict:
    if not arg:
        return {}
    return {k: float(v) for k, v in (item.split('=') for item in arg.split(','))}

def parse_sigfile(path: str) -> list:
    if os.path.isfile(path) and path.endswith('.gmt'):
        genes = []
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                genes.extend(parts[2:])
        return list(dict.fromkeys(genes))
    return [g for g in path.split(',') if g]

def main():
    p = argparse.ArgumentParser(description="EPIC: deconvolute bulk expression like R EPIC().")
    p.add_argument('-i', '--input', required=True,
                   help="Path to the bulk expression matrix CSV file (genes as rows, samples as columns)")
    p.add_argument('--reference', choices=['TRef','BRef','both'], default='TRef',
                   help="Reference dataset to use for deconvolution: TRef, BRef, or both")
    p.add_argument('--mRNA_cell_sub', default=None,
                   help="Optional mRNA per cell substitutions as comma-separated key=value pairs (e.g. 'Tcells=0.5,Monocytes=1.2')")
    p.add_argument('--sigGenes', default=None,
                   help="Comma-separated list of signature genes or path to a GMT file for signatures")
    p.add_argument('--no-scaleExprs', dest='scale_exprs', action='store_false',
                   help="Skip CPM scaling of expression values (use raw counts)")
    p.add_argument('--withoutOtherCells', dest='with_other_cells', action='store_false',
                   help="Exclude an 'otherCells' component in the output cell fractions")
    p.add_argument('--no-constrainedSum', dest='constrained_sum', action='store_false',
                   help="Remove the constraint that the sum of proportions ≤ 1")
    p.add_argument('--rangeBasedOptim', dest='range_based_optim', action='store_true',
                   help="Use range-based optimization accounting for reference variance bounds")
    p.add_argument('--unlog', action='store_true', default=False,
                   help="Treat bulk input as log2(expr) values and convert back via 2**x (default ON)")
    p.add_argument('--solver', choices=['SLSQP','trust-constr'], default='trust-constr',
                   help="Optimization solver to use: 'SLSQP' or 'trust-constr'")
    p.add_argument('--jitter', type=float, default=0.0,
                   help="Relative jitter magnitude for initial proportion estimates (e.g. 1e-6)")
    p.add_argument('--seed', type=int, default=None,
                   help="Random seed for reproducible jitter initialization")
    p.add_argument('-o', '--output', dest='out_file', required=True,
                   help="Path to output file for cellFractions (csv/tsv/txt)")
    args = p.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    sep_in = infer_sep(args.input)
    bulk = pd.read_csv(args.input, sep=sep_in, index_col=0)
    if args.unlog:
        print("DEBUG: applying unlog -> 2**(log2_expr) to bulk matrix")
        bulk = bulk.applymap(lambda x: 2**x)

    ref_pkg = files('iobrpy').joinpath('resources', 'epic_TRef_BRef.pkl')
    with ref_pkg.open('rb') as f:
        ref_data = pickle.load(f)

    def _to_df(raw, meta):
        if isinstance(raw, pd.DataFrame): return raw
        if isinstance(raw, np.ndarray):
            r, c = raw.shape
            idx = meta.get('rownames', meta.get('row_names', meta.get('genes')))
            cols = meta.get('colnames', meta.get('col_names', meta.get('cellTypes')))
            if idx is None or cols is None:
                warnings.warn("ndarray missing metadata; using integer indices.")
                idx, cols = range(r), range(c)
            return pd.DataFrame(raw, index=idx, columns=cols)
        raise TypeError("Raw reference must be DataFrame or ndarray.")

    refs, profs, vars_, flags, sgs = [], [], [], [], []
    if args.reference in ('TRef','both'): refs.append('TRef')
    if args.reference in ('BRef','both'): refs.append('BRef')
    for key in refs:
        dd = ref_data[key]
        profs.append(_to_df(dd['refProfiles'], dd))
        varr = dd.get('refProfiles.var')
        if varr is not None:
            flags.append(True)
            vars_.append(_to_df(varr, dd))
        else:
            flags.append(False)
        sgs.extend(dd.get('sigGenes', []))
    ref_profiles = pd.concat(profs, axis=1)
    ref_profiles = merge_duplicates(ref_profiles, "reference profiles").loc[:, ~ref_profiles.columns.duplicated()]
    if any(flags):
        full_vars = []
        for present, vdf, prof in zip(flags, vars_, profs):
            full_vars.append(vdf if present else pd.DataFrame(0, index=prof.index, columns=prof.columns))
        ref_vars = pd.concat(full_vars, axis=1).loc[:, ref_profiles.columns]
        var_present = True
    else:
        ref_vars = None
        var_present = False
    sig_ref = [g for g in dict.fromkeys(sgs) if g in ref_profiles.index]
    reference = {
        'refProfiles': ref_profiles,
        'refProfiles.var': ref_vars,
        'sigGenes': sig_ref,
        'mRNA_cell': mRNA_cell_default,
        'var_present': var_present
    }

    if not set(bulk.index).intersection(sig_ref):
        warnings.warn("Detected genes in columns; transposing bulk.")
        bulk = bulk.T

    mRNA_sub = parse_keyval(args.mRNA_cell_sub) if args.mRNA_cell_sub else {}
    sig_file = parse_sigfile(args.sigGenes) if args.sigGenes else None

    res = EPIC(
        bulk, reference,
        mRNA_cell=None,
        mRNA_cell_sub=mRNA_sub,
        sig_genes=sig_file,
        scale_exprs=args.scale_exprs,
        with_other_cells=args.with_other_cells,
        constrained_sum=args.constrained_sum,
        range_based_optim=args.range_based_optim,
        solver=args.solver,
        init_jitter=args.jitter
    )

    ext = os.path.splitext(args.out_file)[1].lower()
    sep_out = ',' if ext == '.csv' else '\t'
    res['cellFractions'].columns = [f"{col}_EPIC" for col in res['cellFractions'].columns]
    res['cellFractions'].to_csv(args.out_file, sep=sep_out, index=True)
    print(f"Saved cellFractions ➜ {args.out_file}")
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