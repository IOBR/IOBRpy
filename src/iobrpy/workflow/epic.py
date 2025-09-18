#!/usr/bin/env python3
"""
EPIC (Racle et al. 2017, eLife) — accelerated Python implementation

Goal: keep CLI and outputs compatible with the existing iobrpy entrypoint,
but make the core deconvolution loop much faster and keep NKcells close to R.
Key changes (internal only, no new flags):
  • Use weighted NNLS (Lawson–Hanson) for the default objective instead of
    per-sample SLSQP/trust-constr when range_based_optim is OFF.
  • Vectorized preprocessing (duplicate merge, scaling) and reduced copies.
  • Simplex projection to enforce the "sum ≤ 1" constraint efficiently.
  • NK-preserving floor (1e-12) under the simplex constraint so tiny NK signals
    are not hard-zeroed by NNLS sparsity (closer to R EPIC outcomes).

If you need the exact legacy solver behaviour, pass --rangeBasedOptim to
force the original minimize()-based path (still available).
"""
import argparse
import os
import pickle
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import nnls, minimize
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm
from importlib.resources import files
from iobrpy.utils.print_colorful_message import print_colorful_message

# ---------------- DEFAULT mRNA PER CELL ----------------
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

# Minimal positive floor for NKcells inside the simplex (keeps tiny signals)
_NK_FLOOR = 1e-12

# ---------------- INTERNAL UTILITIES -------------------
def infer_sep(filepath: str) -> str:
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.csv':
        return ','
    if ext in ('.tsv', '.txt'):
        return '\t'
    return None


def merge_duplicates(mat: pd.DataFrame, in_type=None) -> pd.DataFrame:
    """Merge duplicated gene rows by median (fast groupby)."""
    if mat.index.has_duplicates:
        n_dup = mat.index.duplicated().sum()
        warnings.warn(
            f"There are {n_dup} duplicated gene names"
            + (f" in the {in_type}" if in_type else "")
            + "; using median."
        )
        mat = mat.groupby(level=0, sort=False).median(numeric_only=True)
    return mat


def scale_counts(counts: pd.DataFrame,
                 sig_genes=None,
                 renorm_genes=None,
                 norm_fact=None):
    """CPM-like scaling to 1e6, returning (scaled_counts, norm_factor)."""
    counts = counts.loc[:, ~counts.columns.duplicated()]
    if sig_genes is None:
        sig_genes = counts.index.tolist()
    if norm_fact is None:
        if renorm_genes is None:
            renorm_genes = counts.index.tolist()
        # protect against zeros
        norm_fact = counts.loc[renorm_genes].sum(axis=0).replace(0, 1.0)
    sub = counts.loc[sig_genes]
    scaled = sub.div(norm_fact, axis=1) * 1e6
    return scaled, norm_fact


def _to_df(raw, meta):
    if isinstance(raw, pd.DataFrame):
        return raw
    if isinstance(raw, np.ndarray):
        r, c = raw.shape
        idx = meta.get('rownames', meta.get('row_names', meta.get('genes')))
        cols = meta.get('colnames', meta.get('col_names', meta.get('cellTypes')))
        if idx is None or cols is None:
            warnings.warn("ndarray missing metadata; using integer indices.")
            idx, cols = range(r), range(c)
        return pd.DataFrame(raw, index=idx, columns=cols)
    raise TypeError("Raw reference must be DataFrame or ndarray.")


def _parse_keyval(arg: str) -> dict:
    if not arg:
        return {}
    return {k: float(v) for k, v in (item.split('=') for item in arg.split(','))}


def _parse_sigfile(path: str) -> list:
    if os.path.isfile(path) and path.endswith('.gmt'):
        genes = []
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                genes.extend(parts[2:])
        return list(dict.fromkeys(genes))
    return [g for g in path.split(',') if g]


def _project_to_simplex_leq(x: np.ndarray, R: float = 1.0) -> np.ndarray:
    """Project x onto the nonnegative simplex of radius R (sum(x) ≤ R)."""
    x = np.maximum(x, 0.0)
    s = x.sum()
    if s <= R:
        return x
    # Euclidean projection onto simplex (Duchi et al., 2008)
    u = np.sort(x)[::-1]
    cssv = np.cumsum(u)
    rho = np.nonzero(u * np.arange(1, len(u)+1) > (cssv - R))[0][-1]
    theta = (cssv[rho] - R) / (rho + 1.0)
    return np.maximum(x - theta, 0.0)


def _find_nk_index(cell_types) -> int:
    """Best-effort index of NK column (handles naming variants)."""
    for i, name in enumerate(cell_types):
        s = str(name).lower().replace(' ', '').replace('_', '')
        if s.startswith('nk') or 'nkcell' in s or 'naturalkiller' in s:
            return i
    return -1


# ---------------- CORE: EPIC ---------------------------
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
         init_jitter=0.0,
         unlog_bulk=False):
    # 0) optional unlog (compat: main() exposes --unlog; pass via there)
    if unlog_bulk:
        bulk = bulk.applymap(lambda x: 2**x)

    # 1) remove all-NA genes and merge duplicates
    na_all = bulk.isna().all(axis=1)
    if na_all.any():
        warnings.warn(f"{na_all.sum()} genes are NA in all bulk samples; removing.")
        bulk = bulk.loc[~na_all]
    bulk = merge_duplicates(bulk, "bulk samples")

    refP = merge_duplicates(reference['refProfiles'], "reference profiles").loc[:, ~reference['refProfiles'].columns.duplicated()]
    refV = (merge_duplicates(reference['refProfiles.var'], "reference profiles var").loc[:, refP.columns]
            if reference.get('var_present', False) and reference.get('refProfiles.var') is not None else None)

    # 2) signature genes
    common = bulk.index.intersection(refP.index)
    if sig_genes is None:
        sig = [g for g in reference['sigGenes'] if g in common]
    else:
        sig = [g for g in sig_genes if g in common]
    sig = list(dict.fromkeys(sig))
    if len(sig) < refP.shape[1]:
        raise ValueError(f"Only {len(sig)} signature genes < {refP.shape[1]} cell types.")

    # 3) scaling (CPM 1e6) — same factor set for bulk/ref to preserve ratios
    if scale_exprs:
        bulk_s, norm_bulk = scale_counts(bulk, sig, common)
        ref_s,  norm_ref  = scale_counts(refP, sig, common)
        if refV is not None:
            refV_s, _ = scale_counts(refV, sig, common, norm_ref)
        else:
            refV_s = None
    else:
        bulk_s = bulk.loc[sig]
        ref_s  = refP.loc[sig]
        refV_s = refV.loc[sig] if refV is not None else None

    # 4) mRNA per cell
    if mRNA_cell is None:
        mRNA_cell = reference.get('mRNA_cell', mRNA_cell_default).copy()
    if mRNA_cell_sub:
        mRNA_cell.update(mRNA_cell_sub)

    # 5) compute per-gene weights (as in original EPIC)
    if refV_s is not None and not range_based_optim:
        w = (ref_s.div(refV_s + 1e-12)).sum(axis=1).to_numpy()
        med_w = np.median(w[w > 0]) if np.any(w > 0) else 1.0
        w = np.minimum(w, 100.0 * med_w)
    else:
        w = np.ones(len(sig), dtype=float)

    # 6) precompute matrices for fast path
    A = ref_s.to_numpy(dtype=float, copy=False)
    sqrtw = np.sqrt(w, dtype=float)
    Aw = (A.T * sqrtw).T  # weight rows by sqrt(w) once
    gene_order = ref_s.index
    B = bulk_s.loc[gene_order].to_numpy(dtype=float, copy=False)

    n_samples = B.shape[1]
    nC = A.shape[1]

    cell_types = list(ref_s.columns)
    nk_idx = _find_nk_index(cell_types)

    mprops = np.empty((n_samples, nC), dtype=float)
    gof_list = []

    # choose path: fast (NNLS) vs legacy minimize
    use_fast = not range_based_optim

    with tqdm(total=n_samples, desc="EPIC deconvolution", unit="sample", leave=False) as pbar:
        for i in range(n_samples):
            pbar.update(1)
            b = B[:, i]
            if use_fast:
                # Weighted NNLS: minimize ||sqrt(W)(A x - b)||^2  s.t. x ≥ 0
                bw = b * sqrtw
                x, _ = nnls(Aw, bw)

                # ---- NK-preserving projection under constraints ----
                if constrained_sum:
                    if with_other_cells:
                        if nk_idx >= 0:
                            # Reserve a tiny quota for NK, project the rest to radius (1 - eps)
                            eps = _NK_FLOOR
                            # Project to simplex of radius (1 - eps)
                            x = _project_to_simplex_leq(x, R=1.0 - eps)
                            # Give the reserved quota to NK (ensures NK ≥ eps)
                            x[nk_idx] += eps
                        else:
                            x = _project_to_simplex_leq(x, R=1.0)
                    else:
                        s = x.sum()
                        if s > 0:
                            x = x / s
                        else:
                            x[:] = 0.0
                # ----------------------------------------------------

            else:
                # Legacy minimize() path (range-based objective or if user insists)
                if refV_s is None:
                    # standard weighted LS objective
                    def fun(z):
                        r = A.dot(z) - b
                        return np.nansum(w * (r * r))
                else:
                    # range-based optimization as in original code
                    AV = refV_s.to_numpy()
                    def fun(z):
                        vmax = (A + AV).dot(z) - b
                        vmin = (A - AV).dot(z) - b
                        err = np.zeros_like(b)
                        mask = np.sign(vmax) * np.sign(vmin) == 1
                        err[mask] = np.minimum(np.abs(vmax[mask]), np.abs(vmin[mask]))
                        return np.sum(err)

                base = (1 - 1e-5) / nC
                if init_jitter > 0:
                    jitter = init_jitter * (np.random.rand(nC) - 0.5)
                    x0 = np.clip(base * (1 + jitter), 0, None)
                else:
                    x0 = np.full(nC, base)

                # start inside the feasible region and keep tiny NK positive
                if nk_idx >= 0:
                    x0[nk_idx] = max(x0[nk_idx], _NK_FLOOR)

                def make_constraints():
                    cons = [{'type': 'ineq', 'fun': lambda z: z}]
                    if constrained_sum:
                        cMin = 0 if with_other_cells else 0.99
                        cons.append({'type': 'ineq', 'fun': lambda z: np.sum(z) - cMin})
                        cons.append({'type': 'ineq', 'fun': lambda z: 1 - np.sum(z)})
                    return cons

                method = 'trust-constr' if solver == 'trust-constr' else 'SLSQP'
                res = minimize(fun, x0, method=method, constraints=make_constraints(), options={'maxiter': 200})
                x = res.x
                if not with_other_cells and constrained_sum:
                    sx = x.sum()
                    if sx > 0:
                        x = x / sx

            mprops[i, :] = x

            # GOF metrics (same as before)
            b_est = A.dot(x)
            sp = spearmanr(b, b_est)
            pe = pearsonr(b, b_est)
            try:
                a, b0 = np.polyfit(b, b_est, 1)
            except np.linalg.LinAlgError:
                a, b0 = np.nan, np.nan
            a0 = (np.sum(b * b_est) / np.sum(b * b)) if np.sum(b * b) else np.nan
            # weighted RMSE consistent with objective (divide by sqrt(n_genes) for scale)
            resid = (b_est - b) * sqrtw
            rmse = np.sqrt(np.mean(resid * resid))
            resid0 = (-b) * sqrtw
            rmse0 = np.sqrt(np.mean(resid0 * resid0))

            gof_list.append({
                'convergeCode': 0 if use_fast else getattr(res, 'status', np.nan),
                'convergeMessage': 'nnls' if use_fast else str(getattr(res, 'message', '')),
                'RMSE_weighted': rmse,
                'Root_mean_squared_geneExpr_weighted': rmse0,
                'spearmanR': sp.correlation, 'spearmanP': sp.pvalue,
                'pearsonR': pe.statistic, 'pearsonP': pe.pvalue,
                'regline_a_x': a,
                'regline_b': b0,
                'regline_a_x_through0': a0,
                'sum_mRNAProportions': np.sum(x)
            })

    mRNA_df = pd.DataFrame(mprops, index=bulk_s.columns, columns=cell_types)
    if with_other_cells:
        mRNA_df['otherCells'] = 1 - mRNA_df.sum(axis=1)
    # Convert mRNA to cell fractions using per-cell mRNA content
    denom = [mRNA_cell.get(c, mRNA_cell.get('default', 1.0)) for c in mRNA_df.columns]
    cf_raw = mRNA_df.div(denom, axis=1)
    cf = cf_raw.div(cf_raw.sum(axis=1), axis=0)

    gof_df = pd.DataFrame(gof_list, index=bulk_s.columns)
    return {'mRNAProportions': mRNA_df, 'cellFractions': cf, 'fit_gof': gof_df}


# ---------------- MAIN & I/O ----------------
def main():
    p = argparse.ArgumentParser(description="EPIC: deconvolute bulk expression like R EPIC().")
    p.add_argument('-i', '--input', required=True,
                   help="Path to the bulk expression matrix CSV/TSV (genes as rows, samples as columns)")
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
                   help="Treat bulk input as log2(expr) values and convert back via 2**x")
    p.add_argument('--solver', choices=['SLSQP','trust-constr'], default='trust-constr',
                   help="Optimization solver to use for the legacy path")
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

    # Load reference package (TRef/BRef) from iobrpy resources
    ref_pkg = files('iobrpy').joinpath('resources', 'epic_TRef_BRef.pkl')
    with ref_pkg.open('rb') as f:
        ref_data = pickle.load(f)

    # Compose the reference set(s)
    refs = []
    if args.reference in ('TRef','both'): refs.append('TRef')
    if args.reference in ('BRef','both'): refs.append('BRef')

    profs, vars_, flags, sgs = [], [], [], []
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

    # Heuristic: if no overlap, likely genes are in columns → transpose
    if not set(bulk.index).intersection(sig_ref):
        warnings.warn("Detected genes in columns; transposing bulk.")
        bulk = bulk.T

    mRNA_sub = _parse_keyval(args.mRNA_cell_sub) if args.mRNA_cell_sub else {}
    sig_file = _parse_sigfile(args.sigGenes) if args.sigGenes else None

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
        init_jitter=args.jitter,
        unlog_bulk=args.unlog
    )

    # Save cell fractions
    ext = os.path.splitext(args.out_file)[1].lower()
    sep_out = ',' if ext == '.csv' else '\t'
    out_df = res['cellFractions'].copy()
    out_df.columns = [f"{col}_EPIC" for col in out_df.columns]
    out_df.to_csv(args.out_file, sep=sep_out, index=True)

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
