#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from importlib.resources import files
from joblib import Parallel, delayed

try:
    import gseapy as gp
except ImportError:
    gp = None

# Optional pure-python ssGSEA via decoupler (no R/gseapy dependency)
try:
    import decoupler as dc
except ImportError:
    dc = None

try:
    from tqdm.auto import tqdm
except Exception:
    # fallback: if tqdm is not installed, tqdm(...) just returns the iterable unchanged
    def tqdm(x, **kwargs):
        return x

# Supported methods
signature_score_calculation_methods = {
    "pca": "pca",
    "zscore": "zscore",
    "ssgsea": "ssgsea",
    "integration": "integration"
}

# Debug flag
default_debug = False

def _merge_signature_groups(all_sigs, signature_names):
    """
    Merge multiple top-level signature GROUP dicts (name->genes) into one dict.
    - signature_names: list of group names, may contain comma-separated tokens
    - supports special token 'all' (case-insensitive) to use all dict-valued groups
    - if a signature name collides with different gene sets, rename to '<name>__<group>'
    """
    expanded = []
    for tok in signature_names:
        expanded.extend([t.strip() for t in str(tok).split(',') if t.strip()])

    if any(t.lower() == 'all' for t in expanded):
        selected_groups = [k for k, v in all_sigs.items() if isinstance(v, dict)]
    else:
        selected_groups = expanded

    if default_debug:
        print(f"DEBUG: Selected groups -> {selected_groups}")

    combined = {}
    for grp in selected_groups:
        sig_dict = all_sigs.get(grp)
        if not isinstance(sig_dict, dict):
            raise KeyError(f"Signature group '{grp}' must map to dict(name->genes) in the pickle.")
        for sig_name, genes in sig_dict.items():
            if sig_name in combined:
                if set(combined[sig_name]) == set(genes):
                    continue
                new_name = f"{sig_name}__{grp}"
                if default_debug:
                    print(f"DEBUG: name collision on '{sig_name}', renamed to '{new_name}'")
                combined[new_name] = list(genes)
            else:
                combined[sig_name] = list(genes)
    if default_debug:
        print(f"DEBUG: Combined signatures: {len(combined)}")
    return combined

def load_signatures(pkl_path):
    """Load signature collections from a pickle file."""
    return pd.read_pickle(pkl_path)

def _log2eset_like(eset: pd.DataFrame) -> pd.DataFrame:
    """Emulate R's log2eset heuristic (only log when distribution suggests it)."""
    q = np.quantile(eset.values, [0.0, 0.25, 0.5, 0.75, 0.99, 1.0])
    log_judge = (
        (q[4] > 100)
        or (q[5] - q[0] > 50 and q[1] > 0)
        or (q[1] > 0 and q[1] < 1 and q[3] > 1 and q[3] < 2)
    )
    if log_judge:
        eset = eset.copy()
        eset[eset < 0] = 0
        eset = np.log2(eset + 1)
    return eset


def _feature_manipulation_like(eset: pd.DataFrame) -> pd.DataFrame:
    """Filter genes with NA/Inf/non‑numeric/zero‑sd like R's feature_manipulation."""
    eset = eset.copy()
    numeric_cols = [c for c in eset.columns if np.issubdtype(eset[c].dtype, np.number)]
    eset = eset[numeric_cols]

    finite = np.isfinite(eset).all(axis=1)
    eset = eset.loc[finite]

    sd = eset.std(axis=1, ddof=1)
    eset = eset.loc[sd > 0]
    return eset


def preprocess_eset(eset, adjust_eset):
    """Log2-transform conditionally and optionally drop problematic genes."""
    if default_debug:
        print(f"DEBUG: Preprocess: shape={eset.shape}, adjust={adjust_eset}")
    eset_proc = eset.copy()
    if eset_proc.shape[1] < 5000:
        eset_proc = _log2eset_like(eset_proc)
    if adjust_eset:
        eset_proc = _feature_manipulation_like(eset_proc)
    if default_debug:
        print(f"DEBUG: After preprocess: shape={eset_proc.shape}")
    return eset_proc

def filter_signatures(sig_dict, eset, min_genes):
    """Keep only signatures with at least min_genes present in eset."""
    out = {}
    for name, genes in sig_dict.items():
        present = [g for g in genes if g in eset.index]
        if default_debug:
            print(f"DEBUG: {name}: {len(present)} genes in eset")
        if len(present) >= min_genes:
            out[name] = present
    if default_debug:
        print(f"DEBUG: {len(out)} signatures retained (min_genes={min_genes})")
    return out


def _ssgsea_numpy(eset: pd.DataFrame, sigs: dict, min_size: int, parallel_size: int):
    """Lightweight ssGSEA implementation using numpy only (samples as columns)."""
    samples = list(eset.columns)
    sig_names = list(sigs.keys())

    def _score_sample(sample):
        expr = eset[sample]
        order = expr.sort_values(ascending=False)
        genes_order = order.index.to_numpy()
        values_order = order.to_numpy()
        total = len(genes_order)

        sample_scores = {}
        for name, genes in sigs.items():
            hit_mask = np.isin(genes_order, genes)
            Nh = int(hit_mask.sum())
            if Nh < min_size or Nh == 0 or Nh == total:
                sample_scores[name] = 0.0
                continue
            miss = total - Nh
            # weighting exponent matches GSVA/ssGSEA default (0.25)
            hit_weights = np.power(np.abs(values_order), 0.25) * hit_mask
            sum_hit = hit_weights.sum()
            if sum_hit == 0:
                sample_scores[name] = 0.0
                continue
            Phit = np.cumsum(hit_weights / sum_hit)
            Pmiss = np.cumsum(~hit_mask / miss)
            walk = Phit - Pmiss
            sample_scores[name] = walk.sum()
        return sample, sample_scores

    if parallel_size and parallel_size > 1:
        scored = Parallel(n_jobs=int(parallel_size), prefer="threads")(
            delayed(_score_sample)(s) for s in tqdm(samples, desc="ssGSEA samples", unit="sample")
        )
    else:
        scored = [_score_sample(s) for s in tqdm(samples, desc="ssGSEA samples", unit="sample")]

    mat = pd.DataFrame(index=samples, columns=sig_names, dtype=float)
    for sample, d in scored:
        for k, v in d.items():
            mat.at[sample, k] = v
    mat.reset_index(inplace=True)
    mat.rename(columns={'index': 'ID'}, inplace=True)
    return mat

def sig_score_pca(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size=1):
    """
    Parallel PCA implementation: one PCA (PC1) per signature.
    Preserves original semantics (center+scale per gene, PCA(1), flip by corr with mean expression).
    """
    pdata = pd.DataFrame({'ID': eset.columns})
    eset2 = preprocess_eset(eset, adjust_eset)

    min_size = max(mini_gene_count, 2)
    sigs = filter_signatures(sig_dict, eset2, min_size)

    items = list(sigs.items())

    def _one(name, genes):
        valid = sorted(set(genes) & set(eset2.index))
        if len(valid) < 2:
            # Fallback: all zeros if not enough genes
            return name, np.zeros(len(eset2.columns), dtype=float)
        tmp = eset2.loc[valid]               # genes × samples
        mat = tmp.T                          # samples × genes
        # z-score by gene
        mat = mat.sub(mat.mean(axis=0), axis=1)
        mat = mat.div(mat.std(axis=0, ddof=1).replace(0, np.nan), axis=1).fillna(0.0)

        pca = PCA(n_components=1, svd_solver='full', random_state=0)
        pc1 = pca.fit_transform(mat.values)[:, 0]   # length = n_samples

        mean_expr = tmp.mean(axis=0).values         # length = n_samples
        # Spearman was heavy; Pearson on standardized data equals Spearman rank corr approx. Keep Pearson as original.
        corr = np.corrcoef(pc1, mean_expr)[0, 1]
        direction = np.sign(corr) if not np.isnan(corr) else 1.0
        return name, (pc1 * direction)

    if parallel_size and parallel_size > 1:
        results = Parallel(n_jobs=int(parallel_size), prefer="processes")(
            delayed(_one)(name, genes) for name, genes in tqdm(items, desc="PCA signatures", unit="sig")
        )
    else:
        results = [ _one(name, genes) for name, genes in tqdm(items, desc="PCA signatures", unit="sig") ]

    for name, vec in results:
        pdata[name] = vec

    # TME contrasts (unchanged)
    if {'TMEscoreA_CIR','TMEscoreB_CIR'}.issubset(sigs):
        pdata['TMEscore_CIR'] = pdata['TMEscoreA_CIR'] - pdata['TMEscoreB_CIR']
    if {'TMEscoreA_plus','TMEscoreB_plus'}.issubset(sigs):
        pdata['TMEscore_plus'] = pdata['TMEscoreA_plus'] - pdata['TMEscoreB_plus']

    return pdata

def sig_score_zscore(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size=1):
    """
    Parallel zscore implementation: for each signature take mean across genes (colMeans).
    Preserves original semantics.
    """
    pdata = pd.DataFrame({'ID': eset.columns})
    eset2 = preprocess_eset(eset, adjust_eset)

    min_size = max(mini_gene_count, 2)
    sigs = filter_signatures(sig_dict, eset2, min_size)

    items = list(sigs.items())

    def _one(name, genes):
        valid = sorted(set(genes) & set(eset2.index))
        if len(valid) == 0:
            return name, np.zeros(len(eset2.columns), dtype=float)
        mat = eset2.loc[valid]               # genes × samples
        return name, mat.mean(axis=0).values

    if parallel_size and parallel_size > 1:
        results = Parallel(n_jobs=int(parallel_size), prefer="processes")(
            delayed(_one)(name, genes) for name, genes in tqdm(items, desc="zscore signatures", unit="sig")
        )
    else:
        results = [ _one(name, genes) for name, genes in tqdm(items, desc="zscore signatures", unit="sig") ]

    for name, vec in results:
        pdata[name] = vec

    # TME contrasts (unchanged)
    if {'TMEscoreA_CIR','TMEscoreB_CIR'}.issubset(sigs):
        pdata['TMEscore_CIR'] = pdata['TMEscoreA_CIR'] - pdata['TMEscoreB_CIR']
    if {'TMEscoreA_plus','TMEscoreB_plus'}.issubset(sigs):
        pdata['TMEscore_plus'] = pdata['TMEscoreA_plus'] - pdata['TMEscoreB_plus']

    return pdata

def sig_score_ssgsea(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size):
    # Preprocess like R
    eset2 = preprocess_eset(eset, adjust_eset)
    # First filter with original threshold
    sigs = filter_signatures(sig_dict, eset2, mini_gene_count)
    # Then enforce min_size >= 5
    min_size = max(mini_gene_count, 5)
    sigs = filter_signatures(sig_dict, eset2, min_size)

    # Try decoupler's pure python implementation first
    if dc is not None:
        print("Running ssGSEA via decoupler (pure Python)...")
        # decoupler expects samples as rows
        df = eset2.T
        es_mat = dc.ssgsea(
            mat=df,
            net=[(g, name, 1.0) for name, genes in sigs.items() for g in genes],
            source="source",
            target="target",
            weight="weight",
            min_n=min_size,
            max_n=eset2.shape[0],
            exp=0.25,
            norm=True,
            n_proc=parallel_size
        )
        es = es_mat.reset_index()
        es.rename(columns={es.columns[0]: 'ID'}, inplace=True)
    elif gp is not None:
        print("Running ssGSEA via gseapy (this may take a while)...")
        ss = gp.ssgsea(
            data=eset2,
            gene_sets=sigs,
            outdir=None,
            sample_norm_method='rank',    # rank-based kernel = Gaussian
            weight=0.25,                  # match GSVA default tau
            permutation_num=0,
            no_plot=True,
            threads=parallel_size,
            min_size=min_size,
            max_size=eset2.shape[0],
            ssgsea_norm=True
        )

        # Pivot to samples × terms, keep Term as columns; use raw ES to mirror GSVA output
        es = ss.res2d.pivot(index='Term', columns='Name', values='ES').T.reset_index()
        es.rename(columns={'Name': 'ID'}, inplace=True)
    else:
        # Lightweight numpy fallback
        print("Running ssGSEA via numpy fallback (no decoupler/gseapy detected)...")
        es = _ssgsea_numpy(eset2, sigs, min_size, parallel_size)

    if 'TMEscoreA_CIR' in es.columns and 'TMEscoreB_CIR' in es.columns:
        es['TMEscore_CIR'] = es['TMEscoreA_CIR'] - es['TMEscoreB_CIR']
    if 'TMEscoreA_plus' in es.columns and 'TMEscoreB_plus' in es.columns:
        es['TMEscore_plus'] = es['TMEscoreA_plus'] - es['TMEscoreB_plus']
    return es

def sig_score_integration(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size):
    filtered_sigs = {
        name: [g for g in genes if g in eset.index]
        for name, genes in sig_dict.items()
        if len([g for g in genes if g in eset.index]) >= mini_gene_count
    }

    p = sig_score_pca(eset, filtered_sigs, mini_gene_count, adjust_eset)
    p = p.set_index('ID').add_suffix('_PCA')

    z = sig_score_zscore(eset, filtered_sigs, mini_gene_count, adjust_eset)
    z = z.set_index('ID').add_suffix('_zscore')

    eset2 = preprocess_eset(eset, adjust_eset)

    if dc is not None:
        print("Running ssGSEA via decoupler (pure Python)...")
        df = eset2.T
        es_mat = dc.ssgsea(
            mat=df,
            net=[(g, name, 1.0) for name, genes in filtered_sigs.items() for g in genes],
            source="source",
            target="target",
            weight="weight",
            min_n=mini_gene_count,
            max_n=eset2.shape[0],
            exp=0.25,
            norm=True,
            n_proc=parallel_size
        )
        es = es_mat.reset_index()
        es.rename(columns={es.columns[0]: 'ID'}, inplace=True)
    elif gp is not None:
        print("Running ssGSEA via gseapy (this may take a while)...")
        ss = gp.ssgsea(
            data=eset2,
            gene_sets=filtered_sigs,
            outdir=None,
            sample_norm_method='rank',
            weight=0.25,
            permutation_num=0,
            no_plot=True,
            threads=parallel_size,
            min_size=mini_gene_count,
            max_size=eset2.shape[0],
            ssgsea_norm=True
        )
        es = ss.res2d.pivot(index='Term', columns='Name', values='ES').T.reset_index()
        es.rename(columns={'Name': 'ID'}, inplace=True)
    else:
        print("Running ssGSEA via numpy fallback (no decoupler/gseapy detected)...")
        es = _ssgsea_numpy(eset2, filtered_sigs, mini_gene_count, parallel_size)

    if 'TMEscoreA_CIR' in es.columns and 'TMEscoreB_CIR' in es.columns:
        es['TMEscore_CIR'] = es['TMEscoreA_CIR'] - es['TMEscoreB_CIR']
    if 'TMEscoreA_plus' in es.columns and 'TMEscoreB_plus' in es.columns:
        es['TMEscore_plus'] = es['TMEscoreA_plus'] - es['TMEscoreB_plus']
    s = es.set_index('ID').add_suffix('_ssGSEA')

    return pd.concat([p, z, s], axis=1).reset_index()

def calculate_sig_score(eset, signature_names, method, mini_gene_count, adjust_eset, parallel_size):
    resource_pkg = 'iobrpy.resources'
    resource_path = files(resource_pkg).joinpath('calculate_data.pkl')
    all_sigs = pd.read_pickle(resource_path)
    sig_dict = _merge_signature_groups(all_sigs, signature_names)
    if not isinstance(sig_dict, dict) or len(sig_dict) == 0:
        raise KeyError(f"No valid signatures found from groups: {signature_names}")
    m = method.lower()
    if m == 'pca': return sig_score_pca(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size)
    if m == 'zscore': return sig_score_zscore(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size)
    if m == 'ssgsea': return sig_score_ssgsea(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size)
    if m == 'integration': return sig_score_integration(eset, sig_dict, mini_gene_count, adjust_eset, parallel_size)
    raise ValueError("Unknown method")

def main():
    p = argparse.ArgumentParser(description="Calculate signature scores (PCA, z-score, ssGSEA, integration).")
    p.add_argument(
        '-i', '--input',
        dest='input_matrix',
        required=True,
        help='Path to input expression matrix (CSV/TSV, genes×samples)'
    )
    p.add_argument(
        '--signature',
        required=True,
        nargs='+',
        help=('One or more signature GROUP names to use. '
              'Examples: signature_collection signature_tme  (space-separated), '
              'or signature_collection,signature_tme (comma-separated). '
              'Supported groups include: go_bp, go_cc, go_mf, '
              'signature_collection, signature_tme, signature_sc, signature_tumor, '
              'signature_metabolism, kegg, hallmark, reactome, or "all" to use all groups.')
    )
    p.add_argument(
        '--method',
        default='pca',
        choices=list(signature_score_calculation_methods.values()),
        help='Scoring method to apply: "pca", "zscore", "ssgsea" or "integration"'
    )
    p.add_argument(
        '--mini_gene_count',
        type=int,
        default=3,
        help='Minimum number of genes required in a signature to be scored'
    )
    p.add_argument(
        '--adjust_eset',
        action='store_true',
        help='Whether to perform additional Inf/zero‐sd filtering after log2 transform'
    )
    p.add_argument(
        '--parallel_size',
        type=int,
        default=1,
        help='Threads for scoring (PCA/zscore/ssGSEA)'
    )
    p.add_argument(
        '-o', '--output',
        dest='output_matrix',
        required=True,
        help='Path to save the resulting scores matrix'
    )
    p.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug output'
    )
    args = p.parse_args()

    if args.debug:
        global default_debug
        default_debug = True
        print("DEBUG args:", vars(args))

    # load expression
    eset = pd.read_csv(args.input_matrix, sep=None, engine='python', index_col=0)
    if args.debug:
        print("DEBUG eset shape:", eset.shape)

    # calculate
    res = calculate_sig_score(eset, args.signature, args.method,
                              args.mini_gene_count, args.adjust_eset,
                              args.parallel_size)

    # final save (unchanged)
    # for ssgsea we want index (ID column); others index=False
    if args.method.lower() == 'ssgsea':
        res.to_csv(args.output_matrix, index=False)
    else:
        res.to_csv(args.output_matrix, index=False)

    print(f"Saved scores to {args.output_matrix}")

if __name__ == '__main__':
    main()