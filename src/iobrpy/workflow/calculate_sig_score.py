#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import inspect
from sklearn.decomposition import PCA
from importlib.resources import files
from joblib import Parallel, delayed

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


def _run_single_sample_gsea_backend(eset: pd.DataFrame, sigs: dict, min_size: int, parallel_size: int):
    """Try to execute ssGSEA using the ``single-sample-gsea`` package if available.

    The package exposes different entrypoints across releases (e.g., top-level
    ``ssgsea``/``run_ssgsea`` or ``single_sample_gsea.ssgsea.ssgsea``). We scan
    common modules/attributes before falling back to the pure-Python scorer.

    Returns a scored DataFrame or ``None`` when the backend is missing or fails.
    """

    module_candidates = [
        'single_sample_gsea',
        'single_sample_gsea.ssgsea',
        'single_sample_gsea.ssgsea_py',
        'single_sample_gsea.core',
    ]
    mod = None
    for name in module_candidates:  # pragma: no cover - optional dependency
        try:
            mod = __import__(name, fromlist=['*'])
            break
        except Exception as e:
            if default_debug:
                print(f"DEBUG: import {name} failed: {e}")
            continue
    if mod is None:
        return None

    eset = eset.copy()
    eset.index = eset.index.astype(str)
    eset.columns = eset.columns.astype(str)
    gene_sets = {k: list({g for g in v if g in eset.index}) for k, v in sigs.items() if len(v) >= min_size}
    gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= min_size}
    if not gene_sets:
        return None

    fn = None
    attr_candidates = (
        'ssgsea', 'ssGSEA', 'run_ssgsea', 'single_sample_gsea', 'compute_ssgsea',
        'ssgsea_py', 'ssgsea_cpu', 'ssgsea_gpu'
    )
    for attr in attr_candidates:
        cand = getattr(mod, attr, None)
        if callable(cand):
            fn = cand
            break
    if fn is None and hasattr(mod, '__dict__'):
        # Some versions wrap the function inside the module attributes
        for obj in mod.__dict__.values():
            if callable(obj) and getattr(obj, '__name__', '').lower() in {a.lower() for a in attr_candidates}:
                fn = obj
                break
    if fn is None:
        if default_debug:
            print("DEBUG: single_sample_gsea module had no callable ssGSEA entrypoint")
        return None

    def _build_kwargs():
        try:
            sig = inspect.signature(fn)
            params = sig.parameters
        except Exception:
            params = {}

        kwargs = {}
        for name in params:
            lname = name.lower()
            if lname in {'data', 'expression_data', 'expression', 'expression_df', 'exp_df', 'expr'}:
                kwargs[name] = eset
            elif lname in {'gene_sets', 'gs', 'gene_set', 'gene_sets_dict', 'genesets'}:
                kwargs[name] = gene_sets
            elif lname in {'min_size', 'min_sz', 'min_size_expressed'}:
                kwargs[name] = min_size
            elif parallel_size and parallel_size > 1 and lname in {'n_jobs', 'threads', 'n_threads', 'n_workers', 'workers'}:
                kwargs[name] = int(parallel_size)
            elif lname in {'sample_norm', 'normalize', 'norm'}:
                kwargs[name] = True
        if not kwargs:
            kwargs = {'expression_data': eset, 'gene_sets': gene_sets}
        return kwargs

    kwargs = _build_kwargs()
    try:  # pragma: no cover - optional dependency
        res = fn(**kwargs)
    except Exception as e:
        if default_debug:
            print(f"DEBUG: single_sample_gsea call failed: {e}")
        return None

    if isinstance(res, pd.DataFrame):
        df = res.copy()
    elif isinstance(res, dict):
        df = pd.DataFrame(res)
    elif isinstance(res, np.ndarray):
        df = pd.DataFrame(res)
    else:
        if default_debug:
            print(f"DEBUG: single_sample_gsea returned unsupported type: {type(res)}")
        return None

    # Align output to samples × signatures
    if set(df.index) == set(eset.columns):
        mat = df.loc[eset.columns]
    elif set(df.columns) == set(eset.columns):
        mat = df[eset.columns].T
    elif set(df.index) == set(gene_sets.keys()):
        mat = df.loc[gene_sets.keys()].T
    elif set(df.columns) == set(gene_sets.keys()):
        mat = df[gene_sets.keys()]
        if set(mat.index) != set(eset.columns):
            mat = mat.T
    else:
        if default_debug:
            print("DEBUG: Unable to align single_sample_gsea output; falling back")
        return None

    mat = mat.copy()
    for missing in set(gene_sets.keys()) - set(mat.columns):
        mat[missing] = 0.0
    mat = mat.loc[:, list(gene_sets.keys())]
    mat.reset_index(inplace=True)
    mat.rename(columns={'index': 'ID'}, inplace=True)
    return mat


def _ssgsea_pure_python(eset: pd.DataFrame, sigs: dict, min_size: int, parallel_size: int):
    """Pure Python ssGSEA (no gseapy/rpy2/decoupler/numpy dependence in the scorer)."""

    eset = eset.copy()
    eset.index = eset.index.astype(str)
    eset.columns = eset.columns.astype(str)

    finite_rows = np.isfinite(eset).all(axis=1)
    eset = eset.loc[finite_rows]

    sigs = {k: [g for g in v if g in eset.index] for k, v in sigs.items()}
    sigs = {k: v for k, v in sigs.items() if len(v) >= min_size}

    samples = list(eset.columns)
    sig_names = list(sigs.keys())
    gene_index = list(eset.index)
    sig_sets = {k: set(v) for k, v in sigs.items()}

    def _rankdata_desc(values):
        """Average ranks in descending order (higher expression gets rank 1)."""
        pairs = sorted(enumerate(values), key=lambda x: (-x[1], x[0]))
        ranks = [0.0] * len(values)
        i = 0
        while i < len(pairs):
            j = i + 1
            while j < len(pairs) and pairs[j][1] == pairs[i][1]:
                j += 1
            avg_rank = (i + j - 1) / 2.0 + 1.0
            for k in range(i, j):
                ranks[pairs[k][0]] = avg_rank
            i = j
        order = [idx for idx, _ in pairs]
        return ranks, order

    def _score_sample(sample):
        expr = eset[sample].tolist()
        ranks, order = _rankdata_desc(expr)

        ordered_genes = [gene_index[i] for i in order]
        ordered_expr = [expr[i] for i in order]
        ordered_ranks = [ranks[i] for i in order]
        total = len(ordered_genes)

        weights_expr = [pow(abs(v), 0.75) if np.isfinite(v) else 0.0 for v in ordered_expr]
        use_rank_weights = sum(weights_expr) == 0
        base_weights = weights_expr if not use_rank_weights else [pow(r, 0.25) for r in ordered_ranks]

        sample_scores = {}
        for name, genes in sig_sets.items():
            hit_mask = [g in genes for g in ordered_genes]
            Nh = sum(1 for h in hit_mask if h)
            if Nh < min_size or Nh == 0 or Nh == total:
                sample_scores[name] = 0.0
                continue

            miss = total - Nh
            hit_weights = [bw if hit else 0.0 for bw, hit in zip(base_weights, hit_mask)]
            norm_hits = sum(hit_weights)
            if norm_hits == 0:
                sample_scores[name] = 0.0
                continue

            hit_cdf = []
            miss_cdf = []
            h_acc = 0.0
            m_acc = 0.0
            for hw, hit in zip(hit_weights, hit_mask):
                if hit:
                    h_acc += hw / norm_hits
                else:
                    m_acc += 1.0 / miss
                hit_cdf.append(h_acc)
                miss_cdf.append(m_acc)

            running = [h - m for h, m in zip(hit_cdf, miss_cdf)]
            area = 0.0
            for i in range(1, total):
                area += (running[i - 1] + running[i]) / 2.0
            es = area / float(total)
            sample_scores[name] = es

        vals = list(sample_scores.values())
        norm_const = float(np.mean(np.abs(vals))) if vals else 0.0
        if norm_const > 0:
            sample_scores = {k: v / norm_const for k, v in sample_scores.items()}

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


def _ssgsea(eset: pd.DataFrame, sigs: dict, min_size: int, parallel_size: int):
    """Run ssGSEA preferring the single-sample-gsea backend, else fallback to pure Python."""
    mat = _run_single_sample_gsea_backend(eset, sigs, min_size, parallel_size)
    if mat is not None:
        print("Running ssGSEA via single-sample-gsea backend...")
        return mat

    print("Running ssGSEA via pure Python fallback (single-sample-gsea unavailable)...")
    return _ssgsea_pure_python(eset, sigs, min_size, parallel_size)

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

    es = _ssgsea(eset2, sigs, min_size, parallel_size)

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

    es = _ssgsea(eset2, filtered_sigs, mini_gene_count, parallel_size)

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