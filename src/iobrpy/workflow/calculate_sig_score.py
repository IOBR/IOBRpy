#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import zscore, spearmanr
from importlib.resources import files
from joblib import Parallel, delayed

try:
    import gseapy as gp
except ImportError:
    gp = None

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

def preprocess_eset(eset, adjust_eset):
    """Log2-transform and drop Inf/zero‑sd genes like R's log2eset+check_eset."""
    if default_debug:
        print(f"DEBUG: Preprocess: shape={eset.shape}, adjust={adjust_eset}")
    if eset.shape[1] < 5000:
        eset = np.log2(eset + 1)
        # drop Inf rows
        finite = np.isfinite(eset).all(axis=1)
        eset = eset.loc[finite]
    if adjust_eset:
        # drop Inf and zero‑sd again
        finite = np.isfinite(eset).all(axis=1)
        eset = eset.loc[finite]
        eset = eset.loc[eset.std(axis=1) > 0]
    if default_debug:
        print(f"DEBUG: After preprocess: shape={eset.shape}")
    return eset

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
    if gp is None:
        raise ImportError("gseapy required for ssGSEA")
    # Preprocess like R
    eset2 = preprocess_eset(eset, adjust_eset)
    # First filter with original threshold
    sigs = filter_signatures(sig_dict, eset2, mini_gene_count)
    # Then enforce min_size >= 5
    min_size = max(mini_gene_count, 5)
    sigs = filter_signatures(sig_dict, eset2, min_size)

    # Run ssGSEA with R-aligned parameters
    print("Running ssGSEA (this may take a while)...")
    ss = gp.ssgsea(
        data=eset2,
        gene_sets=sigs,
        outdir=None,
        sample_norm_method='rank',    # rank-based kernel = Gaussian
        permutation_num=0,
        no_plot=True,
        threads=parallel_size,
        min_size=min_size,
        ssgsea_norm=True
    )

    # Pivot to samples × terms, keep Term as columns
    nes = ss.res2d.pivot(index='Term', columns='Name', values='NES').T.reset_index()
    nes.rename(columns={'Name': 'ID'}, inplace=True)

    if 'TMEscoreA_CIR' in nes.columns and 'TMEscoreB_CIR' in nes.columns:
        nes['TMEscore_CIR'] = nes['TMEscoreA_CIR'] - nes['TMEscoreB_CIR']
    if 'TMEscoreA_plus' in nes.columns and 'TMEscoreB_plus' in nes.columns:
        nes['TMEscore_plus'] = nes['TMEscoreA_plus'] - nes['TMEscoreB_plus']
    return nes

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

    if gp is None:
        raise ImportError("gseapy required for ssGSEA")

    eset2 = preprocess_eset(eset, adjust_eset)
    
    print("Running ssGSEA (this may take a while)...")
    ss = gp.ssgsea(
        data=eset2,
        gene_sets=filtered_sigs,
        outdir=None,
        sample_norm_method='rank',
        permutation_num=0,
        no_plot=True,
        threads=parallel_size,
        min_size=mini_gene_count,
        ssgsea_norm=True
    )
    nes = ss.res2d.pivot(index='Term', columns='Name', values='NES').T.reset_index()
    nes.rename(columns={'Name': 'ID'}, inplace=True)

    if 'TMEscoreA_CIR' in nes.columns and 'TMEscoreB_CIR' in nes.columns:
        nes['TMEscore_CIR'] = nes['TMEscoreA_CIR'] - nes['TMEscoreB_CIR']
    if 'TMEscoreA_plus' in nes.columns and 'TMEscoreB_plus' in nes.columns:
        nes['TMEscore_plus'] = nes['TMEscoreA_plus'] - nes['TMEscoreB_plus']
    s = nes.set_index('ID').add_suffix('_ssGSEA')

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