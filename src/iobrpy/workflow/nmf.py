#!/usr/bin/env python3
"""
nmf_cluster.py

Command-line script: perform NMF-based clustering on a samples x cell-types matrix and automatically
select the best k within a given range based on silhouette score — excluding k=2.

Input: CSV/TSV with the first column as sample names (index).
Output: saves `clusters.csv`, `pca_plot.png` and `top_features_per_cluster.csv` (top features per cluster).

Example:
python nmf_cluster.py --input matrix.csv --output outdir --kmin 2 --kmax 8 --log1p --features 2-10 --max-iter 5000

Notes:
- --features selects which columns (cell types) from the input to use. It uses 1-based indexing and
  includes both endpoints, e.g. "2-4" picks the 2nd, 3rd and 4th data columns (not the sample column).
- --max-iter controls the maximum number of iterations passed to sklearn.decomposition.NMF (default 1000).
- Clustering method: NMF-only. For each k, samples are assigned to the component with largest W value (argmax).
  k=2 is skipped when searching for the best k.
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF, PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
from iobrpy.utils.print_colorful_message import print_colorful_message

def read_matrix(path):
    # Try comma first, then tab
    try:
        df = pd.read_csv(path, index_col=0)
        return df
    except Exception:
        df = pd.read_csv(path, sep='    ', index_col=0)
        return df


def parse_features_arg(feat_str, n_cols):
    """Parse --features argument. Support 'm-n' or 'm:n' or single 'm' (1-based, inclusive).
       Returns start_idx (inclusive, 0-based), end_idx (exclusive)
    """
    if feat_str is None:
        return 0, n_cols
    feat_str = feat_str.strip()
    sep = None
    if '-' in feat_str:
        sep = '-'
    elif ':' in feat_str:
        sep = ':'

    try:
        if sep:
            parts = feat_str.split(sep)
            if len(parts) != 2:
                raise ValueError
            start = int(parts[0])
            end = int(parts[1])
        else:
            start = int(feat_str)
            end = start
    except Exception:
        raise argparse.ArgumentTypeError("--features format error: use 'm-n' or 'm:n' or single number (1-based)")

    if start < 1 or end < 1 or start > end:
        raise argparse.ArgumentTypeError(f"--features invalid range: {feat_str}")
    start_idx = start - 1
    end_idx = end

    if start_idx >= n_cols:
        raise argparse.ArgumentTypeError(f"--features start column ({start}) exceeds number of columns ({n_cols})")
    if end_idx > n_cols:
        end_idx = n_cols
    return start_idx, end_idx


def ensure_nonneg(X, shift=None):
    minv = X.min()
    if minv < 0:
        if shift is None:
            raise ValueError(f"Input contains negative values (min={minv:.4g}); NMF requires non-negative data. Use --shift to add a constant or preprocess the data.")
        else:
            X = X + abs(minv) + shift
    return X


def fit_nmf(X, n_components, random_state=42, max_iter=1000):
    model = NMF(n_components=n_components, init='nndsvda', random_state=random_state, max_iter=max_iter)
    W = model.fit_transform(X)
    H = model.components_
    rec_err = getattr(model, 'reconstruction_err_', None)
    return model, W, H, rec_err


def try_k_range_nmf_argmax(X, kmin, kmax, random_state=42, max_iter=1000, skip_k_2=False):
    """For each k (optionally skipping k==2 when skip_k_2=True): fit NMF(n_components=k),
    assign labels = argmax(W), compute silhouette score on W.
    Return dict: k -> (silhouette, reconstruction_err)
    """
    n_samples = X.shape[0]
    results = {}
    for k in range(kmin, kmax + 1):
        if k == 2 and skip_k_2:
            print(f"skip k={k} (skip_k_2=True)")
            continue
        if k < 2:
            print(f"skip k={k} (k must be >=2)")
            continue
        if k >= n_samples:
            print(f"skip k={k} (k must be less than number of samples {n_samples})")
            continue
        try:
            model, W, _, rec_err = fit_nmf(X, n_components=k, random_state=random_state, max_iter=max_iter)
            labels = W.argmax(axis=1)
            unique_count = len(np.unique(labels))
            if unique_count <= 1:
                sil = float('-inf')
                print(f"k={k}: only {unique_count} cluster(s) found by argmax -> silhouette not defined")
            else:
                sil = silhouette_score(W, labels)
                if rec_err is None:
                    print(f"k={k}: silhouette={sil:.4f}")
                else:
                    print(f"k={k}: silhouette={sil:.4f}, rec_err={rec_err:.4f}")
            results[k] = (sil, rec_err)
        except Exception as e:
            print(f"k={k} training failed: {e}")
    return results


def save_outputs(outdir, sample_index, cluster_labels, pca_coords=None):
    os.makedirs(outdir, exist_ok=True)
    # Save clusters.csv (sample as index)
    labels_df = pd.DataFrame({'sample': sample_index, 'cluster': cluster_labels})
    labels_df.set_index('sample', inplace=True)
    labels_df.to_csv(os.path.join(outdir, 'clusters.csv'))

    # Save PCA plot
    if pca_coords is not None:
        cat = pd.Categorical(cluster_labels)
        codes = cat.codes
        unique_labels = list(cat.categories)

        fig, ax = plt.subplots(figsize=(7, 5))
        scatter = ax.scatter(pca_coords[:, 0], pca_coords[:, 1], c=codes, s=40)
        for i, s in enumerate(sample_index):
            ax.text(pca_coords[i, 0], pca_coords[i, 1], s, fontsize=8, alpha=0.7)
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title('Samples in NMF-component space (PCA of W)')
        # Create legend mapping colors to cluster names
        handles = []
        for code, name in enumerate(unique_labels):
            handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='k', label=name))
        ax.legend(handles=handles, title='cluster')
        plt.tight_layout()
        fig_path = os.path.join(outdir, 'pca_plot.png')
        fig.savefig(fig_path, dpi=150)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Cluster samples using NMF (argmax on W) and auto-select best k by silhouette score (k!=2).')
    parser.add_argument('--input', '-i', required=True, help='Input matrix file (CSV or TSV). First column should be sample names (index).')
    parser.add_argument('--output', '-o', required=True, help='Output directory where results will be saved.')
    parser.add_argument('--kmin', type=int, default=2, help='Minimum k (inclusive). Default: 2')
    parser.add_argument('--kmax', type=int, default=8, help='Maximum k (inclusive). Default: 8')
    parser.add_argument('--features', type=str, default=None, help="Columns (cell types) to use, e.g. '2-10' or '1:5'. 1-based inclusive.")
    parser.add_argument('--log1p', action='store_true', help='Apply log1p transform to the input (useful for counts).')
    parser.add_argument('--normalize', action='store_true', help='Apply L1 row normalization (each sample sums to 1).')
    parser.add_argument('--shift', type=float, default=None, help='If input contains negatives, add a constant shift to make values non-negative.')
    parser.add_argument('--random-state', type=int, default=42)
    parser.add_argument('--max-iter', type=int, default=1000, help='Maximum iterations for NMF (default: 1000)')
    parser.add_argument('--skip_k_2', action='store_true',help='Skip k=2 when searching for the best k (default: do not skip)')
    args = parser.parse_args()

    # 1) Read input
    print('Reading matrix:', args.input)
    df = read_matrix(args.input)
    print('Input shape (n_samples x total_columns):', df.shape)

    # Parse features
    n_cols = df.shape[1]
    try:
        start_idx, end_idx = parse_features_arg(args.features, n_cols)
    except Exception as e:
        print('ERROR parsing --features:', e)
        sys.exit(1)

    if not (start_idx == 0 and end_idx == n_cols):
        selected_df = df.iloc[:, start_idx:end_idx]
        print(f"Using column range: {start_idx+1} to {end_idx} (total {selected_df.shape[1]} columns)")
    else:
        selected_df = df
        print('Using all columns from input')

    # 2) Preprocessing
    X = selected_df.values.astype(float)
    if args.log1p:
        print('Applying log1p transform')
        X = np.log1p(X)
    if args.normalize:
        print('Applying L1 row normalization')
        X = normalize(X, norm='l1', axis=1)

    # Ensure non-negative
    try:
        X = ensure_nonneg(X, shift=args.shift)
    except ValueError as e:
        print('ERROR:', e)
        sys.exit(1)

    n_samples = X.shape[0]
    if args.kmax >= n_samples:
        print(f"Warning: kmax ({args.kmax}) >= number of samples ({n_samples}). Setting kmax to {n_samples-1}.")
        args.kmax = n_samples - 1

    # 3) Try k range (NMF argmax clustering). k=2 explicitly skipped
    print('Trying k range:', args.kmin, '->', args.kmax)
    results = try_k_range_nmf_argmax(X, args.kmin, args.kmax, random_state=args.random_state, max_iter=args.max_iter, skip_k_2=args.skip_k_2)
    # Filter out ks that failed or gave -inf silhouette
    valid_results = {k: v for k, v in results.items() if v[0] is not None and v[0] != float('-inf')}
    if not valid_results:
        print('No successful results in the provided k range{}.'.format(
            ' (excluding k=2)' if args.skip_k_2 else ''
        ))
        sys.exit(1)

    # 4) Choose best k (max silhouette). If tie, prefer smaller k.
    best_k = max(sorted(valid_results.keys()), key=lambda kk: (valid_results[kk][0], -kk))
    best_sil, best_rec_err = valid_results[best_k]
    print(f"Selected best k = {best_k}, silhouette={best_sil:.4f}, rec_err={best_rec_err:.4f}")

    # 5) Retrain with best_k and produce clusters using argmax on W
    model, W, H, rec_err = fit_nmf(X, n_components=best_k, random_state=args.random_state, max_iter=args.max_iter)
    numeric_labels = W.argmax(axis=1)
    cluster_labels = [f"cluster{int(l)+1}" for l in numeric_labels]

    # 5.1) Save top features per cluster (use H matrix)
    try:
        feature_names = list(selected_df.columns)
        n_feats = len(feature_names)
        top_n = min(15, n_feats)
        rows = []
        for comp_idx in range(H.shape[0]):
            top_idx = np.argsort(-H[comp_idx])[:top_n]
            top_names = [feature_names[i] for i in top_idx]
            # pad if needed (not necessary here since top_n <= n_feats)
            rows.append(top_names)
        col_names = [f'top_{i+1}' for i in range(top_n)]
        top_df = pd.DataFrame(rows, index=[f'cluster{i+1}' for i in range(H.shape[0])], columns=col_names)
        top_df.index.name = 'cluster'
        top_df.to_csv(os.path.join(args.output, 'top_features_per_cluster.csv'))
    except Exception as e:
        print('Failed to save top features per cluster:', e)

    # 6) PCA for visualization (use W space)
    try:
        pca = PCA(n_components=2, random_state=args.random_state)
        W2 = pca.fit_transform(W)
    except Exception:
        W2 = None

    # 7) Save outputs (clusters.csv, pca_plot.png)
    save_outputs(args.output, selected_df.index.tolist(), cluster_labels, pca_coords=W2)

    print('Outputs saved to:', os.path.abspath(args.output))
    print('Files: clusters.csv, pca_plot.png , top_features_per_cluster.csv')

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
