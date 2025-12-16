import argparse
import os
import re
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from iobrpy.utils.print_colorful_message import print_colorful_message
from iobrpy.workflow.count2tpm import feature_manipulation

def parse_args():
    p = argparse.ArgumentParser(description="TME clustering with Hartigan-Wong kmeans and KL index matching R NbClust")
    p.add_argument('-i', '--input',    required=True,
                   help="Input file path (CSV/TSV/TXT)")
    p.add_argument('-o', '--output',   required=True,
                   help="Output file path (CSV/TSV/TXT)")
    p.add_argument('--features', default=None,
                   help="Feature columns to use, e.g. '2:23' (1-based)")
    p.add_argument('--pattern',  default=None,
                   help="Regex to select feature columns by name")
    p.add_argument('--id',       default=None,
                   help="Column name for sample IDs (default: first column)")
    p.add_argument('--scale',    action='store_true', default=True,
                   help="Enable z-score scaling (default: True)")
    p.add_argument('--no-scale', action='store_false', dest='scale',
                   help="Disable scaling")
    p.add_argument('--min_nc',   type=int, default=2,
                   help="Min number of clusters (k for kmeans)")
    p.add_argument('--max_nc',   type=int, default=6,
                   help="Max number of clusters (ignored for kmeans)")
    p.add_argument('--nstart',   type=int, default=10,
                   help="Number of random starts (NbClust uses 10 by default)")
    p.add_argument('--max_iter', type=int, default=10,
                   help="Maximum iterations for Hartigan–Wong k-means (default: 10)")
    p.add_argument('--tol', type=float, default=1e-4,
                   help="Convergence tolerance (default: 1e-4)")
    p.add_argument('--seed', type=int, default=123,
                   help="Seed for reproducible initialization (default: 123, matches R NbClust)")
    p.add_argument('--print_result', action='store_true',
                   help="Print intermediate info")
    p.add_argument('--input_sep',  default=None,
                   help="Field separator for input (auto-detect if not set)")
    p.add_argument('--output_sep', default=None,
                   help="Field separator for output (auto-detect if not set)")
    return p.parse_args()

def detect_sep(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in ['.tsv', '.txt']: return '\t'
    if ext == '.csv': return ','
    return ','

def hartigan_wong(data, k, max_iter=10, tol=1e-4, rng=None):
    n, p = data.shape
    # initialize centers
    rng = rng or np.random
    indices = rng.choice(n, size=k, replace=False)
    centers = data[indices].astype(float).copy()
    # initial assignment
    labels = np.argmin(((data[:, None] - centers[None, :])**2).sum(axis=2), axis=1)
    # handle empty clusters: reassign farthest point
    def fix_empty(labels, centers):
        unique, counts = np.unique(labels, return_counts=True)
        missing = set(range(k)) - set(unique)
        for j in missing:
            # compute distances from each point to its assigned center
            dists = ((data - centers[labels])**2).sum(axis=1)
            idx = np.argmax(dists)
            labels[idx] = j
        return labels
    labels = fix_empty(labels, centers)
    for _ in range(max_iter):
        sums = np.zeros((k, p)); counts = np.zeros(k, int)
        for i, lbl in enumerate(labels):
            sums[lbl] += data[i]; counts[lbl] += 1
        # update centers
        for j in range(k):
            if counts[j] > 0:
                centers[j] = sums[j] / counts[j]
        moved = False
        # reassign points with Hartigan criterion
        for i in range(n):
            lbl = labels[i]; x = data[i]
            curr_cost = ((x - centers[lbl])**2).sum()
            best_delta, best_lbl = tol, lbl
            for j in range(k):
                if j == lbl: continue
                if counts[j] == 0:
                    new_cost = 0
                else:
                    new_center = (centers[j]*counts[j] + x) / (counts[j] + 1)
                    new_cost = ((x - new_center)**2).sum()
                delta = curr_cost - new_cost
                if delta > best_delta:
                    best_delta, best_lbl = delta, j
            if best_lbl != lbl:
                moved = True
                counts[lbl] -= 1; sums[lbl] -= x
                counts[best_lbl] += 1; sums[best_lbl] += x
                labels[i] = best_lbl
                centers[lbl] = sums[lbl]/counts[lbl] if counts[lbl] > 0 else centers[lbl]
                centers[best_lbl] = sums[best_lbl]/counts[best_lbl]
        # after moves, fix any empty
        labels = fix_empty(labels, centers)
        if not moved:
            break
    return labels, centers

def compute_withinss(data, labels):
    W = 0.0
    for lbl in np.unique(labels):
        pts = data[labels == lbl]; center = pts.mean(axis=0)
        W += ((pts - center)**2).sum()
    return W

def best_hartigan_run(data, k, nstart, max_iter, tol, rng):
    """Mimic NbClust's multiple random starts and keep the best withinss run."""
    best_labels = None
    best_w = np.inf
    best_centers = None
    for _ in range(max(nstart, 1)):
        labels, centers = hartigan_wong(data, k, max_iter=max_iter, tol=tol, rng=rng)
        w = compute_withinss(data, labels)
        if w < best_w:
            best_w = w
            best_labels = labels
            best_centers = centers
    return best_labels, best_centers, best_w

def main():
    args = parse_args()
    sep_in = args.input_sep or detect_sep(args.input)
    sep_out = args.output_sep or detect_sep(args.output)
    df = pd.read_csv(args.input, sep=sep_in)
    if args.id and args.id in df.columns:
        ids = df[args.id].astype(str)
        df = df.drop(columns=[args.id])
    else:
        ids = df.iloc[:, 0].astype(str)
        df = df.drop(df.columns[0], axis=1)

    if args.features:
        m = re.match(r'^(\d+):(\d+)$', args.features)
        cols = df.columns[int(m.group(1)) - 1:int(m.group(2))]
    elif args.pattern:
        cols = [c for c in df.columns if re.search(args.pattern, c)]
    else:
        cols = df.columns

    data = df[cols].apply(pd.to_numeric, errors='coerce')
    if args.scale:
        data = (data - data.mean()) / data.std(ddof=1)
        valid_features = feature_manipulation(data, feature=list(data.columns), print_result=args.print_result)
        data = data[valid_features]
    if data.empty:
        raise ValueError("No valid features remain after filtering. Check input data and parameters.")

    X = data.values
    n, p = X.shape

    rng = np.random.RandomState(args.seed)
    kl_scores = {}
    labels_by_k = {}
    withinss_by_k = {}

    main_k_values = list(range(args.min_nc, args.max_nc + 1))

    for k in tqdm(main_k_values, desc="kmeans runs"):
        labels_k, _, w_k = best_hartigan_run(X, k, args.nstart, args.max_iter, args.tol, rng)
        labels_by_k[k] = labels_k
        withinss_by_k[k] = w_k

    # Ensure we have neighboring W values so KL can be computed for min_nc and max_nc as well
    needed_neighbors = set()
    for k in main_k_values:
        if k - 1 >= 1:
            needed_neighbors.add(k - 1)
        needed_neighbors.add(k + 1)
    neighbor_ks = sorted(needed_neighbors - withinss_by_k.keys())
    for k in neighbor_ks:
        rng_neighbor = np.random.RandomState(args.seed + 9973 * k)
        _, _, w_k = best_hartigan_run(X, k, args.nstart, args.max_iter, args.tol, rng_neighbor)
        withinss_by_k[k] = w_k

    # Then compute KL index using neighboring k values (including min_nc and max_nc when possible)
    for k in main_k_values:
        prev_k, next_k = k - 1, k + 1
        if prev_k < 1 or next_k not in withinss_by_k:
            continue
        W_prev = withinss_by_k[prev_k]
        W_k    = withinss_by_k[k]
        W_next = withinss_by_k[next_k]
        num = abs((k-1)**(2/p) * W_prev - k**(2/p) * W_k)
        denom = abs(k**(2/p) * W_k - (k+1)**(2/p) * W_next)
        kl_scores[k] = num / (denom if denom > 0 else 1e-8)

    if args.print_result and kl_scores:
        for k in sorted(kl_scores):
            print(f"k={k} KL={kl_scores[k]:.4f}")

    # Fallback to min_nc if KL cannot be computed
    best_k = max(kl_scores, key=kl_scores.get) if kl_scores else args.min_nc
    if args.print_result:
        print(f"Best k by KL: {best_k}")

    labels = labels_by_k[best_k]
    clusters = [f"TME{l+1}" for l in labels]

    out = pd.DataFrame({'ID': ids, 'cluster': clusters})
    out = pd.concat([out, data.reset_index(drop=True)], axis=1)
    out.to_csv(args.output, sep=sep_out, index=False)
    # Print absolute output file path and its directory for clarity
    abs_out = os.path.abspath(args.output)
    if args.print_result:
        counts = out['cluster'].value_counts().reset_index()
        counts.columns = ['cluster', 'count']
        print(counts.to_string(index=False))
    print(f"tme_cluster results saved to：{abs_out}")
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

if __name__=='__main__':
    main()