import argparse
import os
import re
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

def parse_args():
    p = argparse.ArgumentParser(description="TME clustering with Hartigan-Wong kmeans and KL index, matching R seed per kmeans call")
    p.add_argument('--input',    required=True,
                   help="Input file path (CSV/TSV/TXT)")
    p.add_argument('--output',   required=True,
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
    p.add_argument('--max_iter', type=int, default=10,
                   help="Maximum number of iterations for the Hartigan–Wong k-means algorithm (default: 10)")
    p.add_argument('--tol', type=float, default=1e-4,
                   help="Convergence tolerance for center updates in the clustering loop (default: 1e-4)")
    p.add_argument('--print_result', action='store_true', default=False,
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


def hartigan_wong(data, k, max_iter=10, tol=1e-4):
    n, p = data.shape
    indices = np.random.choice(n, k, replace=False)
    centers = data[indices].copy()
    labels = np.argmin(((data[:, None] - centers[None, :])**2).sum(axis=2), axis=1)
    for _ in range(max_iter):
        moved = False
        sums = np.zeros((k, p)); counts = np.zeros(k, int)
        for i, lbl in enumerate(labels):
            sums[lbl] += data[i]; counts[lbl] += 1
        for j in range(k):
            if counts[j] > 0:
                centers[j] = sums[j] / counts[j]
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
                counts[lbl] -= 1; sums[lbl] -= x
                counts[best_lbl] += 1; sums[best_lbl] += x
                centers[lbl] = sums[lbl]/counts[lbl] if counts[lbl] > 0 else centers[lbl]
                centers[best_lbl] = sums[best_lbl]/counts[best_lbl]
                labels[i] = best_lbl; moved = True
        if not moved: break
    return labels, centers


def compute_withinss(data, labels):
    W = 0.0
    for lbl in np.unique(labels):
        pts = data[labels == lbl]; center = pts.mean(axis=0)
        W += ((pts - center)**2).sum()
    return W


def kl_index(data, lbl_k, lbl_k1, lbl_k2):
    p = data.shape[1]; k = len(np.unique(lbl_k))
    Wk = compute_withinss(data, lbl_k)
    Wk1 = compute_withinss(data, lbl_k1)
    Wk2 = compute_withinss(data, lbl_k2)
    num = abs((k-1)**(2/p)*Wk - k**(2/p)*Wk1)
    denom = abs(k**(2/p)*Wk1 - (k+1)**(2/p)*Wk2)
    return num / (denom if denom > 0 else 1e-8)


def main():
    args = parse_args()
    sep_in = args.input_sep or detect_sep(args.input)
    sep_out = args.output_sep or detect_sep(args.output)
    df = pd.read_csv(args.input, sep=sep_in)
    if args.id and args.id in df.columns:
        ids = df[args.id].astype(str); df = df.drop(columns=[args.id])
    else:
        ids = df.iloc[:, 0].astype(str); df = df.drop(df.columns[0], axis=1)
    if args.features:
        m = re.match(r'^(\d+):(\d+)$', args.features)
        cols = df.columns[int(m.group(1)) - 1:int(m.group(2))]
    elif args.pattern:
        cols = [c for c in df.columns if re.search(args.pattern, c)]
    else:
        cols = df.columns
    data = df[cols].apply(pd.to_numeric, errors='coerce').dropna(axis=1)
    data = data.loc[:, data.std(axis=0) > 0]
    if args.scale:
        data[data.columns] = StandardScaler().fit_transform(data)
    X = data.values; n, p = X.shape

    kl_scores = {}
    for k in tqdm(range(args.min_nc, args.max_nc + 1), desc="KL scoring"):
        np.random.seed(1)
        lbl_k, _ = hartigan_wong(X, k, max_iter=args.max_iter, tol=args.tol)
        np.random.seed(1)
        lbl_k1, _ = hartigan_wong(X, k + 1, max_iter=args.max_iter, tol=args.tol)
        np.random.seed(1)
        lbl_k2, _ = hartigan_wong(X, k + 2, max_iter=args.max_iter, tol=args.tol)
        kl_scores[k] = kl_index(X, lbl_k, lbl_k1, lbl_k2)
        if args.print_result:
            print(f"k={k} KL={kl_scores[k]:.4f}")
    best_k = max(kl_scores, key=kl_scores.get)
    if args.print_result:
        print(f"Best k by KL: {best_k}")

    np.random.seed(1)
    labels, centers = tqdm(hartigan_wong(X, best_k, max_iter=args.max_iter, tol=args.tol),desc=f"Final clustering (k={best_k})")

    sums = centers.sum(axis=1); order = np.argsort(sums)
    mapping = {old: new for new, old in enumerate(order)}
    labels = np.array([mapping[l] for l in labels])
    clusters = [f"TME{l+1}" for l in labels]

    out = pd.DataFrame({'ID': ids, 'cluster': clusters})
    out = pd.concat([out, data.reset_index(drop=True)], axis=1)
    out.to_csv(args.output, sep=sep_out, index=False)
    if args.print_result:
        print(out['cluster'].value_counts())

if __name__=='__main__': 
    main()
