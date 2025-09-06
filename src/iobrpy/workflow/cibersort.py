import numpy as np
import pandas as pd
from sklearn.svm import NuSVR
from scipy.stats import pearsonr
from sklearn.preprocessing import quantile_transform
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse
from importlib.resources import files
import csv

def quantile_normalize(df):
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

def make_unique(index):
    counts = {}
    new_index = []
    for name in index:
        if name not in counts:
            counts[name] = 1
            new_index.append(name)
        else:
            counts[name] += 1
            new_index.append(f"{name}.{counts[name]}")
    return new_index

def core_alg(X, y, absolute=False, abs_method='sig.score', n_jobs=1):
    nu_values = [0.25, 0.5, 0.75]

    def train_model(nu_val):
        model = NuSVR(kernel='linear', nu=nu_val)
        model.fit(X, y)
        return model

    models = Parallel(n_jobs=n_jobs)(delayed(train_model)(nu) for nu in nu_values)

    rmses, corrs = [], []

    for model in models:
        coef = model.dual_coef_
        SV = model.support_vectors_
        weights = np.dot(coef, SV)
        weights[weights < 0] = 0
        w_norm = weights / np.sum(weights)
        k = np.sum(X * w_norm, axis=1)
        rmse = np.sqrt(np.mean((k - y) ** 2))
        corr = pearsonr(k, y)[0]
        rmses.append(rmse)
        corrs.append(corr)

    best_idx = np.argmin(rmses)
    best_model = models[best_idx]
    coef = best_model.dual_coef_
    SV = best_model.support_vectors_
    q = np.dot(coef, SV)
    q[q < 0] = 0

    if not absolute:
        w = q / np.sum(q)
    else:
        w = q

    return {
        "w": w.flatten(),
        "mix_rmse": rmses[best_idx],
        "mix_r": corrs[best_idx]
    }

def do_perm(perm, X, Y, absolute, abs_method, n_jobs=1):
    Y_flat = Y.values.flatten()
    dist = []

    for _ in tqdm(range(perm), desc="Permutation Test Progress", leave=False):
        yr = np.random.choice(Y_flat, size=X.shape[0], replace=True)
        yr = (yr - np.mean(yr)) / np.std(yr)
        result = core_alg(X, yr, absolute, abs_method, n_jobs=n_jobs)
        dist.append(result["mix_r"])

    return np.sort(dist)

def cibersort(input_path, perm=100, QN=True, absolute=False, abs_method='sig.score', n_jobs=1):
    resource_pkg = 'iobrpy.resources'
    lm22_path = files(resource_pkg).joinpath('lm22.txt')
    sig_df = pd.read_csv(lm22_path, sep=' ', index_col=0)
    mixture_df = pd.read_csv(
            input_path,
            sep=None,
            engine='python',
            index_col=0
        )

    mixture_df.index = make_unique(mixture_df.index)

    common_genes = sig_df.index.intersection(mixture_df.index)
    if len(common_genes) == 0:
        raise ValueError("No overlapping genes found between signature and mixture matrices.")

    sig_df = sig_df.loc[common_genes].sort_index()
    mixture_df = mixture_df.loc[common_genes].sort_index()

    X = sig_df.to_numpy()
    Y = mixture_df.to_numpy()

    if np.max(Y) < 50:
        Y = 2 ** Y

    if QN:
        Y = quantile_normalize(pd.DataFrame(Y, index=common_genes, columns=mixture_df.columns)).values

    Y_orig = Y.copy()
    Y_median = max(np.median(Y_orig), 1)

    X = (X - np.mean(X)) / np.std(X)

    nulldist = do_perm(perm, X, pd.DataFrame(Y), absolute, abs_method, n_jobs=n_jobs) if perm > 0 else None

    results = []
    for i in tqdm(range(Y.shape[1]), desc="Mixture Processing Progress", leave=False):
        y = Y[:, i]
        y = (y - np.mean(y)) / np.std(y)
        result = core_alg(X, y, absolute, abs_method, n_jobs=n_jobs)
        w = result["w"]

        pval = 9999
        if nulldist is not None:
            pval = 1 - (np.searchsorted(nulldist, result["mix_r"]) / len(nulldist))

        row = list(w) + [pval, result["mix_r"], result["mix_rmse"]]
        if absolute:
            row.append(np.sum(w))
        results.append(row)

    colnames = list(sig_df.columns) + ['P-value', 'Correlation', 'RMSE']
    if absolute:
        safe_method = abs_method.replace('.', '_')
        colnames.append(f'Absolute_score_({safe_method})')

    result_df = pd.DataFrame(results, columns=colnames, index=mixture_df.columns)
    return result_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, dest="input_path", help="Path to mixture file")
    parser.add_argument("--perm", type=int, default=100, help="Number of permutations")
    parser.add_argument("--QN", type=lambda x: x.lower() == "true", default=True, help="Quantile normalization (True/False)")
    parser.add_argument("--absolute", type=lambda x: x.lower() == "true", default=False, help="Absolute mode (True/False)")
    parser.add_argument("--abs_method", default="sig.score", choices=['sig.score', 'no.sumto1'], help="Absolute scoring method")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads (default=1)")
    parser.add_argument("--output", required=True, help="Output file path")

    args = parser.parse_args()

    result_df = cibersort(
        input_path=args.input_path,
        perm=args.perm,
        QN=args.QN,
        absolute=args.absolute,
        abs_method=args.abs_method,
        n_jobs=args.threads
    )

    if result_df is not None:
        result_df.columns = [col + '_CIBERSORT' for col in result_df.columns]
        result_df.index.name = 'ID'
        delim = ',' if args.output.lower().endswith('.csv') else '\t'
        result_df.to_csv(args.output, sep=delim)
        print(f"[Done] Output written to {args.output} ")
    else:
        print("[Error] No results generated")