import numpy as np
import pandas as pd
from sklearn.svm import NuSVR
from importlib.resources import files
from joblib import Parallel, delayed
import joblib  # needed for custom tqdm_joblib fallback
import argparse

# -------- progress bar (tqdm) with robust fallbacks (Py3.9+) --------
# We prefer tqdm.auto.tqdm; if tqdm.contrib.tqdm_joblib is missing, we provide a
# small compatible shim that works with joblib's BatchCompletionCallBack.
try:
    from tqdm.auto import tqdm  # UI-friendly tqdm
except Exception:
    # Minimal no-op fallback if tqdm is not installed
    class _NoTqdm:
        def __init__(self, *a, **k): pass
        def update(self, *a, **k): pass
        def close(self): pass
        def __enter__(self): return self
        def __exit__(self, exc_type, exc, tb): return False
    def tqdm(*a, **k): return _NoTqdm()

# Try official contrib first; if unavailable, use a local shim.
try:
    from tqdm.contrib import tqdm_joblib  # may be missing on older tqdm
except Exception:
    class tqdm_joblib:
        """
        A lightweight drop-in replacement of tqdm.contrib.tqdm_joblib.
        It updates `tqdm_bar` every time a joblib batch completes.
        """
        def __init__(self, tqdm_bar):
            self.tqdm_bar = tqdm_bar
            self._orig = None

        def __enter__(self):
            # Save original callback
            self._orig = joblib.parallel.BatchCompletionCallBack
            tqdm_bar = self.tqdm_bar

            class TqdmBatchCompletionCallback(self._orig):
                def __call__(self, *args, **kwargs):
                    try:
                        tqdm_bar.update(n=self.batch_size)
                    finally:
                        return super().__call__(*args, **kwargs)

            joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
            return self

        def __exit__(self, exc_type, exc, tb):
            # Restore and close the bar
            joblib.parallel.BatchCompletionCallBack = self._orig
            self.tqdm_bar.close()
            return False

# -------------------- small, fast utilities --------------------
def zscore1d(a: np.ndarray) -> np.ndarray:
    """
    Safe 1D z-score for float32 arrays. Returns zeros if variance is zero.
    """
    a = a.astype(np.float32, copy=False)
    m = a.mean(dtype=np.float64)
    v = a.var(dtype=np.float64)
    if v == 0.0:
        return np.zeros_like(a, dtype=np.float32)
    return (a - m) / np.sqrt(v)

def corr_pearson_fast(a: np.ndarray, b: np.ndarray) -> float:
    """
    Fast Pearson correlation for float32 arrays using a single dot product.
    """
    a = a.astype(np.float32, copy=False)
    b = b.astype(np.float32, copy=False)
    am = a.mean(dtype=np.float64)
    bm = b.mean(dtype=np.float64)
    av = a - am
    bv = b - bm
    denom = np.sqrt((av * av).sum(dtype=np.float64) * (bv * bv).sum(dtype=np.float64))
    if denom == 0.0:
        return 0.0
    return float((av * bv).sum(dtype=np.float64) / denom)

def quantile_normalize_fast(Y: np.ndarray) -> np.ndarray:
    """
    Standard quantile normalization (columns = samples).
    Input: G x N (genes x samples). Output has the same shape.
    """
    order = np.argsort(Y, axis=0)                                      # G x N
    sorted_Y = np.take_along_axis(Y, order, axis=0)                    # G x N
    mean_sorted = sorted_Y.mean(axis=1, keepdims=True)                 # G x 1
    inv_order = np.empty_like(order)
    np.put_along_axis(inv_order, order, np.arange(Y.shape[0])[:, None], axis=0)
    return np.take_along_axis(mean_sorted, inv_order, axis=0)          # G x N

def make_unique(index):
    """
    Make index labels unique by appending .2, .3, ... for duplicates.
    """
    counts = {}
    out = []
    for name in index:
        c = counts.get(name, 0) + 1
        counts[name] = c
        out.append(name if c == 1 else f"{name}.{c}")
    return out

# -------------------- core CIBERSORT solver --------------------
def core_alg(X: np.ndarray, y: np.ndarray, absolute=False, abs_method='sig.score'):
    """
    Solve cell-type weights for a single sample using linear NuSVR with
    three nu values; pick the best by RMSE (tie-break by correlation).
    X: G x C (genes x cell types), standardized.
    y: G (z-scored mixture).
    Returns dict with weights, RMSE, and correlation.
    """
    nu_values = (0.25, 0.5, 0.75)
    best = None
    best_rmse = np.inf
    best_corr = -2.0

    for nu_val in nu_values:
        model = NuSVR(kernel='linear', nu=nu_val)
        model.fit(X, y)  # X: genes x celltypes, y: genes
        coef = model.dual_coef_               # (1, nSV)
        SV = model.support_vectors_           # (nSV, n_features=celltypes)
        w = np.dot(coef, SV).ravel().astype(np.float32)  # celltypes
        w[w < 0] = 0.0                        # clamp negatives

        if absolute:
            w_use = w
        else:
            s = w.sum()
            if s <= 0:
                # fallback to uniform weights to avoid NaNs
                w_use = np.full_like(w, 1.0 / max(len(w), 1), dtype=np.float32)
            else:
                w_use = w / s

        # Reconstruct the mixture and evaluate fit
        k = X @ w_use                          # G
        rmse = float(np.sqrt(np.mean((k - y) ** 2)))
        r = corr_pearson_fast(k, y)

        if (rmse < best_rmse) or (rmse == best_rmse and r > best_corr):
            best_rmse = rmse
            best_corr = r
            best = (w_use, rmse, r)

    return {"w": best[0], "mix_rmse": best[1], "mix_r": best[2]}

# -------------------- permutation helpers --------------------
def _one_perm_mixr(X: np.ndarray, Y_flat: np.ndarray, absolute: bool, abs_method: str, seed: int):
    """
    One permutation draw: sample a gene-length vector with replacement
    from the flattened mixture, z-score it, solve, and return correlation.
    """
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, Y_flat.size, size=X.shape[0], dtype=np.int64)
    yr = zscore1d(Y_flat[idx])
    return core_alg(X, yr, absolute, abs_method)["mix_r"]

def do_perm(perm: int, X: np.ndarray, Y: np.ndarray, absolute: bool, abs_method: str, n_jobs: int):
    """
    Build a null distribution of correlations using permutations.
    Parallelizes across permutations. Returns a sorted 1D array.
    Shows a progress bar if tqdm is available.
    """
    if perm <= 0:
        return None
    Y_flat = Y.ravel().astype(np.float32, copy=False)

    from numpy.random import SeedSequence
    child_seqs = SeedSequence().spawn(int(perm))
    seeds = [int(cs.generate_state(1)[0]) for cs in child_seqs]

    # Progress bar integrated with joblib (works with or without tqdm.contrib)
    with tqdm_joblib(tqdm(total=perm, desc="Permutations", unit="perm",
                          dynamic_ncols=True, mininterval=0.2)):
        mixr_dist = Parallel(n_jobs=max(1, int(n_jobs)), prefer="threads")(
            delayed(_one_perm_mixr)(X, Y_flat, absolute, abs_method, s) for s in seeds
        )
    return np.sort(np.asarray(mixr_dist, dtype=np.float32))

# -------------------- public API --------------------
def cibersort(input_path, perm=100, QN=True, absolute=False, abs_method='sig.score', n_jobs=1):
    """
    Run CIBERSORT-like deconvolution.
    - input_path: mixture matrix file (rows=genes, cols=samples)
    - perm: number of permutations for p-value estimation
    - QN: apply quantile normalization (typical for microarrays; often False for RNA-seq)
    - absolute / abs_method: absolute mode options (compatible naming only)
    - n_jobs: parallel workers for permutations and per-sample solving
    Returns a DataFrame indexed by sample with per-celltype weights and stats.
    """
    # 1) Load signature and mixture
    resource_pkg = 'iobrpy.resources'
    lm22_path = files(resource_pkg).joinpath('lm22.txt')
    # Robust whitespace parsing for LM22 (space-delimited text)
    sig_df = pd.read_csv(lm22_path, sep=r'\s+', engine='python', index_col=0)
    # Auto-detect CSV/TSV for mixture; supports gz
    mix_df = pd.read_csv(input_path, sep=None, engine='python', index_col=0)
    print(f"Loading mixture matrix from: {input_path} (this may take a while)...", flush=True)
    print(f"Mixture matrix loaded: {mix_df.shape[0]} genes × {mix_df.shape[1]} samples.", flush=True)

    # Ensure unique gene names to avoid accidental collapse later
    mix_df.index = make_unique(mix_df.index)

    # 2) Intersect genes
    common = sig_df.index.intersection(mix_df.index)
    if len(common) == 0:
        raise ValueError("No overlapping genes found between signature and mixture matrices.")
    sig_df = sig_df.loc[common].sort_index()
    mix_df = mix_df.loc[common].sort_index()

    # 3) To ndarray (float32)
    X = sig_df.to_numpy(dtype=np.float32, copy=True)   # G x C
    Y = mix_df.to_numpy(dtype=np.float32, copy=True)   # G x N

    # 4) Back-transform if mixture looks log2-like (heuristic)
    if np.max(Y) < 50:  # typical log2 TPM values are < ~20
        np.exp2(Y, out=Y)

    # 5) Optional quantile normalization (columns = samples)
    if QN:
        Y = quantile_normalize_fast(Y)

    # 6) Standardize signature globally (mean/std over all entries)
    X_mean = X.mean(dtype=np.float64)
    X_std = X.std(dtype=np.float64)
    if X_std == 0.0:
        raise ValueError("Signature matrix has zero variance.")
    X = (X - X_mean) / X_std
    X = X.astype(np.float32, copy=False)

    # 7) Z-score mixture per sample (column-wise)
    Yz = (Y - Y.mean(axis=0, keepdims=True)) / (Y.std(axis=0, keepdims=True) + 1e-12)
    Yz = Yz.astype(np.float32, copy=False)

    # 8) Build null distribution via permutations (parallel with progress)
    nulldist = do_perm(perm, X, Y, absolute, abs_method, n_jobs=n_jobs)

    # 9) Solve each sample in parallel (with progress)
    _, N = Yz.shape

    def _solve_one(i):
        res = core_alg(X, Yz[:, i], absolute, abs_method)
        return res["w"], res["mix_r"], res["mix_rmse"]

    with tqdm_joblib(tqdm(total=N, desc="Deconvolving samples", unit="sample",
                          dynamic_ncols=True, mininterval=0.2)):
        outs = Parallel(n_jobs=max(1, int(n_jobs)), prefer="threads")(
            delayed(_solve_one)(i) for i in range(N)
        )

    # 10) Assemble results
    weights = [o[0] for o in outs]
    rs = np.array([o[1] for o in outs], dtype=np.float32)
    rmses = np.array([o[2] for o in outs], dtype=np.float32)

    if nulldist is not None:
        # One-sided p-value: proportion of null correlations >= observed r
        ranks = np.searchsorted(nulldist, rs, side='right')
        pvals = 1.0 - (ranks / len(nulldist))
    else:
        pvals = np.full(N, 9999.0, dtype=np.float32)

    rows = []
    for i in range(N):
        w = weights[i]
        row = list(w) + [float(pvals[i]), float(rs[i]), float(rmses[i])]
        if absolute:
            row.append(float(np.sum(w)))
        rows.append(row)

    colnames = list(sig_df.columns) + ['P-value', 'Correlation', 'RMSE']
    if absolute:
        safe_method = abs_method.replace('.', '_')
        colnames.append(f'Absolute_score_({safe_method})')

    result_df = pd.DataFrame(rows, columns=colnames, index=mix_df.columns)
    return result_df

# -------------------- CLI entry (kept compatible) --------------------
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

    df = cibersort(
        input_path=args.input_path,
        perm=args.perm,
        QN=args.QN,
        absolute=args.absolute,
        abs_method=args.abs_method,
        n_jobs=args.threads
    )
    if df is not None:
        # Keep column suffix compatible with your pipeline
        df.columns = [c + '_CIBERSORT' for c in df.columns]
        df.index.name = 'ID'
        delim = ',' if args.output.lower().endswith('.csv') else '\t'
        df.to_csv(args.output, sep=delim)
        print(f"[Done] Output written to {args.output}")
    else:
        print("[Error] No results generated")