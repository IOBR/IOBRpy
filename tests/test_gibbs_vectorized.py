import time

import numpy as np
import pandas as pd

from iobrpy.bayesprism.gibbs import GibbsSampler, multinomial_rvs


def test_vectorized_sampling_shapes():
    phi = pd.DataFrame(np.array([[0.6, 0.4], [0.4, 0.6]]))
    X_n = np.array([30, 20])
    alpha = np.array([1.0, 1.0])
    gibbs_idx = np.array([0, 1, 2])

    results = GibbsSampler.sample_Z_theta_n(X_n=X_n, phi=phi, alpha=alpha, gibbs_idx=gibbs_idx, seed=0)

    assert results["Z_n"].shape == (phi.shape[1], phi.shape[0])
    assert results["theta_n"].shape == (phi.shape[0],)
    assert results["theta.cv_n"].shape == (phi.shape[0],)


def _legacy_multinomial_sampling(X_n, prob_mat):
    G = prob_mat.shape[1]
    Z = np.empty((G, prob_mat.shape[0]), dtype=int)
    for g in range(G):
        pvals = prob_mat[:, g] / np.sum(prob_mat[:, g])
        Z[g, :] = np.random.multinomial(n=X_n[g], pvals=pvals)
    return Z


def test_vectorized_sampling_performance_regression():
    K, G = 8, 300
    phi = pd.DataFrame(np.abs(np.random.default_rng(0).normal(size=(K, G))))
    theta = np.repeat(1 / K, K)
    X_n = np.random.default_rng(1).integers(low=50, high=100, size=G)
    prob_mat = phi.to_numpy() * theta[:, np.newaxis]

    iterations = 75

    np.random.seed(2)
    start_old = time.perf_counter()
    for _ in range(iterations):
        _legacy_multinomial_sampling(X_n, prob_mat)
    legacy_time = time.perf_counter() - start_old

    np.random.seed(2)
    start_new = time.perf_counter()
    for _ in range(iterations):
        pvals = prob_mat.T / prob_mat.sum(axis=0)[:, np.newaxis]
        multinomial_rvs(X_n, pvals)
    vectorized_time = time.perf_counter() - start_new

    assert vectorized_time <= legacy_time * 0.9

