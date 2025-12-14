import numpy as np
import scipy
import time
from datetime import datetime, timedelta
import multiprocessing
from itertools import repeat
import pandas as pd

from . import references
from . import joint_post
from . import theta_post


#https://stackoverflow.com/questions/55818845/fast-vectorized-multinomial-in-python
def multinomial_rvs(count, p, rng=None, method="binomial"):
    """Vectorized multinomial sampling.

    Parameters
    ----------
    count : array-like
        Total counts for each independent draw. Should broadcast to ``p``
        without relying on the last axis.
    p : ndarray
        Probability matrix where the last axis enumerates categories.
    rng : np.random.RandomState or np.random.Generator, optional
        RNG to use; defaults to a new Generator.
    method : {"binomial", "sequential"}
        "binomial" uses a fast chain of binomial draws (distributionally
        equivalent but consumes RNG state differently from ``Generator.multinomial``).
        "sequential" mirrors the per-vector ``rng.multinomial`` calls for
        deterministic compatibility with the legacy sampler.

    Returns
    -------
    ndarray
        Samples with the same shape as ``p``.
    """
    if rng is None:
        rng = np.random.default_rng()

    p = np.asarray(p)
    count = np.array(count, copy=True)

    if method == "sequential":
        flat_samples = np.array(
            [rng.multinomial(int(n), prob) for n, prob in zip(count.reshape(-1), p.reshape(-1, p.shape[-1]))],
            dtype=int,
        )
        return flat_samples.reshape(p.shape)

    if method != "binomial":
        raise ValueError("Unsupported multinomial sampling method")

    out = np.zeros(p.shape, dtype=int)
    ps = p.cumsum(axis=-1)
    # Conditional probabilities
    with np.errstate(divide='ignore', invalid='ignore'):
        condp = p / ps
    condp[np.isnan(condp)] = 0.0
    for i in range(p.shape[-1]-1, 0, -1):
        binsample = rng.binomial(count, condp[..., i])
        out[..., i] = binsample
        count -= binsample
    out[..., 0] = count
    return out


class GibbsSampler:
    def __init__(self, reference, X, gibbs_control):
        self.reference = reference
        self.X = X
        self.gibbs_control = gibbs_control


    def get_gibbs_idx(gibbs_control):
        chain_length = int(gibbs_control['chain.length'])
        burn_in = int(gibbs_control['burn.in'])
        thinning = int(gibbs_control['thinning'])
        all_idx = np.arange(1, chain_length + 1)
        burned_idx = all_idx[burn_in:]
        thinned_idx = burned_idx[np.arange(0, len(burned_idx), thinning)]
        return thinned_idx


    def rdirichlet(alpha, rng=None, backend="generator"):
        """Dirichlet sampling using the provided RNG for determinism."""

        if rng is None:
            rng = GibbsSampler._make_rng(None, backend=backend)
        x = rng.gamma(alpha, size=len(alpha))
        return x / np.sum(x)


    def _make_rng(seed, backend="generator"):
        """Return a dedicated RNG for Gibbs sampling.

        The backend can be ``generator`` (default MT19937-based Generator) or
        ``randomstate`` to mirror legacy NumPy behaviour. The seed can be a
        SeedSequence, Generator, RandomState, int, or None.
        """

        backend = backend.lower()
        if backend not in {"generator", "randomstate"}:
            raise ValueError("Unsupported RNG backend; use 'generator' or 'randomstate'")

        if backend == "randomstate":
            if isinstance(seed, np.random.RandomState):
                return seed
            if isinstance(seed, np.random.SeedSequence):
                seed = int(seed.generate_state(1, dtype=np.uint32)[0])
            if isinstance(seed, np.random.Generator):
                seed = seed.integers(0, np.iinfo(np.uint32).max, dtype=np.uint32)
            return np.random.RandomState(seed)

        # Generator backend
        if isinstance(seed, np.random.Generator):
            return seed
        if isinstance(seed, np.random.RandomState):
            seed = seed.randint(0, np.iinfo(np.uint32).max, dtype=np.uint32)
        if isinstance(seed, np.random.SeedSequence):
            return np.random.Generator(np.random.MT19937(seed))
        if seed is None:
            return np.random.Generator(np.random.MT19937())
        return np.random.Generator(np.random.MT19937(seed))

    def _spawn_seeds(seed, n_children, backend="generator"):
        """Spawn child seeds by advancing a base RNG once per task."""

        if seed is None:
            return [None] * n_children

        backend = backend.lower()
        base_rng = GibbsSampler._make_rng(seed, backend=backend)

        if backend == "generator":
            return base_rng.integers(0, np.iinfo(np.uint32).max, size=n_children, dtype=np.uint32).tolist()

        if backend == "randomstate":
            return base_rng.randint(0, np.iinfo(np.uint32).max, size=n_children, dtype=np.uint32).tolist()

        raise ValueError("Unsupported RNG backend; use 'generator' or 'randomstate'")

    def sample_Z_theta_n(
        X_n,
        phi,
        alpha,
        gibbs_idx,
        chain_length,
        seed=None,
        rng_backend="generator",
        compute_elbo=False,
        fast_multinomial=False,
    ):

        rng = GibbsSampler._make_rng(seed, backend=rng_backend)

        phi = np.asarray(phi)
        G = phi.shape[1]
        K = phi.shape[0]

        theta_n_i = np.repeat(1 / K, K)
        Z_n_i = np.empty((G, K))

        Z_n_sum = np.zeros((G, K))
        theta_n_sum = np.zeros(K)
        theta_n2_sum = np.zeros(K)

        multinom_coef = 0

        gibbs_idx_arr = np.asarray(gibbs_idx, dtype=int)
        gibbs_idx_set = set(gibbs_idx_arr.tolist())
        iterations = int(np.max(gibbs_idx_arr)) if gibbs_idx_set else int(chain_length)

        for i in range(1, iterations + 1):
            prob_mat = phi * theta_n_i[:, np.newaxis]

            if fast_multinomial and not isinstance(rng, np.random.RandomState):
                prob_mat_norm = prob_mat / prob_mat.sum(axis=0, keepdims=True)
                Z_n_i = multinomial_rvs(
                    count=X_n,
                    p=prob_mat_norm.T,
                    rng=rng,
                    method="binomial",
                )
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    prob_mat_norm = prob_mat / prob_mat.sum(axis=0, keepdims=True)
                prob_mat_norm[np.isnan(prob_mat_norm)] = 0.0

                Z_n_i = np.empty((G, K), dtype=int)
                prob_mat_T = prob_mat_norm.T
                for g_idx, (count_g, prob_g) in enumerate(zip(X_n, prob_mat_T)):
                    Z_n_i[g_idx, :] = rng.multinomial(int(count_g), prob_g)

            Z_nk_i = np.sum(Z_n_i, axis=0)
            theta_n_i = GibbsSampler.rdirichlet(alpha=Z_nk_i + alpha, rng=rng, backend=rng_backend)

            if i in gibbs_idx_set:
                Z_n_sum += Z_n_i
                theta_n_sum += theta_n_i
                theta_n2_sum += theta_n_i**2
                if compute_elbo:
                    multinom_coef += np.sum(np.log(scipy.special.factorial(Z_nk_i))) - np.sum(
                        np.log(scipy.special.factorial(Z_n_i))
                    )

        samples_size = len(gibbs_idx)
        Z_n = Z_n_sum / samples_size
        theta_n = theta_n_sum / samples_size
        theta_cv_n = np.sqrt(theta_n2_sum / samples_size - (theta_n ** 2)) / theta_n
        gibbs_constant = multinom_coef / samples_size

        return {
            'Z_n': Z_n,
            'theta_n': theta_n,
            'theta.cv_n': theta_cv_n,
            'gibbs.constant': gibbs_constant,
        }

    def sample_theta_n(X_n, phi, alpha, gibbs_idx, chain_length, seed=None, rng_backend="generator", fast_multinomial=False):

        rng = GibbsSampler._make_rng(seed, backend=rng_backend)

        phi = np.asarray(phi)
        G = phi.shape[1]
        K = phi.shape[0]

        theta_n_i = np.repeat(1 / K, K)
        Z_n_i = np.empty((G, K))

        theta_n_sum = np.zeros(K)
        theta_n2_sum = np.zeros(K)

        gibbs_idx_arr = np.asarray(gibbs_idx, dtype=int)
        gibbs_idx_set = set(gibbs_idx_arr.tolist())
        iterations = int(np.max(gibbs_idx_arr)) if gibbs_idx_set else int(chain_length)

        for i in range(1, iterations + 1):
            prob_mat = phi * theta_n_i[:, np.newaxis]

            if fast_multinomial and not isinstance(rng, np.random.RandomState):
                prob_mat_norm = prob_mat / prob_mat.sum(axis=0, keepdims=True)
                Z_n_i = multinomial_rvs(
                    count=X_n,
                    p=prob_mat_norm.T,
                    rng=rng,
                    method="binomial",
                )
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    prob_mat_norm = prob_mat / prob_mat.sum(axis=0, keepdims=True)
                prob_mat_norm[np.isnan(prob_mat_norm)] = 0.0

                prob_mat_T = prob_mat_norm.T
                for g_idx, (count_g, prob_g) in enumerate(zip(X_n, prob_mat_T)):
                    Z_n_i[g_idx, :] = rng.multinomial(int(count_g), prob_g)

            theta_n_i = GibbsSampler.rdirichlet(alpha=np.sum(Z_n_i, axis=0) + alpha, rng=rng, backend=rng_backend)

            if i in gibbs_idx_set:
                theta_n_sum += theta_n_i
                theta_n2_sum += theta_n_i**2

        samples_size = len(gibbs_idx)
        theta_n = theta_n_sum / samples_size
        theta_cv_n = np.sqrt(theta_n2_sum / samples_size - (theta_n**2)) / theta_n

        return {'theta_n': theta_n, 'theta.cv_n': theta_cv_n}


    def my_seconds_to_period(x):
        days = round(x // (60 * 60 * 24))
        hours = round((x - days * 60 * 60 * 24) // (60 * 60))
        minutes = round((x - days * 60 * 60 * 24 - hours * 60 * 60) // 60) + 1
        days_str = '' if days == 0 else str(days) + 'days '
        hours_str = '' if (hours == 0 and days == 0) else str(hours) + 'hrs '
        minutes_str = '' if (minutes == 0 and days == 0 and hours == 0) else str(minutes) + 'mins'
        final_str = days_str + hours_str + minutes_str
        return final_str


    def estimate_gibbs_time(self, final, chain_length = 50):
        ref = self.reference
        X = self.X.to_numpy()
        gibbs_control = self.gibbs_control
        fast_mult = gibbs_control.get('fast.multinomial', False)
        rng_backend = gibbs_control.get('rng.backend', 'generator')
        ptm = time.process_time()
        
        if not final:
            assert isinstance(ref, references.RefPhi), "Gibbs is not final but ref is not refPhi"
            GibbsSampler.sample_Z_theta_n(
                X_n = X[0, :],
                phi = ref.phi,
                alpha = gibbs_control['alpha'],
                gibbs_idx = GibbsSampler.get_gibbs_idx(
                    {'chain.length' : chain_length,
                     'burn.in' : chain_length * gibbs_control['burn.in'] / gibbs_control['chain.length'],
                     'thinning' : gibbs_control['thinning']}),
                chain_length=chain_length,
                seed = gibbs_control['seed'],
                rng_backend = rng_backend,
                compute_elbo = False,
                fast_multinomial = fast_mult,
            )
        else:
            if isinstance(ref, references.RefPhi):
                GibbsSampler.sample_theta_n(
                    X_n = X[0, :],
                    phi = ref.phi,
                    alpha = gibbs_control['alpha'],
                    gibbs_idx = GibbsSampler.get_gibbs_idx(
                        {'chain.length' : chain_length,
                         'burn.in' : chain_length * gibbs_control['burn.in'] / gibbs_control['chain.length'],
                         'thinning' : gibbs_control['thinning']}),
                    chain_length=chain_length,
                    seed = gibbs_control['seed'],
                    rng_backend = rng_backend,
                    fast_multinomial = fast_mult,
                )
            if isinstance(ref, references.RefTumor):
                phi_1 = pd.concat([pd.DataFrame(ref.psi_mal.iloc[0, :]).T, ref.psi_env])
                nonzero_idx = np.max(phi_1, axis = 0) > 0
                GibbsSampler.sample_theta_n(
                    X_n = X[0, nonzero_idx],
                    phi = phi_1.loc[:, nonzero_idx],
                    alpha = gibbs_control['alpha'],
                    gibbs_idx = GibbsSampler.get_gibbs_idx(
                        {'chain.length' : chain_length,
                         'burn.in' : chain_length*gibbs_control['burn.in'] / gibbs_control['chain.length'],
                         'thinning' : gibbs_control['thinning']}),
                    chain_length=chain_length,
                    seed = gibbs_control['seed'],
                    rng_backend = rng_backend,
                    fast_multinomial = fast_mult,
                )
        
        total_time = time.process_time() - ptm
        estimated_time = gibbs_control['chain.length'] / chain_length * total_time * np.ceil(X.shape[0] / gibbs_control['n.cores']) * 2
        current_time = datetime.now()
        print("Current time: ", current_time)
        print("Estimated time to complete: ", GibbsSampler.my_seconds_to_period(estimated_time))
        print("Estimated finishing time: ", current_time + timedelta(seconds = estimated_time))


    def run_gibbs_refPhi(self, final, compute_elbo):

        assert isinstance(self.reference, references.RefPhi)
        phi_df = self.reference.phi
        phi = phi_df
        X = self.X.to_numpy()
        gibbs_control = self.gibbs_control
        alpha = gibbs_control['alpha']
        fast_mult = gibbs_control.get('fast.multinomial', False)
        rng_backend = gibbs_control.get('rng.backend', 'generator')
        gibbs_idx = GibbsSampler.get_gibbs_idx(gibbs_control)
        chain_length = gibbs_control['chain.length']
        seed = gibbs_control['seed']
        print("Start run...")

        pool_size = GibbsSampler._resolve_pool_size(gibbs_control['n.cores'])
        X_input = [X[i, :] for i in np.arange(X.shape[0])]
        phi_arr = phi_df.to_numpy()
        seeds = GibbsSampler._spawn_seeds(seed, X.shape[0], backend=rng_backend)

        if not final:
            star_input = list(
                zip(
                    X_input,
                    repeat(phi_arr),
                    repeat(alpha),
                    repeat(gibbs_idx),
                    repeat(chain_length),
                    seeds,
                    repeat(rng_backend),
                    repeat(compute_elbo),
                    repeat(fast_mult),
                )
            )
            gibbs_list = GibbsSampler._starmap_in_pool(
                GibbsSampler.sample_Z_theta_n, star_input, pool_size
            )
            return joint_post.JointPost.new(self.X.index, self.X.columns, phi_df.index, gibbs_list)
        else:
            star_input = list(
                zip(
                    X_input,
                    repeat(phi_arr),
                    repeat(alpha),
                    repeat(gibbs_idx),
                    repeat(chain_length),
                    seeds,
                    repeat(rng_backend),
                    repeat(fast_mult),
                )
            )
            gibbs_list = GibbsSampler._starmap_in_pool(
                GibbsSampler.sample_theta_n, star_input, pool_size
            )
            return theta_post.ThetaPost.new(self.X.index, self.X.columns, gibbs_list)


    def run_gibbs_refTumor(self):

        assert isinstance(self.reference, references.RefTumor)
        psi_mal = self.reference.psi_mal
        psi_env = self.reference.psi_env
        key = self.reference.key
        X = self.X.to_numpy()
        gibbs_control = self.gibbs_control
        alpha = gibbs_control['alpha']
        fast_mult = gibbs_control.get('fast.multinomial', False)
        rng_backend = gibbs_control.get('rng.backend', 'generator')
        gibbs_idx = GibbsSampler.get_gibbs_idx(gibbs_control)
        chain_length = gibbs_control['chain.length']
        seed = gibbs_control['seed']
        print("Start run...")

        psi_env_arr = psi_env.to_numpy()
        seeds = GibbsSampler._spawn_seeds(seed, X.shape[0], backend=rng_backend)
        star_input = []
        for i in range(X.shape[0]):
            psi_mal_n = psi_mal.iloc[i, :].to_numpy()
            phi_n = np.vstack([psi_mal_n, psi_env_arr])
            nonzero_idx = np.max(phi_n, axis = 0) > 0
            child_seed = seeds[i]
            star_input.append((
                X[i, nonzero_idx],
                phi_n[:, nonzero_idx],
                alpha,
                gibbs_idx,
                chain_length,
                child_seed,
                rng_backend,
                fast_mult,
            ))

        pool_size = GibbsSampler._resolve_pool_size(gibbs_control['n.cores'])
        gibbs_list = GibbsSampler._starmap_in_pool(
            GibbsSampler.sample_theta_n, star_input, pool_size
        )

        print("BayesPrism finished.")
        
        return theta_post.ThetaPost.new(self.X.index, [key] + list(psi_env.index), gibbs_list)


    def run(self, final, if_estimate = True, compute_elbo = False):
        if final:
            print("Run Gibbs sampling using updated reference ...")
        else:
            print("Run Gibbs sampling...")
        
        if if_estimate:
            self.estimate_gibbs_time(final = final)
        if isinstance(self.reference, references.RefPhi):
            return GibbsSampler.run_gibbs_refPhi(self, final = final, compute_elbo = compute_elbo)
        if isinstance(self.reference, references.RefTumor):
            return GibbsSampler.run_gibbs_refTumor(self)

    @staticmethod
    def _starmap_in_pool(func, star_input, pool_size):
        task_count = len(star_input)
        chunksize = GibbsSampler._compute_chunksize(task_count, pool_size)
        with multiprocessing.Pool(processes=pool_size) as pool:
            return pool.starmap(func, star_input, chunksize=chunksize)

    @staticmethod
    def _resolve_pool_size(requested):
        if requested <= 0:
            raise ValueError("n.cores must be a positive integer")
        return requested

    @staticmethod
    def _compute_chunksize(task_count, pool_size):
        if task_count <= 0:
            return 1
        if pool_size <= 1:
            return task_count
        return max(1, task_count // (pool_size * 2))
