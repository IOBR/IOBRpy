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
        chain_length = gibbs_control['chain.length']
        burn_in = gibbs_control['burn.in']
        thinning = gibbs_control['thinning']
        all_idx = np.arange(0, chain_length)
        burned_idx = all_idx[int(burn_in):]
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
        """Create child seeds per worker consistent with the chosen backend."""

        if seed is None:
            return [None] * n_children

        backend = backend.lower()
        if backend == "generator":
            seed_seq = seed if isinstance(seed, np.random.SeedSequence) else np.random.SeedSequence(seed)
            return seed_seq.spawn(n_children)

        if backend == "randomstate":
            base_rng = GibbsSampler._make_rng(seed, backend=backend)
            return [base_rng.randint(0, np.iinfo(np.uint32).max, dtype=np.uint32) for _ in range(n_children)]

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

        phi = phi.to_numpy()
        G = phi.shape[1]
        K = phi.shape[0]

        theta_n_i = np.repeat(1 / K, K)
        Z_n_i = np.empty((G, K))

        Z_n_sum = np.zeros((G, K))
        theta_n_sum = np.zeros(K)
        theta_n2_sum = np.zeros(K)

        multinom_coef = 0

        iterations = chain_length if chain_length is not None else (np.max(gibbs_idx) + 1)

        for i in range(iterations):
            prob_mat = phi * theta_n_i[:, np.newaxis]
            prob_mat /= prob_mat.sum(axis=0, keepdims=True)

            Z_n_i = multinomial_rvs(
                count=X_n,
                p=prob_mat.T,
                rng=rng,
                method="binomial" if fast_multinomial else "sequential",
            )

            Z_nk_i = np.sum(Z_n_i, axis=0)
            theta_n_i = GibbsSampler.rdirichlet(alpha=Z_nk_i + alpha, rng=rng, backend=rng_backend)

            if i in gibbs_idx:
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

        phi = phi.to_numpy()
        G = phi.shape[1]
        K = phi.shape[0]

        theta_n_i = np.repeat(1 / K, K)
        Z_n_i = np.empty((G, K))

        theta_n_sum = np.zeros(K)
        theta_n2_sum = np.zeros(K)

        iterations = chain_length if chain_length is not None else (np.max(gibbs_idx) + 1)

        for i in range(iterations):
            prob_mat = phi * theta_n_i[:, np.newaxis]
            prob_mat /= prob_mat.sum(axis=0, keepdims=True)
            Z_n_i = multinomial_rvs(
                count=X_n,
                p=prob_mat.T,
                rng=rng,
                method="binomial" if fast_multinomial else "sequential",
            )

            theta_n_i = GibbsSampler.rdirichlet(alpha=np.sum(Z_n_i, axis=0) + alpha, rng=rng, backend=rng_backend)

            if i in gibbs_idx:
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
        phi = self.reference.phi
        X = self.X.to_numpy()
        gibbs_control = self.gibbs_control
        alpha = gibbs_control['alpha']
        fast_mult = gibbs_control.get('fast.multinomial', False)
        rng_backend = gibbs_control.get('rng.backend', 'generator')
        gibbs_idx = GibbsSampler.get_gibbs_idx(gibbs_control)
        chain_length = gibbs_control['chain.length']
        seed = gibbs_control['seed']
        print("Start run...")

        ctx = multiprocessing.get_context("fork")
        chunk_size = max(1, int(np.ceil(X.shape[0] / (gibbs_control['n.cores'] * 4))))
        if not final:
            seeds = GibbsSampler._spawn_seeds(seed, X.shape[0], backend=rng_backend)
            with ctx.Pool(processes=gibbs_control['n.cores']) as pool:
                X_input = [X[i, :] for i in np.arange(X.shape[0])]
                star_input = zip(
                    X_input,
                    repeat(phi),
                    repeat(alpha),
                    repeat(gibbs_idx),
                    repeat(chain_length),
                    seeds,
                    repeat(rng_backend),
                    repeat(compute_elbo),
                    repeat(fast_mult),
                )
                gibbs_list = pool.starmap(GibbsSampler.sample_Z_theta_n, star_input, chunksize=chunk_size)
            return joint_post.JointPost.new(self.X.index, self.X.columns, phi.index, gibbs_list)
        else:
            seeds = GibbsSampler._spawn_seeds(seed, X.shape[0], backend=rng_backend)
            with ctx.Pool(processes=gibbs_control['n.cores']) as pool:
                X_input = [X[i, :] for i in np.arange(X.shape[0])]
                star_input = zip(
                    X_input,
                    repeat(phi),
                    repeat(alpha),
                    repeat(gibbs_idx),
                    repeat(chain_length),
                    seeds,
                    repeat(rng_backend),
                    repeat(fast_mult),
                )
                gibbs_list = pool.starmap(GibbsSampler.sample_theta_n , star_input, chunksize=chunk_size)
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

        star_input = []
        for i in range(X.shape[0]):
            psi_mal_n = pd.DataFrame(psi_mal.iloc[i, :]).T
            phi_n = pd.concat([psi_mal_n, psi_env])
            nonzero_idx = np.max(phi_n, axis = 0) > 0
            child_seed = GibbsSampler._spawn_seeds(seed, 1, backend=rng_backend)[0]
            star_input.append((
                X[i, nonzero_idx],
                phi_n.loc[:, nonzero_idx],
                alpha,
                gibbs_idx,
                chain_length,
                child_seed,
                rng_backend,
                fast_mult,
            ))

        ctx = multiprocessing.get_context("fork")
        chunk_size = max(1, int(np.ceil(len(star_input) / (gibbs_control['n.cores'] * 4))))
        with ctx.Pool(processes=gibbs_control['n.cores']) as pool:
            gibbs_list = pool.starmap(GibbsSampler.sample_theta_n , star_input, chunksize=chunk_size)
        
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
