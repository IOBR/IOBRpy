import numpy as np
import pandas as pd
import pytest
from importlib.resources import files
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from iobrpy.workflow import cibersort as cibersort_module


def _make_mixture(sig_index, values):
    return pd.DataFrame(
        values,
        index=sig_index[: values.shape[0]],
        columns=["s1", "s2"],
        dtype=np.float32,
    )


def _stub_core_alg(X, y, absolute=False, abs_method="sig.score"):
    n_cell_types = X.shape[1]
    w = np.arange(1, n_cell_types + 1, dtype=np.float32)
    return {"w": w, "mix_r": 0.5, "mix_rmse": 0.1}


@pytest.fixture
def sig_df():
    lm22_path = files("iobrpy.resources").joinpath("lm22.txt")
    return pd.read_csv(lm22_path, sep=r"\s+", engine="python", index_col=0)


def test_absolute_sigscore_scales_by_sample_median(tmp_path, monkeypatch, sig_df):
    mix_df = _make_mixture(sig_df.index, np.array([[1.0, 3.0], [5.0, 7.0]], dtype=np.float32))
    mix_path = tmp_path / "mix.tsv"
    mix_df.to_csv(mix_path, sep="\t")

    monkeypatch.setattr(cibersort_module, "core_alg", _stub_core_alg)

    result = cibersort_module.cibersort(
        input_path=mix_path,
        perm=0,
        QN=False,
        absolute=True,
        abs_method="sig.score",
        n_jobs=1,
    )

    mix_arr = mix_df.to_numpy(dtype=np.float32, copy=True)
    if np.max(mix_arr) < 50:
        np.exp2(mix_arr, out=mix_arr)
    y_median_global = max(float(np.median(mix_arr)), 1.0)
    expected_w_s1 = np.arange(1, sig_df.shape[1] + 1, dtype=np.float32) * (
        np.median(mix_arr[:, 0]) / y_median_global
    )
    expected_w_s2 = np.arange(1, sig_df.shape[1] + 1, dtype=np.float32) * (
        np.median(mix_arr[:, 1]) / y_median_global
    )

    abs_col = "Absolute_score_(sig_score)"
    np.testing.assert_allclose(result.loc["s1", sig_df.columns], expected_w_s1)
    np.testing.assert_allclose(result.loc["s2", sig_df.columns], expected_w_s2)
    assert np.isclose(result.loc["s1", abs_col], expected_w_s1.sum())
    assert np.isclose(result.loc["s2", abs_col], expected_w_s2.sum())


def test_absolute_sigscore_uses_floor_for_global_median(tmp_path, monkeypatch, sig_df):
    mix_df = _make_mixture(sig_df.index, np.array([[0.01, 0.02], [0.03, 0.04]], dtype=np.float32))
    mix_path = tmp_path / "mix_floor.tsv"
    mix_df.to_csv(mix_path, sep="\t")

    monkeypatch.setattr(cibersort_module, "core_alg", _stub_core_alg)

    result = cibersort_module.cibersort(
        input_path=mix_path,
        perm=0,
        QN=False,
        absolute=True,
        abs_method="sig.score",
        n_jobs=1,
    )

    mix_arr = mix_df.to_numpy(dtype=np.float32, copy=True)
    if np.max(mix_arr) < 50:
        np.exp2(mix_arr, out=mix_arr)
    y_median_global = max(float(np.median(mix_arr)), 1.0)
    expected_w_s1 = np.arange(1, sig_df.shape[1] + 1, dtype=np.float32) * (
        np.median(mix_arr[:, 0]) / y_median_global
    )
    expected_w_s2 = np.arange(1, sig_df.shape[1] + 1, dtype=np.float32) * (
        np.median(mix_arr[:, 1]) / y_median_global
    )

    abs_col = "Absolute_score_(sig_score)"
    np.testing.assert_allclose(result.loc["s1", sig_df.columns], expected_w_s1)
    np.testing.assert_allclose(result.loc["s2", sig_df.columns], expected_w_s2)
    assert np.isclose(result.loc["s1", abs_col], expected_w_s1.sum())
    assert np.isclose(result.loc["s2", abs_col], expected_w_s2.sum())


def test_pvalue_uses_r_formula_with_ge_counts(tmp_path, monkeypatch, sig_df):
    mix_df = _make_mixture(sig_df.index, np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32))
    mix_path = tmp_path / "mix_pval.tsv"
    mix_df.to_csv(mix_path, sep="\t")

    rs_values = np.array([0.5, 0.9], dtype=np.float32)
    r_iter = iter(rs_values.tolist())

    def _core_alg_with_r(X, y, absolute=False, abs_method="sig.score"):
        return {
            "w": np.ones(X.shape[1], dtype=np.float32),
            "mix_r": next(r_iter),
            "mix_rmse": 0.0,
        }

    nulldist = np.array([0.1, 0.2, 0.5, 0.5, 0.9, 1.2], dtype=np.float32)

    monkeypatch.setattr(cibersort_module, "core_alg", _core_alg_with_r)
    monkeypatch.setattr(
        cibersort_module,
        "do_perm",
        lambda perm, X, Y, absolute, abs_method, n_jobs: nulldist,
    )

    result = cibersort_module.cibersort(
        input_path=mix_path,
        perm=5,
        QN=False,
        absolute=False,
        abs_method="sig.score",
        n_jobs=1,
    )

    expected_ranks = np.searchsorted(nulldist, rs_values, side="left")
    expected_p = (len(nulldist) - expected_ranks + 1) / (len(nulldist) + 1)

    np.testing.assert_allclose(
        result["P-value"].to_numpy(dtype=np.float32), expected_p.astype(np.float32)
    )


def test_zscore1d_uses_sample_std():
    arr = np.array([1.0, 3.0, 5.0], dtype=np.float32)
    # Sample variance of [1,3,5] is 4 -> std = 2
    expected = np.array([-1.0, 0.0, 1.0], dtype=np.float32)
    out = cibersort_module.zscore1d(arr)
    np.testing.assert_allclose(out, expected)
