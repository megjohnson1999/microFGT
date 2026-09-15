"""Analysis layer — thin scikit-bio wrappers, results written onto the modality."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microfgt import analysis


def _toy_composition():
    counts = np.array(
        [[40, 1, 0, 2], [2, 35, 1, 0], [0, 1, 38, 3], [30, 2, 0, 25],
         [1, 40, 2, 0], [0, 0, 33, 5]],
        dtype=float,
    )
    obs = pd.DataFrame(
        {"group": ["A", "B", "C", "A", "B", "C"]},
        index=[f"s{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=["t1", "t2", "t3", "t4"])
    a = ad.AnnData(X=counts.astype(np.float32), obs=obs, var=var)
    a.layers["counts"] = counts
    return a


def test_relative_abundance_rows_sum_to_one():
    a = analysis.relative_abundance(_toy_composition())
    np.testing.assert_allclose(a.layers["relabund"].sum(axis=1), np.ones(6), atol=1e-6)


def test_clr_handles_zeros_and_is_centered():
    a = analysis.clr_transform(_toy_composition())
    clr = a.layers["clr"]
    assert np.isfinite(clr).all()                       # zeros handled, no -inf
    np.testing.assert_allclose(clr.sum(axis=1), np.zeros(6), atol=1e-6)  # CLR rows sum to 0


def test_alpha_diversity_written_to_obs():
    a = analysis.alpha_diversity(_toy_composition(), metric="shannon")
    assert "alpha_shannon" in a.obs
    assert a.obs["alpha_shannon"].notna().all()


def test_beta_and_ordinate_written_to_obsp_obsm():
    a = _toy_composition()
    dm = analysis.beta_diversity(a, metric="braycurtis")
    assert a.obsp["beta_braycurtis"].shape == (6, 6)
    assert dm.shape == (6, 6)
    analysis.ordinate(a, metric="braycurtis")
    assert a.obsm["X_pcoa"].shape[0] == 6
    assert "proportion_explained" in a.uns["X_pcoa"]


def test_differential_abundance_returns_per_taxon_result():
    a = _toy_composition()
    res = analysis.differential_abundance(a, group_key="group", method="ancom")
    assert list(res.index) == list(a.var_names)         # one row per taxon
    assert "W" in res.columns


def test_diffabund_needs_two_groups():
    a = _toy_composition()
    a.obs["group"] = "A"                                 # collapse to one group
    with pytest.raises(ValueError, match=">=2 groups"):
        analysis.differential_abundance(a, group_key="group")
