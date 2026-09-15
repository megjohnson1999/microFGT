"""``associate`` — two-variable association, dispatching on dtype.

In-memory AnnData with a categorical CST-like variable, a clinical categorical strongly
associated with it, a continuous variable that varies by group, and two correlated continuous
variables. Deterministic (no permutation), so p-values and effect sizes are stable.
"""

import numpy as np
import pandas as pd
import pytest

import anndata as ad

from microfgt.analysis import associate
from microfgt.analysis.results import AnalysisResult


@pytest.fixture
def data():
    cst = ["A"] * 10 + ["B"] * 10 + ["C"] * 10
    # bv strongly tracks cst: A->mostly pos, C->mostly neg, B->mixed.
    bv = (["pos"] * 9 + ["neg"]) + (["pos"] * 5 + ["neg"] * 5) + (["neg"] * 9 + ["pos"])
    ph = np.concatenate([                      # A low, B mid, C high
        4.0 + 0.1 * np.arange(10), 4.6 + 0.1 * np.arange(10), 5.2 + 0.1 * np.arange(10)])
    v1 = np.arange(30, dtype=float)
    v2 = 2.0 * v1 + (np.arange(30) % 3)        # strongly (not perfectly) correlated with v1
    n = 30
    adata = ad.AnnData(
        X=np.ones((n, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {"cst": cst, "bv": bv, "ph": ph, "v1": v1, "v2": v2},
            index=pd.Index([f"s{i:02d}" for i in range(n)], name="sample"),
        ),
        var=pd.DataFrame(index=["t0", "t1"]),
    )
    adata.layers["counts"] = np.ones((n, 2), dtype=np.int64)
    return adata


# --- categorical x categorical --------------------------------------------------------------
def test_cat_cat_chi_square_with_cramers_v(data):
    r = associate(data, "cst", "bv")
    assert isinstance(r, AnalysisResult)
    assert r.stats["test"] == "chi-square"
    assert r.stats["effect"] == "cramers_v"
    assert r.pvalue < 0.05                        # bv tracks cst by construction
    assert r.table.shape == (3, 2)                # contingency table (3 CST x 2 bv)
    assert r.spec["x_type"] == "categorical" and r.spec["y_type"] == "categorical"


def test_cat_cat_2x2_uses_fisher_with_odds_ratio(data):
    r = associate(data, "cst", "bv", subset={"cst": ["A", "C"]})   # -> 2x2
    assert r.stats["test"] == "Fisher exact"
    assert "odds_ratio" in r.stats
    assert r.pvalue < 0.05


def test_association_is_order_invariant(data):
    a = associate(data, "cst", "bv")
    b = associate(data, "bv", "cst")
    assert a.pvalue == pytest.approx(b.pvalue)


# --- categorical x continuous ---------------------------------------------------------------
def test_cat_cont_kruskal_across_three_groups(data):
    r = associate(data, "cst", "ph")             # 3 groups -> Kruskal
    assert r.stats["test"] == "Kruskal–Wallis"
    assert r.pvalue < 0.05
    assert list(r.table.index) == ["A", "B", "C"]
    assert r.table.loc["A", "median"] < r.table.loc["C", "median"]


def test_cat_cont_two_groups_uses_mannwhitney(data):
    r = associate(data, "cst", "ph", subset={"cst": ["A", "C"]})
    assert r.stats["test"] == "Mann–Whitney U"


# --- continuous x continuous ----------------------------------------------------------------
def test_cont_cont_spearman_default(data):
    r = associate(data, "v1", "v2")
    assert r.stats["test"] == "Spearman rho"
    assert r.stats["statistic"] > 0.9            # strongly correlated
    assert r.pvalue < 1e-6


def test_cont_cont_pearson_when_forced(data):
    r = associate(data, "v1", "v2", method="pearson")
    assert r.stats["test"] == "Pearson r"


# --- guards ---------------------------------------------------------------------------------
def test_single_level_after_subset_raises(data):
    with pytest.raises(ValueError, match=">=2 levels"):
        associate(data, "cst", "bv", subset={"cst": "A"})
