"""Viz helpers — smoke tests (Agg backend; matplotlib is in the dev extra)."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

import anndata as ad  # noqa: E402

from microfgt import analysis  # noqa: E402
from microfgt.viz import alpha_boxplot, ordination_scatter  # noqa: E402


def _adata():
    counts = np.array([[40, 1, 0, 2], [2, 35, 1, 0], [0, 1, 38, 3], [30, 2, 0, 25]], dtype=float)
    a = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame({"group": ["A", "B", "A", "B"]}, index=[f"s{i}" for i in range(4)]),
        var=pd.DataFrame(index=["t1", "t2", "t3", "t4"]),
    )
    a.layers["counts"] = counts
    analysis.alpha_diversity(a, metric="shannon")
    analysis.ordinate(a, metric="braycurtis")
    return a


def test_ordination_scatter_returns_axes():
    ax = ordination_scatter(_adata(), color="group")
    assert ax.has_data()


def test_alpha_boxplot_returns_axes():
    ax = alpha_boxplot(_adata(), alpha_key="alpha_shannon", group_key="group")
    assert ax.get_ylabel() == "alpha_shannon"
