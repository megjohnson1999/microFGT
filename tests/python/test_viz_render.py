"""viz.render — one renderer that draws any AnalysisResult from its own plot spec + data.

Builds each verb's result on a small in-memory MuData and checks render() returns a populated
Axes with the expected labels. Agg backend; no display.
"""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402

from microfgt.analysis import associate, compare_abundance, compare_alpha, compare_beta  # noqa: E402
from microfgt.viz import render  # noqa: E402


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    n_per, n_taxa = 12, 6
    low = rng.multinomial(400, [0.7, 0.1, 0.05, 0.05, 0.05, 0.05], size=n_per)
    high = rng.multinomial(400, [0.2, 0.2, 0.15, 0.15, 0.15, 0.15], size=n_per)
    counts = np.vstack([low, high]).astype(np.int64)
    samples = [f"s{i:02d}" for i in range(2 * n_per)]
    taxon = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=[f"taxon{j}" for j in range(n_taxa)]),
    )
    taxon.layers["counts"] = counts
    comp = taxon.copy()
    comp.obs["group"] = ["A"] * n_per + ["B"] * n_per
    comp.obs["bv"] = (["pos"] * 9 + ["neg"] * 3) + (["neg"] * 9 + ["pos"] * 3)
    comp.obs["ph"] = np.r_[4.0 + 0.1 * np.arange(n_per), 5.0 + 0.1 * np.arange(n_per)]
    comp.obs["score"] = np.linspace(0.1, 0.9, 2 * n_per)
    with md.set_options(pull_on_update=False):
        return md.MuData({"composition": comp, "composition_taxon": taxon})


def test_render_box_from_compare_alpha(data):
    ax = render(compare_alpha(data, "group"))
    assert ax.has_data()
    assert ax.get_xlabel() == "group" and ax.get_ylabel() == "alpha_shannon"


def test_render_ordination_from_compare_beta(data):
    r = compare_beta(data, "group", permutations=99)
    assert r.data is not None and {"PC1", "PC2"} <= set(r.data.columns)
    ax = render(r)
    assert ax.has_data()
    assert "PC1" in ax.get_xlabel()          # axis label carries proportion explained


def test_render_heatmap_from_associate_cat_cat(data):
    ax = render(associate(data, "group", "bv"))
    assert ax.has_data()                     # contingency drawn as an image


def test_render_box_from_associate_cat_cont(data):
    ax = render(associate(data, "group", "ph"))
    assert ax.get_ylabel() == "ph"


def test_render_scatter_from_associate_cont_cont(data):
    ax = render(associate(data, "score", "ph"))
    assert ax.has_data()
    assert ax.get_xlabel() == "score" and ax.get_ylabel() == "ph"


def test_render_volcano_from_compare_abundance(data):
    ax = render(compare_abundance(data, "group"))
    assert ax.has_data()
    assert ax.get_ylabel() == "-log10(q)"


def test_result_data_serializes(data):
    d = compare_alpha(data, "group").to_dict()
    assert d["data"] is not None and "alpha_shannon" in d["data"]


def test_render_unknown_kind_raises(data):
    r = compare_alpha(data, "group")
    r.plot["kind"] = "sankey"
    with pytest.raises(ValueError, match="No renderer"):
        render(r)
