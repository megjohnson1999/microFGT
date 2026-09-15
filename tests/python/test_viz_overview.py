"""FGT overview plots for the dashboard Home tab (Agg backend, no display)."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402

from microfgt.viz.overview import (  # noqa: E402
    community_stack,
    cst_bar,
    dominance_hist,
    ordination,
    overview_figures,
)


@pytest.fixture
def mdata():
    rng = np.random.default_rng(0)
    n, t = 40, 8
    counts = rng.multinomial(300, [0.5, 0.2, 0.1, 0.05, 0.05, 0.04, 0.03, 0.03], size=n).astype(np.int64)
    samples = [f"s{i:02d}" for i in range(n)]
    taxon = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(
            {"CST": ["I"] * 20 + ["IV"] * 20,
             "effective_taxa": rng.uniform(1, 8, n),
             "dominance_pct": rng.uniform(20, 95, n),
             "read_count": counts.sum(1)},
            index=pd.Index(samples, name="sample"),
        ),
        var=pd.DataFrame(index=[f"taxon{j}" for j in range(t)]),
    )
    taxon.layers["counts"] = counts
    with md.set_options(pull_on_update=False):
        return md.MuData({"composition_taxon": taxon})


def test_overview_figures_builds_the_expected_set(mdata):
    figs = overview_figures(mdata)
    assert {"cst", "community", "diversity", "dominance", "ordination"} <= set(figs)
    for fig in figs.values():
        assert fig.axes and fig.axes[0].has_data()


def test_individual_plots_return_populated_axes(mdata):
    assert cst_bar(mdata).has_data()
    assert community_stack(mdata, top_n=4).get_ylabel() == "relative abundance"
    assert dominance_hist(mdata).has_data()
    assert "PCoA1" in ordination(mdata, max_samples=30).get_xlabel()


def test_overview_skips_missing_columns():
    # No CST -> cst/diversity plots skipped; community + dominance + ordination still build.
    rng = np.random.default_rng(1)
    counts = rng.multinomial(200, [0.6, 0.3, 0.1], size=12).astype(np.int64)
    taxon = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame({"dominance_pct": rng.uniform(20, 90, 12)},
                         index=pd.Index([f"s{i}" for i in range(12)], name="sample")),
        var=pd.DataFrame(index=["a", "b", "c"]),
    )
    taxon.layers["counts"] = counts
    with md.set_options(pull_on_update=False):
        m = md.MuData({"composition_taxon": taxon})
    figs = overview_figures(m)
    assert "cst" not in figs and "diversity" not in figs
    assert "community" in figs and "dominance" in figs
