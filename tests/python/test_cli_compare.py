"""`microfgt compare` — run a hypothesis-test verb on a .h5mu from the terminal."""

import numpy as np
import pandas as pd
import pytest

import anndata as ad
import mudata as md

from microfgt.cli import main


@pytest.fixture
def h5mu(tmp_path):
    rng = np.random.default_rng(0)
    n_per, n_taxa = 10, 5
    low = rng.multinomial(300, [0.7, 0.1, 0.1, 0.05, 0.05], size=n_per)
    high = rng.multinomial(300, [0.2, 0.2, 0.2, 0.2, 0.2], size=n_per)
    counts = np.vstack([low, high]).astype(np.int64)
    samples = [f"s{i:02d}" for i in range(2 * n_per)]
    taxon = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame({"group": ["A"] * n_per + ["B"] * n_per},
                         index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=[f"taxon{j}" for j in range(n_taxa)]),
    )
    taxon.layers["counts"] = counts
    with md.set_options(pull_on_update=False):
        mdata = md.MuData({"composition_taxon": taxon})
    mdata.obs["group"] = ["A"] * n_per + ["B"] * n_per
    path = tmp_path / "obj.h5mu"
    mdata.write(path)
    return path


def test_compare_alpha_writes_table_and_plot(h5mu, tmp_path, capsys):
    out, fig = tmp_path / "res.csv", tmp_path / "fig.png"
    rc = main(["compare", "-i", str(h5mu), "--verb", "alpha", "--predictors", "group",
               "-o", str(out), "--plot", str(fig)])
    assert rc == 0
    assert "compare_alpha" in capsys.readouterr().out
    assert out.exists() and fig.exists()
    table = pd.read_csv(out, index_col=0)
    assert set(table.index) >= {"A", "B"}


def test_compare_beta_runs(h5mu, capsys):
    rc = main(["compare", "-i", str(h5mu), "--verb", "beta", "--predictors", "group"])
    assert rc == 0
    assert "PERMANOVA" in capsys.readouterr().out


def test_compare_requires_predictors(h5mu):
    with pytest.raises(SystemExit):
        main(["compare", "-i", str(h5mu), "--verb", "alpha"])


def test_compare_associate_needs_x_and_y(h5mu):
    with pytest.raises(SystemExit):
        main(["compare", "-i", str(h5mu), "--verb", "associate"])
