"""``compare_abundance`` — covariate-adjusted, FDR-corrected differential abundance.

In-memory MuData with two taxa planted to move (taxon0 up in group B, taxon1 down) against a
flat background, plus a covariate and paired subjects. Assertions target the planted signal,
covariate handling, and the result contract — not exact p-values. Pure scikit-bio (ANCOM-BC /
Dirichlet-multinomial LME); no R.
"""

import numpy as np
import pandas as pd
import pytest

import anndata as ad
import mudata as md

from microfgt.analysis import compare_abundance
from microfgt.analysis.results import AnalysisResult


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    n_per, n_taxa = 15, 6
    base = [0.25, 0.25, 0.125, 0.125, 0.125, 0.125]
    a = rng.multinomial(1000, base, size=n_per).astype(np.int64)
    b_props = [0.45, 0.05, 0.125, 0.125, 0.125, 0.125]          # taxon0 up, taxon1 down in B
    b = rng.multinomial(1000, b_props, size=n_per).astype(np.int64)
    counts = np.vstack([a, b])
    # a genuinely absent feature (all zero) to exercise the zero-feature drop
    counts = np.hstack([counts, np.zeros((2 * n_per, 1), dtype=np.int64)])
    samples = [f"s{i:02d}" for i in range(2 * n_per)]
    taxa = [f"taxon{j}" for j in range(n_taxa)] + ["absent"]

    taxon = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=taxa),
    )
    taxon.layers["counts"] = counts
    comp = taxon.copy()
    comp.obs["group"] = ["A"] * n_per + ["B"] * n_per
    comp.obs["subject"] = [f"subj{i % n_per}" for i in range(2 * n_per)]
    comp.obs["score"] = rng.normal(0, 1, 2 * n_per)
    with md.set_options(pull_on_update=False):
        mdata = md.MuData({"composition": comp, "composition_taxon": taxon})
    return mdata


def test_ancombc_recovers_planted_signal(data):
    r = compare_abundance(data, "group")
    assert isinstance(r, AnalysisResult)
    assert r.stats["test"] == "ANCOM-BC"
    # only the predictor-of-interest contrast rows survive (no Intercept)
    assert set(r.table["contrast"]) == {"group[T.B]"}
    assert "taxon0" in r.table.index and "taxon1" in r.table.index
    assert r.table.loc["taxon0", "log2fc"] > 0          # up in B
    assert r.table.loc["taxon1", "log2fc"] < 0          # down in B
    assert bool(r.table.loc["taxon0", "signif"])
    # the all-zero feature was dropped from testing, honestly counted
    assert "absent" not in r.table.index
    assert r.notes["n_features_dropped_zero"] == 1


def test_covariate_enters_formula_but_only_primary_reported(data):
    r = compare_abundance(data, ["group", "score"])
    assert r.stats["formula"] == "group + score"
    assert r.spec["covariates"] == ["score"]
    assert set(r.table["contrast"]) == {"group[T.B]"}   # score is adjusted-for, not reported


def test_table_sorted_by_qvalue(data):
    r = compare_abundance(data, "group")
    q = r.table["qvalue"].to_numpy()
    assert np.all(np.diff(q) >= 0)


def test_dirmult_lme_longitudinal_runs(data):
    r = compare_abundance(data, "group", method="dirmult_lme", subject="subject",
                          draws=16, seed=0)
    assert r.stats["test"] == "Dirichlet-multinomial LME"
    assert "taxon0" in r.table.index
    assert r.spec["subject"] == "subject"


def test_dirmult_lme_requires_subject(data):
    with pytest.raises(ValueError, match="subject"):
        compare_abundance(data, "group", method="dirmult_lme")


def test_unknown_method_raises(data):
    with pytest.raises(ValueError, match="Unknown method"):
        compare_abundance(data, "group", method="deseq2")
