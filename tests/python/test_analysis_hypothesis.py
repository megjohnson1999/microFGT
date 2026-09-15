"""Hypothesis-test verbs: alpha (group test / OLS / mixed model) and beta (PERMANOVA).

Uses an in-memory MuData with a deliberately separated structure (no R, no network): a
``lowdiv`` group dominated by one taxon vs a ``highdiv`` even group, so alpha diversity and
community composition both differ by construction. Assertions target the direction and the
result contract, not exact permutation p-values.
"""

import numpy as np
import pandas as pd
import pytest

import anndata as ad
import mudata as md

from microfgt.analysis import compare_alpha, compare_beta
from microfgt.analysis.results import AnalysisResult


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    n_per, n_taxa = 12, 5
    low = rng.multinomial(200, [0.86, 0.05, 0.04, 0.03, 0.02], size=n_per)   # concentrated
    high = rng.multinomial(200, [0.2, 0.2, 0.2, 0.2, 0.2], size=n_per)       # even
    counts = np.vstack([low, high]).astype(np.int64)
    samples = [f"s{i:02d}" for i in range(2 * n_per)]
    group = ["lowdiv"] * n_per + ["highdiv"] * n_per
    taxon = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=[f"taxon{j}" for j in range(n_taxa)]),
    )
    taxon.layers["counts"] = counts
    # Clinical metadata rides on a SEPARATE modality's obs (like phyloseq sample_data on
    # composition.obs) to exercise cross-modality predictor resolution.
    comp = taxon.copy()
    comp.obs["group"] = group
    comp.obs["subject"] = [f"subj{i % n_per}" for i in range(2 * n_per)]   # paired across groups
    comp.obs["score"] = np.linspace(0.1, 0.9, 2 * n_per)
    with md.set_options(pull_on_update=False):
        mdata = md.MuData({"composition": comp, "composition_taxon": taxon})
    return mdata


# --- alpha ----------------------------------------------------------------------------------
def test_compare_alpha_auto_picks_mannwhitney_for_two_groups(data):
    r = compare_alpha(data, "group")                    # predictor lives on composition.obs
    assert isinstance(r, AnalysisResult)
    assert r.stats["test"] == "Mann–Whitney U"       # 2 groups -> Mann–Whitney
    assert r.pvalue < 0.01                            # groups differ strongly by construction
    # lowdiv really is lower-diversity than highdiv.
    assert r.table.loc["lowdiv", "median"] < r.table.loc["highdiv", "median"]
    assert r.notes["n_used"] == 24


def test_compare_alpha_kruskal_when_forced(data):
    r = compare_alpha(data, "group", test="kruskal")
    assert r.stats["test"] == "Kruskal–Wallis"
    assert "epsilon_squared" == r.stats["effect"]


def test_compare_alpha_ols_with_covariate(data):
    r = compare_alpha(data, ["group", "score"], test="ols")
    assert r.stats["test"] == "OLS"
    assert "R2" in r.stats
    # a categorical group term and the continuous covariate both appear
    assert any("group" in t for t in r.table.index)
    assert any("score" in t for t in r.table.index)


def test_compare_alpha_mixedlm_with_subject(data):
    r = compare_alpha(data, "group", subject="subject")
    assert r.stats["test"] == "LinearMixedModel"
    assert r.spec["subject"] == "subject"


def test_compare_alpha_subset_restricts_samples(data):
    r = compare_alpha(data, "group", subset={"group": ["lowdiv", "highdiv"]})
    assert r.notes["n_used"] == 24
    r2 = compare_alpha(data, "score", subset="group == 'lowdiv'", test="ols")
    assert r2.notes["n_used"] == 12


# --- beta -----------------------------------------------------------------------------------
def test_compare_beta_permanova_separates_groups(data):
    r = compare_beta(data, "group", permutations=199)
    assert r.stats["test"] == "PERMANOVA"
    assert r.pvalue <= 0.05
    assert r.stats["R2"] > 0.3                        # strong between-group separation
    assert "dispersion_pvalue" in r.stats
    assert r.notes["n_used"] == 24


def test_compare_beta_flags_extra_predictors_as_single_factor(data):
    r = compare_beta(data, ["group", "score"], permutations=99)
    assert r.notes["ignored_predictors"] == ["score"]
    assert "adonis2" in r.notes["warning"]


def test_compare_beta_needs_two_populated_groups(data):
    with pytest.raises(ValueError, match=">=2 groups"):
        compare_beta(data, "group", subset={"group": "lowdiv"}, permutations=99)


# --- contract -------------------------------------------------------------------------------
def test_result_summary_and_serialization(data):
    r = compare_alpha(data, "group")
    assert "compare_alpha" in r.summary() and "p=" in r.summary()
    d = r.to_dict()
    assert d["verb"] == "compare_alpha" and "table" in d and d["stats"]["test"] == "Mann–Whitney U"
