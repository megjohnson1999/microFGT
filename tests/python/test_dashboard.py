"""Dashboard core: variable introspection (catalog) + the shared verb dispatcher.

The Streamlit UI itself isn't unit-tested (it's a thin shell of widgets); the logic it stands
on — typing obs columns for the dropdowns, and dispatching a spec to the right verb — is pure
and tested here.
"""

import numpy as np
import pandas as pd
import pytest

import anndata as ad
import mudata as md

from microfgt.analysis import VERBS, run_verb
from microfgt.analysis.results import AnalysisResult
from microfgt.dashboard import as_frame, continuous, groupable, variable_catalog


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    n_per, n_taxa = 10, 5
    counts = np.vstack([
        rng.multinomial(300, [0.7, 0.1, 0.1, 0.05, 0.05], size=n_per),
        rng.multinomial(300, [0.2, 0.2, 0.2, 0.2, 0.2], size=n_per),
    ]).astype(np.int64)
    samples = [f"s{i:02d}" for i in range(2 * n_per)]
    taxon = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=[f"taxon{j}" for j in range(n_taxa)]),
    )
    taxon.layers["counts"] = counts
    comp = taxon.copy()
    comp.obs["CST"] = ["I"] * n_per + ["IV"] * n_per          # categorical, on a modality obs
    comp.obs["score"] = rng.random(2 * n_per)                 # continuous
    comp.obs["age"] = rng.normal(30, 5, 2 * n_per)            # continuous
    comp.obs["note"] = ["x"] + [None] * (2 * n_per - 1)       # mostly missing categorical
    with md.set_options(pull_on_update=False):
        return md.MuData({"composition": comp, "composition_taxon": taxon})


# --- catalog --------------------------------------------------------------------------------
def test_catalog_types_variables_from_any_modality(data):
    cat = {v.name: v for v in variable_catalog(data)}
    assert cat["CST"].kind == "categorical" and cat["CST"].n_levels == 2
    assert set(cat["CST"].levels) == {"I", "IV"}
    assert cat["score"].kind == "continuous"
    assert cat["note"].n_missing == 19 and cat["note"].n_present == 1


def test_groupable_and_continuous_helpers(data):
    cat = variable_catalog(data)
    assert "CST" in groupable(cat)          # 2 levels -> usable grouping
    assert "score" not in groupable(cat)    # continuous
    assert "score" in continuous(cat)


def test_as_frame_overview(data):
    frame = as_frame(variable_catalog(data))
    assert set(frame.columns) == {"kind", "levels", "present", "missing"}
    assert "CST" in frame.index


def test_catalog_dedupes_mudata_modality_prefixes(data, tmp_path):
    # On write/read, MuData copies modality obs into global .obs under a 'composition:' prefix.
    # The catalog must present each variable ONCE under its clean name, not twice.
    p = tmp_path / "obj.h5mu"
    data.write(p)
    reloaded = md.read(p)
    names = [v.name for v in variable_catalog(reloaded)]
    assert not any(":" in n for n in names)          # no 'composition:...' leakage
    assert names.count("CST") == 1
    assert names.count("score") == 1


# --- dispatch -------------------------------------------------------------------------------
def test_run_verb_dispatches_each_verb(data):
    assert set(VERBS) == {"alpha", "beta", "associate", "abundance"}
    assert run_verb(data, "alpha", predictors="CST").verb == "compare_alpha"
    assert run_verb(data, "beta", predictors="CST", metric="braycurtis").verb == "compare_beta"
    assert run_verb(data, "associate", x="CST", y="score").verb == "associate"
    r = run_verb(data, "abundance", predictors="CST")
    assert isinstance(r, AnalysisResult) and r.verb == "compare_abundance"


def test_run_verb_validates_inputs(data):
    with pytest.raises(ValueError, match="needs predictors"):
        run_verb(data, "alpha")
    with pytest.raises(ValueError, match="needs both x and y"):
        run_verb(data, "associate", x="CST")
    with pytest.raises(ValueError, match="Unknown verb"):
        run_verb(data, "sankey", predictors="CST")


def test_run_verb_passes_subset_through(data):
    # cont x cont association tolerates a single-group subset (no >=2-group requirement).
    r = run_verb(data, "associate", x="score", y="age", subset={"CST": ["I"]})
    assert r.notes["n_used"] == 10
