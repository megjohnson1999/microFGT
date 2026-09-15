"""Augment descriptors — intrinsic (stored) + the parameterized on-demand threshold count."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microfgt.characterize import describe_composition, taxa_over_threshold


def _taxon_assay(counts_by_sample):
    """counts_by_sample: {sample: {taxon: count}} -> taxon x sample AnnData."""
    df = pd.DataFrame(counts_by_sample).T.fillna(0.0)  # samples x taxa
    adata = ad.AnnData(
        X=df.to_numpy().astype(np.float32),
        obs=pd.DataFrame(index=pd.Index(df.index, name="sample")),
        var=pd.DataFrame(index=pd.Index(df.columns, name="taxon")),
    )
    adata.layers["counts"] = df.to_numpy().astype(np.int64)
    return adata


# --- intrinsic descriptors (stored on the object) -------------------------------------------
def test_dominance_and_effective_taxa():
    comp = _taxon_assay(
        {
            "dominated": {"Lactobacillus_iners": 90, "Gardnerella_vaginalis": 10},
            "diffuse": {"a": 25, "b": 25, "c": 25, "d": 25},   # four even taxa
        }
    )
    out = describe_composition(comp)

    assert out.loc["dominated", "dominant_taxon"] == "Lactobacillus_iners"
    assert out.loc["dominated", "dominance_pct"] == 90.0
    # effective number of taxa is cutoff-free: ~1.4 for the dominated sample, exactly 4 for the
    # four-even diffuse one (exp(Shannon)).
    assert out.loc["dominated", "effective_taxa"] == pytest.approx(1.38, abs=0.05)
    assert out.loc["diffuse", "effective_taxa"] == pytest.approx(4.0, abs=1e-6)
    assert "n_taxa_over_10pct" not in out.columns          # no baked-in cutoff column


def test_unclassified_bucket_counts_as_a_taxon():
    comp = _taxon_assay({"s": {"Unclassified": 80, "Lactobacillus_iners": 20}})
    out = describe_composition(comp)
    assert out.loc["s", "dominant_taxon"] == "Unclassified"
    assert out.loc["s", "dominance_pct"] == 80.0


def test_zero_count_sample_is_undefined_not_forced():
    comp = _taxon_assay({"empty": {"a": 0, "b": 0}, "ok": {"a": 10, "b": 0}})
    out = describe_composition(comp)
    assert pd.isna(out.loc["empty", "dominant_taxon"])
    assert out.loc["empty", "dominance_pct"] == 0.0
    assert pd.isna(out.loc["empty", "effective_taxa"])     # undefined, not 1
    assert out.loc["ok", "effective_taxa"] == pytest.approx(1.0)


# --- parameterized descriptor (computed on demand at a chosen cutoff) ------------------------
def test_taxa_over_threshold_is_adjustable():
    comp = _taxon_assay(
        {
            "dominated": {"Lactobacillus_iners": 90, "Gardnerella_vaginalis": 10},
            "diffuse": {"a": 25, "b": 25, "c": 25, "d": 25},
        }
    )
    # the cutoff is a knob, not a fixed constant — the count responds to it
    at10 = taxa_over_threshold(comp, 0.10)
    assert at10.loc["dominated"] == 1 and at10.loc["diffuse"] == 4
    at30 = taxa_over_threshold(comp, 0.30)
    assert at30.loc["dominated"] == 1 and at30.loc["diffuse"] == 0   # 0.25 no longer clears 0.30
    assert at10.name == "taxa_over_threshold"


def test_taxa_over_threshold_zero_count_is_zero():
    comp = _taxon_assay({"empty": {"a": 0, "b": 0}, "ok": {"a": 10, "b": 5}})
    assert taxa_over_threshold(comp, 0.10).loc["empty"] == 0
