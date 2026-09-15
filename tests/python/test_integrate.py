"""Integrated object (MuData) assembly + honest sample reconciliation."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microfgt.io import build_mudata, import_valencia, import_virgo
from microfgt.io.speciateit import import_speciateit


def _asv_composition(samples, asv_class):
    """ASV-grain composition (samples x ASVs) with a classification per ASV, as the real
    importer produces — so build_mudata materialises composition_taxon + descriptors."""
    asvs = list(asv_class)
    n = len(samples) * len(asvs)
    X = (np.arange(n, dtype=np.float32) + 1).reshape(len(samples), len(asvs))
    var = pd.DataFrame(
        {
            "classification": [asv_class[a] for a in asvs],
            "genus": [asv_class[a].split("_")[0] for a in asvs],
        },
        index=pd.Index(asvs, name="asv"),
    )
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(samples, name="sample")), var=var)
    adata.layers["counts"] = X.astype(np.int64)
    return adata


def _genes(samples, genes):
    X = np.arange(len(samples) * len(genes), dtype=np.float32).reshape(len(samples), len(genes))
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )


def _taxon_shotgun(samples, taxa):
    """Shotgun taxon x sample profile, as import_virgo2_taxonomy will produce: counts with a
    'counts' layer so describe_composition can read it."""
    X = (np.arange(len(samples) * len(taxa), dtype=np.float32) + 1).reshape(
        len(samples), len(taxa)
    )
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=pd.Index(taxa, name="taxon")),
    )
    adata.layers["counts"] = X.astype(np.int64)
    return adata


def test_build_mudata_materialises_taxon_assay_and_descriptors():
    comp = _asv_composition(
        ["s1", "s2", "s3"],
        {"ASV1": "Lactobacillus_iners", "ASV2": "Gardnerella_vaginalis"},
    )
    func = _genes(["s1", "s2", "s3"], ["geneA", "geneB", "geneC"])
    cst = pd.DataFrame(
        {"CST": ["I", "IV-A", "III"], "subCST": ["I-B", "IV-A", "III-A"], "score": [0.99, 0.7, 0.8]},
        index=pd.Index(["s1", "s2", "s3"], name="sample"),
    )

    mdata = build_mudata(composition=comp, function=func, cst=cst)

    # ASV-grain composition + materialised taxon roll-up + function.
    assert set(mdata.mod.keys()) == {"composition", "composition_taxon", "function"}
    assert mdata["composition"].n_vars == 2  # ASV grain preserved
    assert "Lactobacillus_iners" in list(mdata["composition_taxon"].var_names)

    # CST attached, and the augment descriptors alongside it.
    assert list(mdata.obs["CST"]) == ["I", "IV-A", "III"]
    for col in ("dominant_taxon", "dominance_pct", "effective_taxa"):
        assert col in mdata.obs

    recon = mdata.uns["reconciliation"]
    assert recon["shared_all"] == 3
    assert recon["cst_matched"] == 3
    assert recon["cst_unmatched"] == []


def test_reconciliation_reports_mismatch_without_dropping():
    # Honest reconciliation: a CST sample with no assay is reported, not silently dropped;
    # an assay sample with no CST gets NaN, not removed.
    comp = _asv_composition(["s1", "s2"], {"ASV1": "Lactobacillus_iners"})
    cst = pd.DataFrame(
        # s2 has an assay (matches); s_ghost has CST but no assay (unmatched).
        {"CST": ["IV-B", "I"], "subCST": ["IV-B", "I-A"], "score": [0.6, 0.9]},
        index=pd.Index(["s2", "s_ghost"], name="sample"),
    )

    mdata = build_mudata(composition=comp, cst=cst)

    assert mdata["composition"].n_obs == 2          # s1 kept despite no CST
    assert pd.isna(mdata.obs.loc["s1", "CST"])       # NaN, not dropped
    assert mdata.obs.loc["s2", "CST"] == "IV-B"
    recon = mdata.uns["reconciliation"]
    assert recon["cst_matched"] == 1
    assert recon["cst_unmatched"] == ["s_ghost"]     # surfaced, not hidden


def test_cst_sim_columns_routed_to_obsm_not_obs():
    # The 13 <subCST>_sim vectors clutter the sample frame; route them to the taxon assay's
    # obsm and keep only the labels + descriptors in .obs.
    comp = _asv_composition(
        ["s1", "s2"], {"ASV1": "Lactobacillus_iners", "ASV2": "Gardnerella_vaginalis"}
    )
    cst = pd.DataFrame(
        {
            "I-A_sim": [0.9, 0.1], "V_sim": [0.1, 0.8],
            "subCST": ["I-A", "V"], "score": [0.9, 0.8], "CST": ["I", "V"],
        },
        index=pd.Index(["s1", "s2"], name="sample"),
    )

    mdata = build_mudata(composition=comp, cst=cst)

    assert "CST" in mdata.obs and "subCST" in mdata.obs
    assert not any(str(c).endswith("_sim") for c in mdata.obs.columns)
    taxon = mdata["composition_taxon"]
    assert taxon.obsm["cst_sim"].shape == (2, 2)
    assert taxon.uns["cst_sim_columns"] == ["I-A_sim", "V_sim"]


def test_shotgun_only_object_is_complete_and_self_contained():
    # The shotgun arm writes a valid object with no 16S present: function + shotgun taxon
    # modalities, mgCST in .obs, shotgun_-prefixed descriptors, sane reconciliation.
    func = _genes(["m1", "m2"], ["geneA", "geneB", "geneC"])
    taxon_sg = _taxon_shotgun(
        ["m1", "m2"], ["Lactobacillus_crispatus", "Gardnerella_vaginalis"]
    )
    mgcst = pd.DataFrame(
        {"mgCST": ["mgCST 1", "mgCST 12"], "mgCST_score": [0.9, 0.55]},
        index=pd.Index(["m1", "m2"], name="sample"),
    )

    mdata = build_mudata(function=func, composition_taxon_shotgun=taxon_sg, mgcst=mgcst)

    assert set(mdata.mod.keys()) == {"function", "composition_taxon_shotgun"}
    # mgCST labels on the global frame; no 16S CST present.
    assert list(mdata.obs["mgCST"]) == ["mgCST 1", "mgCST 12"]
    assert list(mdata.obs["mgCST_score"]) == [0.9, 0.55]
    assert "CST" not in mdata.obs
    # Shotgun descriptors present, prefixed so they never collide with 16S descriptors.
    for col in ("shotgun_dominant_taxon", "shotgun_dominance_pct", "shotgun_effective_taxa"):
        assert col in mdata.obs
    assert "dominant_taxon" not in mdata.obs   # unprefixed = 16S only, absent here


def test_mgcst_and_cst_coexist_as_separate_obs_columns():
    # 16S + shotgun both present: CST and mgCST are distinct .obs columns, never merged. CST
    # still routes its per-centroid sims to the 16S taxon assay; mgCST has no sim block (VISTA
    # emits only best-match θ), so nothing lands on the shotgun taxon assay's obsm.
    comp = _asv_composition(
        ["s1", "s2"], {"ASV1": "Lactobacillus_iners", "ASV2": "Gardnerella_vaginalis"}
    )
    taxon_sg = _taxon_shotgun(["s1", "s2"], ["Prevotella_bivia", "Sneathia_amnii"])
    cst = pd.DataFrame(
        {"I-A_sim": [0.9, 0.1], "V_sim": [0.1, 0.8],
         "subCST": ["I-A", "V"], "score": [0.9, 0.8], "CST": ["I", "V"]},
        index=pd.Index(["s1", "s2"], name="sample"),
    )
    mgcst = pd.DataFrame(
        {"mgCST": ["mgCST 1", "mgCST 23"], "mgCST_score": [0.7, 0.8]},
        index=pd.Index(["s1", "s2"], name="sample"),
    )

    mdata = build_mudata(
        composition=comp, composition_taxon_shotgun=taxon_sg, cst=cst, mgcst=mgcst
    )

    # Both community-type calls live as separate labels; neither is merged into the other.
    assert list(mdata.obs["CST"]) == ["I", "V"]
    assert list(mdata.obs["mgCST"]) == ["mgCST 1", "mgCST 23"]
    # No sim vectors clutter the sample frame.
    assert not any(str(c).endswith("_sim") for c in mdata.obs.columns)
    # cst_sim on the 16S taxon assay; the shotgun taxon assay carries no mgcst_sim.
    assert mdata["composition_taxon"].obsm["cst_sim"].shape == (2, 2)
    assert "mgcst_sim" not in mdata["composition_taxon_shotgun"].obsm


def test_build_mudata_requires_a_modality():
    with pytest.raises(ValueError):
        build_mudata()


def test_real_fixtures_are_not_a_co_assayed_cohort(real_fixtures, test_data):
    # Documents the honest-scope point: the three real fixtures come from different repos
    # and their sample ids DO NOT overlap, so a true cross-modality MuData can't be built
    # from them. The importers still each validate; integration needs a co-assayed dataset.
    comp = import_speciateit(
        test_data / "speciateit_MC_order7_results.synthetic.txt",
        real_fixtures / "speciateit_test_count_table.csv",
    )
    import tempfile, pathlib
    d = pathlib.Path(tempfile.mkdtemp())
    for name in ("virgo_sub1.out", "virgo_sub2.out"):
        (d / name).write_bytes((real_fixtures / name).read_bytes())
    func = import_virgo(d)
    cst = import_valencia(real_fixtures / "valencia_genuine_output_head.csv")

    comp_ids, func_ids, cst_ids = set(comp.obs_names), set(func.obs_names), set(cst.index)
    assert comp_ids.isdisjoint(func_ids)   # sample1.. vs sub1/sub2
    assert comp_ids.isdisjoint(cst_ids)    # sample1.. vs 1,2,..

    mdata = build_mudata(composition=comp, function=func, cst=cst)
    assert "composition_taxon" in mdata.mod          # roll-up materialised
    assert mdata.uns["reconciliation"]["shared_all"] == 0  # reported, honestly
