"""VIRGO2 importers, grounded in the real ENA/PRJEB34536 fixtures.

Fixtures (public ENA cohort, safe to ship):
  * virgo2_compiled.summary.NR.slice.txt  — the compiled gene x sample matrix
  * virgo2_taxon_annotation.slice.txt     — AnnotationTables/1.VIRGO2.taxon.txt (join on Gene)
  * virgo2_kegg_annotation.slice.txt      — a functional table, ~2/3 of genes annotated
"""

import numpy as np
import pandas as pd
import pytest

from microfgt.io import (
    collapse_virgo2_to_taxon,
    import_virgo2,
    import_virgo2_taxonomy,
)

SUMMARY = "virgo2_compiled.summary.NR.slice.txt"
TAXON = "virgo2_taxon_annotation.slice.txt"
KEGG = "virgo2_kegg_annotation.slice.txt"


def test_import_virgo2_orients_gene_matrix_to_samples_x_genes(real_fixtures):
    raw = pd.read_csv(real_fixtures / SUMMARY, sep="\t", index_col=0)  # genes x samples
    adata = import_virgo2(real_fixtures / SUMMARY)

    # MuData convention: obs = samples, var = genes (the transpose of the on-disk matrix).
    assert adata.n_obs == raw.shape[1]
    assert adata.n_vars == raw.shape[0]
    assert list(adata.obs_names) == [str(c) for c in raw.columns]
    assert "ERR4421550" in adata.obs_names

    # Counts survive round-trip (fractional floats, not ints); a known top cell.
    assert "counts" in adata.layers
    val = adata[["ERR4421550"], ["CIG00063_0001_787"]].layers["counts"][0, 0]
    assert val == pytest.approx(5609.0)


def test_import_virgo2_joins_taxon_and_tolerates_partial_functional_annotation(real_fixtures):
    adata = import_virgo2(
        real_fixtures / SUMMARY,
        taxon_annotation=real_fixtures / TAXON,
        annotations={"kegg": real_fixtures / KEGG},
    )
    # Taxon join on Gene: a gene present in both files gets its Taxa; genes absent -> Unannotated.
    assert "taxon" in adata.var
    assert adata.var.loc["CIG00063_0001_787", "taxon"] == "Lactobacillus_iners"

    # Functional table joined under a name prefix, tolerating genes it does not cover.
    assert "kegg_KEGG" in adata.var
    kegg = pd.read_csv(real_fixtures / KEGG, sep="\t", dtype=str).set_index("Gene")["KEGG"]
    for gene in adata.var_names:
        expected = kegg.get(gene)
        got = adata.var.loc[gene, "kegg_KEGG"]
        if expected is None or pd.isna(expected):
            assert pd.isna(got)                    # missing annotation -> NaN, not a crash
        else:
            assert got == expected
    assert adata.var["kegg_KEGG"].isna().any()     # the slice genuinely has unannotated genes


def test_collapse_virgo2_to_taxon_conserves_counts(real_fixtures):
    func = import_virgo2(real_fixtures / SUMMARY, taxon_annotation=real_fixtures / TAXON)
    taxon = collapse_virgo2_to_taxon(func)

    # Same samples; counts conserved per sample (no gene dropped or double-counted).
    assert list(taxon.obs_names) == list(func.obs_names)
    np.testing.assert_allclose(
        taxon.layers["counts"].sum(axis=1), func.layers["counts"].sum(axis=1)
    )
    # A real taxon lands, rank-comparable with 16S (genus = first token).
    assert "Lactobacillus_iners" in list(taxon.var_names)
    assert taxon.var.loc["Lactobacillus_iners", "genus"] == "Lactobacillus"


def test_import_virgo2_taxonomy_is_the_read_plus_collapse_convenience(real_fixtures):
    one_shot = import_virgo2_taxonomy(real_fixtures / SUMMARY, real_fixtures / TAXON)
    manual = collapse_virgo2_to_taxon(
        import_virgo2(real_fixtures / SUMMARY, taxon_annotation=real_fixtures / TAXON)
    )
    assert list(one_shot.var_names) == list(manual.var_names)
    np.testing.assert_array_equal(one_shot.layers["counts"], manual.layers["counts"])


def test_collapse_requires_taxon_annotation(real_fixtures):
    func = import_virgo2(real_fixtures / SUMMARY)  # no taxon_annotation
    with pytest.raises(ValueError, match="var\\['taxon'\\]"):
        collapse_virgo2_to_taxon(func)
