"""VIRGO importer — validated against the GENUINE fixtures virgo_sub1/sub2.out."""

import numpy as np
import pandas as pd

from microfgt.io import import_virgo


def test_virgo_stacks_real_per_sample_files(real_fixtures, tmp_path):
    # FORMATS.md: one <sample>.out per sample. Stage the two genuine files in a dir.
    for name in ("virgo_sub1.out", "virgo_sub2.out"):
        (tmp_path / name).write_bytes((real_fixtures / name).read_bytes())

    adata = import_virgo(tmp_path)

    # obs = samples (filenames minus .out), var = genes.
    assert list(adata.obs_names) == ["virgo_sub1", "virgo_sub2"]

    sub1 = pd.read_csv(real_fixtures / "virgo_sub1.out", sep="\t", header=None)
    sub2 = pd.read_csv(real_fixtures / "virgo_sub2.out", sep="\t", header=None)
    union = set(sub1[0]) | set(sub2[0])
    assert adata.n_vars == len(union)            # union of genes
    assert adata.n_vars < len(sub1) + len(sub2)  # the two samples share genes

    # Counts round-trip for a known gene present in sub1 (FORMATS.md: V1593031 = 1417).
    gi = list(adata.var_names).index("V1593031")
    si = list(adata.obs_names).index("virgo_sub1")
    assert adata.layers["counts"][si, gi] == 1417
    assert adata.var.loc["V1593031", "gene_length"] == 3663


def test_virgo_zero_fills_genes_absent_from_a_sample(real_fixtures, tmp_path):
    for name in ("virgo_sub1.out", "virgo_sub2.out"):
        (tmp_path / name).write_bytes((real_fixtures / name).read_bytes())
    adata = import_virgo(tmp_path)

    sub1_genes = set(pd.read_csv(real_fixtures / "virgo_sub1.out", sep="\t", header=None)[0])
    sub2_genes = set(pd.read_csv(real_fixtures / "virgo_sub2.out", sep="\t", header=None)[0])
    only_sub1 = (sub1_genes - sub2_genes).pop()

    si2 = list(adata.obs_names).index("virgo_sub2")
    gi = list(adata.var_names).index(only_sub1)
    assert adata.layers["counts"][si2, gi] == 0  # zero-filled, not dropped
    # No double counting: total equals the sum of both files' read counts.
    total = (
        pd.read_csv(real_fixtures / "virgo_sub1.out", sep="\t", header=None)[1].sum()
        + pd.read_csv(real_fixtures / "virgo_sub2.out", sep="\t", header=None)[1].sum()
    )
    assert adata.layers["counts"].sum() == total


def test_virgo_optional_catalog_taxon_join(real_fixtures, tmp_path):
    # The genuine per-sample outputs and the TRIMMED catalog fixture have zero gene
    # overlap (catalog = V1000001.., sample genes = V1593031..), so to validate the
    # join itself we stage a per-sample file using real catalog gene ids.
    catalog = pd.read_csv(real_fixtures / "virgo_taxon.tbl.txt", sep="\t", header=None)
    g0, taxon0 = catalog.iloc[0, 1], catalog.iloc[0, 2]   # e.g. V1000001 -> L. iners
    (tmp_path / "sampleA.out").write_text(f"{g0}\t10\t396\nV9999999\t5\t100\n")

    adata = import_virgo(tmp_path, taxon_table=real_fixtures / "virgo_taxon.tbl.txt")
    assert "taxon" in adata.var.columns
    assert adata.var.loc[g0, "taxon"] == taxon0           # catalog gene annotated
    assert adata.var.loc["V9999999", "taxon"] == "Unannotated"  # gene not in catalog
