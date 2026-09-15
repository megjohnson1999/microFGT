"""phyloseq importer (ASV grain) — the phyloseq entry point.

Exercises the real bundled-R-glue path (Rscript + phyloseq_export.R) against a tiny fixture
phyloseq built in R (``data/phyloseq_tiny.rds``): 3 ASVs x 4 samples, taxa_names ARE the
sequences, tax_table ranks Domain..Genus_Species (one taxon with a missing Genus_Species to
exercise the fallback), sample_data carrying CST/subCST/score. Skips cleanly where R or
phyloseq is unavailable, so the suite still runs without the R toolchain.
"""

import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from microfgt.io import (
    build_mudata,
    collapse_to_taxon,
    existing_cst,
    import_phyloseq,
)

FIXTURE = "phyloseq_tiny.rds"


def _has_phyloseq() -> bool:
    if shutil.which("Rscript") is None:
        return False
    try:
        rc = subprocess.run(
            ["Rscript", "-e", "suppressMessages(library(phyloseq))"],
            capture_output=True, timeout=120,
        ).returncode
    except (OSError, subprocess.SubprocessError):
        return False
    return rc == 0


pytestmark = pytest.mark.skipif(
    not _has_phyloseq(), reason="Rscript + phyloseq not available"
)


@pytest.fixture
def composition(test_data):
    return import_phyloseq(test_data / FIXTURE)


# --- import: ASV grain, sequences + classification retained ---------------------------------
def test_import_is_asv_grain_samples_x_asvs(composition):
    # 3 ASVs -> features; 4 samples -> obs; one feature per ASV (no read-time collapse).
    assert composition.n_vars == 3
    assert composition.n_obs == 4
    assert list(composition.obs_names) == ["s1", "s2", "s3", "s4"]
    # counts oriented samples x ASVs and conserved (ASV1 column = 90,5,40,10).
    np.testing.assert_array_equal(
        composition.layers["counts"][:, 0], np.array([90, 5, 40, 10])
    )
    assert composition.layers["counts"].sum() == 90 + 5 + 40 + 10 + 5 + 85 + 30 + 10 + 5 + 10 + 30 + 80


def test_sequences_are_the_taxa_names(composition):
    # No refseq slot: var['sequence'] comes straight from taxa_names.
    assert composition.var.iloc[0]["sequence"] == "ACGTACGTACGTAAAA"
    assert composition.var["sequence"].notna().all()


def test_classification_from_genus_species_with_fallback(composition):
    cls = composition.var["classification"].tolist()
    # ASV1/ASV2 take Genus_Species directly.
    assert cls[0] == "Lactobacillus_crispatus"
    assert cls[1] == "Lactobacillus_iners"
    # ASV3 has a missing Genus_Species AND Species -> falls back to Genus.
    assert cls[2] == "Gardnerella"
    assert composition.var["genus"].tolist() == ["Lactobacillus", "Lactobacillus", "Gardnerella"]


def test_obs_carries_all_sample_data(composition):
    for col in ("sampleID", "sample_type", "run_label", "CST", "subCST", "score"):
        assert col in composition.obs.columns
    assert composition.obs.loc["s2", "CST"] == "III"
    # read_count = per-sample total.
    assert composition.obs.loc["s1", "read_count"] == 100


# --- CST surfacing + downstream parity ------------------------------------------------------
def test_existing_cst_is_valencia_shaped(composition):
    cst = existing_cst(composition)
    assert list(cst.columns) == ["CST", "subCST", "score"]
    assert cst.loc["s3", "subCST"] == "IV-B0"
    assert list(cst.index) == ["s1", "s2", "s3", "s4"]


def test_collapse_and_describe_run_downstream(composition):
    taxon = collapse_to_taxon(composition)
    # Three distinct classifications -> three taxa; counts conserved.
    assert set(taxon.var_names) == {"Lactobacillus_crispatus", "Lactobacillus_iners", "Gardnerella"}
    assert taxon.layers["counts"].sum() == composition.layers["counts"].sum()

    from microfgt.characterize import describe_composition

    desc = describe_composition(taxon)
    # s1 is 90% ASV1 (crispatus); s4 is 80% Gardnerella.
    assert desc.loc["s1", "dominant_taxon"] == "Lactobacillus_crispatus"
    assert desc.loc["s1", "dominance_pct"] == pytest.approx(90.0)
    assert desc.loc["s4", "dominant_taxon"] == "Gardnerella"


def test_build_mudata_attaches_existing_cst(composition):
    cst = existing_cst(composition)
    mdata = build_mudata(composition=composition, cst=cst)
    # composition_taxon auto-materialised; CST + descriptors on the global obs.
    assert "composition_taxon" in mdata.mod
    assert mdata.obs.loc["s2", "CST"] == "III"
    assert "dominant_taxon" in mdata.obs.columns
    assert mdata.uns["reconciliation"]["cst_matched"] == 4


# --- wired as a real entry point: in-memory workflow + the file-artifact CLI ladder ---------
def test_run_workflow_enters_at_phyloseq(test_data):
    from microfgt.workflow import run_workflow

    mdata = run_workflow({
        "composition": {"phyloseq": str(test_data / FIXTURE)},
        "analysis": {"transforms": ["relabund"]},
    })
    assert set(mdata.mod) >= {"composition", "composition_taxon"}
    # The object's own CST is trusted by default (no method/valencia set).
    assert mdata.obs.loc["s2", "CST"] == "III"
    assert "dominant_taxon" in mdata.obs.columns
    assert "relabund" in mdata["composition"].layers


def test_cli_run_ladder_from_phyloseq(test_data, tmp_path):
    import mudata as md
    import yaml

    from microfgt.cli import main

    out = tmp_path / "out.h5mu"
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(yaml.safe_dump({
        "composition": {"phyloseq": str(test_data / FIXTURE)},
        "output": str(out),
    }))

    assert main(["run", "-c", str(cfg), "--workdir", str(tmp_path / "wd")]) == 0
    assert out.exists()
    # The file-artifact backbone ran the resolved entry-point plan.
    assert (tmp_path / "wd" / "composition.h5ad").exists()
    assert (tmp_path / "wd" / "cst.csv").exists()

    m = md.read(out)
    assert m["composition"].n_vars == 3
    assert "sequence" in m["composition"].var
    assert m.obs.loc["s4", "CST"] == "IV-B"
    assert "dominant_taxon" in m.obs.columns


def test_cli_check_reports_phyloseq_requirements(test_data, tmp_path):
    import yaml

    from microfgt.stages import check

    results = check(yaml.safe_load(
        yaml.safe_dump({"composition": {"phyloseq": str(test_data / FIXTURE)}})
    ))
    labels = " ".join(r.message for r in results)
    assert "phyloseq" in labels  # the phyloseq R package is checked for this entry point
    assert all(r.ok for r in results)
