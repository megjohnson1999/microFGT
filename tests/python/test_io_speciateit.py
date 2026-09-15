"""speciateIT importer (ASV grain) + the taxon roll-up.

NOT real-output validated: no genuine MC_order7_results.txt ships with the fixtures
(speciateIT inputs only). These tests exercise the join against a SYNTHETIC classification
file whose VALUES are the genuine ASV1-ASV10 assignments from the speciateIT README, joined
to the genuine ASV count table (169 samples x 1514 ASVs) and the trimmed genuine FASTA
(ASV1-ASV10). Real-output validation is an IOU for the first real HTCF run.
"""

import numpy as np
import pandas as pd
import pytest

from microfgt.io import collapse_to_taxon, import_speciateit
from microfgt.io.speciateit import _genus_of, make_genus_resolver

RESULTS = "speciateit_MC_order7_results.synthetic.txt"


# --- Tier 1 Fix 1: genus resolution (rank rule + binomial fallback + tree for polyphyly) ----
@pytest.mark.parametrize("label,expected", [
    ("Lactobacillus_iners", "Lactobacillus"),      # ordinary binomial
    ("Prevotella_bivia", "Prevotella"),
    ("Gardnerella_vaginalis", "Gardnerella"),      # the CST IV-B driver — must NOT become NA
    ("Ca_Lachnocurva_vaginae", "Ca_Lachnocurva"),  # Candidatus prefix kept (old bug -> "Ca")
    ("Dorea_A_longicatena", "Dorea_A"),            # genus-level polyphyly suffix kept
    ("g_Prevotella", "Prevotella"),                # genus-rank backoff -> the genus
    ("d_Bacteria", None),                          # domain backoff -> NA (old bug -> "d")
    ("o_Acetivibrionales", None),                  # order backoff -> NA
    ("f_Lachnospiraceae", None),                   # family backoff -> NA
    ("Bacteria", "Bacteria"),                      # single token -> itself
])
def test_genus_treeless_rule(label, expected):
    assert _genus_of(label) == expected


def test_no_genus_is_a_single_letter():
    # The whole point: no rank code leaks through as a fake genus.
    for label in ["d_Bacteria", "g_Prevotella", "o_X", "p_Y", "c_Z", "f_W", "s_V"]:
        g = _genus_of(label)
        assert g is None or len(g) > 1


def _write_tree(tmp_path):
    # Minimal Newick with the two polyphyly shapes that a string rule cannot tell apart:
    #  - Dorea_A_longicatena sits under g_Dorea_A  (the _A belongs to the GENUS)
    #  - Aerococcus_urinae_A sits under g_Aerococcus (the _A belongs to the SPECIES)
    # Plus Gardnerella_vaginalis as a bare leaf (NO g_ ancestor) to prove the binomial fallback.
    db = tmp_path / "vSpeciateIT_TEST"
    db.mkdir()
    (db / "model.tree").write_text(
        "((Dorea_A_longicatena:0.1,Dorea_A_formicigenerans:0.1)g_Dorea_A:0.1,"
        "(Aerococcus_urinae_A:0.1,Aerococcus_christensenii:0.1)g_Aerococcus:0.1,"
        "Gardnerella_vaginalis:0.1)root;\n"
    )
    return db


def test_genus_tree_resolves_polyphyly(tmp_path):
    resolve = make_genus_resolver(_write_tree(tmp_path))
    # The tree disambiguates the two _A shapes — a string rule can't.
    assert resolve("Dorea_A_longicatena") == "Dorea_A"      # _A is the genus
    assert resolve("Aerococcus_urinae_A") == "Aerococcus"   # _A is the species (tree overrides rsplit)
    # A leaf with no g_ ancestor falls back to the binomial, NOT NA.
    assert resolve("Gardnerella_vaginalis") == "Gardnerella"
    # Backoff labels still follow the rank rule even with a tree present.
    assert resolve("g_Prevotella") == "Prevotella"
    assert resolve("d_Bacteria") is None


def test_genus_resolver_without_tree_is_treeless(tmp_path):
    # Missing/absent db -> tree-less rule, never a crash.
    resolve = make_genus_resolver(tmp_path / "does_not_exist")
    assert resolve("Aerococcus_urinae_A") == "Aerococcus_urinae"  # no tree -> binomial rsplit
    assert make_genus_resolver(None)("Lactobacillus_iners") == "Lactobacillus"


def test_tree_override_is_load_bearing(tmp_path):
    # The whole reason the tree exists: for a species-level polyphyly suffix, the string rule and
    # the tree give DIFFERENT genera, and only the tree is right. Assert BOTH answers on the SAME
    # input so a refactor that quietly drops the tree lookup cannot pass green — with-tree MUST
    # differ from without-tree here. (Guards tree-necessity regardless of what any real dataset,
    # e.g. Balle, happens to contain.)
    label = "Aerococcus_urinae_A"                      # _A is on the SPECIES, not the genus
    with_tree = make_genus_resolver(_write_tree(tmp_path))(label)
    without_tree = make_genus_resolver(None)(label)
    assert with_tree == "Aerococcus"                   # tree: nearest g_ ancestor is g_Aerococcus
    assert without_tree == "Aerococcus_urinae"         # string rsplit: strips only the last token
    assert with_tree != without_tree                   # the override is real and load-bearing


# The real vSpeciateDB model (installed by `microfgt setup`), if present on this machine. Pins
# that the synthetic tree above reflects the REAL DB's structure — not just an assumption I baked
# into the fixture. Skips cleanly where the ~2.4 GB model isn't installed (e.g. CI).
import os  # noqa: E402

_REAL_V4V4 = os.path.expanduser(
    "~/Projects/microfgt-refdata/speciateIT/vSpeciateDB_models/vSpeciateIT_V4V4"
)


@pytest.mark.skipif(
    not os.path.exists(os.path.join(_REAL_V4V4, "model.tree")),
    reason="vSpeciateDB V4V4 model not installed",
)
def test_real_vspeciatedb_tree_has_the_polyphyly_structure():
    resolve = make_genus_resolver(_REAL_V4V4)
    assert resolve("Aerococcus_urinae_A") == "Aerococcus"      # tree override present in the REAL DB
    assert resolve("Dorea_A_longicatena") == "Dorea_A"         # genus-level suffix kept
    assert resolve("Gardnerella_vaginalis") == "Gardnerella"   # bare leaf -> binomial, not NA


def _count_table(real_fixtures):
    return pd.read_csv(real_fixtures / "speciateit_test_count_table.csv", index_col=0)


# --- import: ASV grain, sequences retained --------------------------------------------------
def test_import_is_asv_grain_and_conserves_counts(real_fixtures, test_data):
    ct = _count_table(real_fixtures)
    adata = import_speciateit(test_data / RESULTS, real_fixtures / "speciateit_test_count_table.csv")

    # One feature per ASV (no read-time collapse); samples preserved; counts conserved.
    assert list(adata.var_names) == [str(c) for c in ct.columns]
    assert adata.n_vars == ct.shape[1]  # 1514 ASVs
    assert list(adata.obs_names) == [str(s) for s in ct.index]
    assert adata.layers["counts"].sum() == int(ct.to_numpy().sum())


def test_import_keeps_classification_and_genus_per_asv(real_fixtures, test_data):
    adata = import_speciateit(test_data / RESULTS, real_fixtures / "speciateit_test_count_table.csv")

    assert adata.var.loc["ASV2", "classification"] == "Lactobacillus_crispatus"
    assert adata.var.loc["ASV2", "genus"] == "Lactobacillus"
    # Unclassified ASVs keep their identity as features, with a missing classification.
    assert pd.isna(adata.var.loc["ASV11", "classification"])


def test_import_retains_sequences_when_fasta_given(real_fixtures, test_data):
    adata = import_speciateit(
        test_data / RESULTS,
        real_fixtures / "speciateit_test_count_table.csv",
        fasta=real_fixtures / "speciateit_test.fasta",
    )
    # ASV1-ASV10 are in the trimmed FASTA; the rest are not (missing sequence, honestly).
    assert adata.var.loc["ASV1", "sequence"].startswith("TAGGGAAT")
    assert pd.isna(adata.var.loc["ASV11", "sequence"])


def test_import_without_fasta_has_no_sequence_column(real_fixtures, test_data):
    adata = import_speciateit(test_data / RESULTS, real_fixtures / "speciateit_test_count_table.csv")
    assert "sequence" not in adata.var


# --- roll-up: ASV -> taxon ------------------------------------------------------------------
def test_collapse_aggregates_asvs_sharing_a_classification(real_fixtures, test_data):
    ct = _count_table(real_fixtures)
    comp = import_speciateit(test_data / RESULTS, real_fixtures / "speciateit_test_count_table.csv")
    taxon = collapse_to_taxon(comp)

    # ASV2 + ASV6 are both Lactobacillus_crispatus -> their columns aggregate.
    expected = (ct["ASV2"] + ct["ASV6"]).to_numpy()
    ti = list(taxon.var_names).index("Lactobacillus_crispatus")
    np.testing.assert_array_equal(taxon.layers["counts"][:, ti], expected)
    assert taxon.var.loc["Lactobacillus_crispatus", "genus"] == "Lactobacillus"
    # Roll-up conserves total counts.
    assert taxon.layers["counts"].sum() == comp.layers["counts"].sum()


def test_collapse_buckets_unclassified_by_default(real_fixtures, test_data):
    ct = _count_table(real_fixtures)
    comp = import_speciateit(test_data / RESULTS, real_fixtures / "speciateit_test_count_table.csv")
    taxon = collapse_to_taxon(comp)

    assert "Unclassified" in taxon.var_names
    classified_asvs = [f"ASV{i}" for i in range(1, 11)]
    expected_unclassified = ct.drop(columns=classified_asvs).to_numpy().sum()
    ti = list(taxon.var_names).index("Unclassified")
    assert taxon.layers["counts"][:, ti].sum() == expected_unclassified


def test_collapse_keeps_unclassified_per_asv_when_disabled(real_fixtures, test_data):
    # bucket_unclassified=False matches speciateIT's own count_table.py behaviour.
    comp = import_speciateit(test_data / RESULTS, real_fixtures / "speciateit_test_count_table.csv")
    taxon = collapse_to_taxon(comp, bucket_unclassified=False)

    assert "Unclassified" not in taxon.var_names
    assert "ASV11" in taxon.var_names  # unclassified ASV kept as its own "taxon"


def test_collapse_requires_asv_grain():
    import anndata as ad

    bare = ad.AnnData(
        X=np.ones((2, 2), dtype=np.float32),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["taxonA", "taxonB"]),  # no 'classification'
    )
    with pytest.raises(ValueError, match="classification"):
        collapse_to_taxon(bare)


def test_header_autodetect_matches_headerless(real_fixtures, test_data):
    # The headed and headerless variants must produce identical ASV-grain results.
    ct = real_fixtures / "speciateit_test_count_table.csv"
    a = import_speciateit(test_data / RESULTS, ct)
    b = import_speciateit(test_data / "speciateit_MC_order7_results.headed.txt", ct)
    assert list(a.var_names) == list(b.var_names)
    np.testing.assert_array_equal(a.layers["counts"], b.layers["counts"])
    pd.testing.assert_series_equal(a.var["classification"], b.var["classification"])
