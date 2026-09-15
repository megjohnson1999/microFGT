"""Taxonomy → classification mapping in the phyloseq importer (pure Python, no R needed).

Covers the Genus_Species → Species → Genus fallback and GTDB rank-prefix stripping (the HVTN
convention), which the R-gated importer tests don't reach.
"""

import pandas as pd

from microfgt.io.phyloseq import _classification_from_tax


def test_prefers_genus_species_then_species_then_genus():
    tax = pd.DataFrame(
        {
            "Genus_Species": ["Lactobacillus_iners", pd.NA, pd.NA],
            "Species": ["ignored", "Gardnerella_vaginalis", pd.NA],
            "Genus": ["ignored", "ignored", "Prevotella"],
        },
        index=["A", "B", "C"],
    )
    cls = _classification_from_tax(tax)
    assert list(cls) == ["Lactobacillus_iners", "Gardnerella_vaginalis", "Prevotella"]


def test_strips_gtdb_rank_prefixes():
    # HVTN-style: no Genus_Species rank; Species clean binomials; Genus GTDB-prefixed.
    tax = pd.DataFrame(
        {
            "Species": ["Lactobacillus_iners", pd.NA],
            "Genus": ["g_Lactobacillus", "g_Gardnerella"],
            "Domain": ["d_Bacteria", "d__Bacteria"],
        },
        index=["ASV1", "ASV2"],
    )
    cls = _classification_from_tax(tax)
    assert cls.loc["ASV1"] == "Lactobacillus_iners"     # species, already clean
    assert cls.loc["ASV2"] == "Gardnerella"             # genus fallback, g_ prefix stripped


def test_preserves_normal_and_candidatus_names():
    tax = pd.DataFrame(
        {"Genus_Species": ["Ca_Lachnocurva_vaginae", "Lactobacillus_crispatus", "crispatus"]},
        index=["A", "B", "C"],
    )
    cls = _classification_from_tax(tax)
    assert cls.loc["A"] == "Ca_Lachnocurva_vaginae"     # Candidatus prefix is NOT a rank prefix
    assert cls.loc["B"] == "Lactobacillus_crispatus"
    assert cls.loc["C"] == "crispatus"                  # bare epithet untouched


def test_missing_everywhere_is_na():
    tax = pd.DataFrame({"Genus": ["", "NA"], "Species": [pd.NA, pd.NA]}, index=["A", "B"])
    cls = _classification_from_tax(tax)
    assert cls.isna().all()


def test_gtdb_accession_and_doubled_genus_are_cleaned():
    # FRESH-style GTDB tax_table: Species carries an accession suffix, Genus_Species repeats the
    # genus with a space, and unclassified rows are placeholder text. (Regression for the Run 1
    # validation blowup where 'Sneathia_vaginalis(RS_GCF...' leaked into the classification.)
    tax = pd.DataFrame(
        {
            "Genus_Species": [
                "Lactobacillus Lactobacillus_iners(RS_GCF_000160875_1",
                "Fannyhessea Fannyhessea_vaginae(RS_GCF_000159235_2",
                "Bacteria Domain Bacteria Domain",
            ],
            "Species": [
                "Lactobacillus_iners(RS_GCF_000160875_1",
                "Fannyhessea_vaginae(RS_GCF_000159235_2",
                "Bacteria Domain",
            ],
            "Genus": ["Lactobacillus", "Fannyhessea", "Bacteria Domain"],
        },
        index=["A", "B", "C"],
    )
    cls = _classification_from_tax(tax)
    assert cls.loc["A"] == "Lactobacillus_iners"       # accession + doubled genus stripped
    assert cls.loc["B"] == "Fannyhessea_vaginae"
    assert pd.isna(cls.loc["C"])                        # placeholder -> unclassified, not a fake taxon
