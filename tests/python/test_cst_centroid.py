"""Centroid CST — validated against GENUINE Valencia.py output + the swappable seam.

Two validation strengths:
* Exact reproduction of genuine Valencia.py output on the 6-sample head fixture
  (always-on; uses the bundled centroids). This is real-output validation.
* The full 99.9%-vs-VALENCIA gate on the 13k published samples lives in
  test_cst_validation_gate.py (skips if the large dataset isn't staged locally).
"""

import numpy as np
import pandas as pd
import pytest

from microfgt.cst import CST_ORDER, available_methods, classify_cst, register_method
from microfgt.cst.centroid import classify_centroid

_APPENDED = {"subCST", "score", "CST"}


def _reconstruct_input(genuine: pd.DataFrame):
    """Recover (taxa counts, read_count) from a genuine VALENCIA output CSV."""
    taxa_cols = [
        c for c in genuine.columns
        if c not in ("sampleID", "read_count") and not c.endswith("_sim") and c not in _APPENDED
    ]
    counts = genuine[taxa_cols].copy()
    counts.index = genuine["sampleID"].astype(str)
    rc = pd.Series(genuine["read_count"].to_numpy(), index=counts.index)
    return counts, rc


def test_centroid_reproduces_genuine_valencia_output(real_fixtures):
    genuine = pd.read_csv(real_fixtures / "valencia_genuine_output_head.csv")
    counts, rc = _reconstruct_input(genuine)

    # Genuine Valencia.py head fixture is 2020-named, so pin the 2020 centroids (default is 2024).
    out = classify_centroid(counts, read_count=rc, reference="2020")

    # Labels reproduced exactly.
    gi = genuine.set_index(genuine["sampleID"].astype(str))
    assert list(out["subCST"]) == list(gi["subCST"])
    assert list(out["CST"]) == list(gi["CST"])
    np.testing.assert_allclose(out["score"].to_numpy(), gi["score"].to_numpy(), atol=1e-6)

    # Every one of the 13 similarity columns reproduced to numerical precision.
    for cst in CST_ORDER:
        col = f"{cst}_sim"
        np.testing.assert_allclose(
            out[col].to_numpy(), gi[col].to_numpy(), atol=1e-6,
            err_msg=f"{col} mismatch vs genuine Valencia.py output",
        )


def test_centroid_warns_on_total_taxon_name_mismatch():
    # Taxa that share NO names with the reference centroids (e.g. GTDB-with-accession names) must
    # not silently produce confident-looking calls: theta is 0 everywhere and it collapses to the
    # first subCST. That has to surface as a warning, not pass quietly. (Run 1 finding.)
    counts = pd.DataFrame(
        {"Gardnerella_vaginalis(RS_GCF_x)": [100, 50], "Foo_bar": [0, 50]},
        index=["s1", "s2"],
    )
    with pytest.warns(UserWarning, match="reference centroid names"):
        out = classify_centroid(counts)
    assert (out["score"] == 0).all()
    assert out["subCST"].nunique() == 1          # all collapsed to the same (first) subCST


def test_centroid_no_mismatch_warning_on_normal_input(real_fixtures):
    genuine = pd.read_csv(real_fixtures / "valencia_genuine_output_head.csv")
    counts, rc = _reconstruct_input(genuine)
    import warnings as _w

    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        classify_centroid(counts, read_count=rc, reference="2020")
    assert not [w for w in caught if "centroid names" in str(w.message) or "scored 0" in str(w.message)]


def test_classify_cst_dispatches_to_centroid(real_fixtures):
    genuine = pd.read_csv(real_fixtures / "valencia_genuine_output_head.csv")
    counts, rc = _reconstruct_input(genuine)
    via_iface = classify_cst(counts, method="centroid", read_count=rc)
    direct = classify_centroid(counts, read_count=rc)
    pd.testing.assert_frame_equal(via_iface, direct)


def test_seam_is_swappable():
    # The interface is the point: a new method registers and is callable by name.
    assert "centroid" in available_methods()

    def _dummy(composition, **kw):
        return pd.DataFrame({"CST": ["IV-C"]}, index=["s1"])

    register_method("dummy", _dummy)
    assert "dummy" in available_methods()
    assert classify_cst(pd.DataFrame(index=["s1"]), method="dummy").loc["s1", "CST"] == "IV-C"


def test_unknown_method_raises():
    with pytest.raises(ValueError, match="Unknown CST method"):
        classify_cst(pd.DataFrame(index=["s1"]), method="nope")


def test_default_centroids_are_2024_and_match_speciateit_naming():
    """Default centroids are the 2024 set, which carries the modern taxon names speciateIT v6
    emits (the M4 fix). The 2020 set silently lacks these core BV taxa."""
    from microfgt.cst.centroid import load_reference_centroids

    default = load_reference_centroids()
    v2024 = load_reference_centroids("2024")
    v2020 = load_reference_centroids("2020")

    assert list(default.index) == CST_ORDER          # 13 subCSTs, right order
    assert default.equals(v2024)                      # default IS the 2024 set
    assert v2024.shape[1] > v2020.shape[1]            # 2024 has the larger, modern vocabulary
    # taxa speciateIT v6 emits that the 2020 centroids miss but 2024 includes:
    for taxon in ("Fannyhessea_vaginae", "Ca_Lachnocurva_vaginae", "Lactobacillus_mulieris"):
        assert taxon in v2024.columns, f"{taxon} missing from 2024 centroids"
        assert taxon not in v2020.columns, f"{taxon} unexpectedly in 2020 centroids"
