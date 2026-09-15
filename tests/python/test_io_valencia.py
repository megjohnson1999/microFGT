"""VALENCIA importer — validated against GENUINE valencia_genuine_output_head.csv."""

import math

import pandas as pd

from microfgt.io import import_valencia


def test_valencia_pulls_labels_from_trailing_columns(real_fixtures):
    out = import_valencia(real_fixtures / "valencia_genuine_output_head.csv")

    assert list(out.columns) == ["CST", "subCST", "score"]
    assert out.index.name == "sample"

    # Ground truth read directly from the genuine output (row for sampleID "1").
    assert out.loc["1", "CST"] == "I"
    assert out.loc["1", "subCST"] == "I-B"
    assert math.isclose(out.loc["1", "score"], 0.9886471215573692, rel_tol=1e-12)


def test_valencia_keys_by_sampleid_not_sample(real_fixtures):
    # FORMATS.md divergence: real uses 'sampleID' (the mock used 'Sample').
    raw = pd.read_csv(real_fixtures / "valencia_genuine_output_head.csv")
    out = import_valencia(real_fixtures / "valencia_genuine_output_head.csv")
    assert list(out.index) == [str(s) for s in raw["sampleID"]]


def test_valencia_missing_label_column_is_nan_not_error(tmp_path):
    # A partial VALENCIA output (only subCST) should still import.
    p = tmp_path / "partial.csv"
    p.write_text("sampleID,read_count,subCST\nS1,100,III-A\n")
    out = import_valencia(p)
    assert out.loc["S1", "subCST"] == "III-A"
    assert pd.isna(out.loc["S1", "CST"]) and pd.isna(out.loc["S1", "score"])
