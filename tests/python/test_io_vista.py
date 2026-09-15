"""VISTA mgCST importer, grounded in the real ENA/PRJEB34536 fixture vista_mgCSTs.csv."""

import pandas as pd
import pytest

from microfgt.io import import_mgcst

MGCSTS = "vista_mgCSTs.csv"


def test_import_mgcst_reads_call_and_theta(real_fixtures):
    df = import_mgcst(real_fixtures / MGCSTS)

    assert list(df.columns) == ["mgCST", "mgCST_score"]
    assert df.index.name == "sample"
    # The label is kept verbatim as VISTA writes it ("mgCST <n>"), not silently reformatted.
    assert df.loc["ERR4421570", "mgCST"] == "mgCST 1"
    assert df.loc["ERR4421570", "mgCST_score"] == pytest.approx(0.986889578943112)
    # A poorly-fitting sample keeps its low best-match theta (surfaced, not hidden).
    assert df.loc["ERR4421590", "mgCST"] == "mgCST 25"
    assert df.loc["ERR4421590", "mgCST_score"] < 0.1


def test_import_mgcst_covers_every_sample(real_fixtures):
    raw = pd.read_csv(real_fixtures / MGCSTS, index_col=0)
    df = import_mgcst(real_fixtures / MGCSTS)
    assert list(df.index) == [str(i) for i in raw.index]


def test_import_mgcst_rejects_a_non_vista_csv(tmp_path):
    bad = tmp_path / "not_vista.csv"
    pd.DataFrame({"foo": [1]}, index=["s1"]).to_csv(bad)
    with pytest.raises(ValueError, match="VISTA mgCST column"):
        import_mgcst(bad)
