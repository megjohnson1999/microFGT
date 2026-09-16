"""Sample-sheet reader + validator."""

import pytest

from microfgt.io import read_samplesheet, validate_samplesheet


def _touch(d, *names):
    """Create empty files under d and return their names (as the sheet would reference them)."""
    for n in names:
        (d / n).write_text("")
    return names


def _write_sheet(path, rows, header="sample_id,16s_R1,16s_R2,shotgun_R1,shotgun_R2"):
    lines = [header] + [",".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_read_normalises_headers_and_blanks(tmp_path):
    sheet = tmp_path / "s.csv"
    sheet.write_text("Sample_ID, 16S_R1 ,group\nPT01,a.fq,BV\nPT02,,\n")
    df = read_samplesheet(sheet)
    assert list(df.columns) == ["sample_id", "16s_r1", "group"]   # lowercased + stripped
    assert df.loc[1, "16s_r1"] is not df.loc[0, "16s_r1"]         # sanity
    assert df["group"].isna().sum() == 1                          # blank -> NA


def test_valid_both_arms_ok_and_counts(tmp_path):
    _touch(tmp_path, "a1.fq", "a2.fq", "m1.fq", "m2.fq", "b1.fq", "b2.fq")
    sheet = _write_sheet(tmp_path / "s.csv", [
        ["PT01", "a1.fq", "a2.fq", "m1.fq", "m2.fq"],   # both arms
        ["PT02", "b1.fq", "b2.fq", "", ""],             # 16S only
    ])
    rpt = validate_samplesheet(read_samplesheet(sheet), base_dir=tmp_path)
    assert rpt.ok
    assert (rpt.n_samples, rpt.n_16s, rpt.n_shotgun, rpt.n_both) == (2, 2, 1, 1)


def test_duplicate_ids_are_a_problem(tmp_path):
    _touch(tmp_path, "a.fq")
    sheet = _write_sheet(tmp_path / "s.csv",
                         [["PT01", "a.fq", "", "", ""], ["PT01", "a.fq", "", "", ""]])
    rpt = validate_samplesheet(read_samplesheet(sheet), base_dir=tmp_path)
    assert not rpt.ok
    assert any("duplicate" in p for p in rpt.problems)


def test_missing_id_column_is_a_problem(tmp_path):
    sheet = tmp_path / "s.csv"
    sheet.write_text("16s_R1,group\na.fq,BV\n")
    rpt = validate_samplesheet(read_samplesheet(sheet), base_dir=tmp_path)
    assert not rpt.ok and any("sample_id" in p for p in rpt.problems)


def test_missing_file_is_a_problem(tmp_path):
    sheet = _write_sheet(tmp_path / "s.csv", [["PT01", "nope.fq", "", "", ""]])
    rpt = validate_samplesheet(read_samplesheet(sheet), base_dir=tmp_path)
    assert not rpt.ok and any("missing file" in p for p in rpt.problems)


def test_r2_without_r1_is_a_problem(tmp_path):
    _touch(tmp_path, "r2.fq")
    sheet = _write_sheet(tmp_path / "s.csv", [["PT01", "", "r2.fq", "", ""]])
    rpt = validate_samplesheet(read_samplesheet(sheet), base_dir=tmp_path)
    assert not rpt.ok and any("16s_r1" in p for p in rpt.problems)


def test_both_arms_but_none_shared_warns_but_is_usable(tmp_path):
    _touch(tmp_path, "a.fq", "m.fq")
    sheet = _write_sheet(tmp_path / "s.csv",
                         [["PT01", "a.fq", "", "", ""], ["PT02", "", "", "m.fq", ""]])
    rpt = validate_samplesheet(read_samplesheet(sheet), base_dir=tmp_path)
    assert rpt.ok                                        # no blocking problems
    assert rpt.n_both == 0
    assert any("nothing will integrate" in w for w in rpt.warnings)


def test_cli_check_samplesheet(tmp_path, capsys):
    from microfgt.cli import main

    _touch(tmp_path, "a.fq")
    good = _write_sheet(tmp_path / "good.csv", [["PT01", "a.fq", "", "", ""]])
    assert main(["check-samplesheet", "-s", str(good)]) == 0
    assert "sample sheet OK" in capsys.readouterr().out

    bad = _write_sheet(tmp_path / "bad.csv", [["PT01", "missing.fq", "", "", ""]])
    with pytest.raises(SystemExit):
        main(["check-samplesheet", "-s", str(bad)])
