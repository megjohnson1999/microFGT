"""`microfgt setup` — install the 16S prerequisites (speciateIT binary + a vSpeciateDB model).

Network steps (git clone, figshare API, download) are monkeypatched; everything else — region
canonicalisation, OS binary selection, the figshare filename match (incl. the dated V1V9 name),
checksum policy, the written config being `check`-clean, and the ground-truth comparison — runs
for real against local fixtures.
"""

import io
import json
import stat

import pytest

from microfgt import setup16s
from microfgt.stages.registry import canonical_region


def _make_exe(path, body="#!/usr/bin/env bash\nexit 0\n"):
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


# --- region canonicalisation ---------------------------------------------------------------
@pytest.mark.parametrize("raw,expected", [
    ("V3V4", "V3V4"), ("v3v4", "V3V4"), ("V4", "V4V4"), ("v4", "V4V4"),
    ("V1V3", "V1V3"), ("V1V9", "V1V9"), ("V4V4", "V4V4"),
])
def test_canonical_region_accepts_names_and_v4_alias(raw, expected):
    assert canonical_region(raw) == expected


@pytest.mark.parametrize("bogus", ["V2", "16S", "", "V3-V4x"])
def test_canonical_region_rejects_unknown(bogus):
    with pytest.raises(ValueError):
        canonical_region(bogus)


# --- OS binary selection -------------------------------------------------------------------
def test_select_binary_per_os(tmp_path):
    clone = tmp_path / "speciateIT"
    (clone / "bin" / "macosx").mkdir(parents=True)
    (clone / "bin" / "linux").mkdir(parents=True)
    (clone / "bin" / "macosx" / "classify").write_text("mac")
    (clone / "bin" / "linux" / "classify").write_text("linux")
    assert setup16s.select_binary(clone, system="Darwin").name == "classify"
    assert "macosx" in str(setup16s.select_binary(clone, system="Darwin"))
    assert "linux" in str(setup16s.select_binary(clone, system="Linux"))


def test_select_binary_unsupported_os(tmp_path):
    with pytest.raises(setup16s.SetupError):
        setup16s.select_binary(tmp_path, system="Windows")


def test_select_binary_missing_file(tmp_path):
    with pytest.raises(setup16s.SetupError):
        setup16s.select_binary(tmp_path, system="Darwin")


# --- figshare file resolution (the real gotcha: the dated V1V9 name) -----------------------
_FIGSHARE_META = {
    "version": 6,
    "files": [
        {"name": "vSpeciateIT_V3V4.zip", "size": 123, "download_url": "u/v3v4", "computed_md5": "a"},
        {"name": "vSpeciateIT_V4V4.zip", "size": 88, "download_url": "u/v4v4", "computed_md5": "b"},
        {"name": "vSpeciateIT_V1V3.zip", "size": 141, "download_url": "u/v1v3", "computed_md5": "c"},
        {"name": "vSpeciateIT_V1V9_14May2026.zip", "size": 320, "download_url": "u/v1v9",
         "computed_md5": "d"},
    ],
}


@pytest.fixture
def _fake_figshare(monkeypatch):
    def fake_urlopen(url, timeout=0):
        return io.BytesIO(json.dumps(_FIGSHARE_META).encode())
    monkeypatch.setattr(setup16s.urllib.request, "urlopen", fake_urlopen)


@pytest.mark.parametrize("region,name", [
    ("V3V4", "vSpeciateIT_V3V4.zip"),
    ("V4V4", "vSpeciateIT_V4V4.zip"),
    ("V1V3", "vSpeciateIT_V1V3.zip"),
    ("V1V9", "vSpeciateIT_V1V9_14May2026.zip"),   # dated name must still match
])
def test_resolve_figshare_file_matches(_fake_figshare, region, name):
    assert setup16s.resolve_figshare_file(region)["name"] == name


def test_resolve_figshare_v1v3_does_not_swallow_v1v9(_fake_figshare):
    # A naive prefix match on "V1V" would be ambiguous; the token boundary keeps them distinct.
    assert setup16s.resolve_figshare_file("V1V3")["name"] == "vSpeciateIT_V1V3.zip"


def test_resolve_figshare_missing_region(monkeypatch):
    meta = {"files": [{"name": "vSpeciateIT_V3V4.zip", "download_url": "u"}]}
    monkeypatch.setattr(setup16s.urllib.request, "urlopen",
                        lambda url, timeout=0: io.BytesIO(json.dumps(meta).encode()))
    with pytest.raises(setup16s.SetupError):
        setup16s.resolve_figshare_file("V4V4")


# --- checksum policy -----------------------------------------------------------------------
def test_verify_zip_pinned_ok(tmp_path):
    z = tmp_path / "z.zip"
    z.write_bytes(b"payload")
    import hashlib
    sha = hashlib.sha256(b"payload").hexdigest()
    setup16s.PINNED_SHA256["_TEST"] = sha
    try:
        assert setup16s.verify_zip("_TEST", z, require_pinned=True) == sha
    finally:
        del setup16s.PINNED_SHA256["_TEST"]


def test_verify_zip_pinned_mismatch_raises(tmp_path):
    z = tmp_path / "z.zip"
    z.write_bytes(b"payload")
    setup16s.PINNED_SHA256["_TEST"] = "deadbeef"
    try:
        with pytest.raises(setup16s.SetupError):
            setup16s.verify_zip("_TEST", z, require_pinned=False)
    finally:
        del setup16s.PINNED_SHA256["_TEST"]


def test_verify_zip_unpinned_tofu_records_and_warns(tmp_path, capsys):
    z = tmp_path / "z.zip"
    z.write_bytes(b"payload")
    sha = setup16s.verify_zip("V1V9", z, require_pinned=False)   # V1V9 is unpinned
    assert len(sha) == 64
    assert "UNPINNED" in capsys.readouterr().out


def test_verify_zip_unpinned_require_pinned_refuses(tmp_path):
    z = tmp_path / "z.zip"
    z.write_bytes(b"payload")
    with pytest.raises(setup16s.SetupError):
        setup16s.verify_zip("V1V9", z, require_pinned=True)


# --- ground-truth comparison ---------------------------------------------------------------
def test_ground_truth_passes_with_documented_asv9_drift():
    calls = dict(setup16s.GROUND_TRUTH_V3V4)
    calls["ASV9"] = "Sneathia_sanguinegens"           # the documented v6 flip
    ok, report = setup16s.ground_truth_report(calls)
    assert ok
    assert any("9/10" in line for line in report)


def test_ground_truth_fails_on_real_mismatch():
    calls = dict(setup16s.GROUND_TRUTH_V3V4)
    calls["ASV1"] = "Wrong_species"                   # a non-allowed mismatch
    ok, _ = setup16s.ground_truth_report(calls)
    assert not ok


# --- the written config is check-clean + the full run_setup happy path ---------------------
def _fake_clone(tmp_path):
    """A fake speciateIT clone: OS binary dirs, a stub classify, and the 10-ASV test.fasta."""
    clone = tmp_path / "dest" / "speciateIT"
    (clone / "bin" / "macosx").mkdir(parents=True)
    (clone / "bin" / "linux").mkdir(parents=True)
    _make_exe(clone / "bin" / "macosx" / "classify")
    _make_exe(clone / "bin" / "linux" / "classify")
    fasta = "\n".join(f">{a}\nACGT" for a in setup16s.GROUND_TRUTH_V3V4)
    (clone / "test.fasta").write_text(fasta + "\n")
    return clone


def test_written_config_is_check_clean(tmp_path):
    from microfgt.config import load_config
    from microfgt.stages import check

    clone = _fake_clone(tmp_path)
    dest = tmp_path / "dest"
    classify = clone / "bin" / "macosx" / "classify"
    db = dest / "vSpeciateDB_models" / "vSpeciateIT_V3V4"
    db.mkdir(parents=True)

    cfg_path = setup16s.write_config(dest, "V3V4", classify, db, clone / "test.fasta")
    results = check(load_config(str(cfg_path)))
    # The three prereqs the brief names: binary resolves, db path exists, region<->DB match.
    assert all(r.ok for r in results), [r.message for r in results]
    kinds = " ".join(r.message for r in results)
    assert "classify" in kinds and "vSpeciateIT_V3V4" in kinds and "matches" in kinds


def test_run_setup_happy_path(tmp_path, monkeypatch, capsys):
    clone = _fake_clone(tmp_path)
    dest = tmp_path / "dest"
    # Pre-create the model dir so download is skipped (skip_download guards the fetch path).
    (dest / "vSpeciateDB_models" / "vSpeciateIT_V3V4").mkdir(parents=True)

    monkeypatch.setattr(setup16s, "clone_speciateit", lambda d: clone)
    monkeypatch.setattr(setup16s.platform, "system", lambda: "Darwin")

    # Stand in for the real classify run with the documented v6 output (ASV9 flipped).
    def fake_run_classify(classify, db, fasta, outdir):
        from pathlib import Path
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        rows = []
        for asv, sp in setup16s.GROUND_TRUTH_V3V4.items():
            sp = "Sneathia_sanguinegens" if asv == "ASV9" else sp
            rows.append(f"{asv}\t{sp}\t0.97\t46")
        res = outdir / "MC_order7_results.txt"
        res.write_text("\n".join(rows) + "\n")
        return res
    monkeypatch.setattr(setup16s, "run_classify", fake_run_classify)

    rc = setup16s.run_setup("V3V4", str(dest), skip_download=True)
    out = capsys.readouterr().out
    assert rc == 0, out
    assert (dest / "microfgt-16s.yaml").exists()
    assert "all prerequisites satisfied" in out
    assert "9/10" in out                              # the ground-truth line ran


def test_run_setup_warns_on_spaced_dest(tmp_path, monkeypatch, capsys):
    clone = _fake_clone(tmp_path)
    spaced = tmp_path / "a dir with spaces"
    (spaced / "vSpeciateDB_models" / "vSpeciateIT_V3V4").mkdir(parents=True)
    monkeypatch.setattr(setup16s, "clone_speciateit", lambda d: clone)
    monkeypatch.setattr(setup16s.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(setup16s, "run_classify",
                        lambda *a, **k: pytest.skip("not reached"))  # region V3V4 would call it
    # Use V1V3 to skip the ground-truth run entirely and isolate the spaced-dest warning.
    (spaced / "vSpeciateDB_models" / "vSpeciateIT_V1V3").mkdir(parents=True)
    setup16s.run_setup("V1V3", str(spaced), skip_download=True)
    assert "contains a space" in capsys.readouterr().out


def test_run_setup_bad_region_returns_2(tmp_path, capsys):
    assert setup16s.run_setup("V2", str(tmp_path)) == 2
    assert "unknown region" in capsys.readouterr().out
