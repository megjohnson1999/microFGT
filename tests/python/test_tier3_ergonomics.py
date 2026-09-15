"""Tier 3 ergonomics: three newcomer traps from Run 1 turned into clear, actionable messages.

1. space-in-path warning for the speciateIT stage (pipeline, not just `setup`)
2. download retry + short-read detection in `microfgt setup`
3. a conda-env hint when the *missing* prereqs are tools conda provides
"""

import io
import stat

import pytest

from microfgt import setup16s
from microfgt.stages import speciateit_space_warnings


def _make_exe(path, body="#!/usr/bin/env bash\nexit 0\n"):
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


# --- 1. space-in-path warning --------------------------------------------------------------
def test_space_warning_fires_for_spaced_workdir_when_speciateit_runs():
    cfg = {"composition": {
        "reads": {"fastq_dir": "/space free/raw", "region": "V4V4"},
        "speciateit": {"db": "/db/vSpeciateIT_V4V4", "classify": "classify"},
    }}
    warns = speciateit_space_warnings(cfg, workdir="/Users/Megan Johnson/work")
    assert any("workdir" in w and "space" in w for w in warns)


def test_space_warning_flags_spaced_input_fastq_dir():
    cfg = {"composition": {
        "reads": {"fastq_dir": "/Users/Megan Johnson/raw", "region": "V4V4"},
        "speciateit": {"db": "/db/vSpeciateIT_V4V4"},
    }}
    warns = speciateit_space_warnings(cfg, workdir="/tmp/clean")
    assert any("composition.reads.fastq_dir" in w for w in warns)


def test_no_space_warning_when_paths_are_clean():
    cfg = {"composition": {
        "reads": {"fastq_dir": "/clean/raw", "region": "V4V4"},
        "speciateit": {"db": "/db/vSpeciateIT_V4V4"},
    }}
    assert speciateit_space_warnings(cfg, workdir="/clean/work") == []


def test_no_space_warning_when_speciateit_not_in_plan():
    # An existing-outputs entry (no assign stage) never invokes classify, so spaces are harmless.
    cfg = {"composition": {"speciateit": {"results": "/a b/r.txt", "count_table": "/a b/c.csv"}}}
    assert speciateit_space_warnings(cfg, workdir="/a b/work") == []


# --- 2. download retry ---------------------------------------------------------------------
class _Resp(io.BytesIO):
    def __init__(self, data, length):
        super().__init__(data)
        self.headers = {"Content-Length": str(length)}

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_download_retries_then_succeeds(tmp_path, monkeypatch, capsys):
    import urllib.error

    attempts = {"n": 0}

    def flaky_urlopen(url, timeout=0):
        attempts["n"] += 1
        if attempts["n"] == 1:
            raise urllib.error.URLError("transient reset")
        return _Resp(b"hello", 5)

    monkeypatch.setattr(setup16s.urllib.request, "urlopen", flaky_urlopen)
    out = setup16s.download("https://x/y.zip", tmp_path / "y.zip", retries=3)
    assert out.read_bytes() == b"hello"
    assert attempts["n"] == 2                       # failed once, succeeded on the retry
    assert "retrying" in capsys.readouterr().out


def test_download_detects_short_read(tmp_path, monkeypatch):
    # Server claims 100 bytes but delivers 5 -> caught and (after retries) raised, not passed on.
    monkeypatch.setattr(setup16s.urllib.request, "urlopen",
                        lambda url, timeout=0: _Resp(b"short", 100))
    with pytest.raises(setup16s.SetupError):
        setup16s.download("https://x/y.zip", tmp_path / "y.zip", retries=2)


# --- 3. conda-env hint on a MISS -----------------------------------------------------------
def test_conda_env_hint_when_conda_tools_missing(tmp_path, capsys):
    from microfgt.cli import main
    import yaml

    cfg = tmp_path / "c.yaml"
    # Default tool names (cutadapt / Rscript+dada2): in a bare env these MISS on PATH — exactly
    # the newcomer-forgot-to-activate case. The hint keys on those canonical names.
    cfg.write_text(yaml.safe_dump({"composition": {
        "reads": {"fastq_dir": "raw/", "region": "V4V4"},
        "speciateit": {"db": str(tmp_path), "classify": str(_make_exe(tmp_path / "classify"))},
    }}))
    with pytest.raises(SystemExit):           # check signals failure via SystemExit
        main(["check", "-c", str(cfg)])
    out = capsys.readouterr().out
    assert "conda activate" in out            # points the newcomer at the real cause
