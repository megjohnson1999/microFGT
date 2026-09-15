"""`microfgt check` — preflight doctor over the resolved entry point."""

import os
import stat

from microfgt.stages import check


def _make_exe(path, body="#!/usr/bin/env bash\nexit 0\n"):
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_existing_outputs_entry_needs_no_external_tools():
    cfg = {"composition": {"speciateit": {"results": "r.txt", "count_table": "c.csv"}}}
    results = check(cfg)
    assert all(r.ok for r in results)           # those stages have no requirements


def test_missing_cutadapt_reported_for_fastq_entry(tmp_path):
    cfg = {"composition": {"reads": {"fastq_dir": "raw/", "cutadapt": "definitely_missing_xyz"},
                           "speciateit": {"db": str(tmp_path)}}}
    results = check(cfg)
    miss = [r for r in results if not r.ok]
    assert any("definitely_missing_xyz" in r.message for r in miss)


def test_region_db_mismatch_flagged(tmp_path):
    db = tmp_path / "vSpeciateIT_V4V4"           # V4 model...
    db.mkdir()
    classify = _make_exe(tmp_path / "classify")
    cfg = {"composition": {
        "reads": {"fastq_dir": "raw/", "region": "V3V4",   # ...but config says V3V4
                  "cutadapt": str(_make_exe(tmp_path / "cutadapt")),
                  "rscript": str(_make_exe(tmp_path / "Rscript"))},
        "speciateit": {"db": str(db), "classify": str(classify)},
    }}
    results = check(cfg)
    assert any("does not match" in r.message for r in results if not r.ok)


def test_rpackage_check_uses_the_configured_rscript(tmp_path):
    # The doctor must verify R packages with the CONFIGURED Rscript, not a bare 'Rscript' —
    # otherwise it misreports the VISTA layer for anyone whose R is not on PATH. Prove it by
    # flipping the result via two stub Rscripts (one reports every package present, one absent).
    present = _make_exe(tmp_path / "Rscript_present", "#!/usr/bin/env bash\nexit 0\n")
    absent = _make_exe(tmp_path / "Rscript_absent", "#!/usr/bin/env bash\nexit 1\n")

    def randomforestsrc_ok(rscript):
        cfg = {"metagenomics": {"reads": {"fastq_dir": "raw/"}, "rscript": str(rscript)}}
        line = next(r for r in check(cfg) if "randomForestSRC" in r.message)
        return line.ok

    assert randomforestsrc_ok(present) is True     # configured Rscript said "present"
    assert randomforestsrc_ok(absent) is False      # ...and "absent" — so it IS being consulted


def test_virgo2_checksum_guard_catches_version_divergence(tmp_path):
    import hashlib

    virgo2 = tmp_path / "VIRGO2"; virgo2.mkdir()
    script = virgo2 / "VIRGO2.py"; script.write_text("print('the zenodo one')\n")
    good = hashlib.sha256(script.read_bytes()).hexdigest()

    # Match the requirement KIND (2nd token: "OK|MISS <kind> '<name>'"), not a loose substring —
    # pytest's tmp_path is named after the test, so the VIRGO2.py *path* contains "checksum".
    def _is_checksum(r):
        return r.message.split()[1] == "checksum"

    def checksum_ok(sha):
        cfg = {"metagenomics": {"reads": {"fastq_dir": "raw/"},
                                "virgo2_dir": str(virgo2), "virgo2_sha256": sha}}
        return next(r for r in check(cfg) if _is_checksum(r)).ok

    assert checksum_ok(good) is True
    assert checksum_ok("deadbeef") is False          # divergence -> a visible MISS line


def test_no_checksum_check_without_a_pin(tmp_path):
    cfg = {"metagenomics": {"reads": {"fastq_dir": "raw/"}, "virgo2_dir": str(tmp_path)}}
    assert not any(r.message.split()[1] == "checksum" for r in check(cfg))   # opt-in only


def test_region_defaults_use_canonical_model_names_with_v4_alias():
    from microfgt.stages.registry import region_defaults

    # Canonical vSpeciateDB region names all resolve (incl. V4V4 and V1V9, previously missing).
    assert region_defaults("V4V4")["trunc_len"] == [230, 180]
    assert region_defaults("V1V3") == {"trunc_len": [0, 0], "trim_left": [0, 0]}
    assert region_defaults("V1V9") == {"trunc_len": [0, 0], "trim_left": [0, 0]}
    # "V4" (the colloquial name) aliases to V4V4 so users get the defaults either way.
    assert region_defaults("V4") == region_defaults("V4V4")
    # Unknown / missing region -> no defaults, but never errors.
    assert region_defaults("bogus") == {}
    assert region_defaults(None) == {}
