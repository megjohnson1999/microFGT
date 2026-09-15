"""Tier 2: cutadapt drops primer-less reads by default (when primers are set) + extra_args passthrough.

Run 1 (Balle) found 16-27% of reads had no primer yet 100% of pairs were written into DADA2, with
no config knob to change it. These tests pin the new behaviour by capturing the argv cutadapt is
invoked with (the external tool itself is mocked).
"""

import microfgt.orchestrate.cutadapt as cut
import microfgt.stages.registry as reg
from microfgt.stages.model import StageContext
from pathlib import Path


def _make_pair(d):
    (d / "s_R1.fastq.gz").write_bytes(b"")
    (d / "s_R2.fastq.gz").write_bytes(b"")


def _capture_argv(monkeypatch):
    calls = []
    monkeypatch.setattr(cut, "resolve_executable", lambda exe, tool=None: (exe, {}))
    monkeypatch.setattr(cut, "run_command",
                        lambda argv, **kw: calls.append(argv) or type("R", (), {})())
    return calls


def test_discard_untrimmed_flag_added_when_requested(tmp_path, monkeypatch):
    calls = _capture_argv(monkeypatch)
    _make_pair(tmp_path)
    cut.run_cutadapt(tmp_path, tmp_path / "out", fwd_primer="ACGT", rev_primer="TGCA",
                     discard_untrimmed=True)
    assert "--discard-untrimmed" in calls[0]


def test_discard_untrimmed_absent_by_default_in_orchestrator(tmp_path, monkeypatch):
    calls = _capture_argv(monkeypatch)
    _make_pair(tmp_path)
    cut.run_cutadapt(tmp_path, tmp_path / "out", fwd_primer="ACGT")   # default False
    assert "--discard-untrimmed" not in calls[0]


def test_extra_args_appended_verbatim_and_last(tmp_path, monkeypatch):
    calls = _capture_argv(monkeypatch)
    _make_pair(tmp_path)
    cut.run_cutadapt(tmp_path, tmp_path / "out", fwd_primer="ACGT", rev_primer="TGCA",
                     discard_untrimmed=True, extra_args=["--pair-filter=both", "--minimum-length", "50"])
    argv = calls[0]
    # extra_args land after --discard-untrimmed so a user can override the defaults.
    assert argv.index("--pair-filter=both") > argv.index("--discard-untrimmed")
    assert ["--minimum-length", "50"] == argv[argv.index("--minimum-length"):argv.index("--minimum-length") + 2]


# --- the stage plumbing: config -> run_cutadapt kwargs -------------------------------------
def _run_primer_trim_capturing(monkeypatch, tmp_path, reads_cfg):
    captured = {}

    def fake_run_cutadapt(indir, outdir, **kw):
        captured.update(kw)
        return []

    monkeypatch.setattr(reg, "run_cutadapt", fake_run_cutadapt, raising=False)
    monkeypatch.setattr("microfgt.orchestrate.cutadapt.run_cutadapt", fake_run_cutadapt, raising=False)
    artifacts = {"fastq_dir": str(tmp_path / "raw"), "trimmed_reads": str(tmp_path / "trim")}
    ctx = StageContext(Path(tmp_path), artifacts, {"composition": {"reads": reads_cfg}})
    reg._run_primer_trim(ctx)
    return captured


def test_stage_defaults_discard_on_when_primers_configured(monkeypatch, tmp_path):
    cap = _run_primer_trim_capturing(monkeypatch, tmp_path,
                                     {"primers": {"fwd": "GTGYCAGCMGCCGCGGTAA", "rev": "GGACTAC"}})
    assert cap["discard_untrimmed"] is True


def test_stage_defaults_discard_off_when_no_primers(monkeypatch, tmp_path):
    cap = _run_primer_trim_capturing(monkeypatch, tmp_path, {})   # no primers -> would drop all
    assert cap["discard_untrimmed"] is False


def test_stage_discard_and_extra_args_overridable(monkeypatch, tmp_path):
    cap = _run_primer_trim_capturing(monkeypatch, tmp_path, {
        "primers": {"fwd": "GTGYCAGCMGCCGCGGTAA"},
        "discard_untrimmed": False,                       # explicit override wins
        "cutadapt_args": ["--pair-filter=both"],
    })
    assert cap["discard_untrimmed"] is False
    assert cap["extra_args"] == ["--pair-filter=both"]
