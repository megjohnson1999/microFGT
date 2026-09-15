"""Tier 1 Fix 2: the DADA2 trim_left default must honor the primer-trim contract.

cutadapt (the primer_trim stage) already strips primers WHEN they're configured, so applying
the region-default trim_left on top double-trims real sequence. The default must therefore be
0 when primers are configured, the region default when they're not, and an explicit
dada2.trim_left must always win.
"""

from pathlib import Path

import microfgt.stages.registry as reg
from microfgt.stages.model import StageContext


class _FakeRecord:
    def to_dict(self):
        return {"tool": "dada2"}


def _run_denoise_capturing(monkeypatch, tmp_path, reads_cfg):
    captured = {}

    def fake_run_dada2(*args, trunc_len=None, trim_left=None, rscript="Rscript", **kw):
        captured["trunc_len"] = trunc_len
        captured["trim_left"] = trim_left
        return _FakeRecord()

    monkeypatch.setattr(reg, "run_dada2", fake_run_dada2, raising=False)
    monkeypatch.setattr("microfgt.orchestrate.dada2.run_dada2", fake_run_dada2, raising=False)

    workdir = tmp_path
    artifacts = {k: str(tmp_path / f"{k}") for k in
                 ("trimmed_reads", "asv_table", "asv_seqs", "quality_profile")}
    ctx = StageContext(Path(workdir), artifacts, {"composition": {"reads": reads_cfg}})
    reg._run_denoise(ctx)
    return captured


def test_trim_left_zeroed_when_primers_configured(monkeypatch, tmp_path):
    cap = _run_denoise_capturing(monkeypatch, tmp_path, {
        "region": "V4V4",                                   # region default trim_left = [19,20]
        "primers": {"fwd": "GTGYCAGCMGCCGCGGTAA", "rev": "GGACTACNVGGGTWTCTAAT"},
    })
    assert cap["trim_left"] == [0, 0]                        # NOT [19,20]
    assert cap["trunc_len"] == [230, 180]                    # region trunc_len still applies


def test_trim_left_uses_region_default_when_no_primers(monkeypatch, tmp_path):
    cap = _run_denoise_capturing(monkeypatch, tmp_path, {"region": "V4V4"})   # no primers block
    assert cap["trim_left"] == [19, 20]                     # cutadapt was a no-op -> DADA2 trims


def test_explicit_trim_left_always_wins(monkeypatch, tmp_path):
    cap = _run_denoise_capturing(monkeypatch, tmp_path, {
        "region": "V4V4",
        "primers": {"fwd": "GTGYCAGCMGCCGCGGTAA"},
        "dada2": {"trim_left": [5, 6]},
    })
    assert cap["trim_left"] == [5, 6]
