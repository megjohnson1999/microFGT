"""Resolver — the multi-entry ladder: the entry point is just which inputs are present."""

import pytest

from microfgt.stages import StageResolutionError, provided_artifacts, resolve
from microfgt.stages.executors import SnakemakeExecutor


def _plan(config):
    return [s.id for s in resolve("mudata", provided_artifacts(config))]


def test_entry_at_fastqs_runs_full_ladder():
    cfg = {"composition": {"reads": {"fastq_dir": "raw/"}, "speciateit": {"db": "d"}}}
    assert _plan(cfg) == [
        "primer_trim", "denoise", "assign", "import_composition", "cst_classify", "integrate",
    ]


def test_entry_at_asv_table_skips_denoising():
    cfg = {"composition": {"asv_table": "a.csv", "asv_seqs": "a.fasta", "speciateit": {"db": "d"}}}
    assert _plan(cfg) == ["assign", "import_composition", "cst_classify", "integrate"]


def test_entry_at_existing_outputs_skips_tool_running():
    cfg = {"composition": {"speciateit": {"results": "r.txt", "count_table": "c.csv"}}}
    assert _plan(cfg) == ["import_composition", "cst_classify", "integrate"]


def test_valencia_labels_select_the_valencia_cst_producer():
    cfg = {
        "composition": {"speciateit": {"results": "r.txt", "count_table": "c.csv"}},
        "cst": {"valencia": "v.csv"},
    }
    plan = _plan(cfg)
    assert "cst_valencia" in plan and "cst_classify" not in plan


def test_entry_at_phyloseq_trusts_the_objects_cst():
    # A phyloseq .rds enters at import_phyloseq and, by default, uses its own CST.
    cfg = {"composition": {"phyloseq": "ps.rds"}}
    assert _plan(cfg) == ["import_phyloseq", "cst_phyloseq", "integrate"]


def test_phyloseq_entry_reclassifies_when_method_set():
    # An explicit cst.method overrides the object's CST (classify from the composition).
    cfg = {"composition": {"phyloseq": "ps.rds"}, "cst": {"method": "centroid"}}
    plan = _plan(cfg)
    assert plan == ["import_phyloseq", "cst_classify", "integrate"]
    assert "cst_phyloseq" not in plan


def test_phyloseq_entry_honours_valencia_labels():
    cfg = {"composition": {"phyloseq": "ps.rds"}, "cst": {"valencia": "v.csv"}}
    plan = _plan(cfg)
    assert "cst_valencia" in plan and "cst_phyloseq" not in plan


def test_unproducible_target_raises():
    with pytest.raises(StageResolutionError, match="Cannot produce"):
        resolve("mudata", set())            # nothing provided, no way to reach a composition


def test_snakefile_generated_from_same_registry(tmp_path):
    cfg_path = tmp_path / "cfg.yaml"
    import yaml

    config = {"composition": {"reads": {"fastq_dir": "raw/"}, "speciateit": {"db": "d"}},
              "output": str(tmp_path / "out.h5mu")}
    cfg_path.write_text(yaml.safe_dump(config))

    stages = resolve("mudata", provided_artifacts(config))
    text = SnakemakeExecutor().generate(stages, str(cfg_path), str(tmp_path / "wd"))

    assert "rule all:" in text
    assert "rule primer_trim:" in text and "rule integrate:" in text
    # one shell-out per stage, no re-encoding of commands in the Snakefile
    assert text.count('shell: "microfgt _run-stage') == len(stages)
