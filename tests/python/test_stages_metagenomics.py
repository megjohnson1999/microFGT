"""Metagenomics arm — resolver routing, entry-point mapping, and doctor requirements.

Pure graph/config tests (no subprocess): they prove the shotgun chain resolves, that the three
`mudata` producers route correctly by entry shape, and that the audit's walls surface as
`microfgt check` requirements.
"""

from microfgt.stages import provided_artifacts, resolve
from microfgt.stages.check import check


def _ids(target, provided):
    return [s.id for s in resolve(target, provided)]


def test_provided_artifacts_maps_metagenomics_entries():
    cfg = {
        "metagenomics": {"reads": {"fastq_dir": "raw/"}, "virgo2_dir": "/ref/VIRGO2"},
        "mgcst": {"vista_output": "mgCSTs_x.csv"},
    }
    prov = provided_artifacts(cfg)
    assert prov["sg_reads"] == "raw/"
    assert prov["vista_output"] == "mgCSTs_x.csv"
    assert "sg_compiled" not in prov  # no `compiled:` set

    prov2 = provided_artifacts({"metagenomics": {"compiled": "VIRGO2_Compiled.summary.NR.txt"}})
    assert prov2["sg_compiled"] == "VIRGO2_Compiled.summary.NR.txt"


def test_reads_entry_runs_the_full_tool_chain_then_integrate_shotgun():
    ids = _ids("mudata", {"sg_reads"})
    assert ids == [
        "sg_qc", "sg_host_removal", "sg_virgo2_map", "sg_virgo2_compile",
        "import_function", "classify_mgcst", "integrate_shotgun",
    ]


def test_compiled_entry_skips_the_tool_chain():
    ids = _ids("mudata", {"sg_compiled"})
    assert "sg_qc" not in ids and "sg_virgo2_map" not in ids
    assert set(ids) == {"import_function", "classify_mgcst", "integrate_shotgun"}


def test_existing_vista_output_skips_running_vista():
    ids = _ids("mudata", {"sg_compiled", "vista_output"})
    assert "classify_mgcst" not in ids            # VISTA not re-run
    assert "import_mgcst_existing" in ids           # imported instead
    assert ids[-1] == "integrate_shotgun"


def test_16s_entry_is_unchanged_and_picks_plain_integrate():
    ids = _ids("mudata", {"fastq_dir"})
    assert ids[-1] == "integrate"
    assert "integrate_shotgun" not in ids and "integrate_combined" not in ids
    assert ids[0] == "primer_trim"


def test_combined_entry_routes_to_integrate_combined_and_runs_both_chains():
    ids = _ids("mudata", {"fastq_dir", "sg_compiled"})
    assert ids[-1] == "integrate_combined"
    assert {"primer_trim", "denoise", "assign"} <= set(ids)          # 16S chain ran
    assert {"import_function", "classify_mgcst"} <= set(ids)         # shotgun chain ran


def test_snakemake_mode_emits_the_shotgun_rules(tmp_path):
    import yaml

    from microfgt.cli import main

    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(yaml.safe_dump({
        "metagenomics": {"reads": {"fastq_dir": "raw/"}, "virgo2_dir": "/ref/VIRGO2",
                         "host_ref": "/ref/GRCh38.fna.gz", "vista_repo": "/ref/VISTA"},
        "output": str(tmp_path / "o.h5mu"),
    }))
    wd = tmp_path / "wd"
    assert main(["run", "-c", str(cfg), "--workdir", str(wd), "--executor", "snakemake"]) == 0
    snake = (wd / "Snakefile").read_text()
    for rule in ("sg_qc", "sg_host_removal", "sg_virgo2_map", "sg_virgo2_compile",
                 "import_function", "classify_mgcst", "integrate_shotgun"):
        assert f"rule {rule}:" in snake


def test_check_surfaces_the_metagenomics_walls():
    # Point at unset/absent DB paths; the doctor should list the tools + DB slots to satisfy.
    cfg = {"metagenomics": {
        "reads": {"fastq_dir": "raw/"},
        "host_ref": "/nope/GRCh38.fna.gz", "virgo2_dir": "/nope/VIRGO2",
        "vista_repo": "/nope/VISTA",
    }}
    msgs = " ".join(r.message for r in check(cfg))
    for tool in ("fastp", "minimap2", "samtools", "bowtie2", "Rscript"):
        assert tool in msgs
    assert "randomForestSRC" in msgs
    assert "VIRGO2.1.bt2" in msgs and "VISTA_data" in msgs and "GRCh38" in msgs
