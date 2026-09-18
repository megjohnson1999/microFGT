"""Full multi-entry ladder, FASTQs -> .h5mu, via stub tools.

Stubs stand in for cutadapt / Rscript(DADA2) / speciateIT so the WIRING of the whole ladder
(resolve -> local executor -> file artifacts -> import/CST/analysis -> MuData) is validated
end to end. Real cutadapt/DADA2/speciateIT correctness is the deferred HTCF validation IOU
(the same real-run discipline as P3) — stubs prove plumbing, not biology.
"""

import stat

import mudata as md
import pytest
import yaml

pytestmark = pytest.mark.filterwarnings("ignore::FutureWarning")

# Stub cutadapt: copy each input FASTQ to the -o/-p outputs (no real trimming).
STUB_CUTADAPT = '''#!/usr/bin/env python3
import argparse, shutil
p = argparse.ArgumentParser()
p.add_argument("-g"); p.add_argument("-G"); p.add_argument("-o"); p.add_argument("-p")
a, rest = p.parse_known_args()
ins = [x for x in rest if not x.startswith("-")]
shutil.copy(ins[0], a.o)
if a.p and len(ins) > 1:
    shutil.copy(ins[1], a.p)
'''

# Stub Rscript(DADA2): emit an ASV table (samples from trimmed _R1 files) + rep-seqs + qprofile.
STUB_RSCRIPT = '''#!/usr/bin/env python3
import argparse, glob, os
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument("--input"); p.add_argument("--asv-table"); p.add_argument("--asv-seqs")
p.add_argument("--quality-profile"); p.add_argument("--trunc-len"); p.add_argument("--trim-left")
a, _ = p.parse_known_args()
samples = sorted(os.path.basename(f).split("_R1")[0] for f in glob.glob(os.path.join(a.input, "*_R1*")))
with open(a.asv_table, "w") as f:
    f.write("sampleID,ASV1,ASV2\\n")
    for i, s in enumerate(samples):
        f.write(f"{s},{100 + i},{10 + i}\\n")
Path(a.asv_seqs).write_text(">ASV1\\nACGT\\n>ASV2\\nTTTT\\n")
Path(a.quality_profile).write_text("file\\tcycle\\tmean_quality\\n")
'''

# Stub speciateIT classify: emit the documented genuine output, classifying ASV1/ASV2.
STUB_CLASSIFY = '''#!/usr/bin/env python3
import argparse
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument("-d"); p.add_argument("-i"); p.add_argument("-o")
a, _ = p.parse_known_args()
ids = [ln[1:].strip() for ln in Path(a.i).read_text().splitlines() if ln.startswith(">")]
out = Path(a.o); out.mkdir(parents=True, exist_ok=True)
taxa = ["Lactobacillus_iners", "Gardnerella_vaginalis"]
with open(out / "MC_order7_results.txt", "w") as f:
    f.write("Seq\\tClassification\\tpp\\tnDecisions\\n")
    for i, asv in enumerate(ids):
        f.write(f"{asv}\\t{taxa[i % 2]}\\t0.97\\t50\\n")
'''


def _exe(path, body):
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_full_fastq_to_h5mu_ladder(tmp_path):
    from microfgt.cli import main

    # Two samples of (placeholder) paired FASTQs.
    fq = tmp_path / "raw"; fq.mkdir()
    for s in ("sampleA", "sampleB"):
        (fq / f"{s}_R1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (fq / f"{s}_R2.fastq").write_text("@r\nTTTT\n+\nIIII\n")

    cutadapt = _exe(tmp_path / "cutadapt", STUB_CUTADAPT)
    rscript = _exe(tmp_path / "Rscript", STUB_RSCRIPT)
    classify = _exe(tmp_path / "classify", STUB_CLASSIFY)
    db = tmp_path / "vSpeciateIT_V3V4"; db.mkdir()
    out = tmp_path / "result.h5mu"

    config = {
        "composition": {
            "reads": {
                "fastq_dir": str(fq), "region": "V3V4",
                "primers": {"fwd": "AAAA", "rev": "TTTT"},
                "cutadapt": str(cutadapt), "rscript": str(rscript),
            },
            "speciateit": {"db": str(db), "classify": str(classify)},
        },
        "cst": {"method": "centroid"},
        "analysis": {"transforms": ["relabund", "clr"], "alpha": ["shannon"]},
        "output": str(out),
    }
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(yaml.safe_dump(config))

    rc = main(["run", "-c", str(cfg), "--workdir", str(tmp_path / "wd")])
    assert rc == 0
    assert out.exists()

    # Intermediate artifacts landed in the workdir (the file-artifact backbone).
    assert (tmp_path / "wd" / "asv_table.csv").exists()
    assert (tmp_path / "wd" / "quality_profile.tsv").exists()

    m = md.read(out)
    assert m["composition"].n_obs == 2                       # both samples threaded through
    assert set(m["composition"].var_names) == {"ASV1", "ASV2"}   # ASV grain, not collapsed
    assert "sequence" in m["composition"].var               # sequences carried through the ladder
    # Taxon roll-up materialised as its own assay.
    assert set(m["composition_taxon"].var_names) <= {"Lactobacillus_iners", "Gardnerella_vaginalis"}
    assert "CST" in m.obs.columns                            # CST classified end to end
    assert "dominant_taxon" in m.obs.columns                 # augment descriptors present
    assert "clr" in m["composition"].layers
    assert "alpha_shannon" in m["composition"].obs
    # Every 16S tool stage recorded its run under the arm-agnostic provenance key.
    import json

    runs = json.loads(m.uns["tool_runs"])
    assert {"primer_trim", "denoise", "assign"} <= set(runs)
    assert runs["primer_trim"][0]["tool"] == "cutadapt"     # one record per sample
    assert len(runs["primer_trim"]) == 2


def test_full_ladder_from_samplesheet(tmp_path):
    """A run driven by `samples:` — messy source filenames, canonical ids in the sheet — comes
    out keyed by the sheet's sample_ids with the sheet's metadata attached to .obs."""
    from microfgt.cli import main

    # Deliberately messy, non-canonical source names for two samples.
    raw = tmp_path / "raw"; raw.mkdir()
    srcs = {}
    for canon, ugly in (("PT01", "run5_42563_S3"), ("PT02", "run5_42981_S7")):
        r1 = raw / f"{ugly}_R1_001.fastq"; r1.write_text("@r\nACGT\n+\nIIII\n")
        r2 = raw / f"{ugly}_R2_001.fastq"; r2.write_text("@r\nTTTT\n+\nIIII\n")
        srcs[canon] = (r1, r2)

    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "sample_id,16s_R1,16s_R2,group\n"
        + "".join(f"{c},{r1},{r2},{g}\n"
                 for (c, (r1, r2)), g in zip(srcs.items(), ("BV", "Normal")))
    )

    cutadapt = _exe(tmp_path / "cutadapt", STUB_CUTADAPT)
    rscript = _exe(tmp_path / "Rscript", STUB_RSCRIPT)
    classify = _exe(tmp_path / "classify", STUB_CLASSIFY)
    db = tmp_path / "vSpeciateIT_V3V4"; db.mkdir()
    out = tmp_path / "result.h5mu"

    config = {
        "samples": str(sheet),
        "composition": {
            "reads": {
                "region": "V3V4", "primers": {"fwd": "AAAA", "rev": "TTTT"},
                "cutadapt": str(cutadapt), "rscript": str(rscript),
            },
            "speciateit": {"db": str(db), "classify": str(classify)},
        },
        "cst": {"method": "centroid"},
        "output": str(out),
    }
    cfg = tmp_path / "cfg.yaml"; cfg.write_text(yaml.safe_dump(config))

    assert main(["run", "-c", str(cfg), "--workdir", str(tmp_path / "wd")]) == 0
    m = md.read(out)
    assert set(m["composition"].obs_names) == {"PT01", "PT02"}   # canonical ids, not the ugly ones
    assert m.obs.loc["PT01", "group"] == "BV"                    # sheet metadata -> .obs
    assert m.obs.loc["PT02", "group"] == "Normal"


def test_run_emits_snakefile_in_snakemake_mode(tmp_path):
    from microfgt.cli import main

    config = {"composition": {"speciateit": {"results": "r.txt", "count_table": "c.csv"}},
              "output": str(tmp_path / "o.h5mu")}
    cfg = tmp_path / "cfg.yaml"; cfg.write_text(yaml.safe_dump(config))
    wd = tmp_path / "wd"
    assert main(["run", "-c", str(cfg), "--workdir", str(wd), "--executor", "snakemake"]) == 0
    assert (wd / "Snakefile").exists()
    assert "rule integrate:" in (wd / "Snakefile").read_text()


def test_snakemake_mode_supports_sample_sheets(tmp_path):
    """A `samples:` config now works in snakemake mode: the sheet is staged, a staged config is
    written, and the Snakefile's `_run-stage` calls point at that staged config (so they see the
    rewritten FASTQ entry points, not the original `samples:` config)."""
    from microfgt.cli import main

    raw = tmp_path / "raw"; raw.mkdir()
    for s in ("PT01", "PT02"):
        (raw / f"{s}_R1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (raw / f"{s}_R2.fastq").write_text("@r\nTTTT\n+\nIIII\n")
    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "sample_id,16s_R1,16s_R2,group\n"
        f"PT01,{raw/'PT01_R1.fastq'},{raw/'PT01_R2.fastq'},BV\n"
        f"PT02,{raw/'PT02_R1.fastq'},{raw/'PT02_R2.fastq'},Normal\n"
    )
    config = {
        "samples": str(sheet),
        "composition": {"reads": {"region": "V3V4", "primers": {"fwd": "AAAA", "rev": "TTTT"}}},
        "output": str(tmp_path / "o.h5mu"),
    }
    cfg = tmp_path / "cfg.yaml"; cfg.write_text(yaml.safe_dump(config))
    wd = tmp_path / "wd"

    # Without --workdir a sample-sheet snakemake run is refused (staged reads would be temp).
    with pytest.raises(SystemExit):
        main(["run", "-c", str(cfg), "--executor", "snakemake"])

    assert main(["run", "-c", str(cfg), "--workdir", str(wd), "--executor", "snakemake"]) == 0
    staged_cfg = wd / "config.staged.yaml"
    assert staged_cfg.exists()
    staged = yaml.safe_load(staged_cfg.read_text())
    assert "samples" not in staged                                   # sheet resolved away
    assert staged["composition"]["reads"]["fastq_dir"] == str(wd / "staged" / "16s")
    text = (wd / "Snakefile").read_text()
    assert str(staged_cfg) in text                                   # rules read the staged config
    assert str(cfg) not in text                                      # not the original


def test_snakefile_scatters_shotgun_per_sample_and_gathers(tmp_path):
    """The shotgun scatter stages emit one rule per sample (true cluster fan-out), the compile
    stage gathers every per-sample .out, a non-scatter dir artifact stays directory()-wrapped, and
    every interpolated shell path is quoted (survives spaces)."""
    import shlex

    from microfgt.stages.executors import SnakemakeExecutor
    from microfgt.stages.resolve import resolve

    # Two shotgun samples on disk so the generator can discover the scatter set.
    sg = tmp_path / "sg"; sg.mkdir()
    for s in ("SG1", "SG2"):
        (sg / f"{s}_R1.fastq.gz").write_bytes(b"")
        (sg / f"{s}_R2.fastq.gz").write_bytes(b"")

    stages = resolve("mudata", {"fastq_dir", "sg_reads"})
    workdir = tmp_path / "work dir with space"      # bug #1: a spaced path
    config = {
        "composition": {"reads": {"fastq_dir": str(tmp_path / "16s")}},
        "metagenomics": {"reads": {"fastq_dir": str(sg)}},
        "output": str(tmp_path / "o.h5mu"),
    }
    cfg = tmp_path / "cfg.yaml"; cfg.write_text(yaml.safe_dump(config))
    text = SnakemakeExecutor().generate(stages, str(cfg), str(workdir))

    # Per-sample scatter rules, one per (scatter stage, sample), each pinned to --sample.
    for stage in ("sg_qc", "sg_host_removal", "sg_virgo2_map"):
        for s in ("SG1", "SG2"):
            assert f"rule {stage}_{s}:" in text
    assert "--sample SG1" in text and "--sample SG2" in text
    # sg_virgo2_map's per-sample output is the sample's .out (not a wrapped directory).
    assert str(workdir / "mg_virgo2_out" / "SG1.out") in text
    assert "directory(" not in text.split("rule sg_virgo2_map_SG1:")[1].split("rule ")[0]

    # The compile (gather) rule pulls in every per-sample .out.
    compile_block = text.split("rule sg_virgo2_compile:")[1].split("rule ")[0]
    assert str(workdir / "mg_virgo2_out" / "SG1.out") in compile_block
    assert str(workdir / "mg_virgo2_out" / "SG2.out") in compile_block

    # A non-scatter directory artifact (16S primer-trim output) is still directory()-wrapped.
    assert f"directory({str(workdir / 'trimmed')!r})" in text

    # Every shell line quotes the spaced workdir (bug #1) — never bare.
    shell_lines = [ln for ln in text.splitlines() if ln.strip().startswith("shell:")]
    assert shell_lines
    for line in shell_lines:
        assert f"--workdir {shlex.quote(str(workdir))}" in line
        assert f"--workdir {workdir} " not in line


def test_no_file_artifact_nested_under_a_directory_artifact():
    """Invariant: no file artifact's path may live inside a directory artifact's path. Snakemake
    rejects a rule output nested in another rule's directory() output (ChildIOException), which is
    what sank sg_compiled (it used to live inside the sg_virgo2_out map dir)."""
    from pathlib import PurePosixPath

    from microfgt.stages.model import ARTIFACT_FILENAMES, DIRECTORY_ARTIFACTS

    dirs = {k: PurePosixPath(ARTIFACT_FILENAMES[k]) for k in DIRECTORY_ARTIFACTS}
    for key, fn in ARTIFACT_FILENAMES.items():
        if fn is None or key in DIRECTORY_ARTIFACTS:
            continue
        p = PurePosixPath(fn)
        for dkey, dpath in dirs.items():
            assert not p.is_relative_to(dpath), (
                f"artifact {key!r} ({fn}) is nested inside directory artifact {dkey!r} ({dpath}) "
                "— Snakemake will raise ChildIOException; give it its own top-level dir"
            )
