"""VIRGO orchestration — full run->import pipeline via a stub mapping script.

The real mapping needs VIRGO's hosted catalog + bowtie/blast/seqtk, so these tests drive a
stub that emits VIRGO's documented per-sample format (`temp_mapping/<prefix>.out`, 3-col
geneID/read_count/gene_length, no header). Validates the subprocess plumbing, provenance,
and the handoff to import_virgo — not the real mapper.
"""

import stat

import pytest

from microfgt.orchestrate import run_virgo, run_virgo_samples

# Stub runMapping.step1.sh: parse -r/-p/-d, write temp_mapping/<prefix>.out (3-col, no header).
FAKE_STEP1 = r'''#!/usr/bin/env bash
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r) READS="$2"; shift 2;;
    -p) PREFIX="$2"; shift 2;;
    -d) DIR="$2"; shift 2;;
    *) shift;;
  esac
done
mkdir -p temp_mapping
printf 'V1593031\t1417\t3663\nV1607456\t333\t390\n' > "temp_mapping/${PREFIX}.out"
'''


@pytest.fixture
def fake_step1(tmp_path):
    script = tmp_path / "runMapping.step1.sh"
    script.write_text(FAKE_STEP1)
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def test_run_virgo_single_sample(fake_step1, tmp_path):
    workdir = tmp_path / "work"
    out_path, record = run_virgo(
        reads=tmp_path / "sub1.fq",          # stub ignores contents
        prefix="sub1",
        virgo_path=tmp_path / "VIRGO",
        workdir=workdir,
        script=fake_step1,
    )
    assert out_path == workdir / "temp_mapping" / "sub1.out"
    assert out_path.exists()
    assert record.tool == "VIRGO"
    assert record.returncode == 0
    assert "-p" in record.argv and "sub1" in record.argv


def test_run_virgo_samples_stacks_into_anndata(fake_step1, tmp_path):
    workdir = tmp_path / "work"
    adata = run_virgo_samples(
        {"sub1": tmp_path / "sub1.fq", "sub2": tmp_path / "sub2.fq"},
        virgo_path=tmp_path / "VIRGO",
        workdir=workdir,
        script=fake_step1,
    )
    assert adata.n_obs == 2                                  # two samples stacked
    assert set(adata.obs_names) == {"sub1", "sub2"}
    assert "V1593031" in adata.var_names
    assert set(adata.uns["virgo_runs"]) == {"sub1", "sub2"}  # provenance per sample


def test_missing_script_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="VIRGO mapping script"):
        run_virgo(
            reads=tmp_path / "r.fq", prefix="s", virgo_path=tmp_path / "VIRGO",
            workdir=tmp_path / "w", script=tmp_path / "nope.sh",
        )
