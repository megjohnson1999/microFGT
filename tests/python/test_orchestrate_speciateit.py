"""speciateIT orchestration — full run->import pipeline via a stub `classify`.

The real classifier needs ~2.6 GB of vSpeciateDB models, so these tests drive a stub
executable that emits speciateIT's DOCUMENTED genuine output format
(`Seq\\tClassification\\tpp\\tnDecisions`, one row per ASV). This validates the subprocess
plumbing, provenance recording, and the handoff to import_speciateit against the real
header shape — but it is NOT a real-data validation of the classifier (that remains the
P1 IOU, dischargeable by running the real tool on test.fasta).
"""

import os
import stat

import pandas as pd
import pytest

from microfgt.io import collapse_to_taxon
from microfgt.orchestrate import ToolNotFoundError, run_speciateit

# A stub `classify`: parses speciateIT's flags, emits the documented genuine format.
FAKE_CLASSIFY = '''#!/usr/bin/env python3
import argparse
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument("-d"); p.add_argument("-i"); p.add_argument("-o")
p.add_argument("--skip-err-thld", action="store_true")
a, _ = p.parse_known_args()
ids = [ln[1:].strip() for ln in Path(a.i).read_text().splitlines() if ln.startswith(">")]
outdir = Path(a.o); outdir.mkdir(parents=True, exist_ok=True)
taxa = ["Lactobacillus_iners", "Gardnerella_vaginalis"]
with open(outdir / "MC_order7_results.txt", "w") as f:
    f.write("Seq\\tClassification\\tpp\\tnDecisions\\n")   # documented genuine header
    for i, asv in enumerate(ids):
        f.write(f"{asv}\\t{taxa[i % 2]}\\t0.97\\t50\\n")
'''


@pytest.fixture
def fake_classify(tmp_path):
    exe = tmp_path / "classify"
    exe.write_text(FAKE_CLASSIFY)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return exe


def test_run_speciateit_end_to_end(real_fixtures, fake_classify, tmp_path):
    adata = run_speciateit(
        fasta=real_fixtures / "speciateit_test.fasta",
        db=tmp_path / "fake_db",            # stub ignores it
        outdir=tmp_path / "out",
        executable=str(fake_classify),
        count_table=real_fixtures / "speciateit_test_count_table.csv",
    )
    # ASV grain: one feature per ASV; sequences retained from the FASTA (run_speciateit
    # hands its input FASTA to the importer). ASV1-10 classified (alternating), ASV11.. not.
    assert adata.n_vars == 1514        # every ASV column from the count table
    assert adata.n_obs == 169          # samples come from the count table
    assert adata.var.loc["ASV1", "classification"] in {"Lactobacillus_iners", "Gardnerella_vaginalis"}
    assert pd.isna(adata.var.loc["ASV11", "classification"])
    assert adata.var.loc["ASV1", "sequence"].startswith("TAGGGAAT")
    # The taxon roll-up yields the classified taxa + the Unclassified bucket.
    assert set(collapse_to_taxon(adata).var_names) == {
        "Lactobacillus_iners", "Gardnerella_vaginalis", "Unclassified",
    }

    # Provenance recorded (constraint B).
    rec = adata.uns["speciateit_run"]
    assert rec["tool"] == "speciateIT"
    assert rec["returncode"] == 0
    assert rec["argv"][1:] == ["-d", str(tmp_path / "fake_db"),
                               "-i", str(real_fixtures / "speciateit_test.fasta"),
                               "-o", str(tmp_path / "out")]
    assert rec["exe_fingerprint"]["size"] > 0


def test_skip_err_thld_flag_passed(real_fixtures, fake_classify, tmp_path):
    results_path, record = run_speciateit(
        fasta=real_fixtures / "speciateit_test.fasta",
        db=tmp_path / "db",
        outdir=tmp_path / "out",
        executable=str(fake_classify),
        skip_err_thld=True,
    )
    assert results_path.exists()
    assert "--skip-err-thld" in record.argv
    assert record.params["skip_err_thld"] is True


def test_missing_executable_raises_helpful_error(real_fixtures, tmp_path):
    with pytest.raises(ToolNotFoundError, match="speciateIT"):
        run_speciateit(
            fasta=real_fixtures / "speciateit_test.fasta",
            db=tmp_path / "db",
            outdir=tmp_path / "out",
            executable="definitely_not_installed_xyz123",
        )
