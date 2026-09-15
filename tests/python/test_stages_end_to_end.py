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


def test_run_emits_snakefile_in_snakemake_mode(tmp_path):
    from microfgt.cli import main

    config = {"composition": {"speciateit": {"results": "r.txt", "count_table": "c.csv"}},
              "output": str(tmp_path / "o.h5mu")}
    cfg = tmp_path / "cfg.yaml"; cfg.write_text(yaml.safe_dump(config))
    wd = tmp_path / "wd"
    assert main(["run", "-c", str(cfg), "--workdir", str(wd), "--executor", "snakemake"]) == 0
    assert (wd / "Snakefile").exists()
    assert "rule integrate:" in (wd / "Snakefile").read_text()
