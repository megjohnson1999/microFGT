"""Metagenomics arm end-to-end: raw reads -> .h5mu via stub tools.

Stubs stand in for fastp / minimap2+samtools / VIRGO2.py / VISTA(Rscript) so the WIRING of the
whole ladder (resolve -> local executor -> file artifacts -> import/classify -> MuData) is
validated end to end. Real tool correctness is the deferred HTCF validation IOU — stubs prove
plumbing, not biology (same discipline as the 16S ladder test).
"""

import stat
import sys

import mudata as md
import pytest
import yaml

pytestmark = pytest.mark.filterwarnings("ignore::FutureWarning")

# fastp: copy in1/in2 -> out1/out2, write the qc.json.
STUB_FASTP = '''#!/usr/bin/env python3
import argparse, shutil
p = argparse.ArgumentParser()
for f in ("--in1","--in2","--out1","--out2","--json","--html","--thread"): p.add_argument(f)
a, _ = p.parse_known_args()
shutil.copy(a.in1, a.out1); shutil.copy(a.in2, a.out2)
open(a.json, "w").write("{}")
'''

# minimap2: emit a dummy SAM line to stdout (samtools stubs pass it through).
STUB_MINIMAP2 = '''#!/usr/bin/env python3
import sys
sys.stdout.write("@HD\\tVN:1.6\\nread1\\t77\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tACGT\\tIIII\\n")
'''

# samtools: view/sort pass stdin->stdout; fastq writes dummy non-host FASTQs to -1/-2.
STUB_SAMTOOLS = '''#!/usr/bin/env python3
import sys, argparse
sub = sys.argv[1] if len(sys.argv) > 1 else ""
if sub in ("view", "sort"):
    sys.stdin.read()                       # drain the pipe
    sys.stdout.write("passthrough\\n")
elif sub == "fastq":
    p = argparse.ArgumentParser()
    p.add_argument("-1", dest="one"); p.add_argument("-2", dest="two")
    a, _ = p.parse_known_args(sys.argv[2:])
    try: sys.stdin.read()
    except Exception: pass
    for path in (a.one, a.two):
        if path and path != "/dev/null":
            open(path, "w").write("@r\\nACGT\\n+\\nIIII\\n")
'''

# VIRGO2.py: `map` writes <id>.out (gene\\tcount); `compile` stacks *.out -> the wide summary.
STUB_VIRGO2 = '''#!/usr/bin/env python3
import argparse, glob, os
sub = __import__("sys").argv[1]
p = argparse.ArgumentParser()
p.add_argument("-r"); p.add_argument("-o"); p.add_argument("-p"); p.add_argument("-i")
a, _ = p.parse_known_args(__import__("sys").argv[2:])
if sub == "map":
    open(a.o + ".out", "w").write("GENE1\\t100\\nGENE2\\t50\\n")
elif sub == "compile":
    outs = sorted(glob.glob(os.path.join(a.i, "*.out")))
    samples = [os.path.basename(f)[:-4] for f in outs]
    counts = {}
    for f, s in zip(outs, samples):
        for line in open(f):
            g, c = line.split()
            counts.setdefault(g, {})[s] = c
    with open(a.o + ".summary.NR.txt", "w") as fh:
        fh.write("Gene\\t" + "\\t".join(samples) + "\\n")
        for g in sorted(counts):
            fh.write(g + "\\t" + "\\t".join(counts[g].get(s, "0") for s in samples) + "\\n")
'''

# run_VISTA.R: read the compiled header for samples, write mgCSTs_TEST.csv in CWD.
STUB_RUN_VISTA = '''#!/usr/bin/env python3
import sys
compiled = sys.argv[2]        # argv: [run_VISTA.R, compiled, vista_repo]
samples = open(compiled).readline().rstrip("\\n").split("\\t")[1:]
with open("mgCSTs_TEST.csv", "w") as fh:
    fh.write(",mgCST,max_YC_theta\\n")
    for s in samples:
        fh.write(f"{s},mgCST 1,0.95\\n")
'''

TAXON_ANNOTATION = "Cluster\tGene\tTaxa\tCat\n1\tGENE1\tLactobacillus_iners\tMR\n2\tGENE2\tGardnerella_vaginalis\tMR\n"


def _exe(path, body):
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_reads_to_h5mu_via_stub_tools(tmp_path):
    from microfgt.cli import main

    raw = tmp_path / "raw"; raw.mkdir()
    for s in ("sampleA", "sampleB"):
        (raw / f"{s}_R1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (raw / f"{s}_R2.fastq").write_text("@r\nTTTT\n+\nIIII\n")

    fastp = _exe(tmp_path / "fastp", STUB_FASTP)
    minimap2 = _exe(tmp_path / "minimap2", STUB_MINIMAP2)
    samtools = _exe(tmp_path / "samtools", STUB_SAMTOOLS)
    rscript = _exe(tmp_path / "Rscript", STUB_RUN_VISTA)  # stub 'Rscript' == the VISTA runner

    host_ref = tmp_path / "GRCh38.fna.gz"; host_ref.write_text("")
    virgo2_dir = tmp_path / "VIRGO2"; virgo2_dir.mkdir()
    _exe(virgo2_dir / "VIRGO2.py", STUB_VIRGO2)
    (virgo2_dir / "AnnotationTables").mkdir()
    (virgo2_dir / "AnnotationTables" / "1.VIRGO2.taxon.txt").write_text(TAXON_ANNOTATION)
    vista_repo = tmp_path / "VISTA"; vista_repo.mkdir()
    _exe(vista_repo / "run_VISTA.R", STUB_RUN_VISTA)

    out = tmp_path / "result.h5mu"
    config = {
        "metagenomics": {
            "reads": {"fastq_dir": str(raw)},
            "host_ref": str(host_ref), "virgo2_dir": str(virgo2_dir),
            "vista_repo": str(vista_repo), "threads": 2,
            "fastp": str(fastp), "minimap2": str(minimap2), "samtools": str(samtools),
            "python": sys.executable, "rscript": str(rscript),
        },
        "output": str(out),
    }
    cfg = tmp_path / "cfg.yaml"; cfg.write_text(yaml.safe_dump(config))

    assert main(["run", "-c", str(cfg), "--workdir", str(tmp_path / "wd")]) == 0
    assert out.exists()

    m = md.read(out)
    # Both shotgun modalities landed; the tool chain threaded both samples through.
    assert set(m.mod.keys()) == {"function", "composition_taxon_shotgun"}
    assert m["function"].n_obs == 2
    assert set(m["function"].var_names) == {"GENE1", "GENE2"}
    # Genes rolled up to taxa via the annotation join.
    assert set(m["composition_taxon_shotgun"].var_names) == {
        "Lactobacillus_iners", "Gardnerella_vaginalis"
    }
    # mgCST call + score on the global frame; shotgun descriptors beside it.
    assert list(m.obs["mgCST"]) == ["mgCST 1", "mgCST 1"]
    assert "shotgun_dominant_taxon" in m.obs.columns
    # Provenance for every tool stage folded into uns (stored as a JSON string).
    import json

    runs = json.loads(m.uns["tool_runs"])
    assert {
        "sg_qc", "sg_host_removal", "sg_virgo2_map", "sg_virgo2_compile", "classify_mgcst"
    } <= set(runs)
    # VISTA is a real tool invocation, so its run is recorded like the others.
    assert runs["classify_mgcst"][0]["tool"] == "VISTA"


def test_combined_16s_and_metagenomics_from_existing_outputs(tmp_path, real_fixtures, test_data):
    # No subprocess: 16S enters at speciateIT results, metagenomics at a compiled matrix, VISTA
    # via an existing output. The resolver must pick integrate_combined and build both stacks.
    from microfgt.cli import main

    virgo2_dir = tmp_path / "VIRGO2"; (virgo2_dir / "AnnotationTables").mkdir(parents=True)
    (virgo2_dir / "AnnotationTables" / "1.VIRGO2.taxon.txt").write_bytes(
        (real_fixtures / "virgo2_taxon_annotation.slice.txt").read_bytes()
    )
    out = tmp_path / "combined.h5mu"
    config = {
        "composition": {"speciateit": {
            "results": str(test_data / "speciateit_MC_order7_results.synthetic.txt"),
            "count_table": str(real_fixtures / "speciateit_test_count_table.csv"),
        }},
        "cst": {"method": "centroid"},
        "metagenomics": {
            "compiled": str(real_fixtures / "virgo2_compiled.summary.NR.slice.txt"),
            "virgo2_dir": str(virgo2_dir),
        },
        "mgcst": {"vista_output": str(real_fixtures / "vista_mgCSTs.csv")},
        "output": str(out),
    }
    cfg = tmp_path / "cfg.yaml"; cfg.write_text(yaml.safe_dump(config))

    assert main(["run", "-c", str(cfg), "--workdir", str(tmp_path / "wd")]) == 0
    m = md.read(out)
    # Both arms present in one object.
    assert {"composition", "composition_taxon"} <= set(m.mod.keys())       # 16S
    assert {"function", "composition_taxon_shotgun"} <= set(m.mod.keys())  # shotgun
    assert "CST" in m.obs.columns and "mgCST" in m.obs.columns             # both calls, separate
