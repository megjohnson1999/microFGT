"""Metagenomics orchestrators — the fragile bits, in isolation (stub tools).

The end-to-end ladder exercises these through the stage path; here we pin the two behaviors
that are easy to get subtly wrong: the host-removal shell pipeline surviving paths with spaces
(the real dev/HPC environment has them), and VIRGO2's single-end R1+R2 concatenation.
"""

import stat

import pytest

STUB_MINIMAP2 = '''#!/usr/bin/env python3
import sys
sys.stdout.write("read\\n")
'''
STUB_SAMTOOLS = '''#!/usr/bin/env python3
import sys, argparse
sub = sys.argv[1] if len(sys.argv) > 1 else ""
if sub in ("view", "sort"):
    sys.stdin.read(); sys.stdout.write("x\\n")
elif sub == "fastq":
    p = argparse.ArgumentParser(); p.add_argument("-1", dest="one"); p.add_argument("-2", dest="two")
    a, _ = p.parse_known_args(sys.argv[2:])
    try: sys.stdin.read()
    except Exception: pass
    for path in (a.one, a.two):
        if path and path != "/dev/null": open(path, "w").write("@r\\nA\\n+\\nI\\n")
'''
# VIRGO2.py map that reveals its -r input by copying it to <id>.out (so we can inspect the concat).
STUB_VIRGO2_REVEAL = '''#!/usr/bin/env python3
import sys, argparse, shutil
p = argparse.ArgumentParser(); p.add_argument("-r"); p.add_argument("-o"); p.add_argument("-p")
a, _ = p.parse_known_args(sys.argv[2:])
shutil.copy(a.r, a.o + ".out")
'''


def _exe(path, body):
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_host_removal_survives_spaces_in_paths(tmp_path):
    from microfgt.orchestrate import run_host_removal

    # A directory with a space in its name — the pipeline is a shell string, so every path
    # must be shlex.quote-d or this breaks.
    base = tmp_path / "a dir with spaces"; base.mkdir()
    trimmed = base / "trimmed"; trimmed.mkdir()
    (trimmed / "s1_R1.fastq").write_text("@r\nA\n+\nI\n")
    (trimmed / "s1_R2.fastq").write_text("@r\nT\n+\nI\n")
    host_ref = base / "host ref.fna"; host_ref.write_text("")
    nonhost = base / "non host"

    records = run_host_removal(
        trimmed, nonhost, host_ref,
        minimap2=str(_exe(tmp_path / "minimap2", STUB_MINIMAP2)),
        samtools=str(_exe(tmp_path / "samtools", STUB_SAMTOOLS)),
    )
    assert records[0].returncode == 0
    assert (nonhost / "s1_R1.fastq").exists() and (nonhost / "s1_R2.fastq").exists()


def test_virgo2_map_concatenates_r1_then_r2(tmp_path):
    from microfgt.orchestrate import run_virgo2_map

    r1 = tmp_path / "s1_R1.fastq"; r1.write_text("R1DATA\n")
    r2 = tmp_path / "s1_R2.fastq"; r2.write_text("R2DATA\n")
    virgo2_dir = tmp_path / "VIRGO2"; virgo2_dir.mkdir()
    _exe(virgo2_dir / "VIRGO2.py", STUB_VIRGO2_REVEAL)
    outdir = tmp_path / "out"

    out_path, record = run_virgo2_map(r1, r2, "s1", virgo2_dir, outdir, python="python3")
    assert out_path == outdir / "s1.out"
    # VIRGO2 is single-end: the -r input it saw is R1 followed by R2, one file.
    assert out_path.read_text() == "R1DATA\nR2DATA\n"


def test_virgo2_map_missing_script_raises(tmp_path):
    from microfgt.orchestrate import run_virgo2_map

    with pytest.raises(FileNotFoundError, match="VIRGO2.py"):
        run_virgo2_map(
            tmp_path / "s1_R1.fastq", tmp_path / "s1_R2.fastq", "s1",
            tmp_path / "no_virgo2", tmp_path / "out",
        )
