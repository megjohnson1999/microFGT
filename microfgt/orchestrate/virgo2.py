"""Orchestrate VIRGO2 — read mapping + compile (v2), grounded in RECIPE.md.

Two gotchas the recipe pins down and this wrapper encodes:
- **VIRGO2 is single-end only**: ``map`` takes one ``-r`` file, no ``-1/-2``. We concatenate
  R1+R2 into one file first (a gene catalog ignores mate info). The concat is done in Python
  (gzip members and plain FASTQ both concatenate cleanly) — no shell needed.
- **compile is sample-list-agnostic**: it globs ``*.out`` in the dir, so mapping N samples then
  compiling "just works" (this is why the audit's 9-of-11 recovery needed no editing).

VIRGO2.py lives in the configured ``virgo2_dir``; it is run via the resolved ``python3``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from microfgt.orchestrate._run import resolve_executable, run_command

COMPILED_NAME = "VIRGO2_Compiled.summary.NR.txt"


def _virgo2_script(virgo2_dir) -> Path:
    script = Path(virgo2_dir) / "VIRGO2.py"
    if not script.exists():
        raise FileNotFoundError(
            f"VIRGO2.py not found at {script}. Set metagenomics.virgo2_dir to the VIRGO2 "
            "install (contains VIRGO2.py, Index/, AnnotationTables/)."
        )
    return script


def _concat_reads(r1, r2, dest) -> None:
    """Concatenate R1 then R2 into ``dest`` (byte copy; valid for gzip members + plain FASTQ)."""
    with open(dest, "wb") as out:
        for part in (r1, r2):
            with open(part, "rb") as fh:
                shutil.copyfileobj(fh, out)


def run_virgo2_map(
    r1, r2, sample, virgo2_dir, outdir, *,
    threads: int = 4, python: str = "python3", timeout: float | None = None,
):
    """Map one sample's non-host reads against VIRGO2 -> ``outdir/<sample>.out``.

    Returns ``(out_path, RunRecord)``.
    """
    py_exe, py_fp = resolve_executable(python, tool="python3 (VIRGO2)")
    script = _virgo2_script(virgo2_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    combined = outdir / f"{sample}.combined.fq.gz"
    _concat_reads(r1, r2, combined)

    argv = [py_exe, str(script), "map", "-r", str(combined),
            "-o", str(outdir / sample), "-p", str(threads)]
    record = run_command(
        argv, tool="VIRGO2.map", cwd=outdir,
        params={"sample": sample, "virgo2_dir": str(virgo2_dir)},
        exe_fingerprint=py_fp, timeout=timeout,
    )
    out_path = outdir / f"{sample}.out"
    if not out_path.exists():
        raise FileNotFoundError(
            f"VIRGO2 map finished (rc={record.returncode}) but {out_path} was not produced. "
            f"stderr tail:\n{record.stderr_tail}"
        )
    return out_path, record


def run_virgo2_compile(
    outdir, virgo2_dir, *, python: str = "python3", timeout: float | None = None,
):
    """Compile all ``*.out`` in ``outdir`` into one gene x sample matrix.

    Returns ``(compiled_path, RunRecord)`` where compiled_path is
    ``outdir/VIRGO2_Compiled.summary.NR.txt``.
    """
    py_exe, py_fp = resolve_executable(python, tool="python3 (VIRGO2)")
    script = _virgo2_script(virgo2_dir)
    outdir = Path(outdir)

    argv = [py_exe, str(script), "compile", "-i", str(outdir),
            "-o", str(outdir / "VIRGO2_Compiled")]
    record = run_command(
        argv, tool="VIRGO2.compile", cwd=outdir,
        params={"virgo2_dir": str(virgo2_dir)}, exe_fingerprint=py_fp, timeout=timeout,
    )
    compiled = outdir / COMPILED_NAME
    if not compiled.exists():
        raise FileNotFoundError(
            f"VIRGO2 compile finished (rc={record.returncode}) but {compiled} was not "
            f"produced. stderr tail:\n{record.stderr_tail}"
        )
    return compiled, record
