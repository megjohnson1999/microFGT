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
    # Absolute: VIRGO2.py runs with cwd=outdir (run_command below), so any relative path in
    # argv would resolve *inside* outdir and double up (breaks on a relative --workdir).
    outdir = Path(outdir).resolve()
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
    outdir, virgo2_dir, *, compiled_dir=None, python: str = "python3",
    timeout: float | None = None,
):
    """Compile all ``*.out`` in ``outdir`` into one gene x sample matrix.

    The compiled matrix is written to ``compiled_dir`` (default: ``outdir``). The pipeline
    passes a separate dir so the compiled file is NOT nested inside ``outdir`` — a directory
    that is itself a Snakemake rule output (nesting is a ChildIOException).

    Returns ``(compiled_path, RunRecord)`` where compiled_path is
    ``<compiled_dir>/VIRGO2_Compiled.summary.NR.txt``.
    """
    py_exe, py_fp = resolve_executable(python, tool="python3 (VIRGO2)")
    script = _virgo2_script(virgo2_dir)
    outdir = Path(outdir).resolve()  # cwd=outdir below; keep argv paths absolute (see run_virgo2_map)
    dest = Path(compiled_dir).resolve() if compiled_dir else outdir
    dest.mkdir(parents=True, exist_ok=True)

    argv = [py_exe, str(script), "compile", "-i", str(outdir),
            "-o", str(dest / "VIRGO2_Compiled")]
    record = run_command(
        argv, tool="VIRGO2.compile", cwd=outdir,
        params={"virgo2_dir": str(virgo2_dir)}, exe_fingerprint=py_fp, timeout=timeout,
    )
    compiled = dest / COMPILED_NAME
    if not compiled.exists():
        raise FileNotFoundError(
            f"VIRGO2 compile finished (rc={record.returncode}) but {compiled} was not "
            f"produced. stderr tail:\n{record.stderr_tail}"
        )
    return compiled, record
