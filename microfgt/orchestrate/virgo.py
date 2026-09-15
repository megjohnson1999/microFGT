"""Orchestrate VIRGO — run the read mapping, then own the stack-into-matrix.

Command (grounded in VIRGO's ``3_run_VIRGO/runMapping.step1.sh``, not guessed):
    bash runMapping.step1.sh -r <reads.fastq> -p <prefix> -d <VIRGO-path>
which writes ``<workdir>/temp_mapping/<prefix>.out`` — the 3-column (geneID, read_count,
gene_length) per-sample file that :func:`microfgt.io.import_virgo` consumes. VIRGO does NOT
support paired-end mapping; merge reads into one fastq per sample first.

We run step1 per sample into a shared mapping directory, then stack the ``.out`` files into
the gene x sample ``function`` modality with the existing importer.
"""

from __future__ import annotations

from pathlib import Path

from microfgt.io import import_virgo
from microfgt.orchestrate._run import resolve_executable, run_command

STEP1_SCRIPT = "3_run_VIRGO/runMapping.step1.sh"


def run_virgo(
    reads,
    prefix: str,
    virgo_path,
    *,
    workdir=None,
    script=None,
    runner: str = "bash",
    timeout: float | None = None,
):
    """Run VIRGO read-mapping (step1) for one sample.

    Parameters
    ----------
    reads:
        Single merged FASTQ for the sample (no paired-end).
    prefix:
        Sample name; the output is ``temp_mapping/<prefix>.out``.
    virgo_path:
        Path to the VIRGO install (containing ``0_db/`` and ``1_VIRGO/``).
    workdir:
        Directory to run in (``temp_mapping/`` is created under it). Defaults to cwd.
    script:
        Path to ``runMapping.step1.sh`` (defaults to ``<virgo_path>/3_run_VIRGO/...``).
    runner:
        Shell used to run the script (default ``bash``).

    Returns
    -------
    tuple[pathlib.Path, RunRecord]
        Path to ``<prefix>.out`` and the run provenance.
    """
    virgo_path = Path(virgo_path)
    workdir = Path(workdir) if workdir is not None else Path.cwd()
    workdir.mkdir(parents=True, exist_ok=True)
    script_path = Path(script) if script is not None else virgo_path / STEP1_SCRIPT

    runner_exe, runner_fp = resolve_executable(runner, tool="VIRGO (shell)")
    if not script_path.exists():
        raise FileNotFoundError(
            f"VIRGO mapping script not found at {script_path}. Pass script=... or check "
            f"the VIRGO install at {virgo_path}."
        )

    argv = [
        runner_exe, str(script_path),
        "-r", str(reads), "-p", prefix, "-d", str(virgo_path),
    ]
    record = run_command(
        argv,
        tool="VIRGO",
        cwd=workdir,
        params={"reads": str(reads), "prefix": prefix, "virgo_path": str(virgo_path)},
        exe_fingerprint=runner_fp,
        timeout=timeout,
    )

    out_path = workdir / "temp_mapping" / f"{prefix}.out"
    if not out_path.exists():
        raise FileNotFoundError(
            f"VIRGO finished (rc={record.returncode}) but {out_path} was not produced. "
            f"stderr tail:\n{record.stderr_tail}"
        )
    return out_path, record


def run_virgo_samples(reads_by_sample: dict, virgo_path, *, workdir=None, **kwargs):
    """Map several samples, then stack their ``.out`` files into one AnnData.

    Parameters
    ----------
    reads_by_sample:
        ``{sample_prefix: reads_fastq}``.
    virgo_path / workdir / kwargs:
        As :func:`run_virgo`.

    Returns
    -------
    anndata.AnnData
        The ``function`` modality, with per-sample run provenance in ``uns['virgo_runs']``.
    """
    workdir = Path(workdir) if workdir is not None else Path.cwd()
    records = {}
    for prefix, reads in reads_by_sample.items():
        _, record = run_virgo(reads, prefix, virgo_path, workdir=workdir, **kwargs)
        records[prefix] = record.to_dict()

    adata = import_virgo(workdir / "temp_mapping")
    adata.uns["virgo_runs"] = records
    return adata
