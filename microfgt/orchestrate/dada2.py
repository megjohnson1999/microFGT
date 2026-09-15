"""Orchestrate DADA2 — denoising, the second 16S front-end stage.

DADA2 is R/Bioconductor, so we orchestrate it (subprocess ``Rscript``) rather than reinvent
it — the same stance as speciateIT/VIRGO. microFGT owns one glue asset, ``dada2_run.R``,
which denoises primer-trimmed reads into an ASV table + rep-seqs and emits a quality profile.
Region-aware truncation/trimming defaults are chosen by the caller and passed as args; every
value is overridable.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path

from microfgt.orchestrate._run import resolve_executable, run_command


def _bundled_script() -> str:
    return str(resources.files("microfgt.scripts").joinpath("dada2_run.R"))


def run_dada2(
    input_dir,
    asv_table,
    asv_seqs,
    quality_profile,
    *,
    trunc_len=None,
    trim_left=None,
    rscript: str = "Rscript",
    script=None,
    timeout: float | None = None,
):
    """Run the bundled DADA2 script on a directory of primer-trimmed paired FASTQs.

    Parameters
    ----------
    input_dir:
        Directory of primer-trimmed ``*_R1*`` / ``*_R2*`` FASTQs.
    asv_table / asv_seqs / quality_profile:
        Output paths (samples x ASV CSV, rep-seq FASTA, per-position quality TSV).
    trunc_len / trim_left:
        ``(forward, reverse)`` lengths. The soundness knobs — defaulted region-aware by the
        caller, overridable here.

    Returns
    -------
    RunRecord
    """
    exe, fingerprint = resolve_executable(rscript, tool="R (Rscript)")
    script = script or _bundled_script()
    for out in (asv_table, asv_seqs, quality_profile):
        Path(out).parent.mkdir(parents=True, exist_ok=True)

    argv = [exe, str(script), "--input", str(input_dir), "--asv-table", str(asv_table),
            "--asv-seqs", str(asv_seqs), "--quality-profile", str(quality_profile)]
    if trunc_len:
        argv += ["--trunc-len", ",".join(map(str, trunc_len))]
    if trim_left:
        argv += ["--trim-left", ",".join(map(str, trim_left))]

    record = run_command(
        argv, tool="DADA2",
        params={"trunc_len": trunc_len, "trim_left": trim_left, "input": str(input_dir)},
        exe_fingerprint=fingerprint, timeout=timeout,
    )
    for out in (asv_table, asv_seqs):
        if not Path(out).exists():
            raise FileNotFoundError(
                f"DADA2 finished (rc={record.returncode}) but {out} was not produced. "
                f"stderr tail:\n{record.stderr_tail}"
            )
    return record
