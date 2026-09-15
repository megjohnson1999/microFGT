"""Orchestrate fastp — QC (adapter + quality trim), the first metagenomics stage.

Grounded in ``prototype/reference_scripts/RECIPE.md`` (sg_qc): default paired-end fastp with
auto adapter detection + quality trim, one invocation per sample pair, keeping the per-sample
``qc.json`` (read/quality metrics for the object). Like the other wrappers, fastp is located
(PATH / explicit path), run with recorded provenance, and not bundled.
"""

from __future__ import annotations

from pathlib import Path

from microfgt.orchestrate._run import resolve_executable, run_command
from microfgt.orchestrate.cutadapt import discover_pairs


def run_fastp(
    input_dir,
    output_dir,
    *,
    threads: int = 4,
    executable: str = "fastp",
    extra_args=None,
    timeout: float | None = None,
):
    """Trim every FASTQ pair in ``input_dir`` into ``output_dir`` with fastp.

    Trimmed reads keep their original filenames; each sample also gets a ``<sample>.qc.json``
    (and ``.qc.html``) report. Returns one RunRecord per sample.
    """
    exe, fingerprint = resolve_executable(executable, tool="fastp")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs = discover_pairs(input_dir)
    if not pairs:
        raise FileNotFoundError(f"No '*_R1*.fastq*' files found in {input_dir}.")

    records = []
    for sample, r1, r2 in pairs:
        if r2 is None:
            raise FileNotFoundError(
                f"fastp stage expects paired reads; sample {sample!r} has no R2 mate."
            )
        argv = [
            exe,
            "--in1", str(r1), "--in2", str(r2),
            "--out1", str(output_dir / r1.name), "--out2", str(output_dir / r2.name),
            "--thread", str(threads),
            "--json", str(output_dir / f"{sample}.qc.json"),
            "--html", str(output_dir / f"{sample}.qc.html"),
        ]
        argv += list(extra_args or [])
        records.append(
            run_command(
                argv, tool="fastp", params={"sample": sample},
                exe_fingerprint=fingerprint, timeout=timeout,
            )
        )
    return records
