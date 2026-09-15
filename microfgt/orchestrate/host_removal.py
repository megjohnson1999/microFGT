"""Orchestrate host removal — minimap2 vs the host reference, keep non-host pairs.

Grounded in ``prototype/reference_scripts/RECIPE.md`` (sg_host_removal). The logic *is* the
samtools filter: map with minimap2, keep only read pairs where **both mates are unmapped**
(``-f 12``) dropping secondary alignments (``-F 256``), name-sort, then re-emit FASTQ. FGT
samples are very host-heavy (FRESH: ~86% removed), so most depth is spent here.

This is a genuine shell pipeline, so it runs via ``bash -c`` (like the v1 VIRGO wrapper's
``bash`` shell-out) rather than a single argv. Every substituted path is ``shlex.quote``-d —
the pipeline is a shell string, so paths with spaces must be quoted.
"""

from __future__ import annotations

import shlex
from pathlib import Path

from microfgt.orchestrate._run import resolve_executable, run_command
from microfgt.orchestrate.cutadapt import discover_pairs

# The RECIPE.md pipeline, with tool/path/thread placeholders filled by run_host_removal.
_PIPELINE = (
    "{minimap2} -ax sr -t {threads} {host_ref} {r1} {r2} "
    "| {samtools} view -@ {threads} -b -f 12 -F 256 - "
    "| {samtools} sort -n -@ {threads} -m 2G -T {sorttmp} - "
    "| {samtools} fastq -@ {threads} -1 {out1} -2 {out2} -0 /dev/null -s /dev/null -"
)


def run_host_removal(
    input_dir,
    output_dir,
    host_ref,
    *,
    threads: int = 4,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    runner: str = "bash",
    timeout: float | None = None,
):
    """Remove host reads from every trimmed FASTQ pair in ``input_dir`` into ``output_dir``.

    Keeps only both-mates-unmapped pairs (``-f 12``). Non-host reads keep their original
    filenames. Returns one RunRecord per sample.
    """
    mm2_exe, mm2_fp = resolve_executable(minimap2, tool="minimap2")
    st_exe, _ = resolve_executable(samtools, tool="samtools")
    runner_exe, _ = resolve_executable(runner, tool="host-removal (shell)")
    host_ref = Path(host_ref)
    if not host_ref.exists():
        raise FileNotFoundError(
            f"Host reference {host_ref} not found. Set metagenomics.host_ref to the host "
            "genome (e.g. GRCh38.fna.gz)."
        )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs = discover_pairs(input_dir)
    if not pairs:
        raise FileNotFoundError(f"No '*_R1*.fastq*' files found in {input_dir}.")

    records = []
    for sample, r1, r2 in pairs:
        if r2 is None:
            raise FileNotFoundError(
                f"host removal expects paired reads; sample {sample!r} has no R2 mate."
            )
        pipeline = _PIPELINE.format(
            minimap2=shlex.quote(mm2_exe), samtools=shlex.quote(st_exe),
            threads=int(threads), host_ref=shlex.quote(str(host_ref)),
            r1=shlex.quote(str(r1)), r2=shlex.quote(str(r2)),
            sorttmp=shlex.quote(str(output_dir / f"{sample}.sorttmp")),
            out1=shlex.quote(str(output_dir / r1.name)),
            out2=shlex.quote(str(output_dir / r2.name)),
        )
        records.append(
            run_command(
                [runner_exe, "-c", pipeline], tool="host_removal",
                params={"sample": sample, "host_ref": str(host_ref)},
                exe_fingerprint=mm2_fp, timeout=timeout,
            )
        )
    return records
