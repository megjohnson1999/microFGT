"""Orchestrate cutadapt — primer trimming, the first 16S front-end stage.

Removes amplicon primers from paired FASTQs before DADA2. Like the other wrappers, cutadapt
is located (PATH / explicit path), run with recorded provenance, and not bundled.
"""

from __future__ import annotations

import re
from pathlib import Path

from microfgt.orchestrate._run import resolve_executable, run_command

# ENA/SRA mate marker: `_1`/`_2` sitting immediately before the .fastq/.fq[.gz] extension.
# Anchored at the end so a sample name that itself contains `_1` is not mangled.
_ENA_R1 = re.compile(r"_1\.(fastq|fq)(\.gz)?$")


def discover_pairs(fastq_dir):
    """Find (sample, R1, R2) FASTQ pairs by the ``_R1/_R2`` or ``_1/_2`` convention.

    Prefers Illumina ``_R1/_R2`` (also matching ``_R1_001`` etc.); falls back to ENA/SRA
    ``_1/_2`` only when no ``_R1`` files are present, so mixed dirs never double-count.
    """
    fastq_dir = Path(fastq_dir)

    # 1) Illumina _R1/_R2 — unchanged behavior (allows a suffix like _R1_001).
    pairs = []
    for r1 in sorted(fastq_dir.glob("*_R1*.fastq*")):
        r2 = Path(str(r1).replace("_R1", "_R2"))
        sample = r1.name.split("_R1")[0]
        pairs.append((sample, r1, r2 if r2.exists() else None))
    if pairs:
        return pairs

    # 2) Fallback: ENA/SRA _1/_2 (the marker only just before the extension).
    for r1 in sorted(p for p in fastq_dir.glob("*_1.f*") if _ENA_R1.search(p.name)):
        r2 = r1.with_name(_ENA_R1.sub(lambda m: "_2." + m.group(1) + (m.group(2) or ""), r1.name))
        sample = _ENA_R1.sub("", r1.name)
        pairs.append((sample, r1, r2 if r2.exists() else None))
    return pairs


def run_cutadapt(
    input_dir,
    output_dir,
    *,
    fwd_primer: str | None = None,
    rev_primer: str | None = None,
    executable: str = "cutadapt",
    discard_untrimmed: bool = False,
    extra_args=None,
    timeout: float | None = None,
):
    """Trim primers from every FASTQ pair in ``input_dir`` into ``output_dir``.

    Parameters
    ----------
    fwd_primer / rev_primer:
        5' primers passed as cutadapt ``-g`` / ``-G`` (regular 5' adapters — the match plus any
        preceding heterogeneity spacer is removed).
    discard_untrimmed:
        Drop read pairs where the primer was not found (cutadapt ``--discard-untrimmed``). For
        amplicon data this keeps only on-target reads; only meaningful when a primer is given
        (with no adapter, cutadapt would discard everything). See :func:`_run_primer_trim` for the
        primers-configured default.
    extra_args:
        Extra cutadapt args appended verbatim (e.g. ``--pair-filter=both``, ``--minimum-length``),
        so a user can tune behaviour without a code change. Appended last, so they can override
        the flags set above.

    Returns one RunRecord per sample.
    """
    exe, fingerprint = resolve_executable(executable, tool="cutadapt")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs = discover_pairs(input_dir)
    if not pairs:
        raise FileNotFoundError(f"No '*_R1*.fastq*' files found in {input_dir}.")

    records = []
    for sample, r1, r2 in pairs:
        argv = [exe]
        if fwd_primer:
            argv += ["-g", fwd_primer]
        if r2 and rev_primer:
            argv += ["-G", rev_primer]
        if discard_untrimmed:
            argv += ["--discard-untrimmed"]
        argv += ["-o", str(output_dir / r1.name)]
        if r2:
            argv += ["-p", str(output_dir / r2.name)]
        argv += list(extra_args or [])
        argv += [str(r1)] + ([str(r2)] if r2 else [])
        records.append(
            run_command(
                argv, tool="cutadapt", params={"sample": sample},
                exe_fingerprint=fingerprint, timeout=timeout,
            )
        )
    return records
