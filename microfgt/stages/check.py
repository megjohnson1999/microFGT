"""`microfgt check` — preflight dependency doctor driven by the resolved entry point.

Because the resolver already knows which stages will run for the given inputs, check verifies
exactly those stages' requirements (binaries, R packages, DB paths) plus region<->speciateIT
DB consistency — and fails early, helpfully, instead of three stages deep. Turns the
"server-side tools" caveat into a constraint-A feature.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from microfgt.stages.registry import provided_artifacts
from microfgt.stages.resolve import resolve


@dataclass
class CheckResult:
    ok: bool
    message: str


def _r_has_package(pkg: str, rscript: str = "Rscript") -> bool:
    if shutil.which(rscript) is None:
        return False
    try:
        proc = subprocess.run(
            [rscript, "-e", f"quit(status = !requireNamespace('{pkg}', quietly = TRUE))"],
            capture_output=True, timeout=60,
        )
        return proc.returncode == 0
    except Exception:
        return False


def _sha256(path: str) -> str | None:
    import hashlib

    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def _verify(req) -> CheckResult:
    if req.kind == "binary":
        ok = shutil.which(req.name) is not None
    elif req.kind == "path":
        ok = Path(req.name).exists()
    elif req.kind == "rpackage":
        # Verify with the CONFIGURED Rscript, not a bare 'Rscript' — otherwise the doctor
        # misreports the most fragile layer (R/VISTA) for anyone whose R is not on PATH.
        ok = _r_has_package(req.name, req.via or "Rscript")
    elif req.kind == "checksum":
        ok = _sha256(req.name) == req.expected
    else:  # pragma: no cover
        ok = False
    label = f"{req.kind} {req.name!r}"
    return CheckResult(ok, f"OK   {label}" if ok else f"MISS {label} — {req.hint}")


def _region_db_consistency(config, stage_ids) -> list[CheckResult]:
    if "assign" not in stage_ids:
        return []
    reads = (config.get("composition") or {}).get("reads") or {}
    sit = (config.get("composition") or {}).get("speciateit") or {}
    region, db = reads.get("region"), sit.get("db")
    if not region or not db:
        return []
    token = region.replace("-", "").upper()
    ok = token in Path(db).name.replace("-", "").upper()
    msg = (
        f"OK   region {region!r} matches speciateIT db {Path(db).name!r}" if ok
        else f"MISS region {region!r} does not match speciateIT db {Path(db).name!r} "
             f"(check you pointed at the right model)"
    )
    return [CheckResult(ok, msg)]


def speciateit_space_warnings(config: dict, workdir=None, target: str = "mudata") -> list[str]:
    """Advisory: paths feeding speciateIT that contain a space.

    speciateIT's ``classify`` shells out to ``mkdir``/``grep`` unquoted, so a space in a path it
    reads or writes splits the command and the run fails partway (e.g. ``mkdir: /Users/Megan:
    Permission denied``). The offenders are its **input FASTA** and its **output dir / workdir** —
    the ``db`` path is read via ``fopen`` and tolerates spaces. Returns human-readable warning
    lines (empty when the speciateIT ``assign`` stage isn't in the plan, or nothing has a space).
    This is a warning, not a hard failure: a non-bundled ``classify`` may tolerate spaces.
    """
    try:
        stage_ids = {s.id for s in resolve(target, set(provided_artifacts(config)))}
    except Exception:
        return []
    if "assign" not in stage_ids:
        return []
    comp = config.get("composition") or {}
    reads = comp.get("reads") or {}
    candidates = {
        "workdir": workdir,
        "composition.asv_seqs": comp.get("asv_seqs"),
        "composition.reads.fastq_dir": reads.get("fastq_dir"),
    }
    return [
        f"WARN speciateIT's classify breaks on spaces in a path — {name} contains a space "
        f"({value!r}). Use a space-free location or the run may fail mid-pipeline."
        for name, value in candidates.items()
        if value and " " in str(value)
    ]


def check(config: dict, target: str = "mudata") -> list[CheckResult]:
    """Verify the tools/paths needed for the stages that WILL run for these inputs."""
    provided = set(provided_artifacts(config))
    stages = resolve(target, provided)
    stage_ids = {s.id for s in stages}

    results: list[CheckResult] = []
    seen = set()
    for stage in stages:
        for req in stage.requirements(config):
            key = (req.kind, req.name)
            if key in seen:
                continue
            seen.add(key)
            results.append(_verify(req))
    results += _region_db_consistency(config, stage_ids)
    return results
