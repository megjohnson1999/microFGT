"""The resolver — make-semantics over the stage registry.

Given the target artifact and the artifacts already present (the entry point), return the
ordered list of stages to run, skipping any whose output is already supplied. This one
mechanism IS the multi-entry-point ladder: FASTQs, ASV table, existing outputs, or a
composition all "just work" by entering wherever the provided inputs land.
"""

from __future__ import annotations

from microfgt.stages.model import Stage
from microfgt.stages.registry import PRODUCERS, STAGES


class StageResolutionError(Exception):
    """No registered stage can produce a needed artifact from the provided inputs."""


def _resolvable(artifact: str, provided: set[str], stack: frozenset = frozenset()) -> bool:
    """Can ``artifact`` be produced from ``provided`` (directly, or via some producer whose
    inputs are all themselves resolvable)? ``stack`` guards against cycles."""
    if artifact in provided:
        return True
    if artifact in stack:
        return False
    stack = stack | {artifact}
    return any(
        all(_resolvable(i, provided, stack) for i in stage.inputs)
        for stage in PRODUCERS.get(artifact, [])
    )


def _choose_producer(artifact: str, provided: set[str]) -> Stage | None:
    candidates = PRODUCERS.get(artifact, [])
    if not candidates:
        return None
    # 1) Prefer a producer whose input is directly provided (e.g. cst_valencia when a VALENCIA
    # output is supplied, or import_mgcst_existing when a VISTA output is).
    for stage in candidates:
        if any(i in provided for i in stage.inputs):
            return stage
    # 2) Otherwise the *most complete* producer whose every input is resolvable from `provided`.
    # This routes the three `mudata` producers: 16S-only -> integrate, shotgun-only ->
    # integrate_shotgun, both -> integrate_combined (most inputs wins).
    resolvable = [s for s in candidates if all(_resolvable(i, provided) for i in s.inputs)]
    if resolvable:
        return max(resolvable, key=lambda s: len(s.inputs))
    return candidates[0]


def resolve(target: str, provided) -> list[Stage]:
    """Stages needed to produce ``target`` given ``provided`` artifact keys, in run order."""
    provided = set(provided)
    order: list[Stage] = []
    placed: set[str] = set()

    def need(artifact: str) -> None:
        if artifact in provided:
            return
        stage = _choose_producer(artifact, provided)
        if stage is None:
            raise StageResolutionError(
                f"Cannot produce artifact {artifact!r} from the provided inputs "
                f"{sorted(provided)}. Supply it directly or an upstream input."
            )
        if stage.id in placed:
            return
        for inp in stage.inputs:
            need(inp)
        if stage.id not in placed:           # re-check: deps may have placed nothing new
            placed.add(stage.id)
            order.append(stage)

    need(target)
    return order


def all_stage_ids() -> list[str]:
    return [s.id for s in STAGES]
