"""Stage graph: one registry, a resolver (multi-entry ladder), and two executors.

See ``registry.py`` for the stage definitions (the single source of truth).
"""

from microfgt.stages.check import CheckResult, check, speciateit_space_warnings
from microfgt.stages.executors import LocalExecutor, SnakemakeExecutor
from microfgt.stages.model import Requirement, Stage, StageContext
from microfgt.stages.registry import (
    STAGES,
    execute_stage,
    provided_artifacts,
)
from microfgt.stages.resolve import StageResolutionError, resolve

__all__ = [
    "STAGES",
    "Stage",
    "StageContext",
    "Requirement",
    "resolve",
    "StageResolutionError",
    "provided_artifacts",
    "execute_stage",
    "LocalExecutor",
    "SnakemakeExecutor",
    "check",
    "CheckResult",
    "speciateit_space_warnings",
]
