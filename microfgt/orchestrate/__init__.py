"""Orchestration layer (optional) — actually run the external tools.

The user can ingest existing tool outputs (the importers) OR have microFGT run the tools.
microFGT does NOT bundle the tools or their large reference data; it locates an installed
tool (configurable path / PATH), runs it with recorded provenance, and hands the output to
the matching importer.
"""

from microfgt.orchestrate._run import (
    RunRecord,
    ToolNotFoundError,
    resolve_executable,
    run_command,
)
from microfgt.orchestrate.fastp import run_fastp
from microfgt.orchestrate.host_removal import run_host_removal
from microfgt.orchestrate.speciateit import run_speciateit
from microfgt.orchestrate.virgo import run_virgo, run_virgo_samples
from microfgt.orchestrate.virgo2 import run_virgo2_compile, run_virgo2_map
from microfgt.orchestrate.vista import run_vista

__all__ = [
    "RunRecord",
    "ToolNotFoundError",
    "resolve_executable",
    "run_command",
    "run_fastp",
    "run_host_removal",
    "run_speciateit",
    "run_virgo",
    "run_virgo_samples",
    "run_virgo2_compile",
    "run_virgo2_map",
    "run_vista",
]
