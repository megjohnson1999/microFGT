"""Shared subprocess plumbing for the orchestration layer.

The orchestration layer is optional on purpose: a user can ingest existing tool outputs
(the importers) OR have microFGT actually run the tools. Either way the tool's reference
data (speciateIT's ~2.6 GB models, VIRGO's hosted catalog) is far too large to bundle, so
we assume the tool is installed at a configurable path / on PATH and own the rest: build
the command, run it, record provenance (constraint B — reproducibility), surface helpful
errors, and hand the output straight to the matching importer.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path


class ToolNotFoundError(FileNotFoundError):
    """Raised when an external tool executable cannot be located."""


class ToolRunError(subprocess.CalledProcessError):
    """An external tool exited non-zero.

    Subclasses ``CalledProcessError`` (so existing handlers still catch it) but its message
    shows the tool's OWN error output instead of just the exit code — the plain
    ``CalledProcessError`` captures stderr yet never prints it, so the real cause is invisible.
    We do not interpret the message; we relay whatever the tool wrote to stderr, verbatim.
    """

    def __init__(self, tool: str, returncode: int, argv, stdout: str, stderr: str):
        super().__init__(returncode, argv, output=stdout, stderr=stderr)
        self.tool = tool

    def __str__(self) -> str:
        err = (self.stderr or "").strip()
        tail = "\n".join(err.splitlines()[-15:]) if err else "(the tool produced no error output)"
        return f"{self.tool} failed (exit {self.returncode}).\n--- {self.tool} error output ---\n{tail}"


@dataclass
class RunRecord:
    """Provenance for one external-tool invocation (serialize into ``.uns``)."""

    tool: str
    executable: str                     # resolved absolute path
    argv: list[str]
    returncode: int
    cwd: str
    duration_s: float
    started_at: str                     # ISO 8601 UTC
    params: dict = field(default_factory=dict)
    exe_fingerprint: dict = field(default_factory=dict)  # size/mtime of the binary
    stdout_tail: str = ""
    stderr_tail: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def resolve_executable(name_or_path: str, *, tool: str) -> tuple[str, dict]:
    """Resolve an executable to an absolute path, or raise ToolNotFoundError.

    Accepts either a bare command name (looked up on PATH) or an explicit path to the
    binary/script. Returns ``(abs_path, fingerprint)`` where fingerprint records the
    binary's size/mtime for reproducibility (the tools have no reliable ``--version``)."""
    candidate = Path(name_or_path)
    if candidate.exists() or os.sep in name_or_path or (os.altsep and os.altsep in name_or_path):
        resolved = str(candidate.resolve()) if candidate.exists() else None
    else:
        found = shutil.which(name_or_path)
        resolved = str(Path(found).resolve()) if found else None

    if resolved is None or not Path(resolved).exists():
        raise ToolNotFoundError(
            f"Could not find the {tool} executable {name_or_path!r}. Install {tool} and "
            f"either put it on PATH or pass its path explicitly. microFGT does not bundle "
            f"{tool} or its reference data (too large)."
        )
    st = Path(resolved).stat()
    return resolved, {"size": st.st_size, "mtime": st.st_mtime}


def run_command(
    argv: list[str],
    *,
    tool: str,
    cwd: str | os.PathLike | None = None,
    params: dict | None = None,
    exe_fingerprint: dict | None = None,
    timeout: float | None = None,
    check: bool = True,
) -> RunRecord:
    """Run ``argv``, capturing output and provenance into a :class:`RunRecord`.

    Raises :class:`ToolRunError` (a ``CalledProcessError`` subclass whose message shows the
    tool's own stderr) if ``check`` and the command exits non-zero."""
    started = datetime.now(timezone.utc).isoformat()
    t0 = time.monotonic()
    proc = subprocess.run(
        argv, cwd=cwd, capture_output=True, text=True, timeout=timeout,
    )
    duration = time.monotonic() - t0

    record = RunRecord(
        tool=tool,
        executable=argv[0],
        argv=list(argv),
        returncode=proc.returncode,
        cwd=str(cwd or os.getcwd()),
        duration_s=duration,
        started_at=started,
        params=params or {},
        exe_fingerprint=exe_fingerprint or {},
        stdout_tail=proc.stdout[-2000:],
        stderr_tail=proc.stderr[-2000:],
    )
    if check and proc.returncode != 0:
        raise ToolRunError(tool, proc.returncode, argv, proc.stdout, proc.stderr)
    return record
