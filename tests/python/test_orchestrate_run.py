"""run_command surfaces an external tool's own error output on failure (no hardcoded messages)."""

import subprocess

import pytest

from microfgt.orchestrate._run import ToolRunError, run_command


def test_run_command_returns_record_on_success():
    rec = run_command(["sh", "-c", "echo hi; exit 0"], tool="FakeTool")
    assert rec.returncode == 0
    assert rec.tool == "FakeTool"


def test_run_command_surfaces_tool_stderr_on_failure():
    # The tool writes an arbitrary message to stderr and exits non-zero.
    argv = ["sh", "-c", "echo 'the tool broke: bad input file' >&2; exit 3"]
    with pytest.raises(ToolRunError) as excinfo:
        run_command(argv, tool="FakeTool")

    err = excinfo.value
    assert err.returncode == 3
    assert err.tool == "FakeTool"
    # Backward-compatible: existing handlers catching CalledProcessError still work.
    assert isinstance(err, subprocess.CalledProcessError)

    msg = str(err)
    assert "FakeTool" in msg                         # names the tool + exit code
    assert "exit 3" in msg
    assert "the tool broke: bad input file" in msg   # relays the tool's OWN words, verbatim


def test_run_command_no_check_returns_record_even_on_failure():
    rec = run_command(["sh", "-c", "exit 5"], tool="FakeTool", check=False)
    assert rec.returncode == 5
