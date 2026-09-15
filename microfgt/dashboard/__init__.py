"""Layer 3 — the dashboard: a thin Streamlit presenter over the analysis verbs.

The importable pieces here are the *pure* dashboard logic (variable introspection); the
Streamlit UI lives in ``app.py`` and is launched via ``microfgt dashboard`` (needs the
``[app]`` extra), never imported as part of the package.
"""

from microfgt.dashboard.catalog import (
    Variable,
    as_frame,
    continuous,
    groupable,
    variable_catalog,
)

__all__ = ["Variable", "variable_catalog", "groupable", "continuous", "as_frame"]
