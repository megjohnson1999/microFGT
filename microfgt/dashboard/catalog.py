"""Introspect an object's sample variables — what the dashboard's dropdowns offer.

The dashboard is a spec-builder: it reads the object's obs columns, sorts them into
categorical (groupings) vs continuous, and hands them to widgets. That introspection is the
one genuinely-dashboard piece of logic, so it lives here as a *pure, testable* function —
the Streamlit file stays a thin shell of widgets over this + ``run_verb`` + ``render``.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from microfgt.analysis._frame import is_categorical, merged_obs


@dataclass
class Variable:
    """A sample variable the user can pick, with what a widget needs to present it."""

    name: str
    kind: str                       # "categorical" | "continuous"
    n_present: int
    n_missing: int
    n_levels: int | None = None     # categorical only
    levels: list | None = None      # categorical only (sorted, capped)


def variable_catalog(data, *, max_levels: int = 100) -> list[Variable]:
    """Every sample variable (from global obs + all modalities), typed for the widgets."""
    obs = merged_obs(data)
    out: list[Variable] = []
    for col in obs.columns:
        s = obs[col]
        present = int(s.notna().sum())
        missing = int(s.isna().sum())
        if is_categorical(s):
            levels = sorted(s.dropna().astype(str).unique())
            out.append(Variable(col, "categorical", present, missing,
                                n_levels=len(levels), levels=levels[:max_levels]))
        else:
            out.append(Variable(col, "continuous", present, missing))
    return out


def groupable(catalog: list[Variable], *, min_levels: int = 2, max_levels: int = 50) -> list[str]:
    """Categorical variables usable as a grouping/predictor (2..N levels)."""
    return [v.name for v in catalog
            if v.kind == "categorical" and min_levels <= (v.n_levels or 0) <= max_levels]


def continuous(catalog: list[Variable]) -> list[str]:
    return [v.name for v in catalog if v.kind == "continuous"]


def as_frame(catalog: list[Variable]) -> pd.DataFrame:
    """A tidy overview table of the variables (shown on the dashboard landing view)."""
    return pd.DataFrame(
        [{"variable": v.name, "kind": v.kind, "levels": v.n_levels,
          "present": v.n_present, "missing": v.n_missing} for v in catalog]
    ).set_index("variable")
