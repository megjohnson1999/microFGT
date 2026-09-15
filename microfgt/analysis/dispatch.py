"""One dispatcher from a (verb name + selections) spec to the matching analysis verb.

This is the seam the CLI and the dashboard both call, so "run an analysis from a spec of
choices" lives in exactly one place — the dashboard is provably the same path a power user
takes, never a parallel implementation.
"""

from __future__ import annotations

from microfgt.analysis.abundance import compare_abundance
from microfgt.analysis.association import associate
from microfgt.analysis.hypothesis import compare_alpha, compare_beta
from microfgt.analysis.results import AnalysisResult

VERBS = ("alpha", "beta", "associate", "abundance")


def run_verb(
    data,
    verb: str,
    *,
    predictors=None,
    x: str | None = None,
    y: str | None = None,
    subject: str | None = None,
    metric: str | None = None,
    method: str | None = None,
    subset=None,
) -> AnalysisResult:
    """Run one analysis ``verb`` from a normalized selection spec.

    ``predictors`` (str or list) drives ``alpha`` / ``beta`` / ``abundance`` — first is the
    predictor of interest, the rest covariates. ``x`` / ``y`` drive ``associate``. ``metric``
    and ``method`` fall back to each verb's own default when None.
    """
    if verb == "alpha":
        _require(predictors, verb, "predictors")
        return compare_alpha(data, predictors, metric=metric or "shannon",
                             subject=subject, subset=subset)
    if verb == "beta":
        _require(predictors, verb, "predictors")
        return compare_beta(data, predictors, metric=metric or "braycurtis", subset=subset)
    if verb == "abundance":
        _require(predictors, verb, "predictors")
        return compare_abundance(data, predictors, method=method or "ancombc",
                                subject=subject, subset=subset)
    if verb == "associate":
        if not (x and y):
            raise ValueError("verb 'associate' needs both x and y.")
        return associate(data, x, y, method=method or "auto", subset=subset)
    raise ValueError(f"Unknown verb {verb!r}; available: {VERBS}.")


def _require(value, verb, name):
    if not value:
        raise ValueError(f"verb {verb!r} needs {name}.")
