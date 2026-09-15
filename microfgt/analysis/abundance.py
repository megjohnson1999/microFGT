"""``compare_abundance`` — which taxa differ, adjusting for covariates? (FDR-corrected)

The differential-abundance verb, completing the family. Unlike the single-factor ANCOM
*primitive* in :mod:`microfgt.analysis.diffabund`, this is the test a researcher actually
runs: a **covariate-adjusted, FDR-corrected** model of each taxon's abundance against a
predictor of interest plus nuisance covariates, returning the uniform
:class:`~microfgt.analysis.results.AnalysisResult`.

Engine (the hybrid decision, resolved by what's real): scikit-bio 0.7 ships **ANCOM-BC**
(Lin & Peddada 2020 — one of the brand-name adjusted methods) and a **Dirichlet-multinomial
linear mixed model** for repeated measures, both pure-Python and formula-based. So the
covariate-adjusted core needs no R, and validates on real data today. Orchestrating R for
MaAsLin2 / DESeq2 stays a genuine extension point behind the same ``method=`` switch — added
when a project needs those exact tools (neither is installed here, so shipping that glue now
would be untested dead code).

* ``method="ancombc"`` (default) — cross-sectional, covariate-adjusted (ANCOM-BC).
* ``method="dirmult_lme"`` — longitudinal: a per-subject random intercept (needs ``subject``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from microfgt.analysis._frame import analysis_frame, get_assay
from microfgt.analysis.results import AnalysisResult

_CORE = ["contrast", "log2fc", "pvalue", "qvalue", "signif"]


def compare_abundance(
    data,
    predictors,
    *,
    method: str = "ancombc",
    modality=None,
    subject: str | None = None,
    subset=None,
    alpha: float = 0.05,
    p_adjust: str = "fdr_bh",
    pseudocount: float = 1.0,
    layer: str = "counts",
    **kwargs,
) -> AnalysisResult:
    """Covariate-adjusted differential abundance of each taxon against ``predictors``.

    Parameters
    ----------
    data:
        MuData or AnnData; predictors resolve from global obs OR any modality's obs.
    predictors:
        ``obs`` column name(s). The **first** is the predictor of interest (reported); any
        others are covariates the model adjusts for. A single string is fine.
    method:
        ``"ancombc"`` (default, cross-sectional) or ``"dirmult_lme"`` (repeated measures, needs
        ``subject``).
    modality:
        Which assay's features to test (default: the taxon roll-up, else composition).
    subject:
        ``obs`` subject-id column — the random-intercept group for ``dirmult_lme``.
    subset:
        Restrict samples first: a query string or ``{column: value(s)}``.
    alpha, p_adjust:
        Significance level and multiple-testing correction (default Benjamini–Hochberg).

    Returns
    -------
    AnalysisResult
        ``table`` = per-feature results for the predictor of interest (``log2fc`` / ``pvalue``
        / ``qvalue`` / ``signif`` / ``contrast``), sorted by ``qvalue``; ``stats`` carries the
        method, formula, and count of significant features.
    """
    preds = [predictors] if isinstance(predictors, str) else list(predictors)
    primary = preds[0]
    adata, mod = get_assay(data, modality)

    cols = preds + ([subject] if subject else [])
    frame, notes = analysis_frame(data, cols, modality=mod, subset=subset)

    counts = pd.DataFrame(
        np.asarray(adata.layers[layer] if layer in adata.layers else adata.X, dtype=float),
        index=adata.obs_names.astype(str), columns=adata.var_names.astype(str),
    ).reindex(frame.index)
    keep_feat = counts.columns[counts.sum(axis=0) > 0]          # all-zero features are undefined
    counts = counts[keep_feat]
    n_dropped_feat = int(adata.n_vars - counts.shape[1])

    formula = " + ".join(_ref(p) for p in preds)
    tidy = _run(method, counts, frame, formula, subject, alpha, p_adjust, pseudocount, kwargs)
    table = _primary_rows(tidy, primary).sort_values("qvalue")

    n_sig = int(table["signif"].sum())
    stats = {"test": "ANCOM-BC" if method == "ancombc" else "Dirichlet-multinomial LME",
             "formula": formula, "primary": primary, "n_features": int(table.shape[0]),
             "n_significant": n_sig, "alpha": alpha, "p_adjust": p_adjust}
    notes = {**notes, "n_features_tested": int(counts.shape[1]),
             "n_features_dropped_zero": n_dropped_feat}

    return AnalysisResult(
        verb="compare_abundance",
        table=table,
        stats=stats,
        spec={"predictors": preds, "primary": primary, "covariates": preds[1:],
              "method": method, "subject": subject, "subset": subset, "modality": mod},
        plot={"kind": "volcano", "effect": "log2fc", "significance": "qvalue",
              "label": "feature", "note": f"log2 fold-change vs -log10(q) for {primary}"},
        notes=notes,
    )


def _ref(name: str) -> str:
    """A patsy-safe formula reference (bare when it's a valid identifier, else Q('...'))."""
    return name if name.isidentifier() else "Q(%r)" % name


def _run(method, counts, frame, formula, subject, alpha, p_adjust, pseudocount, kwargs):
    if method == "ancombc":
        from skbio.stats.composition import ancombc

        # scikit-bio's ANCOM-BC requires strictly positive input (it does not zero-handle).
        raw = ancombc(counts + pseudocount, frame, formula, alpha=alpha, p_adjust=p_adjust,
                      **kwargs).reset_index()
        return _standardize(raw, effect="Log2(FC)")
    if method == "dirmult_lme":
        from skbio.stats.composition import dirmult_lme

        if not subject:
            raise ValueError("method='dirmult_lme' needs a subject column for the random intercept.")
        raw = dirmult_lme(counts, frame, formula, grouping=subject, pseudocount=pseudocount,
                          p_adjust=p_adjust, **kwargs)
        return _standardize(raw, effect="Log2(FC)")
    raise ValueError(f"Unknown method {method!r} (have: 'ancombc', 'dirmult_lme').")


def _standardize(raw: pd.DataFrame, effect: str) -> pd.DataFrame:
    """Both skbio DA outputs -> a common tidy frame: feature, contrast, log2fc, p, q, signif."""
    return pd.DataFrame({
        "feature": raw["FeatureID"].astype(str),
        "contrast": raw["Covariate"].astype(str),
        "log2fc": raw[effect].astype(float),
        "pvalue": raw["pvalue"].astype(float),
        "qvalue": raw["qvalue"].astype(float),
        "signif": raw["Signif"].astype(bool),
    })


def _primary_rows(tidy: pd.DataFrame, primary: str) -> pd.DataFrame:
    """Keep only the predictor-of-interest's contrast rows (drop Intercept + covariates)."""
    token = _ref(primary)
    contrast = tidy["contrast"]
    mask = (contrast == token) | contrast.str.startswith(token + "[")
    out = tidy[mask].copy()
    if out.empty:                       # continuous or unusual labelling — fall back to non-nuisance
        drop = {"Intercept"}
        out = tidy[~tidy["contrast"].isin(drop)].copy()
    return out.set_index("feature")[_CORE]
