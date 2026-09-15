"""Hypothesis-test verbs — "is this community feature associated with the variable I picked?"

These turn the compute primitives (``diversity``, ``transforms``) into the tests a researcher
actually runs in an R Markdown: alpha diversity across a group (with covariate adjustment and
repeated-measures support), and community-level PERMANOVA on beta diversity. Each returns the
uniform :class:`~microfgt.analysis.results.AnalysisResult`.

Engine stance (the hybrid decision): the portable core is Python — scipy for the
nonparametric group tests, statsmodels for (mixed) linear models, scikit-bio for PERMANOVA /
dispersion. Multi-term, covariate-adjusted PERMANOVA the way reviewers know it (vegan's
``adonis2``) is orchestrated in R as a follow-up method behind the same verb; the Python path
here is single-factor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from microfgt.analysis._frame import analysis_frame, get_assay, is_categorical
from microfgt.analysis.results import AnalysisResult


# --- alpha diversity ~ predictors -----------------------------------------------------------
def compare_alpha(
    data,
    predictors,
    *,
    metric: str = "shannon",
    modality=None,
    subject: str | None = None,
    subset=None,
    test: str = "auto",
    layer: str = "counts",
) -> AnalysisResult:
    """Test per-sample alpha diversity against one or more ``obs`` predictors.

    Parameters
    ----------
    data:
        A MuData (predictors resolved from global obs OR a modality's obs) or an AnnData.
    predictors:
        ``obs`` column name(s) — what to compare/adjust by. A single string is fine.
    metric:
        Alpha metric (scikit-bio; default ``"shannon"``).
    modality:
        Which assay to compute diversity on (default: the taxon roll-up, else composition).
    subject:
        ``obs`` column of subject id for repeated measures -> a linear mixed model with a
        per-subject random intercept.
    subset:
        Restrict samples first: a query string or ``{column: value(s)}``.
    test:
        ``"auto"`` (default) picks: mixed model if ``subject`` given; else a nonparametric
        group test for a single categorical predictor (Mann–Whitney for 2 groups,
        Kruskal–Wallis for >2); else OLS. Force with ``"kruskal"``/``"mannwhitney"``/``"ols"``
        /``"mixedlm"``.

    Returns
    -------
    AnalysisResult
        ``table`` = per-group summary (group test) or model coefficients (OLS/mixed);
        ``stats`` carries the headline test/statistic/p-value/effect size.
    """
    preds = [predictors] if isinstance(predictors, str) else list(predictors)
    adata, mod = get_assay(data, modality)

    cols = preds + ([subject] if subject else [])
    frame, notes = analysis_frame(data, cols, modality=mod, subset=subset)

    # Alpha per sample, aligned onto the (subset, non-missing) frame.
    from skbio.diversity import alpha_diversity as _skbio_alpha

    X = np.asarray(adata.layers[layer] if layer in adata.layers else adata.X, dtype=float)
    ids = list(adata.obs_names.astype(str))
    alpha = _skbio_alpha(metric, X, ids=ids).reindex(frame.index)
    y = pd.to_numeric(alpha, errors="coerce")
    finite = np.isfinite(y.to_numpy())
    dropped_alpha = int((~finite).sum())
    frame = frame[finite]
    frame["alpha"] = y[finite].to_numpy()
    notes = {**notes, "n_used": len(frame), "n_dropped_alpha": dropped_alpha, "metric": metric}

    chosen = _choose_alpha_test(test, preds, subject, frame)
    if chosen in ("kruskal", "mannwhitney"):
        table, stats = _alpha_group_test(frame, preds[0], chosen)
    else:
        table, stats = _alpha_model(frame, preds, subject, chosen)

    y_name = f"alpha_{metric}"
    plot_data = frame[[preds[0]]].copy()
    plot_data[y_name] = frame["alpha"].to_numpy()

    return AnalysisResult(
        verb="compare_alpha",
        table=table,
        stats=stats,
        spec={"outcome": y_name, "predictors": preds, "subject": subject,
              "subset": subset, "modality": mod, "metric": metric},
        plot={"kind": "box", "y": y_name, "x": preds[0],
              "modality": mod, "note": "boxplot of alpha by the first predictor"},
        notes=notes,
        data=plot_data,
    )


def _choose_alpha_test(test, preds, subject, frame) -> str:
    if test != "auto":
        return test
    if subject:
        return "mixedlm"
    if len(preds) == 1 and is_categorical(frame[preds[0]]):
        return "mannwhitney" if frame[preds[0]].astype(str).nunique() == 2 else "kruskal"
    return "ols"


def _alpha_group_test(frame, group_col, kind):
    from scipy.stats import kruskal, mannwhitneyu

    groups = frame[group_col].astype(str)
    levels = sorted(groups.unique())
    samples = [frame.loc[(groups == g).to_numpy(), "alpha"].to_numpy() for g in levels]
    if any(len(s) == 0 for s in samples) or len(levels) < 2:
        raise ValueError(
            f"Need >=2 non-empty groups in {group_col!r}; got {dict(zip(levels, map(len, samples)))}."
        )
    n = int(sum(map(len, samples)))
    if kind == "mannwhitney":
        if len(levels) != 2:
            raise ValueError(f"mannwhitney needs exactly 2 groups; {group_col!r} has {len(levels)}.")
        stat, p = mannwhitneyu(samples[0], samples[1], alternative="two-sided")
        n1, n2 = len(samples[0]), len(samples[1])
        effect = float(1 - 2 * stat / (n1 * n2))            # rank-biserial correlation
        test_name, effect_name = "Mann–Whitney U", "rank_biserial"
    else:
        stat, p = kruskal(*samples)
        k = len(levels)
        effect = float((stat - k + 1) / (n - k)) if n > k else np.nan  # epsilon-squared
        test_name, effect_name = "Kruskal–Wallis", "epsilon_squared"

    table = pd.DataFrame(
        {"n": [len(s) for s in samples],
         "median": [float(np.median(s)) for s in samples],
         "mean": [float(np.mean(s)) for s in samples]},
        index=pd.Index(levels, name=group_col),
    )
    stats = {"test": test_name, "statistic": float(stat), "pvalue": float(p),
             "effect_size": effect, "effect": effect_name}
    return table, stats


def _q(name: str) -> str:
    """Patsy-safe reference to an arbitrary column name (may contain spaces/punctuation)."""
    return "Q(%r)" % name


def _alpha_model(frame, preds, subject, kind):
    import statsmodels.formula.api as smf

    terms = [f"C({_q(p)})" if is_categorical(frame[p]) else _q(p) for p in preds]
    formula = "alpha ~ " + " + ".join(terms)
    data = frame.rename(columns={c: c for c in frame.columns})  # keep names; Q() handles them

    if kind == "mixedlm":
        if not subject:
            raise ValueError("mixedlm requires a subject column for the random intercept.")
        model = smf.mixedlm(formula, data, groups=data[subject])
        res = model.fit()
        table = pd.DataFrame({"coef": res.params, "std_err": res.bse, "pvalue": res.pvalues})
        table.index.name = "term"
        nonint = [t for t in table.index if t not in ("Intercept", "Group Var")]
        stats = {"test": "LinearMixedModel", "formula": formula, "groups": subject,
                 "pvalue": float(table.loc[nonint[0], "pvalue"]) if nonint else None}
        return table, stats

    res = smf.ols(formula, data).fit()
    ci = res.conf_int()
    table = pd.DataFrame({"coef": res.params, "std_err": res.bse, "pvalue": res.pvalues,
                          "ci_low": ci[0], "ci_high": ci[1]})
    table.index.name = "term"
    stats = {"test": "OLS", "formula": formula, "R2": float(res.rsquared),
             "pvalue": float(res.f_pvalue), "F": float(res.fvalue)}
    return table, stats


# --- beta diversity: PERMANOVA --------------------------------------------------------------
def compare_beta(
    data,
    predictors,
    *,
    metric: str = "braycurtis",
    modality=None,
    subset=None,
    permutations: int = 999,
    dispersion: bool = True,
    ordinate: bool = True,
    layer: str = "counts",
) -> AnalysisResult:
    """PERMANOVA: does community composition differ across a categorical ``obs`` predictor?

    Single-factor (scikit-bio). Zero-count samples are excluded from the distance matrix and
    counted. Also runs a dispersion test (PERMDISP) by default, since a PERMANOVA signal can
    be driven by unequal within-group spread rather than a location shift.

    Multi-term / covariate-adjusted PERMANOVA (``adonis2``) is the R-orchestrated follow-up
    behind this same verb; for now extra predictors beyond the first are noted, not fit.
    """
    from skbio.stats.distance import permanova, permdisp

    from microfgt.analysis.diversity import _nonzero_distance_matrix

    preds = [predictors] if isinstance(predictors, str) else list(predictors)
    group_col = preds[0]
    adata, mod = get_assay(data, modality)
    frame, notes = analysis_frame(data, [group_col], modality=mod, subset=subset)

    dm, _idx, _ids, n_skipped = _nonzero_distance_matrix(adata, metric, layer)
    keep = [s for s in dm.ids if s in set(frame.index)]
    dm = dm.filter(keep)
    grouping = frame.loc[keep, group_col].astype(str)

    counts = grouping.value_counts()
    if grouping.nunique() < 2 or (counts < 2).any():
        raise ValueError(
            f"PERMANOVA needs >=2 groups each with >=2 samples in {group_col!r}; "
            f"got {counts.to_dict()} (after excluding {n_skipped} zero-count samples)."
        )

    pa = permanova(dm, grouping, permutations=permutations)
    r2 = _permanova_r2(dm, grouping)
    stats = {"test": "PERMANOVA", "pseudo_F": float(pa["test statistic"]),
             "pvalue": float(pa["p-value"]), "R2": r2, "permutations": permutations,
             "metric": metric}

    table = pd.DataFrame(
        {"n": counts.reindex(sorted(counts.index)).to_numpy()},
        index=pd.Index(sorted(counts.index), name=group_col),
    )
    if dispersion:
        pd_res = permdisp(dm, grouping, permutations=permutations)
        stats["dispersion_test"] = "PERMDISP"
        stats["dispersion_pvalue"] = float(pd_res["p-value"])

    notes = {**notes, "n_used": int(grouping.shape[0]), "n_skipped_zero_count": n_skipped}
    if len(preds) > 1:
        notes["ignored_predictors"] = preds[1:]
        notes["warning"] = ("Python PERMANOVA is single-factor; covariate-adjusted adonis2 "
                            "is the R-orchestrated follow-up.")

    plot_data = None
    if ordinate:
        from skbio.stats.ordination import pcoa

        res = pcoa(dm, number_of_dimensions=2)
        coords = res.samples.iloc[:, :2]
        plot_data = pd.DataFrame(
            {"PC1": coords.iloc[:, 0].to_numpy(), "PC2": coords.iloc[:, 1].to_numpy(),
             group_col: grouping.reindex(coords.index).to_numpy()},
            index=pd.Index(coords.index, name="sample"),
        )
        stats["proportion_explained"] = [float(v) for v in res.proportion_explained.iloc[:2]]

    return AnalysisResult(
        verb="compare_beta",
        table=table,
        stats=stats,
        spec={"predictors": preds, "subset": subset, "modality": mod, "metric": metric},
        plot={"kind": "ordination", "x": "PC1", "y": "PC2", "color": group_col, "modality": mod,
              "note": "PCoA scatter colored by the predictor"},
        notes=notes,
        data=plot_data,
    )


def _permanova_r2(dm, grouping) -> float:
    """R^2 = SSA/SST for a one-way PERMANOVA (the fraction of variation explained)."""
    d2 = np.asarray(dm.data) ** 2
    n = d2.shape[0]
    sst = d2[np.triu_indices(n, k=1)].sum() / n
    labels = np.asarray(list(grouping))
    ssw = 0.0
    for g in np.unique(labels):
        idx = np.where(labels == g)[0]
        ng = len(idx)
        if ng > 1:
            sub = d2[np.ix_(idx, idx)]
            ssw += sub[np.triu_indices(ng, k=1)].sum() / ng
    return float((sst - ssw) / sst) if sst > 0 else float("nan")
