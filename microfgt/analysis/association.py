"""``associate`` — is one sample variable associated with another?

The canonical FGT use is **CST vs a clinical variable**: does community state type track BV
status, a diagnosis, pregnancy outcome (categorical), or pH / Nugent score / age (continuous)?
But nothing here hardcodes CST — it's a general two-variable association verb that dispatches
on the column dtypes, which also gives us continuous-vs-continuous correlation for free (the
seed of cross-modal "taxon vs metabolite" later).

Dispatch (``method="auto"``):

* **categorical × categorical** — contingency table + chi-square (Cramér's V effect size); a
  2×2 table uses Fisher's exact instead (with an odds ratio), and low expected cell counts are
  flagged since chi-square is unreliable there.
* **categorical × continuous** — the continuous values compared across the groups
  (Mann–Whitney for 2, Kruskal–Wallis for >2).
* **continuous × continuous** — Spearman rank correlation (robust default; ``"pearson"`` to
  force linear).

Pure Python (scipy). Returns the uniform :class:`~microfgt.analysis.results.AnalysisResult`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from microfgt.analysis._frame import analysis_frame, is_categorical
from microfgt.analysis.results import AnalysisResult


def associate(
    data,
    x: str,
    y: str,
    *,
    method: str = "auto",
    modality=None,
    subset=None,
) -> AnalysisResult:
    """Test association between two ``obs`` variables ``x`` and ``y``.

    Parameters
    ----------
    data:
        MuData or AnnData; ``x`` / ``y`` resolve from global obs OR any modality's obs.
    x, y:
        The two ``obs`` variable names. Order does not matter.
    method:
        ``"auto"`` (default) dispatches on dtype (see module docstring). Force with
        ``"chi2"`` / ``"fisher"`` (categorical×categorical), ``"kruskal"`` / ``"mannwhitney"``
        (categorical×continuous), or ``"spearman"`` / ``"pearson"`` (continuous×continuous).
    subset:
        Restrict samples first: a query string or ``{column: value(s)}``.

    Returns
    -------
    AnalysisResult
        ``table`` = the contingency table (cat×cat), per-group summary (cat×cont), or a
        one-row correlation summary (cont×cont); ``stats`` carries the headline
        test/statistic/p-value/effect size.
    """
    frame, notes = analysis_frame(data, [x, y], modality=modality, subset=subset)
    if frame.shape[0] < 3:
        raise ValueError(f"Too few samples ({frame.shape[0]}) after subset/missing-drop to test.")

    x_cat, y_cat = is_categorical(frame[x]), is_categorical(frame[y])
    kind = _dispatch(method, x_cat, y_cat)

    plot_data = None
    if kind in ("chi2", "fisher"):
        table, stats = _cat_cat(frame, x, y, kind)
        plot = {"kind": "heatmap", "x": x, "y": y, "note": "contingency table"}
    elif kind in ("kruskal", "mannwhitney"):
        cat, cont = (x, y) if x_cat else (y, x)
        table, stats = _cat_cont(frame, cat, cont, kind)
        plot = {"kind": "box", "x": cat, "y": cont, "note": f"{cont} across {cat}"}
        plot_data = frame[[cat, cont]].copy()
        plot_data[cont] = pd.to_numeric(plot_data[cont], errors="coerce")
    else:
        table, stats = _cont_cont(frame, x, y, kind)
        plot = {"kind": "scatter", "x": x, "y": y, "note": "with correlation"}
        plot_data = frame[[x, y]].apply(pd.to_numeric, errors="coerce")

    return AnalysisResult(
        verb="associate",
        table=table,
        stats=stats,
        spec={"x": x, "y": y, "method": stats["test"], "subset": subset,
              "x_type": "categorical" if x_cat else "continuous",
              "y_type": "categorical" if y_cat else "continuous"},
        plot=plot,
        notes=notes,
        data=plot_data,
    )


def _dispatch(method: str, x_cat: bool, y_cat: bool) -> str:
    if method != "auto":
        return method
    if x_cat and y_cat:
        return "chi2"                       # narrowed to fisher for 2x2 inside _cat_cat
    if x_cat != y_cat:
        return "kruskal"                    # narrowed to mannwhitney for 2 groups inside
    return "spearman"


def _cramers_v(chi2: float, table: np.ndarray) -> float:
    n = table.sum()
    r, c = table.shape
    denom = n * (min(r, c) - 1)
    return float(np.sqrt(chi2 / denom)) if denom > 0 else float("nan")


def _cat_cat(frame, x, y, kind):
    from scipy.stats import chi2_contingency, fisher_exact

    ct = pd.crosstab(frame[x].astype(str), frame[y].astype(str))
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        raise ValueError(
            f"Association needs >=2 levels in each of {x!r} ({ct.shape[0]}) and "
            f"{y!r} ({ct.shape[1]})."
        )
    observed = ct.to_numpy()

    if kind == "fisher" or (kind == "chi2" and observed.shape == (2, 2)):
        odds, p = fisher_exact(observed)
        stats = {"test": "Fisher exact", "pvalue": float(p), "odds_ratio": float(odds),
                 "effect_size": float(odds), "effect": "odds_ratio"}
    else:
        chi2, p, dof, expected = chi2_contingency(observed)
        v = _cramers_v(chi2, observed)
        stats = {"test": "chi-square", "statistic": float(chi2), "dof": int(dof),
                 "pvalue": float(p), "effect_size": v, "effect": "cramers_v"}
        low = int((expected < 5).sum())
        if low:
            stats["low_expected_cells"] = low
            stats["warning"] = (f"{low} cell(s) have expected count <5; chi-square is "
                                "unreliable here — prefer Fisher/exact or collapse levels.")
    return ct, stats


def _cat_cont(frame, cat_col, cont_col, kind):
    from scipy.stats import kruskal, mannwhitneyu

    values = pd.to_numeric(frame[cont_col], errors="coerce")
    groups = frame[cat_col].astype(str)
    levels = sorted(groups.unique())
    samples = [values[(groups == g).to_numpy()].dropna().to_numpy() for g in levels]
    if len(levels) < 2 or any(len(s) < 1 for s in samples):
        raise ValueError(f"Need >=2 non-empty groups in {cat_col!r}; got "
                         f"{dict(zip(levels, map(len, samples)))}.")
    n = int(sum(map(len, samples)))
    if kind == "mannwhitney" or (kind == "kruskal" and len(levels) == 2):
        stat, p = mannwhitneyu(samples[0], samples[1], alternative="two-sided")
        n1, n2 = len(samples[0]), len(samples[1])
        stats = {"test": "Mann–Whitney U", "statistic": float(stat), "pvalue": float(p),
                 "effect_size": float(1 - 2 * stat / (n1 * n2)), "effect": "rank_biserial"}
    else:
        stat, p = kruskal(*samples)
        k = len(levels)
        stats = {"test": "Kruskal–Wallis", "statistic": float(stat), "pvalue": float(p),
                 "effect_size": float((stat - k + 1) / (n - k)) if n > k else np.nan,
                 "effect": "epsilon_squared"}
    table = pd.DataFrame(
        {"n": [len(s) for s in samples],
         "median": [float(np.median(s)) if len(s) else np.nan for s in samples],
         "mean": [float(np.mean(s)) if len(s) else np.nan for s in samples]},
        index=pd.Index(levels, name=cat_col),
    )
    return table, stats


def _cont_cont(frame, x, y, kind):
    from scipy.stats import pearsonr, spearmanr

    a = pd.to_numeric(frame[x], errors="coerce")
    b = pd.to_numeric(frame[y], errors="coerce")
    keep = a.notna() & b.notna()
    a, b = a[keep].to_numpy(), b[keep].to_numpy()
    if len(a) < 3:
        raise ValueError("Need >=3 paired non-missing values for a correlation.")
    if kind == "pearson":
        r, p = pearsonr(a, b)
        name, coef = "Pearson r", "r"
    else:
        r, p = spearmanr(a, b)
        name, coef = "Spearman rho", "rho"
    stats = {"test": name, "statistic": float(r), "pvalue": float(p),
             "effect_size": float(r), "effect": coef}
    table = pd.DataFrame({coef: [float(r)], "pvalue": [float(p)], "n": [len(a)]},
                         index=pd.Index([f"{x}~{y}"], name="pair"))
    return table, stats
