"""Minimal matplotlib helpers. Each returns an Axes so callers can compose/save.

:func:`render` is the single entry the dashboard (and the CLI) use: hand it any
:class:`~microfgt.analysis.results.AnalysisResult` and it draws that result from its own
``plot`` spec + ``data`` — one renderer for every verb, so the surface never re-encodes how a
result is drawn.
"""

from __future__ import annotations

import numpy as np


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover - exercised only without matplotlib
        raise ImportError(
            "Plotting requires matplotlib. Install it with: pip install 'microfgt[viz]'."
        ) from e
    return plt


def render(result, ax=None):
    """Draw an :class:`AnalysisResult` from its declarative plot spec. Returns an Axes.

    Dispatches on ``result.plot['kind']`` (box / ordination / scatter / heatmap / volcano),
    drawing from ``result.data`` when present, else ``result.table``. This is the shared
    rendering path Layer 3 sits on.
    """
    kind = result.plot.get("kind")
    renderer = {
        "box": _render_box,
        "ordination": _render_ordination,
        "scatter": _render_scatter,
        "heatmap": _render_heatmap,
        "volcano": _render_volcano,
    }.get(kind)
    if renderer is None:
        raise ValueError(f"No renderer for plot kind {kind!r} (have: box/ordination/scatter/"
                         "heatmap/volcano).")
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots()
    renderer(result, ax)
    ax.set_title(result.summary(), fontsize="small")
    return ax


def _render_box(result, ax):
    x, y = result.plot["x"], result.plot["y"]
    d = result.data
    groups = d[x].astype(str)
    labels = sorted(groups.unique())
    ax.boxplot([d.loc[(groups == g).to_numpy(), y].dropna().to_numpy() for g in labels])
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlabel(x)
    ax.set_ylabel(y)


def _render_ordination(result, ax):
    d = result.data
    xk, yk, color = result.plot["x"], result.plot["y"], result.plot.get("color")
    prop = result.stats.get("proportion_explained")
    if color and color in d:
        labels = d[color].astype(str).to_numpy()
        for lab in sorted(set(labels)):
            m = labels == lab
            ax.scatter(d[xk].to_numpy()[m], d[yk].to_numpy()[m], label=lab, s=18)
        ax.legend(title=color, fontsize="small")
    else:
        ax.scatter(d[xk].to_numpy(), d[yk].to_numpy(), s=18)
    ax.set_xlabel(f"{xk} ({prop[0]:.1%})" if prop else xk)
    ax.set_ylabel(f"{yk} ({prop[1]:.1%})" if prop else yk)


def _render_scatter(result, ax):
    x, y = result.plot["x"], result.plot["y"]
    d = result.data
    ax.scatter(d[x].to_numpy(), d[y].to_numpy(), s=18)
    ax.set_xlabel(x)
    ax.set_ylabel(y)


def _render_heatmap(result, ax):
    ct = result.table                       # contingency table (index x columns of counts)
    im = ax.imshow(ct.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels(ct.columns, rotation=45, ha="right", fontsize="small")
    ax.set_yticks(range(ct.shape[0]))
    ax.set_yticklabels(ct.index, fontsize="small")
    ax.set_xlabel(result.plot.get("y", ct.columns.name))
    ax.set_ylabel(result.plot.get("x", ct.index.name))
    ax.figure.colorbar(im, ax=ax, label="count")


def _render_volcano(result, ax):
    t = result.table                        # per-feature: log2fc + qvalue + signif
    lfc = t["log2fc"].to_numpy()
    q = np.clip(t["qvalue"].to_numpy(), 1e-300, 1.0)
    nlq = -np.log10(q)
    sig = t["signif"].to_numpy().astype(bool)
    ax.scatter(lfc[~sig], nlq[~sig], s=12, color="0.6", label="ns")
    ax.scatter(lfc[sig], nlq[sig], s=14, color="crimson", label="significant")
    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, color="0.4")
    ax.set_xlabel("log2 fold-change")
    ax.set_ylabel("-log10(q)")
    ax.legend(fontsize="small")


def ordination_scatter(adata, key: str = "X_pcoa", color: str | None = None, ax=None):
    """Scatter of the first two ordination axes, optionally colored by an ``obs`` column."""
    plt = _plt()
    coords = np.asarray(adata.obsm[key])
    if ax is None:
        _, ax = plt.subplots()
    if color is not None and color in adata.obs:
        labels = adata.obs[color].astype(str).to_numpy()
        for lab in sorted(set(labels)):
            m = labels == lab
            ax.scatter(coords[m, 0], coords[m, 1], label=lab, s=18)
        ax.legend(title=color, fontsize="small")
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=18)
    ax.set_xlabel(f"{key}[0]")
    ax.set_ylabel(f"{key}[1]")
    return ax


def alpha_boxplot(adata, alpha_key: str, group_key: str, ax=None):
    """Boxplot of an alpha-diversity ``obs`` column grouped by another ``obs`` column."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots()
    groups = adata.obs[group_key].astype(str)
    vals = adata.obs[alpha_key].to_numpy()
    labels = sorted(set(groups))
    ax.boxplot([vals[(groups == g).to_numpy()] for g in labels])
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)            # matplotlib 3.9+ dropped boxplot(labels=)
    ax.set_xlabel(group_key)
    ax.set_ylabel(alpha_key)
    return ax
