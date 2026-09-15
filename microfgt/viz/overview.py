"""Standard FGT overview plots — the "what is this cohort?" visuals for the dashboard Home tab.

These take the whole object (a MuData) and draw the canonical FGT views: CST makeup, the iconic
per-sample community-composition stack, the descriptor distributions that carry the two-regime
story (dominance, effective_taxa), and a PCoA of cohort structure. Each returns a matplotlib
Axes/Figure; missing columns are skipped gracefully. Prototype-level styling (tab20 palette) —
polish is a follow-up.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from microfgt.analysis._frame import get_assay, merged_obs
from microfgt.viz.plots import _plt


def _taxon_rel(mdata) -> pd.DataFrame:
    """samples x taxa relative-abundance frame from the taxon-grain assay."""
    adata, _ = get_assay(mdata, prefer=("composition_taxon", "composition"))
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    counts = pd.DataFrame(
        np.asarray(X, dtype=float),
        index=adata.obs_names.astype(str), columns=adata.var_names.astype(str),
    )
    return counts.div(counts.sum(axis=1), axis=0).fillna(0.0)


def _short(label: str, n: int = 28) -> str:
    label = str(label)
    return label if len(label) <= n else label[: n - 1] + "…"


def cst_bar(mdata, group: str = "CST", ax=None):
    """Cohort makeup: sample count per CST."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 3))
    counts = merged_obs(mdata)[group].astype(str).value_counts().sort_index()
    ax.bar(counts.index, counts.to_numpy(), color="#4C72B0")
    ax.set_ylabel("samples")
    ax.set_xlabel(group)
    ax.set_title(f"Cohort by {group} (n={int(counts.sum())})", fontsize="medium")
    return ax


def community_stack(mdata, top_n: int = 12, group: str = "CST", ax=None):
    """The iconic per-sample composition stack: top-N taxa relative abundance, sorted by group."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(11, 3.5))
    rel = _taxon_rel(mdata)
    top = rel.mean().sort_values(ascending=False).head(top_n).index.tolist()
    mat = rel[top].copy()
    mat["Other"] = (1.0 - mat[top].sum(axis=1)).clip(lower=0)

    obs = merged_obs(mdata)
    g = obs[group].astype(str).reindex(mat.index) if group in obs.columns else pd.Series("all", index=mat.index)
    # order: group block, then within block by the sample's dominant top-taxon and its abundance
    order = pd.DataFrame(
        {"g": g, "dom": mat[top].idxmax(axis=1), "val": mat[top].max(axis=1)}
    ).sort_values(["g", "dom", "val"], ascending=[True, True, False]).index
    mat = mat.loc[order]
    g = g.loc[order]

    cols = top + ["Other"]
    palette = list(plt.cm.tab20(np.linspace(0, 1, len(top)))) + [(0.8, 0.8, 0.8, 1.0)]
    x = np.arange(len(mat))
    ax.stackplot(x, *[mat[c].to_numpy() for c in cols], colors=palette, labels=[_short(c) for c in cols])
    ax.set_xlim(0, max(1, len(mat)))
    ax.set_ylim(0, 1)
    ax.set_ylabel("relative abundance")
    ax.set_xlabel(f"samples (each column = one sample), sorted by {group}")

    # group-block separators + labels along the top
    if g.nunique() > 1:
        changes = np.where(g.to_numpy()[1:] != g.to_numpy()[:-1])[0] + 1
        for c in changes:
            ax.axvline(c, color="white", lw=0.8)
        bounds = [0, *changes, len(g)]
        for a, b in zip(bounds[:-1], bounds[1:]):
            ax.text((a + b) / 2, 1.02, str(g.iloc[a]), ha="center", va="bottom", fontsize="x-small")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize="x-small", frameon=False)
    return ax


def descriptor_box(mdata, value: str = "effective_taxa", group: str = "CST", ax=None):
    """A descriptor distribution by group (e.g. effective_taxa by CST)."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 3))
    obs = merged_obs(mdata)
    vals = pd.to_numeric(obs[value], errors="coerce")
    g = obs[group].astype(str)
    levels = sorted(g.dropna().unique())
    ax.boxplot([vals[(g == lv).to_numpy()].dropna().to_numpy() for lv in levels], showfliers=False)
    ax.set_xticks(range(1, len(levels) + 1))
    ax.set_xticklabels(levels)
    ax.set_xlabel(group)
    ax.set_ylabel(value)
    ax.set_title(f"{value} by {group}", fontsize="medium")
    return ax


def dominance_hist(mdata, value: str = "dominance_pct", ax=None):
    """Dominance distribution — the bimodal 'two-regime' signal."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 3))
    vals = pd.to_numeric(merged_obs(mdata)[value], errors="coerce").dropna().to_numpy()
    ax.hist(vals, bins=30, color="#55A868")
    ax.axvline(50, ls="--", lw=0.9, color="0.3")
    ax.set_xlabel("dominance (% of the top taxon)")
    ax.set_ylabel("samples")
    ax.set_title("Dominance distribution (bimodal → two regimes)", fontsize="medium")
    return ax


def ordination(mdata, color: str = "CST", max_samples: int = 600, ax=None):
    """PCoA (Bray–Curtis) of cohort structure, colored by group. Subsampled for speed on big cohorts."""
    from skbio.diversity import beta_diversity
    from skbio.stats.ordination import pcoa

    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    adata, _ = get_assay(mdata, prefer=("composition_taxon", "composition"))
    X = np.asarray(adata.layers["counts"] if "counts" in adata.layers else adata.X, dtype=float)
    ids = adata.obs_names.astype(str).tolist()
    idx = np.where(X.sum(axis=1) > 0)[0]
    subsampled = len(idx) > max_samples
    if subsampled:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(idx, size=max_samples, replace=False))
    sub_ids = [ids[i] for i in idx]
    dm = beta_diversity("braycurtis", X[idx], ids=sub_ids)
    res = pcoa(dm, number_of_dimensions=2)
    coords = res.samples.iloc[:, :2].to_numpy()
    prop = res.proportion_explained.iloc[:2].to_numpy()

    g = merged_obs(mdata)[color].astype(str).reindex(sub_ids) if color in merged_obs(mdata).columns else None
    if g is not None:
        for lv in sorted(g.dropna().unique()):
            m = (g == lv).to_numpy()
            ax.scatter(coords[m, 0], coords[m, 1], s=10, label=lv)
        ax.legend(title=color, fontsize="x-small", frameon=False)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=10)
    ax.set_xlabel(f"PCoA1 ({prop[0]:.1%})")
    ax.set_ylabel(f"PCoA2 ({prop[1]:.1%})")
    title = "Cohort structure (Bray–Curtis PCoA)"
    if subsampled:
        title += f" — {max_samples}-sample subsample"
    ax.set_title(title, fontsize="medium")
    return ax


def overview_figures(mdata, top_n: int = 12, group: str = "CST") -> dict:
    """Build all available overview figures as a name -> Figure dict (skips ones needing absent columns)."""
    plt = _plt()
    obs = merged_obs(mdata)
    figs: dict = {}

    def _fig(fn, **kw):
        fig, ax = plt.subplots(figsize=kw.pop("figsize", (5, 3)))
        fn(mdata, ax=ax, **kw)
        fig.tight_layout()
        return fig

    if group in obs.columns:
        figs["cst"] = _fig(cst_bar, group=group)
    figs["community"] = _fig(lambda m, ax: community_stack(m, top_n=top_n, group=group, ax=ax), figsize=(11, 3.8))
    if "effective_taxa" in obs.columns and group in obs.columns:
        figs["diversity"] = _fig(descriptor_box, value="effective_taxa", group=group)
    if "dominance_pct" in obs.columns:
        figs["dominance"] = _fig(dominance_hist)
    figs["ordination"] = _fig(ordination, figsize=(5, 4))
    return figs
