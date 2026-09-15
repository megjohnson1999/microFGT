"""CHARACTERIZE layer — per-sample descriptors that *augment* the CST label.

CST (VALENCIA) is one column that flattens community structure into a class. These descriptors
read out the structure it flattens — without replacing it or computing a second CST. They are
deterministic functions of the taxon-grain composition.

Two kinds, deliberately kept apart:

* **Intrinsic** (one right value, no free parameter) — stored as ``.obs`` columns:
  - ``dominant_taxon`` — the taxon with the highest relative abundance
  - ``dominance_pct``  — that taxon's abundance, as a percentage (0–100)
  - ``effective_taxa`` — the *effective number of taxa* (Hill q=1 = ``exp(Shannon)``): a
    **cutoff-free** read of evenness/diffuseness. ~1 when one organism dominates; large when the
    community is spread across many. This is the boundary-free answer to "how diffuse is it?"

* **Parameterized** (value depends on a knob) — NOT stored, computed on demand via
  :func:`taxa_over_threshold`. "How many taxa exceed X%?" has no privileged X, so baking a
  single cutoff into the object would smuggle back the arbitrary boundary CSTs already impose.
  The dashboard drives the threshold from a slider; the cutoff is an exploration knob, not a
  stored constant.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def _relative_abundance(composition_taxon: ad.AnnData):
    """(samples x taxa relative-abundance DataFrame, per-sample totals). Zero-count -> all zeros."""
    X = (
        composition_taxon.layers["counts"]
        if "counts" in composition_taxon.layers
        else composition_taxon.X
    )
    counts = pd.DataFrame(
        np.asarray(X, dtype=float),
        index=composition_taxon.obs_names.astype(str),
        columns=composition_taxon.var_names.astype(str),
    )
    totals = counts.sum(axis=1)
    rel = counts.div(totals, axis=0).fillna(0.0)
    return rel, totals


def describe_composition(composition_taxon: ad.AnnData) -> pd.DataFrame:
    """Per-sample **intrinsic** (parameter-free) descriptors from a taxon-grain composition.

    Parameters
    ----------
    composition_taxon:
        Taxon x sample AnnData (the ``composition_taxon`` assay from
        :func:`microfgt.io.collapse_to_taxon`).

    Returns
    -------
    pandas.DataFrame
        Indexed by sample, columns ``dominant_taxon`` / ``dominance_pct`` / ``effective_taxa``.
        Computed over *all* taxa present, including any ``Unclassified`` bucket. Zero-count
        samples get ``dominant_taxon = NaN``, ``dominance_pct = 0``, ``effective_taxa = NaN``.
        The cutoff-dependent "# taxa over X%" is intentionally *not* here — see
        :func:`taxa_over_threshold`.
    """
    rel, totals = _relative_abundance(composition_taxon)

    dominance = rel.max(axis=1)
    dominant = rel.idxmax(axis=1).where(dominance > 0, other=np.nan)

    # Effective number of taxa = exp(Shannon entropy, natural log). Cutoff-free evenness.
    p = rel.to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(p > 0, np.log(p), 0.0)
    shannon = -(p * logp).sum(axis=1)                       # natural-log Shannon per sample
    effective = pd.Series(np.exp(shannon), index=rel.index).where(totals > 0, other=np.nan)

    return pd.DataFrame(
        {
            "dominant_taxon": dominant.to_numpy(),
            "dominance_pct": (dominance * 100).to_numpy(),
            "effective_taxa": effective.to_numpy(),
        },
        index=pd.Index(rel.index, name="sample"),
    )


def taxa_over_threshold(composition_taxon: ad.AnnData, threshold: float = 0.10) -> pd.Series:
    """Per-sample count of taxa whose relative abundance exceeds ``threshold`` (a fraction 0–1).

    A **parameterized, on-demand** view — deliberately *not* stored on the object, because there
    is no privileged threshold. Materialize it at whatever cutoff the user chooses (the dashboard
    drives this from a slider). Zero-count samples get 0.
    """
    rel, totals = _relative_abundance(composition_taxon)
    count = (rel > threshold).sum(axis=1).where(totals > 0, other=0).astype(int)
    return pd.Series(count.to_numpy(), index=pd.Index(rel.index, name="sample"),
                     name="taxa_over_threshold")
