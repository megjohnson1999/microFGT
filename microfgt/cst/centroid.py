"""Centroid CST classifier — faithful reimplementation of VALENCIA (``Valencia.py``).

This is the one blessed CST method (constraint B): nearest-centroid by Yue–Clayton theta
to 13 fixed reference subCST centroids, exactly as VALENCIA does it. It sits behind
:func:`microfgt.cst.classify_cst` — the interface exists for genuine variants of this same
standard (e.g. a custom centroid set), not for rival CST classifiers. What VALENCIA can't
say about diffuse/continuum communities is surfaced by *augmenting* the CST label with
interpretable descriptors, not by computing CST a second way.

Fidelity notes vs ``Valencia.py``:
* Yue–Clayton theta = ``sum(p*q) / (sum((p-q)^2) + sum(p*q))`` over the taxon union.
  It is a sum over taxa, so taxon *order* is irrelevant — only the union set (zero-filled)
  and the relative abundances matter. We compute it vectorized (identical arithmetic).
* Relative abundance uses the *given* ``read_count`` (``Valencia.py`` divides by the
  ``read_count`` column, not by the recomputed taxon sum).
* ``subCST`` = argmax of the 13 ``_sim`` columns (pandas ``idxmax`` → first-max tie-break,
  same as VALENCIA); ``score`` = max; ``CST`` = subCST collapsed.
* Reference centroids are already relative abundances (VALENCIA uses them un-normalized).
"""

from __future__ import annotations

import warnings
from importlib import resources

import anndata as ad
import numpy as np
import pandas as pd

# subCST order, matching Valencia.py:86 and the bundled centroids file.
CST_ORDER = [
    "I-A", "I-B", "II", "III-A", "III-B", "IV-A", "IV-B",
    "IV-C0", "IV-C1", "IV-C2", "IV-C3", "IV-C4", "V",
]
# subCST -> CST collapse (Valencia.py:135).
_COLLAPSE = {
    "I-A": "I", "I-B": "I", "III-A": "III", "III-B": "III",
    "IV-C0": "IV-C", "IV-C1": "IV-C", "IV-C2": "IV-C", "IV-C3": "IV-C", "IV-C4": "IV-C",
}
# Bundled centroid sets, selectable by name. VALENCIA's centroids are keyed by taxon NAME, so
# the set must match the taxonomy assigner's naming (see design/method_log.md M4):
#   "2024" — VALENCIA2 centroids (356 taxa, modern names). Matches speciateIT v6 output, so it
#            is the DEFAULT: the 2020 set silently drops core BV taxa (Fannyhessea/Ca. Lachnocurva/
#            L. mulieris/…) that speciateIT v6 emits, degrading CST IV calls.
#   "2020" — original VALENCIA centroids (199 taxa, older names). Kept for reproducing the
#            published paper and for the validation gate (whose gold-standard data is 2020-named).
_CENTROID_SETS = {
    "2024": "VALENCIA2_CST_centroids_19Aug2024.csv",
    "2020": "cst_centroids_012920.csv",
}
_DEFAULT_CENTROIDS = "2024"


def load_reference_centroids(reference=None) -> pd.DataFrame:
    """Load subCST x taxon reference centroids (relative abundances), indexed by subCST.

    ``reference`` may be a bundled set name (``"2024"`` default, or ``"2020"``) or a path to a
    custom centroids CSV. The default is the 2024 set so the centroid method matches speciateIT's
    modern taxonomy out of the box (UX constraint A; see M4). Pass ``reference="2020"`` to
    reproduce the original VALENCIA paper."""
    if reference is None:
        reference = _DEFAULT_CENTROIDS
    if isinstance(reference, str) and reference in _CENTROID_SETS:
        with resources.files("microfgt.data").joinpath(_CENTROID_SETS[reference]).open() as fh:
            df = pd.read_csv(fh)
    else:
        df = pd.read_csv(reference)
    df = df.set_index("sub_CST")
    return df.reindex(CST_ORDER)


def _counts_and_read_count(composition, read_count=None):
    """Normalize input to (counts DataFrame [samples x taxa], read_count Series)."""
    if isinstance(composition, ad.AnnData):
        X = composition.layers["counts"] if "counts" in composition.layers else composition.X
        counts = pd.DataFrame(
            np.asarray(X),
            index=composition.obs_names.astype(str),
            columns=composition.var_names.astype(str),
        )
        if read_count is None and "read_count" in composition.obs:
            read_count = pd.Series(
                composition.obs["read_count"].to_numpy(), index=counts.index
            )
    elif isinstance(composition, pd.DataFrame):
        counts = composition.copy()
        counts.index = counts.index.astype(str)
    else:
        raise TypeError(
            "composition must be an anndata.AnnData (composition modality) or a "
            "samples x taxa pandas.DataFrame."
        )

    if read_count is None:
        read_count = counts.sum(axis=1)
    elif not isinstance(read_count, pd.Series):
        read_count = pd.Series(np.asarray(read_count), index=counts.index)
    else:
        read_count = read_count.copy()
        read_count.index = read_count.index.astype(str)
        read_count = read_count.reindex(counts.index)
    return counts, read_count


def classify_centroid(composition, reference=None, read_count=None) -> pd.DataFrame:
    """Assign each sample to a (sub)CST by nearest reference centroid (Yue–Clayton theta).

    Parameters
    ----------
    composition:
        ``composition`` modality as an AnnData (taxa = ``var``, counts in ``layers['counts']``
        or ``X``; ``obs['read_count']`` used if present) or a samples x taxa DataFrame.
    reference:
        Path to a reference-centroids CSV (default: VALENCIA's bundled centroids).
    read_count:
        Optional per-sample total reads (Series or array). Defaults to ``obs['read_count']``
        if present, else the per-sample taxon sum.

    Returns
    -------
    pandas.DataFrame
        Indexed by sample, with the 13 ``<subCST>_sim`` columns, then ``subCST``, ``score``,
        ``CST`` — the same trailing columns VALENCIA appends.
    """
    centroids = load_reference_centroids(reference)
    counts, rc = _counts_and_read_count(composition, read_count)

    # Taxon union (sample taxa first, then centroid-only taxa); order is irrelevant to theta.
    sample_taxa = list(counts.columns)
    overlap = set(map(str, sample_taxa)) & set(map(str, centroids.columns))
    seen = set(sample_taxa)
    all_taxa = sample_taxa + [t for t in centroids.columns if t not in seen]

    rel = counts.reindex(columns=all_taxa).fillna(0.0).div(rc, axis=0).fillna(0.0)
    C = centroids.reindex(columns=all_taxa).fillna(0.0)

    P = rel.to_numpy(dtype=float)   # n_samples x n_taxa  (relative abundances)
    Q = C.to_numpy(dtype=float)     # 13 x n_taxa         (centroid abundances)

    prod = P @ Q.T                                  # sum(p*q),  n x 13
    diff_sq = (P**2).sum(1)[:, None] + (Q**2).sum(1)[None, :] - 2 * prod  # sum((p-q)^2)
    with np.errstate(divide="ignore", invalid="ignore"):
        sim = np.where((diff_sq + prod) > 0, prod / (diff_sq + prod), 0.0)

    sim_df = pd.DataFrame(
        sim, index=counts.index, columns=[f"{c}_sim" for c in CST_ORDER]
    )
    out = sim_df.copy()
    out["subCST"] = sim_df.idxmax(axis=1).str.replace("_sim", "", regex=False)
    out["score"] = sim_df.max(axis=1)
    out["CST"] = out["subCST"].replace(_COLLAPSE)
    out.index.name = "sample"

    # Guard the silent-garbage case: theta is 0 when a sample shares no taxa with the reference
    # centroids, and idxmax then collapses it to the first subCST with score 0 — a confident-
    # looking wrong call. Surface it rather than return it quietly. The usual cause is a
    # taxon-naming mismatch (e.g. GTDB names like 'Gardnerella_vaginalis(RS_GCF_...)' vs the
    # speciateIT/VALENCIA name 'Gardnerella_vaginalis'), not real biology.
    if not overlap:
        warnings.warn(
            f"CST: none of the {len(sample_taxa)} input taxa match the reference centroid names, "
            f"so every sample scored 0 and collapsed to '{CST_ORDER[0]}'. The CST calls are not "
            "meaningful — this is a taxon-naming mismatch; reconcile taxon names to the VALENCIA "
            "convention (underscored binomials, no accession suffixes).",
            stacklevel=2,
        )
    elif (out["score"] <= 0).any():
        n0 = int((out["score"] <= 0).sum())
        warnings.warn(
            f"CST: {n0} of {len(out)} sample(s) scored 0 (no overlap with reference taxa) and "
            f"defaulted to '{CST_ORDER[0]}' — check those samples' taxonomy.",
            stacklevel=2,
        )
    return out
