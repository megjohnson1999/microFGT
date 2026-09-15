"""Alpha/beta diversity + ordination — BOUGHT from scikit-bio, not built.

Results are written back onto the modality in AnnData-idiomatic places: alpha -> ``obs``,
beta distances -> ``obsp`` + ``uns``, ordination coordinates -> ``obsm`` + ``uns``.

Beta diversity / ordination are undefined for zero-count samples (Bray–Curtis is 0/0), so
those samples are excluded from the computation and left as NaN in the result rather than
silently dropped from the object; the count is recorded in ``uns``.
"""

from __future__ import annotations

import numpy as np
from skbio.diversity import alpha_diversity as _skbio_alpha
from skbio.diversity import beta_diversity as _skbio_beta
from skbio.stats.ordination import pcoa as _pcoa

from microfgt.analysis._util import get_count_matrix


def alpha_diversity(adata, metric: str = "shannon", layer: str = "counts", key_added=None):
    """Per-sample alpha diversity; store in ``obs[key_added]`` (default ``alpha_<metric>``)."""
    ids = list(adata.obs_names)
    series = _skbio_alpha(metric, get_count_matrix(adata, layer), ids=ids)
    key = key_added or f"alpha_{metric}"
    adata.obs[key] = series.reindex(ids).to_numpy()
    return adata


def _nonzero_distance_matrix(adata, metric, layer):
    """Bray–Curtis-style distances on samples that have reads; returns (dm, idx, ids)."""
    X = get_count_matrix(adata, layer)
    ids = list(adata.obs_names)
    mask = X.sum(axis=1) > 0
    idx = np.where(mask)[0]
    sub_ids = [ids[i] for i in idx]
    dm = _skbio_beta(metric, X[mask], ids=sub_ids)
    return dm, idx, ids, int((~mask).sum())


def beta_diversity(adata, metric: str = "braycurtis", layer: str = "counts", key_added=None):
    """Pairwise beta-diversity distances; store the square matrix in ``obsp[key_added]``.

    Zero-count samples are excluded (NaN in the stored matrix). Returns the scikit-bio
    ``DistanceMatrix`` for the samples that had reads."""
    dm, idx, ids, n_skipped = _nonzero_distance_matrix(adata, metric, layer)
    key = key_added or f"beta_{metric}"
    full = np.full((len(ids), len(ids)), np.nan)
    full[np.ix_(idx, idx)] = dm.data
    adata.obsp[key] = full
    adata.uns[key] = {"metric": metric, "ids": ids, "n_skipped_zero_count": n_skipped}
    return dm


def ordinate(adata, metric: str = "braycurtis", layer: str = "counts", key_added: str = "X_pcoa"):
    """PCoA over a beta-diversity distance matrix; coords -> ``obsm[key_added]``.

    Zero-count samples are excluded and left as NaN coordinates. Returns the scikit-bio
    ``OrdinationResults``."""
    dm, idx, ids, n_skipped = _nonzero_distance_matrix(adata, metric, layer)
    res = _pcoa(dm)
    coords = res.samples
    full = np.full((len(ids), coords.shape[1]), np.nan)
    full[idx] = coords.to_numpy()
    adata.obsm[key_added] = full
    adata.uns[key_added] = {
        "metric": metric,
        "proportion_explained": res.proportion_explained.to_numpy(),
        "n_skipped_zero_count": n_skipped,
    }
    return res
