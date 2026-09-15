"""Compositional transforms — BOUGHT from scikit-bio, not built.

Microbiome counts are compositional, so defaults use relative abundance / CLR and never
naive stats on raw counts (constraint B). These are thin wrappers that write the transform
into a new ``layer`` of the modality, leaving the raw counts intact.

Zero-count samples (no reads at all) can't be compositionally transformed; rather than
silently drop them (constraint B — honest reconciliation) we leave their transformed row
zero-filled and record how many in ``uns['analysis_notes']``.
"""

from __future__ import annotations

import numpy as np
from skbio.stats.composition import clr as _clr
from skbio.stats.composition import closure, multi_replace

from microfgt.analysis._util import get_count_matrix


def _note(adata, key, value):
    adata.uns.setdefault("analysis_notes", {})[key] = value


def relative_abundance(adata, layer: str = "counts", key_added: str = "relabund"):
    """Close each sample to sum 1 (relative abundance); store in ``layers[key_added]``."""
    X = get_count_matrix(adata, layer)
    rowsum = X.sum(axis=1, keepdims=True)
    rel = np.divide(X, rowsum, out=np.zeros_like(X, dtype=float), where=rowsum > 0)
    adata.layers[key_added] = rel
    n_zero = int((rowsum.ravel() == 0).sum())
    if n_zero:
        _note(adata, "relabund_zero_count_samples", n_zero)
    return adata


def clr_transform(adata, layer: str = "counts", key_added: str = "clr"):
    """Centered log-ratio transform; store in ``layers[key_added]``.

    Zeros are handled by multiplicative replacement on the closed composition before the
    log-ratio (CLR is undefined at zero). Zero-count samples are left as a zero row."""
    X = get_count_matrix(adata, layer)
    mask = X.sum(axis=1) > 0
    out = np.zeros((X.shape[0], X.shape[1]), dtype=float)
    if mask.any():
        out[mask] = _clr(multi_replace(closure(X[mask])))
    adata.layers[key_added] = out
    n_zero = int((~mask).sum())
    if n_zero:
        _note(adata, "clr_zero_count_samples", n_zero)
    return adata
