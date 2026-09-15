"""Helpers shared by the analysis wrappers."""

from __future__ import annotations

import numpy as np


def get_count_matrix(adata, layer: str | None = "counts") -> np.ndarray:
    """Return a samples x features float matrix from ``layer`` (or ``X`` if absent)."""
    if layer is not None and layer in adata.layers:
        X = adata.layers[layer]
    else:
        X = adata.X
    return np.asarray(X, dtype=float)
