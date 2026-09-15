"""CST layer — the ``classify_cst(composition, ...)`` interface, one blessed method.

CST is a single field standard: **VALENCIA**, ported as ``centroid`` and validated to
99.94% subCST agreement. It is not computed multiple ways and it is not competed against
rival classifiers — the community structure VALENCIA flattens is surfaced by *augmenting*
CST with interpretable descriptors (dominant taxon, % dominant, effective # taxa), not by
swapping in alternative CST methods.

``register_method`` remains as a genuine extension point (e.g. a user's own centroid set or
a re-derivation), but ``centroid`` is the one method microFGT ships and blesses. There is no
namespacing and no per-method provenance stamp: with a single standard method there is
nothing to disambiguate.
"""

from __future__ import annotations

from typing import Callable

from microfgt.cst.centroid import CST_ORDER, classify_centroid, load_reference_centroids

_METHODS: dict[str, Callable] = {}


def register_method(name: str, fn: Callable) -> None:
    """Register a CST method behind the ``classify_cst`` interface.

    An extension point, not an invitation to compute CST several ways: microFGT ships and
    blesses exactly one method (``centroid`` = VALENCIA). Use this to plug in a variant of
    that same standard (e.g. a custom centroid set), not a rival classifier.
    """
    _METHODS[name] = fn


def available_methods() -> list[str]:
    return sorted(_METHODS)


def classify_cst(composition, method: str = "centroid", **kwargs):
    """Assign samples to CSTs from their taxonomic composition.

    Parameters
    ----------
    composition:
        The ``composition`` modality (AnnData) or a samples x taxa DataFrame.
    method:
        CST method name (default ``"centroid"`` — faithful VALENCIA, the one blessed method).
    **kwargs:
        Passed through to the method (e.g. ``reference``, ``read_count`` for centroid).

    Returns
    -------
    pandas.DataFrame indexed by sample with ``<subCST>_sim`` columns, ``subCST``, ``score``,
    ``CST`` — ready to attach as sample-level annotation via ``build_mudata(cst=...)``.
    """
    if method not in _METHODS:
        raise ValueError(
            f"Unknown CST method {method!r}; available: {available_methods()}"
        )
    return _METHODS[method](composition, **kwargs)


register_method("centroid", classify_centroid)

__all__ = [
    "classify_cst",
    "classify_centroid",
    "load_reference_centroids",
    "register_method",
    "available_methods",
    "CST_ORDER",
]
