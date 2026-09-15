"""mgCST layer — the ``classify_mgcst(function, ...)`` seam, mirroring the CST layer.

mgCST is the shotgun-side community-type call: VISTA's random forest over the VIRGO2 gene
(``function``) profile, with a YC-θ score. It is the metagenomic parallel to VALENCIA/CST, and
it is kept **separate** from CST — the two are not number-comparable (different biases,
references, resolutions), so they live in separate ``.obs`` columns and are reconciled only as
an explicit analysis verb, never merged (see ``design/shotgun_arm_design.md``).

This module is the method seam, structurally identical to :mod:`microfgt.cst`. The one blessed
method is **VISTA** (``classify_mgcst_vista`` in the orchestration layer), which shells out to R;
it is registered at import. An unknown method name still raises helpfully, exactly like
:func:`microfgt.cst.classify_cst`.

Note the division of labour, parallel to the 16S side: ``classify_mgcst`` *computes* mgCST
(runs VISTA); :func:`microfgt.io.import_mgcst` *imports* an existing VISTA output (the analogue
of :func:`microfgt.io.import_valencia`).
"""

from __future__ import annotations

from typing import Callable

_METHODS: dict[str, Callable] = {}


def register_method(name: str, fn: Callable) -> None:
    """Register an mgCST method behind the ``classify_mgcst`` interface.

    An extension point, not an invitation to compute mgCST several ways: microFGT blesses one
    method (VISTA). Use this to plug in a variant of that same standard, not a rival classifier.
    The VISTA method registers itself here once the orchestration layer is built.
    """
    _METHODS[name] = fn


def available_methods() -> list[str]:
    return sorted(_METHODS)


def classify_mgcst(function, method: str = "vista", **kwargs):
    """Assign samples to mgCSTs from their VIRGO2 gene (``function``) profile.

    Parameters
    ----------
    function:
        The ``function`` modality (gene x sample AnnData) from
        :func:`microfgt.io.import_virgo2`.
    method:
        mgCST method name (default ``"vista"``). The blessed VISTA method is registered by the
        orchestration layer; until then no method is available and this raises.
    **kwargs:
        Passed through to the method.

    Returns
    -------
    pandas.DataFrame indexed by sample with ``mgCST`` (the community-type label) and
    ``mgCST_score`` (θ, VISTA's best-match ``max_YC_theta``) — same shape as
    :func:`microfgt.io.import_mgcst`, ready for ``build_mudata(mgcst=...)``. VISTA emits no
    per-centroid similarities, so there is no ``mgcst_sim`` block.
    """
    if method not in _METHODS:
        raise ValueError(
            f"Unknown mgCST method {method!r}; available: {available_methods()}. "
            "To use an existing VISTA output instead of running it, import it with "
            "microfgt.io.import_mgcst."
        )
    return _METHODS[method](function, **kwargs)


def _register_builtin_methods() -> None:
    # VISTA is the one blessed method; it shells out to R (orchestration layer).
    from microfgt.orchestrate.vista import classify_mgcst_vista

    register_method("vista", classify_mgcst_vista)


_register_builtin_methods()

__all__ = ["classify_mgcst", "register_method", "available_methods"]
