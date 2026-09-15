"""The classify_mgcst method seam — structurally parallel to the CST seam.

The blessed VISTA method ships with the shotgun orchestration layer (a later increment), so out
of the box no method is registered and the seam must fail helpfully rather than silently.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microfgt import mgcst


def _function(samples, genes):
    X = np.arange(len(samples) * len(genes), dtype=np.float32).reshape(len(samples), len(genes))
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(samples, name="sample")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )


def test_vista_is_the_registered_method():
    # VISTA (the R random forest) is the one blessed method, registered at import.
    assert mgcst.available_methods() == ["vista"]


def test_classify_mgcst_raises_on_unknown_method():
    func = _function(["s1", "s2"], ["g1", "g2", "g3"])
    with pytest.raises(ValueError, match="Unknown mgCST method"):
        mgcst.classify_mgcst(func, method="does-not-exist")


def test_vista_method_needs_compiled_or_function():
    # Dispatch reaches the VISTA method; with neither compiled= nor a function it explains why.
    with pytest.raises(ValueError, match="compiled=|function"):
        mgcst.classify_mgcst(None, vista_repo="/nope", outdir="/tmp/x")


def test_register_method_makes_it_dispatchable():
    func = _function(["s1", "s2"], ["g1", "g2"])
    sentinel = pd.DataFrame({"mgCST": [1, 2]}, index=["s1", "s2"])
    try:
        mgcst.register_method("fake", lambda function, **kw: sentinel)
        assert "fake" in mgcst.available_methods()
        out = mgcst.classify_mgcst(func, method="fake")
        assert list(out["mgCST"]) == [1, 2]
    finally:
        mgcst._METHODS.pop("fake", None)  # keep the global registry clean for other tests
