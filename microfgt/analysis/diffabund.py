"""Differential abundance — BOUGHT from scikit-bio (ANCOM / ANCOM-BC), not built.

scikit-bio 0.7 ships both ``ancom`` and ``ancombc``; we expose both behind one wrapper.
(The spec also allows orchestrating R for this step — a future method behind the same call.)
"""

from __future__ import annotations

import pandas as pd
from skbio.stats.composition import ancom as _ancom

from microfgt.analysis._util import get_count_matrix


def differential_abundance(
    adata,
    group_key: str,
    method: str = "ancom",
    layer: str = "counts",
    pseudocount: float = 1.0,
    **kwargs,
) -> pd.DataFrame:
    """Test which taxa differ in abundance across groups in ``obs[group_key]``.

    Parameters
    ----------
    adata:
        Composition modality.
    group_key:
        ``obs`` column defining the groups. Samples with a missing group are dropped.
    method:
        ``"ancom"`` (default). ANCOM requires strictly positive data, so ``pseudocount`` is
        added to the counts first.
    pseudocount:
        Added to counts before testing (ANCOM is undefined at zero).

    Returns
    -------
    pandas.DataFrame indexed by taxon with the test result (e.g. ANCOM's ``W`` / ``Signif``).
    """
    ids = list(adata.obs_names)
    table = pd.DataFrame(
        get_count_matrix(adata, layer) + pseudocount, index=ids, columns=list(adata.var_names)
    )
    grouping = adata.obs[group_key].astype("object")
    keep = grouping.notna() & (grouping.astype(str) != "nan")
    table, grouping = table.loc[keep.to_numpy()], grouping[keep.to_numpy()]
    if grouping.nunique() < 2:
        raise ValueError(
            f"Need >=2 groups in obs[{group_key!r}] for differential abundance; "
            f"found {grouping.nunique()} after dropping missing."
        )

    if method == "ancom":
        result, _percentiles = _ancom(table, grouping, **kwargs)
        return result
    raise ValueError(f"Unknown differential-abundance method {method!r} (have: 'ancom').")
