"""VISTA importer — the mgCST community-type call (shotgun analogue of import_valencia).

REAL shape (FORMATS.md, validated against ``vista_mgCSTs.csv``): ``run_VISTA.R`` writes six
files; the **authoritative per-sample call** is ``mgCSTs_<...>.csv`` — a CSV whose first
(unnamed) column is the sample id, plus ``mgCST`` (the community-type label, e.g. ``"mgCST 11"``)
and ``max_YC_theta`` (the YC-θ of the *best-matching* mgCST). Two things the audit pinned down
and this importer honours:

* **No per-centroid similarities.** VISTA emits only ``max_YC_theta`` (best match), not θ against
  all 25 centroids — so, unlike CST, there is no ``mgcst_sim`` vector to route to ``.obsm``.
* **No scalar subtype in the call file.** The finer mgSs level lives in
  ``norm_counts_mgSs_mgCST_<...>.csv`` as a *feature matrix* (mgSs x sample), not a per-sample
  label, so mgCST_subtype is not a column here. (An mgSs modality is a later increment.)

We reshape only: label + θ into a sample-keyed frame ready for ``build_mudata(mgcst=...)``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def import_mgcst(mgcsts_csv) -> pd.DataFrame:
    """Import VISTA's ``mgCSTs_<...>.csv`` -> a sample-keyed mgCST call frame.

    Parameters
    ----------
    mgcsts_csv:
        Path to VISTA's ``mgCSTs_<...>.csv`` (first column = sample id, then ``mgCST`` and
        ``max_YC_theta``).

    Returns
    -------
    pandas.DataFrame
        Indexed by sample, columns ``mgCST`` (the label as written) and ``mgCST_score`` (θ, from
        ``max_YC_theta``). Ready for :func:`microfgt.io.build_mudata` via ``mgcst=`` — its labels
        land on the global ``.obs`` beside (never merged with) the 16S ``CST``.
    """
    raw = pd.read_csv(mgcsts_csv, index_col=0)
    raw.index = raw.index.astype(str)
    raw.index.name = "sample"
    missing = {"mgCST", "max_YC_theta"} - set(raw.columns)
    if missing:
        raise ValueError(
            f"{mgcsts_csv} is missing VISTA mgCST column(s) {sorted(missing)} "
            f"(have: {list(raw.columns)}); expected VISTA's mgCSTs_*.csv."
        )
    return pd.DataFrame(
        {"mgCST": raw["mgCST"].astype(str).to_numpy(),
         "mgCST_score": raw["max_YC_theta"].astype(float).to_numpy()},
        index=raw.index,
    )
