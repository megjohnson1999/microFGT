"""VALENCIA importer — community state typing (CST) labels.

REAL shape (FORMATS.md, validated against ``valencia_genuine_output_head.csv``):
VALENCIA writes **one wide CSV** = the input (``sampleID``, ``read_count``, then one
column per taxon) **plus appended columns**: 13 ``<subCST>_sim`` similarity columns,
then ``subCST`` (argmax of the 13), ``score`` (max sim), and ``CST`` (subCST collapsed)
(``Valencia.py:125-135``).

GLUE we own: pull the per-sample labels out of the trailing columns, keyed by
``sampleID``, so the user gets a clean sample-level annotation table to attach to the
integrated object — not a 200-column wide CSV to wrangle.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

_LABEL_COLS = ("CST", "subCST", "score")


def import_valencia(csv_path, sample_col: str = "sampleID") -> pd.DataFrame:
    """Read a VALENCIA output CSV and return its CST/subCST/score, keyed by sample.

    Parameters
    ----------
    csv_path:
        Path to VALENCIA's one wide output CSV.
    sample_col:
        Name of the sample-id column (default ``"sampleID"``).

    Returns
    -------
    pandas.DataFrame
        Indexed by sample (as ``str``), columns ``CST``, ``subCST``, ``score``. A label
        column absent from the file is filled with NaN rather than raising, so partial
        VALENCIA outputs still import. Intended to be attached as sample-level annotation
        (``.obs``) on the integrated MuData.
    """
    d = pd.read_csv(csv_path)
    if sample_col not in d.columns:
        raise KeyError(
            f"Sample column {sample_col!r} not found in {csv_path}. "
            f"Columns start with: {list(d.columns[:5])}..."
        )
    out = pd.DataFrame(index=pd.Index(d[sample_col].astype(str), name="sample"))
    for col in _LABEL_COLS:
        out[col] = d[col].to_numpy() if col in d.columns else np.nan
    return out
