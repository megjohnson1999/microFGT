"""Resolve "which variables, on which samples" — the spec side of the analysis contract.

Every verb needs the same thing first: a clean sample-level frame holding the requested
predictor columns, restricted to a subset, with missing rows honestly dropped and counted.
This is where "the user picks what to look at" (principle #5: free-form metadata, roles
assigned at analysis time) becomes concrete.

Metadata may live in the global ``MuData.obs`` (CST, descriptors, clinical passed to
``build_mudata``) OR only in a modality's own ``obs`` (e.g. a phyloseq's ``sample_data``,
which rides on ``composition.obs``). We resolve against a *merged* view so a predictor is
found wherever it sits — the analyst shouldn't have to know.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def is_mudata(data) -> bool:
    return hasattr(data, "mod")


def pick_modality(data, modality=None, prefer=("composition_taxon", "composition")) -> str | None:
    """Choose which assay to compute on. ``None`` means the input is already an AnnData."""
    if not is_mudata(data):
        return None
    if modality is not None:
        if modality not in data.mod:
            raise KeyError(f"No modality {modality!r}; have {list(data.mod)}.")
        return modality
    for m in prefer:
        if m in data.mod:
            return m
    return list(data.mod)[0]


def get_assay(data, modality=None, prefer=("composition_taxon", "composition")):
    """Return ``(adata, modality_name)`` for the chosen assay (name is ``None`` for AnnData)."""
    mod = pick_modality(data, modality, prefer)
    return (data[mod], mod) if mod is not None else (data, None)


def merged_obs(data, modality=None) -> pd.DataFrame:
    """A sample-level view where a predictor is found wherever it sits, named once.

    Clinical/sample metadata is sample-level, but it may physically ride on any assay's obs
    (e.g. a phyloseq's ``sample_data`` lands on ``composition.obs``) rather than the global
    ``MuData.obs``. So we union global obs with *every* modality's obs — global first, then
    modalities in order — and the first definition of a column name wins.

    MuData stores a modality's obs in the global frame under a ``<modality>:`` prefix (this
    happens automatically on write/read), so the same variable shows up as both
    ``composition:HIV_status`` (global) and ``HIV_status`` (modality). That prefix is a storage
    artifact, not a distinct variable: we strip it and collapse the two to one clean name.
    The ``modality`` argument is accepted for symmetry but does not restrict what is visible.
    """
    frames = []
    prefixes: tuple[str, ...] = ()
    if is_mudata(data):
        frames.append(data.obs)
        prefixes = tuple(f"{m}:" for m in data.mod)
        for mod in data.mod:
            frames.append(data[mod].obs)
    else:
        frames.append(data.obs)
    out = pd.DataFrame(index=pd.Index(frames[0].index.astype(str), name="sample"))
    for fr in frames:
        fr = fr.copy()
        fr.index = fr.index.astype(str)
        for col in fr.columns:
            clean = _strip_modality_prefix(str(col), prefixes)
            if clean not in out.columns:                    # first definition wins
                out[clean] = fr[col].reindex(out.index)
    return out


def _strip_modality_prefix(col: str, prefixes: tuple[str, ...]) -> str:
    """Drop a leading ``<modality>:`` (MuData's storage prefix) to get the clean variable name."""
    for p in prefixes:
        if col.startswith(p):
            return col[len(p):]
    return col


def _apply_subset(obs: pd.DataFrame, subset) -> pd.DataFrame:
    """Restrict rows: a pandas query string, or a dict ``{col: value | [values]}``."""
    if subset is None:
        return obs
    if isinstance(subset, str):
        return obs.query(subset)
    if isinstance(subset, dict):
        mask = pd.Series(True, index=obs.index)
        for col, val in subset.items():
            if col not in obs.columns:
                raise KeyError(f"subset column {col!r} not in obs.")
            allowed = val if isinstance(val, (list, tuple, set)) else [val]
            mask &= obs[col].isin(list(allowed))
        return obs[mask]
    raise TypeError("subset must be a query string or a {column: value(s)} dict.")


def _missing_mask(frame: pd.DataFrame, columns) -> pd.Series:
    """True where any required column is missing (NaN, or the string sentinels 'nan'/'')."""
    mask = pd.Series(False, index=frame.index)
    for col in columns:
        s = frame[col]
        m = s.isna()
        if s.dtype == object or str(s.dtype).startswith("string"):
            m = m | s.astype("string").str.strip().isin(["", "nan", "NA", "None", "<NA>"])
        mask = mask | m
    return mask


def analysis_frame(data, columns, *, modality=None, subset=None, dropna=True):
    """Build the sample-level frame for ``columns``, subset + missing-dropped.

    Returns ``(frame, notes)`` where ``notes`` records ``n_total`` (after subset), ``n_used``,
    and ``n_dropped`` (rows removed for a missing value in any required column).
    """
    columns = list(dict.fromkeys(columns))                  # de-dup, keep order
    obs = merged_obs(data, modality)
    missing = [c for c in columns if c not in obs.columns]
    if missing:
        raise KeyError(
            f"Variable(s) {missing} not found in obs. Available: {list(obs.columns)}"
        )
    obs = _apply_subset(obs, subset)
    frame = obs[columns].copy()
    n_total = len(frame)
    if dropna and columns:
        keep = ~_missing_mask(frame, columns)
        frame = frame[keep]
    return frame, {"n_total": n_total, "n_used": len(frame), "n_dropped": n_total - len(frame)}


def is_categorical(series: pd.Series) -> bool:
    """A predictor is categorical unless it is a numeric dtype (then it's continuous)."""
    if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
        return False
    return True
