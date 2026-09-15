"""Config loading, in one place for both executors.

A YAML section written with no body — ``mgcst:`` with nothing under it — parses to ``None``,
not ``{}``. The key is *present*, so ``config.get("mgcst", {})`` returns that ``None`` and the
next ``.get(...)`` raises ``AttributeError``. Rather than police the ``.get`` idiom at every
call site, we normalize once on load: drop the ``None``-valued keys so an empty section reads
as an *absent* one everywhere downstream.

We drop rather than coerce to ``{}`` deliberately. Some sites use a section's mere presence as
a truthiness signal (an empty ``analysis:`` means "no analysis requested" -> falsy); coercing
to ``{}`` would flip that to truthy. Dropping preserves it: absent stays falsy.
"""

from __future__ import annotations

from pathlib import Path


def normalize_config(obj):
    """Recursively drop ``None``-valued keys from a parsed-config mapping."""
    if isinstance(obj, dict):
        return {k: normalize_config(v) for k, v in obj.items() if v is not None}
    if isinstance(obj, list):
        return [normalize_config(v) for v in obj]
    return obj


def load_config(path) -> dict:
    """Read a YAML config from ``path`` and normalize empty sections away."""
    import yaml

    raw = yaml.safe_load(Path(path).read_text()) or {}
    return normalize_config(raw)
