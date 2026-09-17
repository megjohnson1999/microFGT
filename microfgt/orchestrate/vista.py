"""Orchestrate VISTA — the mgCST classifier (random forest + YC-θ), grounded in RECIPE.md.

Two gotchas the recipe pins down:
- **arg2 is the dir CONTAINING ``VISTA_data/``** (the ``vista_repo``), not ``VISTA_data``
  itself; ``run_VISTA.R`` prepends ``/VISTA_data/...`` internally. The README's example fails.
- **VISTA writes to the current working directory** (timestamped filenames), with no ``-o``
  flag — so we run it with ``cwd=`` a dedicated output dir and then glob the result.

The runner produces VISTA's six files; the authoritative call (``mgCSTs_*.csv``) is handed to
:func:`microfgt.io.import_mgcst`. This is the concrete method behind the ``classify_mgcst`` seam.
"""

from __future__ import annotations

from pathlib import Path

from microfgt.orchestrate._run import resolve_executable, run_command


def run_vista(
    compiled, vista_repo, outdir, *, rscript: str = "Rscript", timeout: float | None = None,
):
    """Run ``run_VISTA.R`` on a VIRGO2 compiled matrix -> VISTA's ``mgCSTs_*.csv``.

    Parameters
    ----------
    compiled:
        Path to ``VIRGO2_Compiled.summary.NR.txt``.
    vista_repo:
        The VISTA repo dir (contains ``run_VISTA.R`` and ``VISTA_data/``).
    outdir:
        Directory VISTA runs in and writes its (timestamped) outputs to.

    Returns
    -------
    tuple[pathlib.Path, RunRecord]
        Path to the ``mgCSTs_*.csv`` call file and the run provenance.
    """
    rs_exe, rs_fp = resolve_executable(rscript, tool="Rscript (VISTA)")
    vista_repo = Path(vista_repo)
    script = vista_repo / "run_VISTA.R"
    if not script.exists():
        raise FileNotFoundError(
            f"run_VISTA.R not found at {script}. Set metagenomics.vista_repo to the VISTA "
            "repo (contains run_VISTA.R and VISTA_data/)."
        )
    compiled = Path(compiled).resolve()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # arg2 = the dir CONTAINING VISTA_data/ (absolute); VISTA writes to CWD, so cwd=outdir.
    argv = [rs_exe, str(script.resolve()), str(compiled), str(vista_repo.resolve())]
    record = run_command(
        argv, tool="VISTA", cwd=outdir,
        params={"compiled": str(compiled), "vista_repo": str(vista_repo)},
        exe_fingerprint=rs_fp, timeout=timeout,
    )
    hits = sorted(outdir.glob("mgCSTs_*.csv"))
    if not hits:
        raise FileNotFoundError(
            f"VISTA finished (rc={record.returncode}) but no mgCSTs_*.csv appeared in {outdir}. "
            f"stderr tail:\n{record.stderr_tail}"
        )
    return hits[-1], record


def _write_compiled_from_function(function, dest) -> Path:
    """Write a ``function`` (gene x sample) AnnData back to the VIRGO2 compiled matrix format
    (``Gene\\t<sample>…``) so VISTA can consume it when only the AnnData is in hand."""
    import numpy as np
    import pandas as pd

    X = function.layers["counts"] if "counts" in function.layers else function.X
    df = pd.DataFrame(
        np.asarray(X).T, index=function.var_names.astype(str), columns=function.obs_names.astype(str)
    )
    df.index.name = "Gene"
    df.to_csv(dest, sep="\t")
    return Path(dest)


def _restore_sample_ids(df, compiled):
    """Undo VISTA's positional ``V1..Vn`` relabeling of sample columns.

    ``run_VISTA.R`` loses the real sample ids and R relabels the compiled matrix's columns
    positionally as ``V1, V2, …`` (``V1`` = the ``Gene`` column, ``V2`` = the first sample, …).
    Map them back by column position, using the compiled matrix's own header. ``Vk`` anchors to
    the k-th column regardless of VISTA dropping/reordering low-count samples. Indices that are
    not ``V<k>`` (e.g. a future VISTA that preserves ids) are left untouched.
    """
    import warnings

    import pandas as pd

    try:
        cols = pd.read_csv(compiled, sep="\t", nrows=0).columns.tolist()
    except Exception as e:  # pragma: no cover - defensive; falls back to VISTA's raw index
        warnings.warn(f"could not read compiled matrix header to restore sample ids: {e}", stacklevel=2)
        return df
    remap = {f"V{k + 1}": str(cols[k]) for k in range(len(cols))}  # V1->Gene, V2->sample1, …
    unmapped = [str(i) for i in df.index if str(i).startswith("V") and str(i) not in remap]
    if unmapped:
        warnings.warn(
            f"VISTA emitted V-labels not in the compiled matrix's columns ({unmapped[:3]}); "
            "left unmapped — check the VIRGO2 compiled matrix and VISTA output align.",
            stacklevel=2,
        )
    out = df.copy()
    out.index = pd.Index([remap.get(str(i), str(i)) for i in df.index], name=df.index.name)
    return out


def classify_mgcst_vista(
    function=None, *, vista_repo, outdir, compiled=None,
    rscript: str = "Rscript", timeout: float | None = None, return_record: bool = False,
):
    """The VISTA method behind :func:`microfgt.mgcst.classify_mgcst`.

    Runs VISTA and parses the call. Pass ``compiled=`` (the VIRGO2 matrix path, the stage path)
    or a ``function`` AnnData (materialized to a compiled matrix first, the Python-API path).
    Returns the sample-keyed mgCST frame from :func:`microfgt.io.import_mgcst`.

    ``return_record=True`` returns ``(frame, RunRecord)`` instead, so the stage layer can record
    the VISTA invocation as provenance (like the other shotgun tool stages). The seam keeps the
    bare-frame default so :func:`microfgt.mgcst.classify_mgcst` stays DataFrame-returning.
    """
    from microfgt.io import import_mgcst

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if compiled is None:
        if function is None:
            raise ValueError(
                "classify_mgcst(method='vista') needs either compiled=<VIRGO2 matrix> or a "
                "function AnnData to derive it from."
            )
        compiled = _write_compiled_from_function(
            function, outdir / "VIRGO2_Compiled.summary.NR.txt"
        )
    mgcsts_csv, record = run_vista(compiled, vista_repo, outdir, rscript=rscript, timeout=timeout)
    df = _restore_sample_ids(import_mgcst(mgcsts_csv), compiled)
    return (df, record) if return_record else df
