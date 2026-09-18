"""Sample sheet — the user-declared map from sample id to that sample's files per arm.

This is the source of truth for sample identity. Rather than the tool *guessing* a sample
id out of arbitrary FASTQ filenames (fragile across labs) or trying to auto-link the 16S and
shotgun arms (unsafe — a wrong guess silently mis-pairs patients), the user supplies a small
CSV that says, for each biological sample: its ``sample_id`` and where its reads live in each
arm. The tool reads and *honours* it. Any extra columns are per-sample metadata.

Format (one row per sample; a CSV)::

    sample_id, 16s_R1,          16s_R2,          shotgun_R1,      shotgun_R2,      group
    PT01,      16s/PT01_R1.fq,  16s/PT01_R2.fq,  mgx/PT01_R1.fq,  mgx/PT01_R2.fq,  BV
    PT02,      16s/PT02_R1.fq,  16s/PT02_R2.fq,  ,                ,                Normal

- ``sample_id`` is required, non-empty, and unique — it is the join key across arms.
- The read columns are optional and recognised case-insensitively. A sample may have one arm
  or both; leave a cell blank when that arm wasn't sequenced for that sample.
- **v1 supports paired-end 16S only.** The DADA2 denoise stage (``scripts/dada2_run.R``)
  requires both ``16s_R1`` and ``16s_R2`` and errors on unpaired input — single-end amplicon
  reads are not yet supported (tracked in ``design/refinement_backlog.md``). A blank
  ``16s_R2`` is accepted by the sheet but the 16S run will fail without it.
- Every other column is carried through as sample metadata (wired to the integrated object's
  ``.obs`` in a later step).

This module only *reads and validates* the sheet; consuming it in the pipeline run is a
separate step.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Recognised read-file columns, matched case-insensitively. R2 is optional (single-end -> R1
# only); each arm's presence for a sample is decided by whether it has any read file listed.
_ARMS = {
    "16s": ("16s_r1", "16s_r2"),
    "shotgun": ("shotgun_r1", "shotgun_r2"),
}
_READ_COLS = {c for cols in _ARMS.values() for c in cols}
ID_COL = "sample_id"


@dataclass
class SampleSheetReport:
    """Outcome of validating a sample sheet."""

    n_samples: int = 0
    n_16s: int = 0
    n_shotgun: int = 0
    n_both: int = 0
    problems: list[str] = field(default_factory=list)   # blocking errors
    warnings: list[str] = field(default_factory=list)   # non-blocking heads-ups

    @property
    def ok(self) -> bool:
        return not self.problems

    def summary(self) -> str:
        return (
            f"{self.n_samples} samples "
            f"(16S={self.n_16s}, shotgun={self.n_shotgun}, both={self.n_both}); "
            f"{len(self.problems)} problem(s), {len(self.warnings)} warning(s)"
        )


def _norm_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Lower/strip header names so recognition is case- and whitespace-insensitive."""
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return df


def read_samplesheet(path) -> pd.DataFrame:
    """Read a sample-sheet CSV into a DataFrame (headers normalised, blank cells -> NA).

    Does not validate — call :func:`validate_samplesheet` for that.
    """
    df = pd.read_csv(path, dtype=str)
    df = _norm_cols(df)
    # Treat empty strings / whitespace-only cells as missing.
    df = df.apply(lambda col: col.map(lambda v: v.strip() if isinstance(v, str) else v))
    df = df.replace("", pd.NA)
    return df


def _arm_present(row: pd.Series, arm: str) -> bool:
    """True if the row lists at least one read file for the given arm."""
    return any(col in row and pd.notna(row[col]) for col in _ARMS[arm])


def validate_samplesheet(df: pd.DataFrame, base_dir=None) -> SampleSheetReport:
    """Validate a sample sheet: unique ids, files exist, and report per-arm overlap.

    Parameters
    ----------
    df:
        A sample sheet as returned by :func:`read_samplesheet`.
    base_dir:
        Directory to resolve relative file paths against (defaults to the current dir).

    Returns
    -------
    SampleSheetReport
        ``.problems`` are blocking (empty => the sheet is usable); ``.warnings`` are
        non-blocking heads-ups (e.g. arms that will never integrate).
    """
    rpt = SampleSheetReport()
    base = Path(base_dir) if base_dir is not None else None

    if ID_COL not in df.columns:
        rpt.problems.append(f"missing required '{ID_COL}' column")
        return rpt

    ids = df[ID_COL]
    if ids.isna().any():
        rpt.problems.append(f"{int(ids.isna().sum())} row(s) have an empty {ID_COL}")
    dupes = [str(v) for v in ids[ids.duplicated(keep=False) & ids.notna()].unique()]
    if dupes:
        rpt.problems.append(f"duplicate {ID_COL}(s): {', '.join(sorted(dupes))}")

    if not (_READ_COLS & set(df.columns)):
        rpt.warnings.append(
            "no recognised read columns (16s_R1/16s_R2/shotgun_R1/shotgun_R2) — "
            "the sheet lists no files to run from"
        )

    rpt.n_samples = int(ids.notna().sum())
    n16 = nsg = nboth = 0
    missing_files: list[str] = []
    for _, row in df.iterrows():
        sid = row[ID_COL]
        if pd.isna(sid):
            continue
        has16, hassg = _arm_present(row, "16s"), _arm_present(row, "shotgun")
        n16 += has16
        nsg += hassg
        nboth += has16 and hassg
        # A listed R2 with no R1 is almost certainly a mistake.
        for arm, (r1, r2) in _ARMS.items():
            if r2 in row and pd.notna(row.get(r2)) and (r1 not in row or pd.isna(row.get(r1))):
                rpt.problems.append(f"{sid}: {arm} has {r2} but no {r1}")
        # Every listed file must exist.
        for col in _READ_COLS & set(df.columns):
            val = row.get(col)
            if pd.notna(val):
                p = Path(val)
                if base is not None and not p.is_absolute():
                    p = base / p
                if not p.exists():
                    missing_files.append(f"{sid}/{col}: {val}")
    rpt.n_16s, rpt.n_shotgun, rpt.n_both = n16, nsg, nboth

    if missing_files:
        shown = "; ".join(missing_files[:10])
        more = "" if len(missing_files) <= 10 else f" (+{len(missing_files) - 10} more)"
        rpt.problems.append(f"{len(missing_files)} missing file(s): {shown}{more}")

    # Declare-time mirror of the build-time cross-arm guard: if both arms are used but no
    # single sample carries both, they can never integrate.
    if n16 and nsg and nboth == 0:
        rpt.warnings.append(
            "both arms are present but no sample has files for both — nothing will integrate "
            "across arms; give a sample both its 16S and shotgun files under one sample_id"
        )

    return rpt


def _resolve(val, base: Path | None) -> Path:
    p = Path(val)
    if base is not None and not p.is_absolute():
        p = base / p
    return p


def _stage_arm(df: pd.DataFrame, arm: str, dest_dir: Path, base: Path | None) -> int:
    """Symlink each sample's reads for one arm into ``dest_dir`` under canonical
    ``{sample_id}_R1``/``_R2`` names, so the existing filename-based discovery picks up the
    sample_id. Returns the number of samples staged for this arm."""
    r1col, r2col = _ARMS[arm]
    n = 0
    for _, row in df.iterrows():
        sid = row[ID_COL]
        if pd.isna(sid) or not _arm_present(row, arm):
            continue
        dest_dir.mkdir(parents=True, exist_ok=True)
        for read, col in (("R1", r1col), ("R2", r2col)):
            val = row.get(col)
            if col in row and pd.notna(val):
                # Normalise the staged extension to .fastq[.gz] (what discovery globs), but keep
                # the real compression so downstream tools read it correctly.
                suffix = ".fastq.gz" if str(val).endswith(".gz") else ".fastq"
                link = dest_dir / f"{sid}_{read}{suffix}"
                if link.is_symlink() or link.exists():
                    link.unlink()
                link.symlink_to(_resolve(val, base).resolve())
        n += 1
    return n


def stage_samplesheet(config: dict, workdir) -> dict:
    """Honour a ``samples: <sheet.csv>`` config key by materialising the sheet into the run.

    Validates the sheet (raising on blocking problems), symlinks each sample's reads into
    per-arm dirs under canonical ``{sample_id}_R1/_R2`` names, points the FASTQ entry points at
    those dirs, and writes the sheet's metadata columns to ``<workdir>/samplesheet_obs.csv`` (the
    integrate stage attaches them to the object's ``.obs``). Returns the (possibly modified)
    config; a no-op when there is no ``samples`` key.
    """
    sheet = config.get("samples")
    if not sheet:
        return config

    sheet = Path(sheet)
    base = sheet.parent
    df = read_samplesheet(sheet)
    rpt = validate_samplesheet(df, base_dir=base)
    if not rpt.ok:
        raise ValueError("sample sheet has problems:\n  - " + "\n  - ".join(rpt.problems))

    workdir = Path(workdir)
    config = dict(config)

    if _stage_arm(df, "16s", workdir / "staged" / "16s", base):
        comp = dict(config.get("composition") or {})
        reads = dict(comp.get("reads") or {})
        reads["fastq_dir"] = str(workdir / "staged" / "16s")
        comp["reads"] = reads
        config["composition"] = comp

    if _stage_arm(df, "shotgun", workdir / "staged" / "shotgun", base):
        mg = dict(config.get("metagenomics") or {})
        reads = dict(mg.get("reads") or {})
        reads["fastq_dir"] = str(workdir / "staged" / "shotgun")
        mg["reads"] = reads
        config["metagenomics"] = mg

    meta_cols = [c for c in df.columns if c not in _READ_COLS and c != ID_COL]
    if meta_cols:
        obs = df.loc[df[ID_COL].notna(), [ID_COL] + meta_cols].set_index(ID_COL)
        obs.to_csv(workdir / "samplesheet_obs.csv")

    return config
