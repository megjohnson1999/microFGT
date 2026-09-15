"""Orchestrate speciateIT — run the classifier, then own the ASV->sample join.

Command (grounded in the speciateIT README, not guessed):
    classify -d <vSpeciateDB dir> -i <fasta> -o <outDir>   [--skip-err-thld]
Output is always ``<outDir>/MC_order7_results.txt``. speciateIT classifies sequences, not
samples, so we hand the result to :func:`microfgt.io.import_speciateit` together with the
ASV count table (and the input FASTA, so the sequences are retained) to produce the
ASV x sample ``composition`` modality.

Note: a successful real run here also discharges the P1 real-output validation IOU on
``import_speciateit`` — run on ``test.fasta``, then validate the importer against the
genuine ``MC_order7_results.txt`` instead of the synthetic fixture.
"""

from __future__ import annotations

from pathlib import Path

from microfgt.io import import_speciateit
from microfgt.orchestrate._run import RunRecord, resolve_executable, run_command

RESULTS_FILENAME = "MC_order7_results.txt"


def run_speciateit(
    fasta,
    db,
    outdir,
    *,
    executable: str = "classify",
    skip_err_thld: bool = False,
    count_table=None,
    timeout: float | None = None,
    **import_kwargs,
):
    """Run speciateIT ``classify`` on a FASTA of ASVs.

    Parameters
    ----------
    fasta:
        Input FASTA (one record per ASV; headers are the ASV ids).
    db:
        Path to the vSpeciateDB model directory (e.g. ``vSpeciateIT_V3V4``).
    outdir:
        Output directory; ``classify`` writes ``MC_order7_results.txt`` here.
    executable:
        ``classify`` binary name (on PATH) or explicit path.
    skip_err_thld:
        Pass ``--skip-err-thld`` to force species-level annotation.
    count_table:
        If given (ASV count table CSV), the genuine results are imported immediately and an
        AnnData is returned with the run provenance in ``uns['speciateit_run']``. Otherwise
        the path to ``MC_order7_results.txt`` and the :class:`RunRecord` are returned.

    Returns
    -------
    anndata.AnnData (if ``count_table`` given) | tuple[pathlib.Path, RunRecord]
    """
    exe, fingerprint = resolve_executable(executable, tool="speciateIT")
    fasta, db, outdir = Path(fasta), Path(db), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    argv = [exe, "-d", str(db), "-i", str(fasta), "-o", str(outdir)]
    if skip_err_thld:
        argv.append("--skip-err-thld")

    record = run_command(
        argv,
        tool="speciateIT",
        params={"db": str(db), "fasta": str(fasta), "skip_err_thld": skip_err_thld},
        exe_fingerprint=fingerprint,
        timeout=timeout,
    )

    results_path = outdir / RESULTS_FILENAME
    if not results_path.exists():
        raise FileNotFoundError(
            f"speciateIT finished (rc={record.returncode}) but {results_path} was not "
            f"produced. stderr tail:\n{record.stderr_tail}"
        )

    if count_table is None:
        return results_path, record

    adata = import_speciateit(results_path, count_table, fasta=fasta, **import_kwargs)
    adata.uns["speciateit_run"] = record.to_dict()
    return adata
