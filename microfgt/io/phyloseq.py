"""phyloseq importer — the phyloseq entry point of the multi-entry design.

A phyloseq ``.rds`` is the standard 16S container and the shape the real FGT datasets
arrive in (a dada2 run wrapped up with its taxonomy and, often, an already-computed CST).
``import_phyloseq`` reads one and produces the **same ASV-grain ``composition`` AnnData**
that :func:`microfgt.io.import_speciateit` produces, so ``collapse_to_taxon`` -> CST ->
descriptors all just run downstream — the caller enters the workflow *late*, trusting the
phyloseq's existing taxonomy and CST rather than re-running speciateIT/VALENCIA.

phyloseq is R/Bioconductor, so (as with dada2) we don't reinvent it: microFGT owns one glue
asset, ``scripts/phyloseq_export.R``, which flattens the object's four slots to CSVs via a
subprocess ``Rscript`` call, and this module reshapes those into the AnnData.

Grain + slot mapping (validated against the FRESH ``ps_with_cst`` object):

* ``otu_table``   -> ``X`` / ``layers['counts']`` (samples x ASVs), ASV ids -> ``var_names``
* ``taxa_names``  -> ``var['sequence']`` (**the sequences ARE the taxa names; no refseq slot**)
* ``tax_table``   -> ``var['classification']`` from ``Genus_Species`` (fallback ``Species`` /
  ``Genus``) — the label :func:`collapse_to_taxon` groups on — plus ``var['genus']``
* ``sample_data`` -> ``obs`` (kept whole); if it carries ``CST`` / ``subCST`` / ``score`` they
  are surfaced by :func:`existing_cst` for ``build_mudata(cst=...)`` to attach.

v1 trusts the phyloseq's existing taxonomy + CST. Re-running speciateIT from the sequences is
a deliberately later option, not this path.
"""

from __future__ import annotations

import re
import tempfile
from contextlib import nullcontext
from importlib import resources
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from microfgt.io.speciateit import _genus_of

_LABEL_COLS = ("CST", "subCST", "score")
_MISSING = {"", "NA", "nan", "None", "<NA>"}
# GTDB-style rank prefixes: a single rank letter (kingdom/domain/phylum/…/species) followed by
# one or two underscores, e.g. ``g_Lactobacillus`` / ``d__Bacteria``. A lone rank letter + ``_``
# never starts a real binomial (``Lactobacillus_iners``) or a ``Ca_`` (Candidatus) name, so this
# strips the toolchain prefix without touching genuine labels.
_GTDB_PREFIX = re.compile(r"^[kdpcofgs]__?")


def _bundled_script() -> str:
    return str(resources.files("microfgt.scripts").joinpath("phyloseq_export.R"))


def _clean(series: pd.Series) -> pd.Series:
    """A tax-rank column as strings, normalised to a single usable taxon label.

    Strips GTDB rank prefixes (``g_Lactobacillus`` -> ``Lactobacillus``), strips a GTDB
    accession suffix (``Lactobacillus_iners(RS_GCF_000160875_1`` -> ``Lactobacillus_iners``),
    and collapses a doubled-genus join (``Lactobacillus Lactobacillus_iners`` -> the FRESH-style
    ``Genus_Species`` column repeats the genus -> ``Lactobacillus_iners``). Blanks / R ``NA``
    sentinels become <NA>; so does any value that *still* contains a space after those repairs
    (a GTDB placeholder like ``Bacteria Domain``), so it falls through to a usable rank rather
    than becoming a fake taxon. Normal binomials and ``Ca_`` names are untouched."""
    s = series.astype("string").str.replace(_GTDB_PREFIX, "", regex=True)
    s = s.str.replace(r"\(.*$", "", regex=True)                       # drop GTDB accession suffix
    s = s.str.replace(r"^(\S+)\s+(\1(?:_.*)?)$", r"\2", regex=True)    # collapse doubled genus
    s = s.str.strip()
    return s.mask(s.isin(_MISSING) | s.str.contains(" ", na=False))


def _classification_from_tax(tax: pd.DataFrame) -> pd.Series:
    """Per-ASV classification: ``Genus_Species``, falling back to ``Species`` then ``Genus``.

    Case-insensitive on rank names. Missing everywhere -> <NA> (that ASV rolls up into the
    ``Unclassified`` bucket downstream, honestly, rather than being dropped)."""
    lut = {str(c).lower(): c for c in tax.columns}

    def rank(name: str) -> pd.Series:
        if name in lut:
            return _clean(tax[lut[name]])
        return pd.Series(pd.NA, index=tax.index, dtype="string")

    out = rank("genus_species")
    out = out.fillna(rank("species"))
    out = out.fillna(rank("genus"))
    return out


def import_phyloseq(
    rds_path,
    *,
    rscript: str = "Rscript",
    script=None,
    export_dir=None,
    timeout: float | None = None,
) -> ad.AnnData:
    """Read a phyloseq ``.rds`` -> ASV-grain ``composition`` AnnData (import_speciateit's shape).

    Parameters
    ----------
    rds_path:
        Path to a phyloseq ``.rds`` (an ``otu_table`` + ``tax_table`` + ``sample_data``; the
        taxa names are the ASV sequences, as dada2 leaves them).
    rscript:
        ``Rscript`` executable — a bare name (looked up on PATH) or an explicit path. microFGT
        does not bundle R/phyloseq.
    script:
        Override the bundled ``phyloseq_export.R`` (testing / patched glue).
    export_dir:
        If given, the flat CSV exports are written here and kept (handy for inspection /
        presentation); otherwise a temp dir is used and cleaned up.
    timeout:
        Seconds to allow the R subprocess.

    Returns
    -------
    anndata.AnnData
        ``obs`` = the whole ``sample_data`` (plus ``read_count``, the per-sample total).
        ``var`` = ASVs, with ``classification`` (from ``Genus_Species`` / ``Species`` /
        ``Genus``), ``genus`` (first token of it), and ``sequence`` (the taxa name). ``X`` /
        ``layers['counts']`` hold ASV counts. This is the ``composition`` modality — roll it
        up with :func:`microfgt.io.collapse_to_taxon`. If ``sample_data`` carried CST columns,
        pull them out with :func:`existing_cst` for ``build_mudata(cst=...)``.
    """
    # Deferred import: the io layer must not pull in the orchestrate package at load time
    # (orchestrate imports back from io -> circular).
    from microfgt.orchestrate._run import resolve_executable, run_command

    exe, fingerprint = resolve_executable(rscript, tool="R (Rscript)")
    script = script or _bundled_script()

    ctx = nullcontext(str(export_dir)) if export_dir is not None else tempfile.TemporaryDirectory()
    with ctx as outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        argv = [exe, str(script), "--rds", str(rds_path), "--outdir", str(outdir)]
        record = run_command(
            argv, tool="phyloseq_export",
            params={"rds": str(rds_path)},
            exe_fingerprint=fingerprint, timeout=timeout,
        )

        out = Path(outdir)
        for name in ("counts.csv", "taxa_names.csv", "tax_table.csv", "sample_data.csv"):
            if not (out / name).exists():
                raise FileNotFoundError(
                    f"phyloseq_export finished (rc={record.returncode}) but {name} was not "
                    f"produced. stderr tail:\n{record.stderr_tail}"
                )
        counts = pd.read_csv(out / "counts.csv", index_col=0)      # samples x ASV ids
        taxa = pd.read_csv(out / "taxa_names.csv", index_col=0)     # ASV id -> sequence
        tax = pd.read_csv(out / "tax_table.csv", index_col=0)       # ASV id -> ranks
        sdata = pd.read_csv(out / "sample_data.csv", index_col=0)   # sample -> variables

    counts.index = counts.index.astype(str)
    asvs = [str(c) for c in counts.columns]
    taxa.index = taxa.index.astype(str)
    tax.index = tax.index.astype(str)

    classification = _classification_from_tax(tax).reindex(asvs)
    sequence = taxa["sequence"].astype("string").reindex(asvs) if "sequence" in taxa else pd.Series(pd.NA, index=asvs)
    genus = [_genus_of(c) if isinstance(c, str) else None for c in classification]
    var = pd.DataFrame(
        {
            "classification": [c if isinstance(c, str) else None for c in classification],
            "genus": genus,
            "sequence": [s if isinstance(s, str) else None for s in sequence],
        },
        index=pd.Index(asvs, name="asv"),
    )

    obs = sdata.copy()
    obs.index = obs.index.astype(str)
    obs = obs.reindex(counts.index)
    obs.index.name = "sample"
    counts_i = counts.to_numpy().astype(np.int64)
    obs.insert(0, "read_count", counts_i.sum(axis=1))

    adata = ad.AnnData(X=counts_i.astype(np.float32), obs=obs, var=var)
    adata.layers["counts"] = counts_i
    adata.uns["phyloseq_run"] = record.to_dict()
    adata.uns["has_cst"] = bool([c for c in _LABEL_COLS if c in sdata.columns])
    return adata


def existing_cst(composition: ad.AnnData) -> pd.DataFrame | None:
    """Pull the phyloseq's existing CST out of ``obs`` as a sample-keyed table.

    Mirrors :func:`microfgt.io.import_valencia`'s output (index = sample, columns
    ``CST`` / ``subCST`` / ``score``) so it drops straight into ``build_mudata(cst=...)``.
    Returns ``None`` if the imported ``sample_data`` carried none of those columns (nothing to
    attach — classify from the composition with :func:`microfgt.cst.classify_cst` instead).
    """
    present = [c for c in _LABEL_COLS if c in composition.obs.columns]
    if not present:
        return None
    out = pd.DataFrame(index=pd.Index(composition.obs_names.astype(str), name="sample"))
    for col in _LABEL_COLS:
        out[col] = composition.obs[col].to_numpy() if col in composition.obs.columns else np.nan
    return out
