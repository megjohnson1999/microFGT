"""VIRGO / VIRGO2 importers — functional profiling (gene counts).

**VIRGO v1** (:func:`import_virgo`), REAL shape (FORMATS.md, validated against ``virgo_sub*.out``):
one file per sample, named ``<sample>.out``, TSV, NO header, 3 columns
``geneID \\t read_count \\t gene_length``. Sample identity is the filename.

**VIRGO2** (:func:`import_virgo2`), REAL shape (FORMATS.md, validated against
``virgo2_compiled.summary.NR.slice.txt``): ``VIRGO2.py compile`` writes exactly **one** wide
matrix, ``VIRGO2_Compiled.summary.NR.txt`` — TSV, header ``Gene\\t<sample>\\t<sample>…``,
one row per gene, float counts. Unlike v1 the compiled matrix carries **no annotation
columns**: taxon/KEGG/EC/… all live in separate ``AnnotationTables/`` files joined on
``Gene`` (never ``Cluster``). Shotgun taxon composition is therefore something microFGT
*derives* (join genes → taxon, sum per taxon), not a file VIRGO2 emits — mirroring the 16S
``import_speciateit`` → :func:`collapse_to_taxon` split.

GLUE we own (so the user doesn't): v1 stacks per-sample files; v2 orients the compiled matrix
to MuData convention (obs = samples, var = genes) and joins the requested annotation tables,
tolerating genes an annotation table doesn't cover.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

_COLS = ["gene_id", "read_count", "gene_length"]


def _genus_of(taxon: str) -> str:
    """First token of a VIRGO2 ``Taxa`` string (e.g. 'Lactobacillus_iners' -> 'Lactobacillus').

    Kept identical to the 16S convention (:func:`microfgt.io.speciateit._genus_of`) so the
    shotgun taxon assay is rank-comparable with the 16S ``composition_taxon`` (the doc's
    'reconcile the format, not the content')."""
    return taxon.split("_")[0].split(" ")[0]


def import_virgo(directory, pattern: str = "*.out", taxon_table=None) -> ad.AnnData:
    """Stack a directory of VIRGO per-sample ``.out`` files into an AnnData.

    Parameters
    ----------
    directory:
        Directory containing one ``<sample>.out`` file per sample.
    pattern:
        Glob for the per-sample files (default ``"*.out"``). The sample name is the
        filename with this suffix stripped.
    taxon_table:
        Optional path to VIRGO's ``1.taxon.tbl.txt`` (``Cluster\\tgeneID\\ttaxon\\tlength``).
        If given, a ``taxon`` column is added to ``var`` (catalog annotation join).

    Returns
    -------
    anndata.AnnData
        ``obs`` = samples, ``var`` = genes (with ``gene_length``; ``taxon`` if a catalog
        was supplied). ``X`` holds read counts; ``layers["counts"]`` is the integer copy.
    """
    directory = Path(directory)
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No VIRGO per-sample files matching {pattern!r} in {directory}. "
            "Expected one '<sample>.out' file per sample (3-column TSV, no header)."
        )

    suffix_len = len(pattern.lstrip("*"))  # "*.out" -> strip ".out"
    counts_by_sample: dict[str, pd.Series] = {}
    gene_length: dict[str, int] = {}
    for f in files:
        sample = f.name[:-suffix_len] if suffix_len else f.name
        df = pd.read_csv(f, sep="\t", header=None, names=_COLS)
        if sample in counts_by_sample:
            raise ValueError(f"Duplicate sample name {sample!r} from file {f}")
        counts_by_sample[sample] = df.set_index("gene_id")["read_count"]
        gene_length.update(zip(df["gene_id"], df["gene_length"]))

    # genes x samples, zero-filled; sorted for deterministic ordering.
    wide = pd.DataFrame(counts_by_sample).fillna(0).astype(np.int64).sort_index()
    samples = list(wide.columns)
    genes = list(wide.index)

    var = pd.DataFrame(
        {"gene_id": genes, "gene_length": [gene_length[g] for g in genes]},
        index=pd.Index(genes, name="gene"),
    )
    if taxon_table is not None:
        tax = pd.read_csv(
            taxon_table, sep="\t", header=None,
            names=["cluster", "gene_id", "taxon", "length"],
        )
        gene2taxon = dict(zip(tax["gene_id"], tax["taxon"]))
        var["taxon"] = [gene2taxon.get(g, "Unannotated") for g in genes]

    obs = pd.DataFrame(index=pd.Index(samples, name="sample"))
    counts = wide.T.to_numpy()  # samples x genes
    adata = ad.AnnData(X=counts.astype(np.float32), obs=obs, var=var)
    adata.layers["counts"] = counts
    return adata


UNANNOTATED = "Unannotated"


def _join_annotation(var: pd.DataFrame, table_path, *, key: str, prefix: str) -> None:
    """Left-join a VIRGO2 ``AnnotationTables/`` file onto ``var`` (keyed on gene id).

    Every non-key column of the table is added to ``var`` under ``<prefix>_<column>``, aligned
    to the gene index; genes the table does not cover get NaN (the join must tolerate partial
    annotation — e.g. the real KEGG table annotates ~2/3 of genes). Mutates ``var`` in place.
    """
    tab = pd.read_csv(table_path, sep="\t", dtype=str)
    if key not in tab.columns:
        raise ValueError(
            f"{table_path} has no {key!r} column (columns: {list(tab.columns)}); VIRGO2 "
            "annotation tables join on the gene id."
        )
    tab = tab.drop_duplicates(subset=key).set_index(key)
    for col in tab.columns:
        var[f"{prefix}_{col}"] = tab[col].reindex(var.index).to_numpy()


def import_virgo2(summary, *, taxon_annotation=None, annotations=None) -> ad.AnnData:
    """Import the VIRGO2 compiled gene matrix -> the **gene x sample** ``function`` modality.

    Parameters
    ----------
    summary:
        Path to ``VIRGO2_Compiled.summary.NR.txt`` (TSV, header ``Gene\\t<sample>…``, one row
        per gene, float counts).
    taxon_annotation:
        Optional path to ``AnnotationTables/1.VIRGO2.taxon.txt`` (cols ``Cluster, Gene, Taxa,
        Cat``). When given, ``var['taxon']`` (from ``Taxa``) is attached, joined on ``Gene``;
        genes the table does not cover get ``"Unannotated"``. This is what
        :func:`collapse_virgo2_to_taxon` reads to derive the shotgun taxon assay.
    annotations:
        Optional ``{name: path}`` of further ``AnnotationTables/`` files (e.g.
        ``{"kegg": ".../3.VIRGO2.kegg.txt"}``). Each table's non-``Gene`` columns are joined
        onto ``var`` prefixed with ``name`` (``kegg_KEGG``, …), tolerating missing genes.

    Returns
    -------
    anndata.AnnData
        ``obs`` = samples, ``var`` = genes (with ``taxon`` / annotation columns when supplied).
        ``X`` and ``layers['counts']`` hold the (fractional) counts. This is the ``function``
        modality VISTA/mgCST classifies and the shotgun taxon assay derives from.
    """
    summary = Path(summary)
    wide = pd.read_csv(summary, sep="\t", index_col=0)  # genes x samples
    wide.index = wide.index.astype(str)
    wide.columns = [str(c) for c in wide.columns]
    if wide.empty:
        raise ValueError(f"{summary} parsed to an empty matrix; expected 'Gene\\t<sample>…'.")

    genes = list(wide.index)
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    if taxon_annotation is not None:
        tax = pd.read_csv(taxon_annotation, sep="\t", dtype=str).drop_duplicates(subset="Gene")
        gene2taxon = dict(zip(tax["Gene"], tax["Taxa"]))
        var["taxon"] = [gene2taxon.get(g, UNANNOTATED) for g in genes]
    for name, path in (annotations or {}).items():
        _join_annotation(var, path, key="Gene", prefix=name)

    counts = wide.T.to_numpy(dtype=np.float64)  # samples x genes
    adata = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(index=pd.Index(wide.columns, name="sample")),
        var=var,
    )
    adata.layers["counts"] = counts
    return adata


def collapse_virgo2_to_taxon(function: ad.AnnData) -> ad.AnnData:
    """Roll a VIRGO2 ``function`` (gene x sample) assay up to a **taxon x sample** assay.

    Aggregates gene counts by ``var['taxon']`` (from :func:`import_virgo2` with a taxon
    annotation). Counts are conserved. This is the shotgun analogue of
    :func:`microfgt.io.collapse_to_taxon`, producing the ``composition_taxon_shotgun`` modality
    — derived by microFGT, not read from VISTA.
    """
    if "taxon" not in function.var:
        raise ValueError(
            "collapse_virgo2_to_taxon needs a function assay with var['taxon'] "
            "(import_virgo2 with taxon_annotation=...)."
        )
    X = function.layers["counts"] if "counts" in function.layers else function.X
    counts = pd.DataFrame(
        np.asarray(X, dtype=np.float64),
        index=function.obs_names.astype(str),
        columns=function.var_names.astype(str),
    )  # samples x genes

    labels = [str(t) if not pd.isna(t) else UNANNOTATED for t in function.var["taxon"]]
    grouped = counts.T.groupby(pd.Index(labels, name="taxon")).sum().T  # samples x taxa
    taxa = list(grouped.columns)

    var = pd.DataFrame(
        {"genus": [t if t == UNANNOTATED else _genus_of(t) for t in taxa]},
        index=pd.Index(taxa, name="taxon"),
    )
    obs = pd.DataFrame(
        {"read_count": grouped.sum(axis=1).to_numpy()},
        index=pd.Index(grouped.index, name="sample"),
    )
    taxon_counts = grouped.to_numpy(dtype=np.float64)
    adata = ad.AnnData(X=taxon_counts.astype(np.float32), obs=obs, var=var)
    adata.layers["counts"] = taxon_counts
    return adata


def import_virgo2_taxonomy(summary, taxon_annotation) -> ad.AnnData:
    """Derive the shotgun **taxon x sample** ``composition_taxon_shotgun`` modality.

    Convenience wiring of :func:`import_virgo2` (join genes -> taxon) then
    :func:`collapse_virgo2_to_taxon` (sum per taxon). ``summary`` is the compiled gene matrix;
    ``taxon_annotation`` is ``AnnotationTables/1.VIRGO2.taxon.txt``.
    """
    function = import_virgo2(summary, taxon_annotation=taxon_annotation)
    return collapse_virgo2_to_taxon(function)
