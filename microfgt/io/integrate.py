"""Assemble per-modality assays into one sample-keyed MuData (layer 1, the core object).

A MuData holds sample-keyed assays in one object: ``composition`` (**ASV x sample**, from
speciateIT — the source of truth, carrying sequences), its materialised taxon roll-up
``composition_taxon`` (taxon x sample, what CST and the descriptors read), and ``function``
(gene x sample, from VIRGO), with CST + augment descriptors + clinical variables as
sample-level annotations. This is the currency every later layer reads/writes.

**Honest sample reconciliation** (constraint B): we never silently drop or double-count
samples across assays. MuData keeps the *union* of sample ids; ``build_mudata`` reports
how the assays and the CST table overlap so mismatches surface instead of hiding.

A caveat the fixtures make concrete: ``prototype/real_fixtures/`` carries outputs from
three *different* tool repos whose sample ids do NOT match (``sample1...`` vs ``sub1/sub2``
vs ``1,2,...``) — they are not a co-assayed cohort. So these fixtures validate each
importer and exercise the assembly/reconciliation mechanics, but a true end-to-end
cross-modality integration needs a real co-assayed dataset (an open question in the spec).
"""

from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import mudata as md
import pandas as pd

from microfgt.characterize import describe_composition
from microfgt.io.speciateit import collapse_to_taxon


@dataclass
class Reconciliation:
    """How sample ids overlapped across the assembled assays + CST table."""

    n_samples: int          # size of the union
    per_modality: dict      # modality -> n samples
    shared_all: int         # samples present in every modality
    cst_matched: int        # CST rows that matched a sample in the union
    cst_unmatched: list     # CST sample ids with no assay

    def summary(self) -> str:
        mods = ", ".join(f"{k}={v}" for k, v in self.per_modality.items())
        s = (
            f"{self.n_samples} samples in union ({mods}); "
            f"{self.shared_all} shared across all assays"
        )
        if self.per_modality:
            s += f"; CST matched {self.cst_matched}, unmatched {len(self.cst_unmatched)}"
        return s


def attach_cst_annotations(mdata: md.MuData, cst: pd.DataFrame) -> None:
    """Attach CST results to ``mdata``, keeping the sample annotation frame clean.

    The CST *label* columns (``CST`` / ``subCST`` / ``score`` and any user columns) go onto
    the global ``.obs`` by sample id. The 13 ``<subCST>_sim`` per-centroid similarity vectors
    are routed to ``composition_taxon.obsm['cst_sim']`` (with the column names in
    ``.uns['cst_sim_columns']``) — retrievable on the object, but not cluttering ``.obs``.
    If there is no ``composition_taxon`` modality, the sims are dropped (they are recomputable
    from :func:`microfgt.cst.classify_cst`).
    """
    cst = cst.copy()
    cst.index = cst.index.astype(str)
    union = list(mdata.obs_names)

    sim_cols = [c for c in cst.columns if str(c).endswith("_sim")]
    label_cols = [c for c in cst.columns if c not in sim_cols]

    aligned = cst.reindex(union)
    for col in label_cols:
        mdata.obs[col] = aligned[col].to_numpy()

    if sim_cols and "composition_taxon" in mdata.mod:
        taxon = mdata.mod["composition_taxon"]
        sims = cst.reindex(taxon.obs_names.astype(str))[sim_cols]
        taxon.obsm["cst_sim"] = sims.to_numpy(dtype=float)
        taxon.uns["cst_sim_columns"] = [str(c) for c in sim_cols]


def attach_mgcst_annotations(mdata: md.MuData, mgcst: pd.DataFrame) -> None:
    """Attach mgCST (shotgun community-type) results to ``mdata``'s global ``.obs``.

    The mgCST columns (``mgCST`` / ``mgCST_score`` and any user columns) go onto the global
    ``.obs`` by sample id, **beside** — never merged with — the 16S ``CST`` columns (they are not
    number-comparable). Unlike CST there is no per-centroid similarity block to route to
    ``.obsm``: VISTA emits only the best-match θ (``max_YC_theta``), not θ against all centroids
    (see :func:`microfgt.io.import_mgcst`), so this is a plain label attach.
    """
    mgcst = mgcst.copy()
    mgcst.index = mgcst.index.astype(str)
    aligned = mgcst.reindex(list(mdata.obs_names))
    for col in mgcst.columns:
        mdata.obs[col] = aligned[col].to_numpy()


def build_mudata(
    composition: ad.AnnData | None = None,
    function: ad.AnnData | None = None,
    cst: pd.DataFrame | None = None,
    obs: pd.DataFrame | None = None,
    *,
    composition_taxon: ad.AnnData | None = None,
    composition_taxon_shotgun: ad.AnnData | None = None,
    mgcst: pd.DataFrame | None = None,
    descriptors: bool = True,
) -> md.MuData:
    """Build a MuData from the available modalities and attach sample-level annotations.

    Parameters
    ----------
    composition:
        **ASV x sample** AnnData from :func:`microfgt.io.import_speciateit` (source of
        truth, carrying sequences).
    function:
        gene x sample AnnData from :func:`microfgt.io.import_virgo`.
    cst:
        Sample-keyed CST/subCST/score from :func:`microfgt.io.import_valencia` or
        :func:`microfgt.cst.classify_cst`. Label columns go to ``.obs``; any
        ``<subCST>_sim`` similarity vectors are routed to
        ``composition_taxon.obsm['cst_sim']`` (see :func:`attach_cst_annotations`).
    obs:
        Optional extra sample-level (clinical) annotations, indexed by sample id.
    composition_taxon:
        Taxon x sample roll-up of ``composition``. If omitted and ``composition`` is
        ASV-grain, it is materialised automatically via :func:`microfgt.io.collapse_to_taxon`
        (the taxon view is always present).
    composition_taxon_shotgun:
        Shotgun taxon x sample profile (from :func:`microfgt.io.import_virgo2_taxonomy`),
        mounted as its own modality parallel to the 16S ``composition_taxon``. Retained
        faithfully and never merged with the 16S taxonomy.
    mgcst:
        Sample-keyed mgCST + score (θ) from :func:`microfgt.io.import_mgcst` or
        :func:`microfgt.mgcst.classify_mgcst`. Its columns go to ``.obs`` beside — not merged
        with — 16S ``CST`` (see :func:`attach_mgcst_annotations`).
    descriptors:
        If True (default), compute the intrinsic augment descriptors (dominant taxon,
        % dominant, effective # taxa) from the taxon roll-up and attach them to ``.obs``.
        When a shotgun taxon modality is present its descriptors are attached too, under a
        ``shotgun_`` prefix.

    Returns
    -------
    mudata.MuData
        Modalities under ``composition`` / ``composition_taxon`` / ``function`` keys; CST,
        descriptors, and clinical variables joined onto the global ``.obs`` by sample id. A
        :class:`Reconciliation` is stored in ``mdata.uns["reconciliation"]`` (and its
        one-line summary in ``...["reconciliation_summary"]``).
    """
    mods: dict[str, ad.AnnData] = {}
    if composition is not None:
        mods["composition"] = composition
        if composition_taxon is None and "classification" in composition.var:
            composition_taxon = collapse_to_taxon(composition)
    if composition_taxon is not None:
        mods["composition_taxon"] = composition_taxon
    if composition_taxon_shotgun is not None:
        mods["composition_taxon_shotgun"] = composition_taxon_shotgun
    if function is not None:
        mods["function"] = function
    if not mods:
        raise ValueError(
            "build_mudata needs at least one modality "
            "(composition / composition_taxon_shotgun / function)."
        )

    # Adopt mudata's forthcoming default: don't auto-pull modality obs/var into the
    # global frame. We attach sample annotations explicitly below, so we don't rely on it.
    with md.set_options(pull_on_update=False):
        mdata = md.MuData(mods)

    # Reconciliation report over the union of sample ids (MuData's global obs).
    union = list(mdata.obs_names)
    union_set = set(union)
    per_modality = {k: a.n_obs for k, a in mods.items()}
    shared_all = len(set.intersection(*(set(a.obs_names) for a in mods.values())))

    cst_matched, cst_unmatched = 0, []
    if cst is not None:
        cst_idx = cst.index.astype(str)
        cst_matched = len([s for s in cst_idx if s in union_set])
        cst_unmatched = [s for s in cst_idx if s not in union_set]
        # Labels -> global .obs; <subCST>_sim vectors -> composition_taxon.obsm['cst_sim'].
        attach_cst_annotations(mdata, cst)

    # mgCST (shotgun) rides beside CST, never merged: labels -> global .obs,
    # mgcst_*_sim -> composition_taxon_shotgun.obsm['mgcst_sim'].
    if mgcst is not None:
        attach_mgcst_annotations(mdata, mgcst)

    if obs is not None:
        obs = obs.copy()
        obs.index = obs.index.astype(str)
        aligned = obs.reindex(union)
        for col in aligned.columns:
            mdata.obs[col] = aligned[col].to_numpy()

    # Augment descriptors: deterministic, always-wanted summaries of the taxon roll-up,
    # attached as plain .obs columns beside CST (they augment it, never replace it).
    if descriptors and composition_taxon is not None:
        desc = describe_composition(composition_taxon).reindex(union)
        for col in desc.columns:
            mdata.obs[col] = desc[col].to_numpy()

    # Shotgun descriptors: the same intrinsic summaries over the shotgun taxon profile,
    # under a shotgun_ prefix so they sit beside (never collide with) the 16S descriptors.
    if descriptors and composition_taxon_shotgun is not None:
        sdesc = describe_composition(composition_taxon_shotgun).reindex(union)
        for col in sdesc.columns:
            mdata.obs[f"shotgun_{col}"] = sdesc[col].to_numpy()

    recon = Reconciliation(
        n_samples=len(union),
        per_modality=per_modality,
        shared_all=shared_all,
        cst_matched=cst_matched,
        cst_unmatched=cst_unmatched,
    )
    mdata.uns["reconciliation_summary"] = recon.summary()
    mdata.uns["reconciliation"] = recon.__dict__
    return mdata
