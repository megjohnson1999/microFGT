"""The 16S stage registry — the one place the workflow DAG is defined.

Both executors and ``microfgt check`` consume this. Stage run functions reuse the existing
importers / CST / analysis / orchestration code, so there is no logic duplicated between the
in-memory API and the file-artifact pipeline — only the orchestration differs.

Metagenomic (QC/host-removal -> VIRGO) front-end mirrors this and is the next pass.
"""

from __future__ import annotations

from microfgt.stages.model import Requirement, Stage, StageContext, artifact_paths

# Region-aware DADA2 defaults, keyed to the vSpeciateDB model regions (the model dirs are named
# vSpeciateIT_<REGION>: V1V3 / V1V9 / V3V4 / V4V4). Starting points only — truncation genuinely
# depends on run quality, so these are overridable and the quality profile is emitted to inform
# overrides. (trim_left trims primers when DADA2 does it; leave it 0 if cutadapt already has.)
REGION_DEFAULTS = {
    "V1V3": {"trunc_len": [0, 0], "trim_left": [0, 0]},
    "V3V4": {"trunc_len": [280, 230], "trim_left": [17, 21]},
    "V4V4": {"trunc_len": [230, 180], "trim_left": [19, 20]},   # 515F/806R
    # V1V9 is full-length 16S (~1.5 kb): typically long-read (PacBio) and entered at an ASV
    # table, not Illumina paired DADA2 — so no useful paired-truncation default.
    "V1V9": {"trunc_len": [0, 0], "trim_left": [0, 0]},
}
# Colloquial region name -> the vSpeciateDB dir name, so a config may say the common "V4" or the
# exact model-region "V4V4" and get the same defaults.
_REGION_ALIASES = {"V4": "V4V4"}


def region_defaults(region: str | None) -> dict:
    """DADA2 defaults for a region, accepting colloquial aliases (e.g. ``"V4"`` -> ``"V4V4"``)."""
    return REGION_DEFAULTS.get(_REGION_ALIASES.get(region, region), {})


def canonical_region(region: str) -> str:
    """Normalise a region string to its canonical vSpeciateDB model-region name.

    Accepts the four model-region names (``V1V3 | V1V9 | V3V4 | V4V4``) case-insensitively
    plus the colloquial aliases (e.g. ``"V4"`` -> ``"V4V4"``). Raises ``ValueError`` for
    anything else, so ``microfgt setup`` never fetches a model for a region that isn't real.
    """
    if region is None:
        raise ValueError("no region given (expected one of V1V3 / V1V9 / V3V4 / V4V4)")
    token = _REGION_ALIASES.get(region.upper(), region.upper())
    if token not in REGION_DEFAULTS:
        raise ValueError(
            f"unknown region {region!r} — expected one of "
            f"{' / '.join(sorted(REGION_DEFAULTS))} (or the alias "
            f"{' / '.join(sorted(_REGION_ALIASES))})"
        )
    return token


# --- config -> provided artifacts (the entry point is just which inputs are present) -------
def provided_artifacts(config: dict) -> dict[str, str]:
    provided: dict[str, str] = {}
    comp = config.get("composition") or {}
    if "reads" in comp and comp["reads"].get("fastq_dir"):
        provided["fastq_dir"] = comp["reads"]["fastq_dir"]
    if comp.get("asv_table"):
        provided["asv_table"] = comp["asv_table"]
    if comp.get("asv_seqs"):
        provided["asv_seqs"] = comp["asv_seqs"]
    if comp.get("composition_h5ad"):
        provided["composition"] = comp["composition_h5ad"]
    sit = comp.get("speciateit") or {}
    if sit.get("results"):
        provided["speciateit_results"] = sit["results"]
    if sit.get("count_table"):
        provided["asv_table"] = sit["count_table"]
    cst_cfg = config.get("cst") or {}
    if comp.get("phyloseq"):
        provided["phyloseq_rds"] = comp["phyloseq"]
        # A phyloseq entry supplies its own CST by default; an explicit cst source
        # (method / valencia) overrides it, so only offer phyloseq as the CST producer
        # when neither is set.
        if not cst_cfg.get("method") and not cst_cfg.get("valencia"):
            provided["phyloseq_cst"] = comp["phyloseq"]
    if cst_cfg.get("valencia"):
        provided["valencia_output"] = cst_cfg["valencia"]

    # --- metagenomics (shotgun) arm entry points ---
    mg = config.get("metagenomics") or {}
    reads = mg.get("reads") or {}
    if reads.get("fastq_dir"):
        provided["sg_reads"] = reads["fastq_dir"]
    if mg.get("compiled"):
        provided["sg_compiled"] = mg["compiled"]
    if (config.get("mgcst") or {}).get("vista_output"):
        provided["vista_output"] = config["mgcst"]["vista_output"]
    return provided


# --- stage run functions (reuse existing code) ---------------------------------------------
def _run_primer_trim(ctx: StageContext) -> None:
    from microfgt.orchestrate.cutadapt import run_cutadapt

    reads = (ctx.config.get("composition") or {}).get("reads") or {}
    primers = reads.get("primers") or {}
    primers_configured = bool(primers.get("fwd") or primers.get("rev"))
    # Default to dropping primer-less (off-target) reads WHEN primers are configured — the Run 1
    # shakedown found ~16-27% of reads had no primer yet all flowed into DADA2. Only safe with a
    # primer set (cutadapt with no adapter would discard everything), and overridable via config.
    discard_untrimmed = reads.get("discard_untrimmed", primers_configured)
    records = run_cutadapt(
        ctx.path("fastq_dir"), ctx.path("trimmed_reads"),
        fwd_primer=primers.get("fwd"), rev_primer=primers.get("rev"),
        executable=reads.get("cutadapt", "cutadapt"),
        discard_untrimmed=discard_untrimmed,
        extra_args=reads.get("cutadapt_args") or None,
    )
    _write_provenance(ctx, "primer_trim", records)


def _run_denoise(ctx: StageContext) -> None:
    from microfgt.orchestrate.dada2 import run_dada2

    reads = (ctx.config.get("composition") or {}).get("reads") or {}
    defaults = region_defaults(reads.get("region"))
    d = reads.get("dada2") or {}
    # trim_left contract: cutadapt (the primer_trim stage) already removes the primers WHEN
    # they're configured, so applying the region-default trim_left on top would double-trim real
    # sequence. Key the default on whether primers were configured — not on stage presence: a
    # cutadapt run with no primers set is a no-op, and there trim_left is the legitimate primer
    # strip. An explicit dada2.trim_left always wins.
    primers = reads.get("primers") or {}
    primers_configured = bool(primers.get("fwd") or primers.get("rev"))
    default_trim_left = [0, 0] if primers_configured else defaults.get("trim_left")
    record = run_dada2(
        ctx.path("trimmed_reads"), ctx.path("asv_table"), ctx.path("asv_seqs"),
        ctx.path("quality_profile"),
        trunc_len=d.get("trunc_len", defaults.get("trunc_len")),
        trim_left=d.get("trim_left", default_trim_left),
        rscript=reads.get("rscript", "Rscript"),
    )
    _write_provenance(ctx, "denoise", [record])


def _run_assign(ctx: StageContext) -> None:
    from microfgt.orchestrate import run_speciateit

    sit = (ctx.config.get("composition") or {}).get("speciateit") or {}
    _, record = run_speciateit(
        fasta=ctx.path("asv_seqs"), db=sit["db"],
        outdir=ctx.path("speciateit_results").parent,
        executable=sit.get("classify", "classify"),
        skip_err_thld=sit.get("skip_err_thld", False),
    )
    _write_provenance(ctx, "assign", [record])


def _run_import_composition(ctx: StageContext) -> None:
    from microfgt.io import import_speciateit

    # Carry ASV sequences when the FASTA is available (dada2-emitted or user-provided).
    asv_seqs = ctx.path("asv_seqs")
    # The vSpeciateDB dir (model.tree) drives authoritative species->genus resolution.
    db = ((ctx.config.get("composition") or {}).get("speciateit") or {}).get("db")
    adata = import_speciateit(
        ctx.path("speciateit_results"), ctx.path("asv_table"),
        fasta=asv_seqs if asv_seqs.exists() else None, db=db,
    )
    adata.write(ctx.path("composition"))


def _run_import_phyloseq(ctx: StageContext) -> None:
    from microfgt.io import import_phyloseq

    rscript = (ctx.config.get("composition") or {}).get("phyloseq_rscript", "Rscript")
    adata = import_phyloseq(ctx.path("phyloseq_rds"), rscript=rscript)
    adata.write(ctx.path("composition"))


def _run_cst_phyloseq(ctx: StageContext) -> None:
    import anndata as ad

    from microfgt.io import existing_cst

    # The composition h5ad written by import_phyloseq already carries CST in its obs;
    # pull it out rather than re-reading the .rds.
    adata = ad.read_h5ad(ctx.path("composition"))
    cst = existing_cst(adata)
    if cst is None:
        raise ValueError(
            "cst source is the phyloseq object, but it carried no CST/subCST/score. "
            "Set cst.method to classify from the composition instead."
        )
    cst.to_csv(ctx.path("cst"))


def _run_cst_classify(ctx: StageContext) -> None:
    import anndata as ad

    from microfgt.cst import classify_cst
    from microfgt.io import collapse_to_taxon

    adata = ad.read_h5ad(ctx.path("composition"))
    # CST reads the taxon roll-up; collapse ASV-grain composition first.
    taxon = collapse_to_taxon(adata) if "classification" in adata.var else adata
    method = (ctx.config.get("cst") or {}).get("method", "centroid")
    classify_cst(taxon, method=method).to_csv(ctx.path("cst"))


def _run_cst_valencia(ctx: StageContext) -> None:
    from microfgt.io import import_valencia

    import_valencia(ctx.path("valencia_output")).to_csv(ctx.path("cst"))


def _run_integrate(ctx: StageContext) -> None:
    """Shared final assembler for all three integrate producers (16S / metagenomics / combined).

    Reads whichever modality artifacts are present in the workdir, derives the shotgun taxon
    roll-up from the function assay, folds in provenance, and assembles one MuData. The three
    Stage entries differ only in their declared inputs (so the resolver runs the right chains);
    the assembly logic is identical."""
    import anndata as ad
    import pandas as pd

    from microfgt.io import build_mudata, collapse_virgo2_to_taxon, import_virgo
    from microfgt.workflow import apply_analysis

    def _read_csv(key):
        df = pd.read_csv(ctx.path(key), index_col=0)
        df.index = df.index.astype(str)
        return df

    comp = ad.read_h5ad(ctx.path("composition")) if ctx.path("composition").exists() else None
    cst = _read_csv("cst") if ctx.path("cst").exists() else None
    mgcst = _read_csv("mgcst") if ctx.path("mgcst").exists() else None

    function = taxon_sg = None
    if ctx.path("function").exists():
        function = ad.read_h5ad(ctx.path("function"))
        if "taxon" in function.var:                      # derive the shotgun taxon modality
            taxon_sg = collapse_virgo2_to_taxon(function)
    else:                                                # backward-compat: v1 function from config
        fcfg = ctx.config.get("function") or {}
        if "virgo" in fcfg:
            function = import_virgo(fcfg["virgo"]["dir"])

    if comp is not None and ctx.config.get("analysis"):
        apply_analysis(comp, ctx.config["analysis"])

    mdata = build_mudata(
        composition=comp, function=function, cst=cst,
        composition_taxon_shotgun=taxon_sg, mgcst=mgcst,
    )
    prov = _collect_provenance(ctx.workdir)
    if prov:
        import json

        # One arm-agnostic key for every tool invocation (16S + shotgun), keyed by stage id;
        # each RunRecord carries its own `tool`. Store as a JSON string: the nested RunRecords
        # (argv lists, param dicts) are not h5-serializable as a raw nested dict, but a scalar
        # string round-trips cleanly.
        mdata.uns["tool_runs"] = json.dumps(prov)
    mdata.write(ctx.path("mudata"))


# --- metagenomics (shotgun) stage run functions (reuse the orchestrators) ------------------
def _mg(ctx: StageContext) -> dict:
    return ctx.config.get("metagenomics") or {}


def _mg_require(ctx: StageContext, key: str):
    val = _mg(ctx).get(key)
    if not val:
        raise ValueError(
            f"metagenomics.{key} is required for this stage — set it in the config "
            "(microfgt check reports it up front)."
        )
    return val


def _write_provenance(ctx: StageContext, stage_id: str, records) -> None:
    import json

    prov_dir = ctx.workdir / "provenance"
    prov_dir.mkdir(parents=True, exist_ok=True)
    payload = [r.to_dict() for r in records]
    (prov_dir / f"{stage_id}.json").write_text(json.dumps(payload, indent=2))


def _collect_provenance(workdir) -> dict:
    import json
    from pathlib import Path

    prov_dir = Path(workdir) / "provenance"
    if not prov_dir.is_dir():
        return {}
    return {p.stem: json.loads(p.read_text()) for p in sorted(prov_dir.glob("*.json"))}


def _run_sg_qc(ctx: StageContext) -> None:
    from microfgt.orchestrate import run_fastp

    records = run_fastp(
        ctx.path("sg_reads"), ctx.path("sg_trimmed"),
        threads=_mg(ctx).get("threads", 4), executable=_mg(ctx).get("fastp", "fastp"),
    )
    _write_provenance(ctx, "sg_qc", records)


def _run_sg_host_removal(ctx: StageContext) -> None:
    from microfgt.orchestrate import run_host_removal

    records = run_host_removal(
        ctx.path("sg_trimmed"), ctx.path("sg_nonhost"), _mg_require(ctx, "host_ref"),
        threads=_mg(ctx).get("threads", 4),
        minimap2=_mg(ctx).get("minimap2", "minimap2"),
        samtools=_mg(ctx).get("samtools", "samtools"),
    )
    _write_provenance(ctx, "sg_host_removal", records)


def _run_sg_virgo2_map(ctx: StageContext) -> None:
    from microfgt.orchestrate import run_virgo2_map
    from microfgt.orchestrate.cutadapt import discover_pairs

    virgo2_dir = _mg_require(ctx, "virgo2_dir")
    outdir = ctx.path("sg_virgo2_out")
    pairs = discover_pairs(ctx.path("sg_nonhost"))
    if not pairs:
        raise FileNotFoundError(f"No host-removed FASTQ pairs found in {ctx.path('sg_nonhost')}.")
    records = []
    for sample, r1, r2 in pairs:
        _, record = run_virgo2_map(
            r1, r2, sample, virgo2_dir, outdir,
            threads=_mg(ctx).get("threads", 4), python=_mg(ctx).get("python", "python3"),
        )
        records.append(record)
    _write_provenance(ctx, "sg_virgo2_map", records)


def _run_sg_virgo2_compile(ctx: StageContext) -> None:
    from microfgt.orchestrate import run_virgo2_compile

    _, record = run_virgo2_compile(
        ctx.path("sg_virgo2_out"), _mg_require(ctx, "virgo2_dir"),
        python=_mg(ctx).get("python", "python3"),
    )
    _write_provenance(ctx, "sg_virgo2_compile", [record])


def _run_import_function(ctx: StageContext) -> None:
    from pathlib import Path

    from microfgt.io import import_virgo2

    # Annotations are joined from the VIRGO2 install's AnnotationTables (taxon is used to derive
    # the shotgun taxon modality). virgo2_dir is optional for a bare compiled entry.
    taxon_ann = None
    virgo2_dir = _mg(ctx).get("virgo2_dir")
    annotations = dict(_mg(ctx).get("annotations", {}))
    if virgo2_dir:
        cand = Path(virgo2_dir) / "AnnotationTables" / "1.VIRGO2.taxon.txt"
        if cand.exists():
            taxon_ann = cand
    adata = import_virgo2(
        ctx.path("sg_compiled"), taxon_annotation=taxon_ann, annotations=annotations or None
    )
    adata.write(ctx.path("function"))


def _run_classify_mgcst(ctx: StageContext) -> None:
    from microfgt.orchestrate.vista import classify_mgcst_vista

    df, record = classify_mgcst_vista(
        compiled=ctx.path("sg_compiled"), vista_repo=_mg_require(ctx, "vista_repo"),
        outdir=ctx.workdir / "vista", rscript=_mg(ctx).get("rscript", "Rscript"),
        return_record=True,
    )
    df.to_csv(ctx.path("mgcst"))
    _write_provenance(ctx, "classify_mgcst", [record])


def _run_import_mgcst_existing(ctx: StageContext) -> None:
    from microfgt.io import import_mgcst

    import_mgcst(ctx.path("vista_output")).to_csv(ctx.path("mgcst"))


# --- requirements (entry-dependent; checked by the doctor) ---------------------------------
def _req_primer_trim(cfg):
    return [Requirement("binary", ((cfg.get("composition") or {}).get("reads") or {}).get("cutadapt", "cutadapt"),
                        "install cutadapt and put it on PATH")]


def _req_denoise(cfg):
    rscript = ((cfg.get("composition") or {}).get("reads") or {}).get("rscript", "Rscript")
    return [Requirement("binary", rscript, "install R"),
            Requirement("rpackage", "dada2", "install the Bioconductor dada2 package",
                        via=rscript)]


def _req_import_phyloseq(cfg):
    rscript = (cfg.get("composition") or {}).get("phyloseq_rscript", "Rscript")
    return [Requirement("binary", rscript, "install R"),
            Requirement("rpackage", "phyloseq", "install the Bioconductor phyloseq package",
                        via=rscript)]


def _req_assign(cfg):
    sit = (cfg.get("composition") or {}).get("speciateit") or {}
    reqs = [Requirement("binary", sit.get("classify", "classify"),
                        "install speciateIT (the classify binary)")]
    if sit.get("db"):
        reqs.append(Requirement("path", sit["db"], "path to the vSpeciateDB model directory"))
    return reqs


# --- metagenomics req_fns (the audit's walls -> actionable check errors; per RECIPE.md) ----
def _req_sg_qc(cfg):
    mg = cfg.get("metagenomics") or {}
    return [Requirement("binary", mg.get("fastp", "fastp"), "install fastp (QC/trim)")]


def _req_sg_host_removal(cfg):
    mg = cfg.get("metagenomics") or {}
    reqs = [
        Requirement("binary", mg.get("minimap2", "minimap2"), "install minimap2 (host removal)"),
        Requirement("binary", mg.get("samtools", "samtools"), "install samtools (host removal)"),
    ]
    if mg.get("host_ref"):
        reqs.append(Requirement("path", mg["host_ref"],
                                "host genome for removal (e.g. GRCh38.fna.gz)"))
    return reqs


def _req_sg_virgo2_map(cfg):
    from pathlib import Path

    mg = cfg.get("metagenomics") or {}
    reqs = [
        Requirement("binary", mg.get("python", "python3"), "python3 to run VIRGO2.py"),
        Requirement("binary", mg.get("bowtie2", "bowtie2"), "install bowtie2 (VIRGO2 mapping)"),
    ]
    if mg.get("virgo2_dir"):
        d = Path(mg["virgo2_dir"])
        reqs.append(Requirement("path", str(d / "VIRGO2.py"), "VIRGO2.py in the VIRGO2 install"))
        reqs.append(Requirement("path", str(d / "Index" / "VIRGO2.1.bt2"),
                                "the VIRGO2 bowtie2 index (build once with VIRGO2.py install)"))
        # Version-identity guard: the GitHub and Zenodo VIRGO2.py diverge (GitHub's has an
        # args.threads crash) and the failure is silent + un-debuggable for a newcomer. Pin the
        # known-good Zenodo checksum via metagenomics.virgo2_sha256 to make divergence a visible
        # MISS line rather than a mystery.
        if mg.get("virgo2_sha256"):
            reqs.append(Requirement(
                "checksum", str(d / "VIRGO2.py"),
                "VIRGO2.py does not match the pinned checksum — use the Zenodo VIRGO2.py "
                "(DOI 10.5281/zenodo.18703182), not the GitHub one",
                expected=mg["virgo2_sha256"],
            ))
    return reqs


def _req_sg_virgo2_compile(cfg):
    from pathlib import Path

    mg = cfg.get("metagenomics") or {}
    reqs = [Requirement("binary", mg.get("python", "python3"), "python3 to run VIRGO2.py")]
    if mg.get("virgo2_dir"):
        reqs.append(Requirement("path", str(Path(mg["virgo2_dir"]) / "VIRGO2.py"),
                                "VIRGO2.py in the VIRGO2 install"))
    return reqs


_VISTA_RPACKAGES = ("randomForestSRC", "pheatmap", "dplyr", "data.table", "R.utils")


def _req_classify_mgcst(cfg):
    from pathlib import Path

    mg = cfg.get("metagenomics") or {}
    rscript = mg.get("rscript", "Rscript")
    reqs = [Requirement("binary", rscript, "install R (VISTA classifier)")]
    reqs += [Requirement("rpackage", pkg, f"install the R package {pkg} (VISTA)", via=rscript)
             for pkg in _VISTA_RPACKAGES]
    if mg.get("vista_repo"):
        d = Path(mg["vista_repo"])
        reqs.append(Requirement("path", str(d / "run_VISTA.R"), "run_VISTA.R in the VISTA repo"))
        reqs.append(Requirement("path", str(d / "VISTA_data" / "volume"),
                                "VISTA_data/ bundle (fetch from figshare)"))
    return reqs


# --- the registry --------------------------------------------------------------------------
STAGES = [
    Stage("primer_trim", ("fastq_dir",), ("trimmed_reads",), _run_primer_trim, _req_primer_trim),
    Stage("denoise", ("trimmed_reads",), ("asv_table", "asv_seqs", "quality_profile"),
          _run_denoise, _req_denoise),
    Stage("assign", ("asv_seqs",), ("speciateit_results",), _run_assign, _req_assign),
    Stage("import_composition", ("speciateit_results", "asv_table"), ("composition",),
          _run_import_composition),
    Stage("import_phyloseq", ("phyloseq_rds",), ("composition",),
          _run_import_phyloseq, _req_import_phyloseq),
    Stage("cst_classify", ("composition",), ("cst",), _run_cst_classify),
    Stage("cst_valencia", ("valencia_output",), ("cst",), _run_cst_valencia),
    # phyloseq CST reads the composition h5ad's obs; phyloseq_cst is the selection gate
    # (present only when the phyloseq object is the chosen CST source).
    Stage("cst_phyloseq", ("composition", "phyloseq_cst"), ("cst",), _run_cst_phyloseq),
    # --- metagenomics (shotgun) arm: reads -> compiled -> {function, mgcst} ---
    Stage("sg_qc", ("sg_reads",), ("sg_trimmed",), _run_sg_qc, _req_sg_qc),
    Stage("sg_host_removal", ("sg_trimmed",), ("sg_nonhost",),
          _run_sg_host_removal, _req_sg_host_removal),
    Stage("sg_virgo2_map", ("sg_nonhost",), ("sg_virgo2_out",),
          _run_sg_virgo2_map, _req_sg_virgo2_map),
    Stage("sg_virgo2_compile", ("sg_virgo2_out",), ("sg_compiled",),
          _run_sg_virgo2_compile, _req_sg_virgo2_compile),
    Stage("import_function", ("sg_compiled",), ("function",), _run_import_function),
    # import-existing is registered BEFORE classify so a provided VISTA output wins over
    # re-running VISTA (both otherwise have a directly-provided input, sg_compiled).
    Stage("import_mgcst_existing", ("vista_output",), ("mgcst",), _run_import_mgcst_existing),
    Stage("classify_mgcst", ("sg_compiled",), ("mgcst",),
          _run_classify_mgcst, _req_classify_mgcst),
    # --- final assembly: three producers of `mudata`, one shared run fn. Registered most-
    # complete first so the resolver's "most-complete resolvable" tie-break routes combined runs.
    Stage("integrate_combined", ("composition", "cst", "function", "mgcst"), ("mudata",),
          _run_integrate),
    Stage("integrate", ("composition", "cst"), ("mudata",), _run_integrate),
    Stage("integrate_shotgun", ("function", "mgcst"), ("mudata",), _run_integrate),
]
STAGE_BY_ID = {s.id: s for s in STAGES}

# artifact key -> stages that can produce it (an artifact may have >1 producer, e.g. cst).
PRODUCERS: dict[str, list[Stage]] = {}
for _s in STAGES:
    for _out in _s.outputs:
        PRODUCERS.setdefault(_out, []).append(_s)


def execute_stage(stage_id: str, workdir, config: dict, output: str | None = None) -> None:
    """Run one stage by id (used by the local executor and the ``_run-stage`` CLI alike)."""
    from pathlib import Path

    stage = STAGE_BY_ID[stage_id]
    paths = artifact_paths(workdir, config, output)
    for out_key in stage.outputs:
        Path(paths[out_key]).parent.mkdir(parents=True, exist_ok=True)
    stage.run(StageContext(Path(workdir), paths, config))
