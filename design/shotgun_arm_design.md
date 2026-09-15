# microFGT — Shotgun (Metagenomic) Arm Design

*Working design doc for the shotgun/metagenomic arm. Complements `design_vision.md`
(the why), `fgt_workflow_map.md` (orchestrate-vs-reimplement), and `implementation_plan.md`
(build decisions). Derived from a hands-on end-to-end audit (2026-08-19/20): building and
running QC → host removal → VIRGO2 → VISTA from scratch on real data (ENA PRJEB34536 +
11 real FRESH shotgun metagenomes), so the decisions below are grounded in what actually
broke and what actually worked, not speculation. Last updated 2026-08-20.*

---

## In one sentence

The shotgun arm ingests metagenomic FASTQs and lands a **`function` (gene) modality**,
a **shotgun taxonomic profile**, and an **mgCST community-type call** into the same MuData
object the 16S arm builds — and its whole reason to exist is that **users should not have
to endure the setup ordeal this arm was reverse-engineered from.**

## Core principle (the thing to hold onto)

**Opinionated about software, unopinionated about data location. Ease = no cryptic
failures + clear guidance, not magic auto-install.**

- **Software dependencies → the tool owns them.** A conda/bioconda recipe installs microFGT
  + all external tools (bowtie2, samtools, fastp, minimap2, R + randomForestSRC + VISTA
  R deps). microFGT stays a **thin orchestrator** (`orchestrate/_run.py` finds tools on
  PATH/config; it does not bundle them). An **Apptainer** image is an *optional*
  reproducibility layer — **Docker is not required** (it's unavailable on most HPC).
- **Databases → configure + verify, do NOT provision or prescribe location.** Environments
  are heterogeneous and unpredictable, and labs frequently already have these DBs centrally
  (e.g. the audit lab had GRCh38 / VIRGO2 / kraken2 in a shared `/ref/...`). So: config
  points at DB paths (existing or user-downloaded); the `microfgt check` doctor verifies and
  emits actionable errors ("VIRGO2 index missing at `<path>` → get/build it like X"); docs
  list which DBs, where to get each, and the workarounds; a `fetch-references` helper is
  **optional**, for greenfield users, and must respect existing installs. Never force a
  provisioning mechanism.

## The shotgun pipeline (stages)

Add these to the stage DAG in `stages/registry.py`, parallel to the 16S
`primer_trim → denoise → assign → import_composition → cst → integrate` chain. The audit
produced a **working reference implementation** of each (see "Reference implementation" below).

| Stage id | Does | In → Out | Notes |
|---|---|---|---|
| `sg_qc` | fastp trim | fastqs → trimmed | standard |
| `sg_host_removal` | minimap2 vs GRCh38, keep both-mates-unmapped | trimmed → non-host | FGT is very host-heavy (audit: 18–94% host); host fraction is unpredictable per sample |
| `sg_virgo2_map` | bowtie2 vs VIRGO2 index, per sample | non-host → per-sample `.out` | **VIRGO2 is single-end**: concatenate R1+R2, map combined (undocumented convention) |
| `sg_virgo2_compile` | `VIRGO2.py compile` | `.out`s → `VIRGO2_Compiled.summary.NR.txt` | |
| `import_function` | compiled matrix → `function` modality | matrix → AnnData(gene × sample) | **`import_virgo2`** added *alongside* the kept v1 `import_virgo`. Reads the single wide `VIRGO2_Compiled.summary.NR.txt`; annotations are joined from separate `AnnotationTables/` files on `Gene` (see below). |
| `classify_mgcst` | VISTA random forest + YC-θ | function → mgCST + θ | **NEW method**, mirror `classify_cst(method=)` seam |
| `integrate` | fold modality into MuData | → `.h5mu` | existing stage; reads `config["function"][...]` |

Each stage's `req_fn` is where the audit's walls become preflight checks (GRCh38 present,
VIRGO2 index built, VISTA env + `VISTA_data` present, spike reference available).

## The integrated object — what the shotgun arm writes

The MuData modality model means the shotgun arm writes the **same self-contained set whether
or not 16S is present**; combination is the reconciliation layer's job, not the feeder's.

**Modalities (feature × sample AnnData):**
- **`function`** — gene × sample (VIRGO2 gene counts). The shotgun-unique payoff 16S can't
  produce. `.var` carries VIRGO2 gene annotations (taxon, KEGG/EC/CAZy/VOG/AMR/phage), each
  **joined from a separate `AnnotationTables/` file on `Gene`** (the compiled matrix itself has
  no annotation columns), tolerating genes a table doesn't cover; functional roll-ups
  (pathway × sample) derive from it.
- **shotgun taxon composition** (`composition_taxon_shotgun`) — taxon × sample, a shotgun view
  parallel to 16S's `composition_taxon`. **There is no `VIRGO2.py taxonomy` output**: microFGT
  *derives* this itself by joining the gene matrix → `AnnotationTables/1.VIRGO2.taxon.txt` on
  `Gene` and summing per taxon (`import_virgo2` → `collapse_virgo2_to_taxon`), exactly mirroring
  the 16S `import_speciateit` → `collapse_to_taxon` split. This deliberately decouples taxonomy
  from running the classifier.

**Global `.obs` (sample-level):** `mgCST`, `mgCST_score` (θ) — stored **separately** from 16S
`CST`; `virgo2_mapping_rate`, `genes_detected`, `host_fraction` (reference-fit); `absolute_load`
(from the spike, see below); shotgun `dominant_taxon` / `effective_taxa` descriptors
(`shotgun_`-prefixed, parallel to the 16S augment descriptors). **No `mgCST_subtype`**: VISTA's
call file (`mgCSTs_*.csv`) carries only the mgCST label + best-match `max_YC_theta`. The finer
mgSs level lives in `norm_counts_mgSs_mgCST_*.csv` as a **feature matrix** (mgSs × sample), not a
per-sample label — a possible future mgSs modality, not an `.obs` column.

**`.obsm`:** *(none for mgCST)* — VISTA emits only the best-match θ, not θ against all 25
centroids, so there is **no `mgcst_sim`** vector (this is where mgCST differs from CST's
`cst_sim`). A low `mgCST_score` is the reference-fit signal instead.

**`.uns`:** `reference_fit` report; provenance (`virgo2_runs` / `vista_runs`); participates
in `reconciliation`.

**Shotgun-only vs shotgun+16S:** the shotgun contribution is identical either way. Shotgun-only
= a complete valid object with just the shotgun modalities. Shotgun+16S = both stacks coexist,
`reconciliation` (union of samples) tracks who has which, and cross-modal analysis runs on the
**shared subset** via the `associate`/`compare` verbs.

## Decision: keep 16S and shotgun taxonomy SEPARATE (v1, probably beyond)

Do **not** merge 16S and shotgun taxonomy into one reconciled matrix. They measure different
things (different biases, resolutions, references); a merge destroys information and bakes in an
un-undoable, possibly-wrong opinion; and reconciliation is a research problem that would stall
v1. Instead:
- **Retain both faithfully** as separate modalities; **CST and mgCST as separate `.obs`
  columns** (they are not number-comparable).
- **Reconcile the *format*, not the content** — store taxonomy in a consistent schema (same
  `.var` rank columns / taxon-string convention) so the two are *comparable at any rank the
  user picks*, without being merged.
- **Reconciliation becomes an optional, transparent analysis verb** — a user who wants to
  compare calls `compare`/`associate` on the shared subset at a chosen rank and sees
  agreement/disagreement (often the *disagreement* is the interesting signal). It is never a
  hidden object-build step.
- **Document a default recommendation without enforcing it:** 16S for community/CST-level
  comparability with the field; shotgun for species/strain/genomespecies resolution + function.

## Reference decision: just use VIRGO2

VIRGO2 is the field standard; building a new FGT gene catalog would be silly (years of
multi-lab effort, curated annotations; a single-cohort catalog would be *more* biased). The
audit showed VIRGO2's apparently-low mapping on the South-African FRESH cohort (~55–80%, vs
~79% on a European cohort) is **mostly benign**, not catalog failure:
- **Spike-in** — the ZymoBIOMICS D6320 spike (*Imtechella* / *Allobacillus* halotolerans)
  is present in the reads and VIRGO2 correctly doesn't map it (it's a vaginal catalog). *Bonus:
  the spike is in the shotgun data → absolute microbial load is recoverable from shotgun, not
  just 16S → the `absolute_load` `.obs` field above.*
- **Non-coding / rRNA** — VIRGO2 is a **gene** catalog (coding only), so intergenic/rRNA reads
  of *well-represented* organisms can't map even though the organism is in the catalog. (A
  fair gene-vs-gene comparison: VMGC's gene catalog maps only ~3 points better than VIRGO2;
  the widely-quoted "83.8% vs 71.7%" gain is VMGC *genomes* vs VIRGO *genes* — a genome-vs-gene
  confound, not comprehensiveness.)
- **Real catalog gaps** = the small remainder.

Two consequences for the tool:
1. **Trust VIRGO2 for the catalog/gene layer.** "Low mapping" ≠ failure. If real gaps ever
   matter for a cohort, the answer is reference-free **assembly**, not a new catalog.
2. **Keep the reference-bias caution on the *classifier* (mgCST/VISTA), which is 100%
   North-American-built** (Holm 2023; the catalog VIRGO2 spans 5 continents). Treat a low-θ
   mgCST call as "this community doesn't match the reference types well," not a confident label.
3. **Report reference-fit, don't assert completeness** — surface mapping rate + θ per sample
   (the `reference_fit` machinery above) so a user *sees* whether the tools fit *their* cohort.

*(Myth to avoid: the "≈98.7% North American" figure is not a real sourced statistic — it
describes VIRGO v1 / mgCST, not VIRGO2, and is unverifiable. Don't cite it.)*

## Dependencies & databases (the real obstacle)

The audit confirmed: the entire difficulty is deps + DBs; the analysis itself is trivial once
set up. Deps are a *solved* class (conda); DBs are the *ongoing* one (large, externally hosted,
awkwardly distributed, version-drift-prone). The provisioning **recipe** (from the audit) that
docs + the optional `fetch-references` should encode:
- **GRCh38** (host removal reference) + a minimap2/bowtie2 index.
- **VIRGO2**: get it from **Zenodo (DOI 10.5281/zenodo.18703182), NOT git-lfs** (git-lfs is the
  primary documented route and it breaks on HPC → 133-byte pointer stubs). Zenodo `Index.tar`
  is an **empty placeholder** → build the bowtie2 index once (`VIRGO2.py install`, ~2.5 GB).
  Pin the **Zenodo `VIRGO2.py`**, not the GitHub one (GitHub's `install` has an `args.threads`
  crash; the two channels ship divergent versions).
- **VISTA**: clone repo + fetch `VISTA_data.tar.gz` from figshare
  (DOI 10.6084/m9.figshare.28684934). Gotcha: `run_VISTA.R` arg2 is the dir *containing*
  `VISTA_data/` (not `VISTA_data` itself), and it writes to CWD.
- **R env**: R ≥ 4.3 + randomForestSRC, pheatmap, dplyr, data.table, R.utils. (VISTA's Python/
  Streamlit deps are for its app, not the classifier — skip them.)

## Implementation to-dos (concrete)

**Increment 1 — DONE** (object/import core, commit `ad6c314`): `import_virgo2` (+ annotation
joins) alongside the kept v1 `import_virgo`; `collapse_virgo2_to_taxon` / `import_virgo2_taxonomy`
(derived shotgun taxon); `import_mgcst` (VISTA `mgCSTs_*.csv`); `build_mudata` generalized
(`composition_taxon_shotgun` + `mgcst` kwargs, `shotgun_` descriptors). Grounded in public
ENA/PRJEB34536 fixtures.

**Increment 2 — DONE** (the running front-end): the `sg_qc → sg_host_removal → sg_virgo2_map →
sg_virgo2_compile → import_function → classify_mgcst → integrate` chain in `stages/registry.py`,
so `microfgt run` takes raw reads → `.h5mu`. Orchestrators (`orchestrate/fastp|host_removal|virgo2|
vista.py`) grounded in `prototype/reference_scripts/RECIPE.md`; VISTA registered behind the
`classify_mgcst` seam; `req_fn`s surface every wall via `microfgt check`; resolver routes 16S /
metagenomics / combined via three `integrate` producers; config block is top-level `metagenomics:`.
Stub-tool e2e proves the plumbing; real-tool correctness is the HTCF-run IOU.

**Remaining:**
1. Reference-fit + spike-derived absolute load as first-class `.obs`/`.uns` outputs.
2. Snakemake/Slurm **resource** scaling (array-job mem/time from the RECIPE resource notes:
   host removal peak ~13 G; VIRGO2 map OOM'd at 32 G / cleared at 64 G on the deepest samples).
3. (Optional) an mgSs feature modality from `norm_counts_mgSs_mgCST_*.csv`.
4. Real-tool validation on HTCF (discharge the stub-only IOU).

## Reference implementation (the working audit pipeline)

Not in this repo, but the scripts that prove each stage work end-to-end live on HTCF:
`/scratch/sahlab/Megan/metagenomics_fresh/scripts/` (02_preproc, 03a_virgo2_map,
03b_virgo2_compile, 03c_vista — array-job versions) and `/scratch/sahlab/Megan/metagenomics_fgt/`
(the original single-run audit + the built VIRGO2 index, `envs/vista`, `resources/VISTA`).

## Open / deferred

- Cross-modal `compare`/`associate` verb behavior on the shared 16S∩shotgun subset (the
  "reconciliation as a verb" — spec the rank-selection + agreement reporting).
- Functional roll-ups (KEGG/pathway modules) as derived modalities vs `.var` annotations.
- Assembly / MAG arm for the genuinely-novel unmapped fraction (only if a cohort needs it).
- Validation dataset for cross-modal reconciliation: **PIN / PRJNA876771** (same-DNA co-assayed
  16S+shotgun, per `candidate_datasets.md`) — FRESH is *not* suitable (its shotgun and 16S are
  different collections, so no same-sample link exists).
