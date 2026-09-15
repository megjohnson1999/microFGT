# microFGT — Implementation Plan & Design Spec

**Read this first, then `design/fgt_workflow_map.md` (the reasoning behind it) and
`prototype/real_fixtures/FORMATS.md` (the authoritative tool-output formats).**
This document is the spec for building microFGT. It is deliberately tight — do not
expand scope, do not generate scaffolding/docs/tests beyond what's specified.

---

## What microFGT is

A flexible, tunable, **one-stop** tool for **female genital tract (FGT) microbiome
analysis** — general tool development for the field, usable on anyone's FGT data.
**Not** glue for any one cohort or project.

It either **runs** the standard FGT tools or **reimplements their concepts**, and lets
method choices across the *whole* workflow be explored and improved — vs. the existing
tools, which are rigid, hard to install, and frustrating to use. (Direction set by
Scott, 2026-06-16: "be open to other methods" at every decision point, not just CST.)

## Two hard constraints (the architecture must satisfy both)

### A. User-friendly — TOP PRIORITY
- **One-command install** (pip/conda + optional container). This is the single biggest
  thing people are frustrated by today. It is a feature, not a footnote.
- **Turnkey CLI that runs the whole workflow**, with a clean Python API underneath for
  power users. Design CLI/config-first; API second.
- **Sensible defaults that just work**, every method choice overridable.
- **Robust input handling** — above all taxa-name reconciliation (VALENCIA's #1 pain
  point). The tool fixes names; the user shouldn't have to.
- **The tool owns the glue** — joins, reshaping, annotation joins happen inside, never
  in the user's scripts. Helpful, specific error messages.

### B. Scientifically sound
- **Validate every method against a ground-truth reference before trusting it.** Centroid
  CST must reproduce VALENCIA (target 99.9% on the published data — we have it). Importers
  must reproduce real tool output shapes. Any alternative method is always diffed against
  the centroid baseline on the same data.
- **Compositional correctness** — microbiome counts are compositional. Defaults use
  relative abundance / CLR; never naive stats on raw counts.
- **Reproducibility is intrinsic** — pin tool versions, record parameters, deterministic
  runs. (Also answers the recurring "track what I did" need.)
- **Honest sample reconciliation** across assays — no silent dropping or double-counting.

## Language & container (decided)
- **Python.** Rationale: easiest install story (the top constraint); richest
  clustering/ML/stats ecosystem for CST method experimentation; VALENCIA is already
  Python so the centroid port is near-copy; the heavy tools are language-neutral
  subprocess calls. Language was evaluated open — Python won on the two constraints above.
- **Integrated object = MuData** (`muon`) — the Python multi-assay container, sample-keyed
  assays in one object. (R/MultiAssayExperiment was the alternative; rejected only because
  install + method-experimentation favor Python, not because of the container.)

---

## Architecture — layered, each concern swappable

1. **Integrated object (core).** A MuData holding sample-keyed assays:
   `composition` (taxon×sample), `function` (gene×sample); CST and clinical/sample
   variables as sample-level annotations. This is the currency every layer reads/writes.
2. **Import layer.** `import_speciateit` (ASV→taxon joined to the ASV count table →
   taxon×sample), `import_virgo` (stack per-sample `<sample>.out` files → gene×sample),
   `import_valencia` (one wide CSV → CST/subCST/score). Map **real** tool formats per
   `FORMATS.md` — do not guess formats; the validated prototype importers + fixtures are
   the starting point.
3. **Orchestration layer (optional).** Wrappers that actually *run* speciateIT / VIRGO /
   etc. via subprocess. Optional on purpose — a user can ingest existing outputs OR have
   microFGT run the tools.
4. **CST method layer.** A `classify_cst(composition, method=...)` interface. **Centroid
   is the first plugin** (reimplemented from `Valencia.py`, validated). Alternatives drop
   in behind the same interface — the seam exists from day one because we know ≥2
   implementations are coming.
5. **Analysis layer.** Normalization/transforms, diversity, differential abundance.
   **Build vs. buy:** the FGT-specific value is layers 2–4; the commodity steps
   (diversity, CLR/transforms, ordination) **lean on existing libraries** — reimplementing
   them would reinvent validated wheels and *hurt* soundness. Differential abundance:
   use Python ports of ANCOM-BC/ALDEx2, or orchestrate R for just that step. Do not build
   these natively.
6. **Viz + interface.** Python API; turnkey CLI over the top (see constraint A).

---

## Per-tool decisions (from the workflow map; grounded in the fixtures)

- **Taxonomic assignment — speciateIT:** *orchestrate.* Trained C++ classifier, not worth
  reinventing. Value added = wrap install + own the ASV→sample join (sample identity lives
  in the count table, not speciateIT's output).
- **CST — VALENCIA:** *reimplement the concept, openly.* The one rich method-development
  target. Nearest-centroid to 13 fixed reference centroids → forces every sample to a
  nearest type, no novelty detection, no calibrated/soft assignment, brittle taxa-name
  matching. v1 = faithful centroid (the baseline + yardstick); later = robust name
  reconciliation, then model-based clustering / classifiers / soft assignment.
- **Function — VIRGO:** *orchestrate.* Catalog-backed; not worth reinventing. Value added =
  wrap the run + own the stack-into-matrix + catalog annotation joins.

---

## Phased build order (walking skeleton first)

- **P0 — Skeleton + install.** Package layout, one-command install (pip/conda), minimal CI.
  Because UX is the priority, the install story is a P0 deliverable, not deferred.
- **P1 — Integrated object + importers.** Build a real MuData from the fixtures.
  *Milestone: real tool outputs go in, integrated object comes out.* Importers validated
  against `FORMATS.md` shapes.
- **P2 — Centroid CST behind the interface.** Reimplement from `Valencia.py`; validate
  ≥99.9% match vs VALENCIA on the published data. *Milestone: validated baseline + the
  swappable seam exist.*
- **P3 — Orchestration wrappers.** Actually run speciateIT/VIRGO; the "run the tools" half.
  - *Fast-follow (logged 2026-06-18):* a **Snakemake/Slurm cluster-execution layer** so the
    wrappers run as batch jobs on HTCF (WashU's cluster), not just local subprocesses —
    finishes the real batch story. (Per-sample VIRGO/speciateIT runs are embarrassingly
    parallel.) The first real HTCF run is also the trigger to discharge the speciateIT
    real-output validation IOU.
- **P4 — Analysis + viz + turnkey CLI.** Diversity, transforms, diff. abundance (bought,
  not built); the end-to-end CLI.
- **P5 — Alternative CST methods.** Compared against the centroid baseline on the same
  data. *The research payoff.*

Each phase is independently demoable. Validation gate at P1 (importer fidelity) and P2
(centroid reproduces VALENCIA) before moving on.

---

## Guardrails (explicit do-NOTs)
- **Do not carry forward the old R code.** `R/`, `R-new/`, the FGTExperiment S4 class, and
  the TESTING_*/INTEGRATION_*/ENHANCEMENT_* doc sprawl are dead weight. Build clean.
- **Do not reimplement commodity stats** (diversity, ordination, diff. abundance) — buy them.
- **Do not hardcode CST to centroid** — it goes behind the `classify_cst` interface from
  the start.
- **Do not scope to any one cohort.** General tool.
- **Do not guess tool formats.** Read `FORMATS.md` and the fixtures.

## Carry-forward assets (real, validated — do not regenerate)
- `prototype/real_fixtures/` — genuine VIRGO/VALENCIA/speciateIT outputs + `FORMATS.md`.
- `prototype/importers.R` — validated importer logic to **port to Python** (logic, not language).
- Validation target: `Valencia.py` reproduces the paper's CST for 99.9% of 13,231 samples.

## Open questions for implementation
- Exact CLI command surface (`microfgt run` / `import` / `classify` …) — design in P0/P4.
- Whether orchestration wrappers ship per-tool installers or assume tools on PATH.
- Diff-abundance: Python port vs. R orchestration — decide when P4 lands.
- Real co-assayed public dataset for end-to-end validation beyond fixtures (e.g. a public
  FGT cohort with both 16S and metagenomics).
