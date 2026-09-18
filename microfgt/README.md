# microFGT

A flexible, tunable, **one-stop** tool for **female genital tract (FGT) microbiome
analysis**. It runs the standard FGT tools (speciateIT, VALENCIA, VIRGO) or
reimplements their concepts, and lets method choices across the whole workflow be
explored — vs. the existing tools, which are rigid, hard to install, and frustrating
to use.

This is the Python rebuild. Two hard constraints drive every decision: **user-friendliness**
(one-command install, a turnkey CLI, and the tool owns all the glue — joins, reshaping,
annotation — so the user never scripts it) and **scientific soundness** (every method is
validated against a ground-truth reference before it's trusted; compositional correctness;
reproducible, provenance-recorded runs; honest sample reconciliation across assays).

## Install

> ⚠️ **Platform: Linux.** v1 targets **Linux**; **macOS is not supported** (the shotgun stack
> has no Apple-Silicon `osx-arm64` build, and speciateIT/conda break on spaces in macOS home
> paths). Run microFGT on a Linux host.

```bash
pip install -e ".[dev]"     # from a clone, with test deps
```

or with conda:

```bash
conda env create -f environment.yml
conda activate microfgt
pip install -e .
```

The package itself is pure Python and pip-installs in one command (the UX-first
constraint, exercised in CI on 3.10 + 3.12). The heavy external tools and their
reference data are **not** bundled — see [External tools & reference data](#external-tools--reference-data).

## What it does today

microFGT covers **two arms** — amplicon **16S** and shotgun **metagenomics** — as stage
graphs defined in exactly one place (`microfgt/stages/registry.py`) and consumed by two
executors (a local runner and a Snakemake/Slurm generator), so the laptop and cluster
pipelines can never drift apart:

- **16S:** `primer-trim → denoise (DADA2) → assign (speciateIT) → import → classify CST → integrate`
- **Metagenomics:** `fastp QC → host removal (minimap2) → VIRGO2 map → VIRGO2 compile → import → classify mgCST (VISTA) → integrate`

You don't pick a stage to start at. **The entry point is just whichever inputs your
config provides** — microFGT resolves the shortest path to a finished object:

| You have… | Enters at | Runs |
|---|---|---|
| **16S** — Paired FASTQs (`composition.reads.fastq_dir`) | the top | cutadapt → DADA2 → speciateIT → CST → analysis |
| **16S** — An ASV table + rep-seqs | `assign` | speciateIT → CST → analysis |
| **16S** — Existing speciateIT / VALENCIA output | `import` | CST → analysis |
| **Shotgun** — Paired FASTQs (`metagenomics.reads.fastq_dir`) | the top | fastp → host removal → VIRGO2 → VISTA/mgCST → analysis |
| **Shotgun** — A `VIRGO2_Compiled` matrix (`metagenomics.compiled`) | `import` | VISTA/mgCST → analysis |
| **Shotgun** — Existing VISTA output (`mgcst.vista_output`) | `import` | analysis |

The output is one `.h5mu` (MuData) holding sample-keyed assays — `composition`
(ASV×sample; the 16S source of truth, carrying each ASV's sequence), its taxon roll-up
`composition_taxon` (taxon×sample), `function` (VIRGO2 gene×sample) with its own taxon
roll-up `composition_taxon_shotgun`, and the `mgcst` call — with the CST label, its
**augment descriptors** (dominant taxon, % dominant, effective # taxa), and analysis results
attached as sample-level annotations. CST stays one method (VALENCIA); the descriptors
surface the community structure a single class flattens, without competing it.

### Command line (turnkey, config-first)

```bash
microfgt check   -c config.yaml                 # preflight: are the tools/paths this entry point needs present?
microfgt run     -c config.yaml -o out.h5mu     # import → CST → analysis → one .h5mu
microfgt run     -c config.yaml --executor snakemake   # emit a Snakefile to submit on a cluster instead
microfgt classify -i out.h5mu -o out.h5mu -m centroid  # (re)classify CST on an existing object
microfgt analyze  -i out.h5mu -o out.h5mu --transform clr --alpha shannon
microfgt compare  -i out.h5mu --verb alpha --predictors CST   # a hypothesis-test verb (alpha/beta/abundance/associate)
microfgt dashboard -i out.h5mu                                # launch the Streamlit dashboard (needs the [app] extra)
```

`microfgt check` is a preflight doctor — it inspects the *resolved* entry point and
reports missing binaries, R packages, DB paths, and region↔DB mismatches up front,
with actionable messages, instead of failing deep in a run. Copy-and-run example configs,
one per mode, are at the repo root: `example_16s_config.yaml` (runnable against the repo
fixtures), `example_metagenomics_config.yaml`, and `example_combined_config.yaml`.

### Sample sheet (how you name samples and combine arms)

microFGT does **not** guess a sample's id from its FASTQ filename, and it never tries to
auto-match your 16S and shotgun files — a wrong guess would silently pair the wrong patients.
Instead you hand it a small CSV that says, per sample, its id and where its reads live:

```csv
sample_id,16s_R1,16s_R2,shotgun_R1,shotgun_R2,group
PT01,16s/PT01_R1.fastq.gz,16s/PT01_R2.fastq.gz,mgx/PT01_R1.fastq.gz,mgx/PT01_R2.fastq.gz,BV
PT02,16s/PT02_R1.fastq.gz,16s/PT02_R2.fastq.gz,,,Normal
```

- `sample_id` is required and unique — it's the key that ties a sample's 16S and shotgun reads
  together, so **give the same sample the same id in both arms**. That declared id is what
  flows through the whole run, so the arms integrate correctly regardless of how the raw
  filenames were named.
- Read columns are optional; leave a cell blank for an arm a sample doesn't have.
  **Note: v1's 16S arm is paired-end only** — a 16S sample needs both `_R1` and `_R2` (the
  DADA2 stage errors on unpaired input; single-end amplicon support is not implemented yet).
- Any extra columns (e.g. `group`) are carried onto the object as per-sample `.obs` metadata.

Point a config at it with a top-level `samples:` key — per-run settings (region, primers, DB
paths) still live in the config; the sheet only supplies the per-sample files and metadata:

```yaml
samples: sample_sheet.csv
composition:
  reads: { region: V3V4, primers: { fwd: ..., rev: ... } }
  speciateit: { db: /path/to/vSpeciateIT_V3V4 }
```

Validate a sheet before a run — it checks that ids are unique and every listed file exists, and
reports how the arms overlap:

```bash
microfgt check-samplesheet -s sample_sheet.csv
```

See `example_sample_sheet.csv` at the repo root. (Sample-sheet runs currently use the local
executor; the Snakemake path is a follow-up.)

### Python API (for power users)

```python
from microfgt.io import import_virgo, import_valencia, import_speciateit, build_mudata

func = import_virgo("path/to/virgo_outputs/")          # dir of <sample>.out files
comp = import_speciateit("MC_order7_results.txt", "count_table.csv", fasta="asvs.fasta")
cst  = import_valencia("valencia_output.csv")          # sample-keyed CST/subCST/score

mdata = build_mudata(composition=comp, function=func, cst=cst)
# build_mudata materialises the taxon roll-up (composition_taxon) from the ASV-grain
# composition and attaches the augment descriptors automatically.
```

`import_speciateit` returns an **ASV×sample** composition (each ASV keeps its
classification, genus, and — with `fasta=` — its sequence); `collapse_to_taxon(comp)` is the
taxon roll-up CST reads.

`from microfgt.cst import classify_cst` and the `microfgt.analysis` module (compositional
transforms, α/β diversity, PCoA, ANCOM differential abundance — all bought from
scikit-bio, none reimplemented) operate on the same objects.

## What's validated vs. wired

This is the honest state, because it's what determines whether a result can be trusted.
**Stubs prove the plumbing; they do not prove the biology.** Anything marked ⚠️ has had
its wiring exercised end-to-end but has never run against the real tool's output.

> **Milestone (2026-09-17):** both running front-ends (16S *and* shotgun) have now been
> validated on **real HTCF runs** against genuine vaginal 16S + shotgun data (PIN, and the
> matched mother–daughter V3V4 cohort PRJNA779415) — moving the rows below from ⚠️ to ✅.
> Honest caveats: small *n* (a handful of samples), correctness is *spot-checked* (biologically
> sensible calls + the V3V4 ground-truth), **not** benchmarked against a labeled gold standard;
> and v1 is **paired-end only** (the DADA2 stage requires R1+R2). The cluster path is still unrun.

| Piece | Status |
|---|---|
| **Centroid CST → VALENCIA** | ✅ **Validated against ground truth (on the 2020 centroids)** — 99.94% subCST agreement on all 13,231 published samples (≥99.9% target), plus exact reproduction of genuine `Valencia.py` output on the head fixture. **Caveat:** that 99.94% is the *2020 paper-reproduction* set; the shipped **default is the 2024 (VALENCIA2) set** — recommended because it matches modern speciateIT taxonomy, but not itself benchmarked against a 2024-named gold standard (none is published). Select the validated set with `cst: {reference: "2020"}` (or `microfgt classify --reference 2020`) for paper-comparable calls. **Re-confirmed 2026-09-17 on pandas-3.0 / scikit-bio-0.7 — 99.94% holds** (no drift from the newer libraries). |
| `import_virgo` | ✅ Real-output validated (`virgo_sub*.out`). |
| `import_valencia` | ✅ Real-output validated (`valencia_genuine_output_head.csv`). |
| `import_virgo2` / `import_mgcst` (VISTA) | ✅ Real-output validated against the public ENA/PRJEB34536 fixtures (compiled gene matrix + taxon/KEGG annotation; genuine `vista_mgCSTs.csv`). |
| `import_speciateit` | ✅ **Real-output validated (2026-09-17)** — parsed genuine `MC_order7_results.txt` on real HTCF runs (PIN V1V3 + the matched mother–daughter V3V4 cohort PRJNA779415); the V3V4 `setup` ground-truth check matches 9/10 README species (+1 documented DB-version drift). |
| Orchestration wrappers (speciateIT, VIRGO2) | ✅ **Run against the real tools on HTCF (2026-09-17)** — speciateIT classifies real ASVs correctly; VIRGO2 maps + compiles real vaginal shotgun reads. |
| Preprocessing ladder (cutadapt, DADA2) | ✅ **Real-run validated (2026-09-17)** — real vaginal 16S denoised end-to-end (PIN V1V3; PRJNA779415 V3V4 with the built-in default `trunc_len`, no tuning), yielding biologically-correct taxa + CSTs. |
| Shotgun running front-end (fastp, minimap2 host removal, VIRGO2 map/compile, VISTA) | ✅ **Run end-to-end on real vaginal shotgun (2026-09-17)** → genuine VISTA mgCST calls, and the **clean both-arms integration validated** on the matched mother–daughter V3V4 cohort (PRJNA779415): 5/5 samples integrated across both arms, with 16S CST and shotgun mgCST **concordant per sample** (e.g. *L. crispatus* → CST I / mgCST 4; *Gardnerella* → CST IV-B). |
| Snakemake/Slurm executor | ⚠️ Snakefile is generated and asserted, but **never yet run on a cluster** — and has known path-quoting / directory-output bugs (see `design/refinement_backlog.md`). Real runs to date use the local executor. |

Run the suite (the CST gate needs VALENCIA's ~8 MB published dataset — stage it
with `python validation/fetch_valencia_published_data.py`, else it skips):

```bash
pytest -q
```

### Open validation IOUs

**Discharged 2026-09-17** by real HTCF runs on genuine vaginal data (PIN V1V3; the matched
mother–daughter V3V4 cohort PRJNA779415):

1. ~~**speciateIT real-output**~~ — `import_speciateit` now parses genuine
   `MC_order7_results.txt`, and the V3V4 `setup` ground-truth check matches 9/10 README species
   (+1 documented DB-version drift).
2. ~~**Real tool run (front-ends)**~~ — both running front-ends validated on real reads: 16S
   (cutadapt → DADA2 → speciateIT) yields biologically-correct taxa/CSTs; shotgun
   (fastp → minimap2 → VIRGO2 → VISTA) runs end-to-end to genuine mgCST calls.

**Still open:**

- **Cluster path** — the Snakemake/Slurm executor has never run on a cluster (the real runs
  used the local executor) and has known path-quoting / directory-output bugs
  (`design/refinement_backlog.md`). Discharge by a real cluster submission after those fixes.
- **Scale + correctness benchmarking** — real runs so far are small-*n* and correctness is
  spot-checked (biologically sensible calls; V3V4 ground-truth), not benchmarked against a
  labeled gold standard.

## External tools & reference data

microFGT ships **no installers and no reference data** (the tools' DBs are too large to
bundle); the orchestration layer locates an installed tool by configured path or PATH.
You only need the tools your **entry point** uses — run `microfgt check -c config.yaml` at
any time to see exactly what's still missing for *your* config, with install hints.

What each step requires, at a glance:

| Step | Tool | Reference data |
|---|---|---|
| Primer trim | cutadapt | none (primers are config strings) |
| Denoise | R + Bioconductor `dada2` | **none** — denoise-only; taxonomy is speciateIT's job |
| Taxonomy | speciateIT (`classify`) | **vSpeciateDB models, ~2.6 GB** (separate download) |
| Shotgun QC | fastp | none |
| Shotgun host removal | minimap2 + samtools | host genome, e.g. GRCh38 (separate download) |
| Shotgun gene profiling | VIRGO2 (`VIRGO2.py`) | VIRGO2 catalog + `Index/`, `AnnotationTables/` (separate download) |
| Shotgun mgCST | VISTA (Rscript) | `VISTA_data/` (ships in the VISTA repo) |

So the entire 16S **front-end** (real cutadapt + real DADA2 on real reads) runs with only
small installs and **no multi-GB download** — the models are needed solely for the
taxonomy-assignment step. Per-tool setup follows.

### cutadapt + DADA2 (16S front-end — only if you start from raw FASTQs)

Both are on bioconda and come with `environment.yml`:

```bash
conda env create -f environment.yml && conda activate microfgt && pip install -e .
```

Neither needs a reference database. (If you enter at an existing ASV table or existing tool
outputs, you don't need these at all.)

### speciateIT (taxonomy — turns ASVs into the `composition` matrix)

A C++ 16S classifier from the Ravel Lab: <https://github.com/ravel-lab/speciateIT>.

1. **Binary.** A precompiled `classify` ships in the repo under `bin/linux` / `bin/macosx`
   — clone the repo and put that directory on your `PATH` (no build needed in the common
   case). <!-- TODO: verify build-from-source steps + whether a bioconda package or
   container exists; web access was unavailable when this was written. -->
2. **Models (~2.6 GB, separate download).** The `vSpeciateDB` reference models are *not*
   bundled with the binary. They come as region-specific directories named
   `vSpeciateIT_<REGION>` for **V1V3**, **V3V4**, or **V4** — download the one matching your
   amplicon. <!-- TODO: confirm the exact figshare/Zenodo URL + the distributed directory
   names before relying on this; not verifiable when written. -->
3. **Point microFGT at it** in your config:

   ```yaml
   composition:
     speciateit:
       classify: classify                 # binary name on PATH, or an explicit path
       db: /path/to/vSpeciateIT_V3V4       # the downloaded model directory
   ```

4. **Verify:** `microfgt check -c config.yaml` confirms the binary resolves, the `db` path
   exists, and the region matches the model directory before you run.

speciateIT runs as `classify -d <db> -i <asvs.fasta> -o <outdir>`, always writing
`MC_order7_results.txt`; microFGT then owns the join to the ASV count table.

### Shotgun metagenomics (fastp → host removal → VIRGO2 → VISTA)

The metagenomics arm runs the whole chain for you; you install the tools and point the
config at them. `microfgt check -c config.yaml` verifies every binary and DB path (and can
checksum `VIRGO2.py`) before a multi-hour run starts.

- **fastp** (read QC) and **minimap2 + samtools** (host-read removal against a host genome
  such as GRCh38) — both on bioconda.
- **VIRGO2** — the vaginal gene catalog + mapping pipeline from the Ravel Lab. Install
  unpacks to a directory holding `VIRGO2.py`, `Index/`, and `AnnotationTables/`. VIRGO2 maps
  **single-end**, so microFGT concatenates each pair's R1+R2 before mapping, then compiles the
  per-sample outputs into the `function` modality (and its taxon roll-up).
  <!-- TODO: confirm the canonical VIRGO2 catalog download location + size and runtime deps
  before relying on this. -->
- **VISTA** (Rscript) classifies the compiled matrix into **mgCSTs**; `VISTA_data/` ships in
  the VISTA repo, so no separate DB download.

Point microFGT at all four in the config:

```yaml
metagenomics:
  reads:
    fastq_dir: raw_shotgun/        # paired *_R1*/*_R2* FASTQs — runs the full chain
  host_ref:   /path/to/GRCh38.fna.gz
  virgo2_dir: /path/to/VIRGO2      # has VIRGO2.py, Index/, AnnotationTables/
  vista_repo: /path/to/VISTA       # contains VISTA_data/
  threads: 8
mgcst:
  # vista_output: mgCSTs_x.csv     # OR import an existing VISTA call instead of running it
```

Enter partway in by giving `metagenomics.compiled` (a `VIRGO2_Compiled` matrix) or
`mgcst.vista_output` (an existing VISTA call) instead of raw reads.

## Build status (phased — walking skeleton first)

- **P0 — Skeleton + install** ✅ one-command install, CI.
- **P1 — Integrated object + importers** ✅ MuData + importers grounded in
  `prototype/real_fixtures/FORMATS.md` (see the validation table above).
- **P2 — Centroid CST behind `classify_cst`** ✅ faithful VALENCIA port behind a swappable
  `classify_cst(composition, method=...)` seam; validated 99.94% vs the paper.
- **P3 — Orchestration wrappers** ✅ `run_speciateit` / `run_virgo`, each locating an
  installed tool, recording provenance, and handing output to the matching importer.
- **P3.5 — Preprocessing front-end + multi-entry workflow** ✅ one stage registry, two
  executors (local resolver + Snakefile generator), the multi-entry ladder, and
  `microfgt check`. DADA2 orchestrated with region-aware, overridable defaults.
- **P4 — Analysis + viz + turnkey CLI** ✅ commodity stats bought from scikit-bio
  (transforms, α/β diversity, PCoA, ANCOM); the `compare` hypothesis-test verbs; minimal
  matplotlib viz behind the `viz` extra; the config-first CLI; and a Streamlit dashboard
  (`microfgt dashboard`, the `[app]` extra). Zero-count samples handled honestly (excluded
  from compositional steps, recorded in `uns`, never silently dropped).
- **P5 — Shotgun metagenomics arm** ✅ VIRGO2 + VISTA/mgCST importers (validated on real
  ENA/PRJEB34536 fixtures) and the running front-end (fastp → minimap2 host removal → VIRGO2
  map/compile → VISTA) added as stages in the *same* registry and executors, sharing the
  multi-entry resolver. Real-tool correctness is part of the real-run IOU (stub-validated end
  to end today).
- **Augment descriptors (not rival CST)** ✅ CST is one blessed method (VALENCIA); the
  diffuse/continuum structure it flattens is read out by *augmenting* the label with
  per-sample descriptors — dominant taxon, % dominant, effective # taxa — attached alongside CST.
  (This replaces the earlier "alternative CST methods" framing: the object carries ASVs and
  their sequences as the source of truth, with a materialised `composition_taxon` roll-up.)
