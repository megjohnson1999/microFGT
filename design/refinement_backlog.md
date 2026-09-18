# 16S Refinement Backlog

*Things to tighten on the **16S path we already have**, before returning to the multi-omics /
shotgun phase (see `candidate_datasets.md`). Started 2026-08-14. Not a rigid plan — a living
list; pick items as they make sense.*

Status key: **TODO** · **NEXT** (queued to do soon) · **DONE** (kept for continuity)

---

## Object / data model

- **NEXT — Lift imported sample metadata to the global obs.** Right now `import_phyloseq` puts
  the phyloseq's `sample_data` on the **composition** modality's obs, while microFGT's own
  additions (CST, descriptors) sit on the **global** `MuData.obs`; `merged_obs` unions them so
  the dashboard shows one clean list. But clinical metadata (`HIV_status`, `PID`, `week`, …)
  describes the **sample**, not the 16S assay — so it should live at the **global** level so it
  applies across *all* modalities. Small refactor in `build_mudata` (lift the modality's
  sample_data to global obs). **Do this before a second modality lands**, when it starts to
  matter. *(Requested 2026-08-14.)*

## Importer generalization (`import_phyloseq` → "accept standard formats", principle #3)

- **TODO — refseq slot support.** Today the importer assumes sequences ARE the `taxa_names`. A
  phyloseq named `ASV1…ASVn` with sequences in `refseq()` would keep the ids and lose the real
  sequences. Read `refseq()` when present.
- **TODO — unrecognized rank *names*.** Classification looks for ranks named
  `Genus_Species`/`Species`/`Genus`. A QIIME2/mothur object with `Rank1…Rank7` imports
  structurally but lands **all-`Unclassified`**. Map/relabel common rank schemes. (GTDB-style
  *prefixes* on recognized ranks are now handled — see Recently done.)
- **TODO — placeholder taxonomy labels.** Both FRESH and HVTN carry "no-resolution" placeholders
  like `Bacteria Domain` / `Lactobacillus Genus`. Optional display-cleanup (trust-but-tidy),
  separate from the GTDB-prefix strip.

## Dashboard depth

- **TODO — per-modality coverage view.** Surface `build_mudata`'s `Reconciliation` in the UI
  ("N have 16S, M have shotgun, K have both"). Low-value now (one modality), high-value the
  moment a second feeder lands — pairs with the metadata-to-global refactor above.
- **TODO — composition / CST landing view.** Stacked-bar composition + CST distribution as a
  landing overview, so the object is legible before running a verb.
- **TODO — subset UX.** Make it easy to drop controls / replicates (`sample_type`, `is_rep`)
  for clean stats; maybe a default "real samples only" toggle.

## Analysis breadth / cleanliness

- **TODO — verbose taxonomy labels.** The FRESH taxonomy carries placeholder labels like
  `"Lactobacillus Lactobacillus Genus"` / `"Bacteria Domain Bacteria Domain"` (faithful to the
  source). Optional: a display-cleanup pass (trust-but-tidy), without re-deriving taxonomy.
- **TODO — more alpha metrics + prevalence filtering.** Expose Simpson/observed/etc. cleanly;
  add low-count / low-prevalence taxon filtering as a pre-step.

## Docs & robustness

- **TODO — end-to-end 16S workflow doc.** Write up phyloseq → CST → verbs → dashboard so it's
  reproducible (README or a design note).
- **TODO — display polish.** Wide contingency tables (`associate` on high-cardinality vars) dump
  awkwardly in the CLI/terminal; better rendering/truncation. Better error messages on bad
  selections.

### General-tool findings — HTCF real-run exercise (2026-09-17)

Surfaced by running microFGT end-to-end on real data (PIN cohort) on HTCF. All are
**general-tool** improvements — they help any user on any dataset, NOT PIN-specific tuning.
Roughly ordered by user impact.

- **BUG — `ancombc` differential-abundance verb crashes on every call.** `analysis/abundance.py`
  `_run()` does `from skbio.stats.composition import ancombc`, but scikit-bio has no `ancombc`
  (only `ancom`, `dirmult_ttest`, `dirmult_lme`). `ancombc` is the DEFAULT in `run_verb` / CLI
  `compare` / dashboard, so the flagship DA verb is dead on arrival. Fix: call a real scikit-bio
  function with a matching signature and pick a valid default. (Working DA lives on the separate
  `differential_abundance`/`ancom` path via `--diffabund-group`.)
- **BUG (FIXED 2026-09-17) — VISTA renames samples → mgCST fails to join → all-NaN → write crash.**
  `run_VISTA.R` loses the compiled matrix's sample ids and R relabels columns positionally as
  `V1..Vn` (`V1`=Gene col, `V2`=first sample, …). microFGT imported those `V2/V3` ids, which
  don't join the real sample ids (`25056/22460`), so the global `.obs['mgCST']` came out all-NaN
  and `.h5mu` write died with "Can't implicitly convert non-string objects to strings" writing
  `/obs` key `mgCST`. Fixed in `orchestrate/vista.py` (`_restore_sample_ids`): remap `Vk` →
  compiled matrix's k-th column, robust to VISTA dropping/reordering low-count samples. Hit live
  on the first combined PIN run (VISTA classified fine — mgCST 18 & 11 — just mislabeled).
- **BUG (FIXED 2026-09-17) — VIRGO2 orchestration doubles paths under a relative `--workdir`.**
  `orchestrate/virgo2.py` ran VIRGO2.py with `cwd=outdir` but passed `-r`/`-o` as paths relative
  to the *parent* cwd, so a relative `--workdir` made them resolve *inside* `outdir` (doubled) →
  bowtie2 couldn't find its input → no `.sam` → VIRGO2 crashed (its own `files2del` typo and
  missing `check=True` masked the cause). Fixed by resolving `outdir` to absolute in
  `run_virgo2_map`/`run_virgo2_compile`. Hit live on the first combined PIN run; workaround was
  an absolute `--workdir`. (Fix is in the canonical repo; HTCF clone still needs the pull.)
- **BUG (FIXED 2026-09-17) — `microfgt check` crashed with a raw traceback** on any config that
  doesn't resolve (empty, mistyped key, or a sample-sheet config). `stages/check.py check()`
  called `resolve()` unguarded. Fixed: `check()` catches `StageResolutionError` → clean "no
  usable entry point" result.
- **GAP (FIXED 2026-09-17) — `microfgt check` couldn't preflight a sample-sheet (`samples:`)
  run.** It didn't stage the sheet, so `provided_artifacts` was empty → resolve failed (→ the
  crash above). Fixed: `cli._cmd_check` now stages the sheet (as `run` does) before checking, so
  the doctor works on the exact config style users are told to use — and validates the sheet
  (unique ids, declared files exist) as part of the preflight. (Both hit live on the first PIN run.)
- **BUG — `check`'s R-package probe gives a false MISS on heavy packages.** `check.py
  _r_has_package` runs `requireNamespace` with a 60 s subprocess timeout; `dada2` (big
  Bioconductor tree) exceeds it on first load, so `check` reports `MISS rpackage 'dada2'` even
  when it's installed and working (the real run used it fine). Also makes `check` slow (>120 s).
  Fix: raise/remove the timeout, warm/caching probe, or a lighter check.
- **GAP — no per-sample low-read-depth guard/warning.** microFGT excludes *zero*-count samples
  honestly, but a sample surviving with very few reads still gets a confident CST and no flag. In
  the PIN shakeout, `trunc_len [250,250]` collapsed sample 25056 to **81 reads**, yet the tool
  emitted CST IV-B for it (vs III-A at 29,431 reads with `[0,0]`) with no warning. Add a
  low-depth guard: flag samples below a read-count/retention threshold in `.uns` + CLI output so
  results built on ~no data are surfaced, not silently trusted. (Related: DADA2 `trunc_len`
  materially changes the CST call, so parameter sensitivity should be visible to the user.)
- **GAP — DADA2 quality profile reports R1 only.** `quality_profile.tsv` carries no R2, so a
  user tuning `trunc_len` on paired data is blind on the reverse read. Add R2 to the profile.
- **GAP — no config validation.** A mistyped key (`fastqdir`, `compositon:`) is silently dropped
  → entry point vanishes → cryptic "cannot produce artifact" error, no hint the key was ignored
  (`config.py`). Add a key allow-list / "ignored keys" warning; optionally stat input paths.
- **BUG — resolver misroutes a direct `composition` entry.** `stages/resolve.py:39-40` step-1
  shortcut greedily picks `integrate_combined` when `composition` is provided directly (prebuilt
  h5ad), then dies on the shotgun `sg_reads` artifact for a pure-16S run. Fix: exclude the
  multi-input mudata producers from the step-1 shortcut. (Common fastq/asv-table paths unaffected.)
- **BUG — Snakemake executor breaks on real paths.** `stages/executors.py` interpolates
  workdir/config/output into the `shell:` string unquoted (breaks on spaces — incl. the repo's
  own path) and emits directory artifacts as file `output:` rather than `directory(...)`. Never
  run on a real cluster. Fix before trusting the cluster/scale path.

- **TODO — single-end 16S support (v1 is paired-end only).** The DADA2 denoise stage
  (`scripts/dada2_run.R`) requires matching `_R1`/`_R2` and errors on unpaired input — there is
  no single-end branch, and `mergePairs` is called unconditionally. The sample-sheet + README
  docs were corrected 2026-09-17 to stop implying single-end works (they said "`_R2` optional").
  Add a single-end path (denoise `_R1` alone, skip merge) so single-end amplicon datasets can
  run — e.g. the PIN cohort's 16S is 68/72 *genuinely* single-end (confirmed via
  `fasterq-dump`: 1 read/spot, R2 never deposited), leaving only ~3 usable paired subjects.
  Decision (2026-09-17): paired-end-only is acceptable for v1; this is the come-back item.

### Analysis-layer coverage (design gaps — the object *structure* supports these; the verbs don't yet)

Surfaced discussing whether the integrated object is well-shaped for common downstream analyses
(2026-09-17). Consensus: the MuData container is a solid, standard, extensible foundation and
serves single-modality analysis well out of the box; the gaps are all "verbs to exploit the
structure," fixable without reshaping the object.

- **Cross-modal analysis verbs (highest priority — it's the tool's whole point).** The object
  sample-aligns the 16S and shotgun modalities + reports the shared subset, but the analysis
  verbs are single-modality. Add joint verbs — e.g. `function ~ CST`, `taxon ↔ gene` — that run
  on the **shared subset** with honest `n_used` / `n_dropped-for-missing-modality` bookkeeping
  (same reconciliation the single-modality verbs already do). Until then multi-omics analysis
  means manual mudata wrangling. (Also flagged in `design/candidate_datasets.md`.)
- **First-class longitudinal support.** FGT is heavily longitudinal (pregnancy cohorts; lab has
  "committed to longitudinal"). Handleable today via `subject`/`visit` obs columns + the compare
  verbs' `--subject`, but no first-class subject×timepoint structure — trajectories, CST/mgCST
  state-transitions, time-series lean on the user.
- **Phylogenetic diversity (UniFrac / Faith's PD).** Needs an ASV tree; `composition` keeps the
  sequences to *build* one, but no tree is attached, so it isn't turnkey. (Mild for FGT
  specifically — low-diversity, species-level — but common in general microbiome work.)
- **More omics + strain resolution.** Container is extensible (add metabolite/cytokine
  modalities — MOMS-PI-style), which is the right bones, but the tool only builds 16S + shotgun
  today, and mgCST subtypes (`mgSs`, strain-ish) are deferred — ironic for a mother–daughter
  shared-strain cohort like PRJNA779415.

---

## Recently done (context)

- **DONE — validated on a 2nd real cohort (HVTN, 4,856 samples × 34,904 ASVs).** Two-regime
  `effective_taxa` signal replicated (CST I ~1.8 → CST IV-C ~10.5); rich clinical metadata
  (STIs, contraceptives, HIV time-to-event). Surfaced + fixed the GTDB-prefix gap below.
- **DONE — strip GTDB rank prefixes in `import_phyloseq`** (`g_Lactobacillus` → `Lactobacillus`).
  HVTN's GTDB taxonomy was leaking `d_`/`g_` prefixes into 71% of labels and breaking genus
  extraction; now clean, without touching normal binomials or `Ca_` (Candidatus) names.
- **DONE — cutoff-free `effective_taxa` descriptor + adjustable `taxa_over_threshold`** (the
  stored-vs-derived split; dashboard slider).
- **DONE — switched demo object to the full merged FRESH dataset** (5,659 samples ×
  21,964 ASVs, rich metadata incl. `HIV_status`/`PID`/`week`) — was mistakenly using one run
  (893, MD1048). See `[[fresh-dataset]]` memory.
- **DONE — collapse MuData `modality:` obs-column prefixes** so the variable list shows each
  variable once (was duplicated `composition:HIV_status` vs `HIV_status`).
- **DONE — dashboard loads MuData via `st.cache_resource`** (was `cache_data`, which pickled and
  failed on real objects).
