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
