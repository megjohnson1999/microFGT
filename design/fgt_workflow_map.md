# FGT Analysis Workflow Map

Working artifact for deciding what microFGT should be. **Not** a build spec yet —
this maps the FGT microbiome analysis workflow, finds where the existing tools are
rigid or disappointing, and from that the structure/language/first-spike choices
fall out.

This is **general tool development** — a tool for FGT microbiome analysis as a field,
usable on anyone's data. Not scoped to any one project or cohort.

## Goal (the new direction, 2026-06-16)

A flexible, tunable, **one-stop** FGT analysis tool — runs the tools (or reimplements
their concepts) and lets method choices across the *whole* workflow be explored and
improved, vs. the rigid take-it-or-leave-it tools people are frustrated by. "Be open
to other methods" applies at every decision point, not just CST. Language/structure
is **open** — let it follow from this map, don't assume R.

## How to fill each step

For every decision point below:
- **What happens here** — the analytical decision the user is making.
- **Current tool** — what exists (speciateIT / VALENCIA / VIRGO / DADA2 / mia / etc.).
- **Where it's rigid / disappointing** — install pain, inflexibility, format lock-in,
  questionable method. Ground against the real fixtures in `prototype/real_fixtures/`
  (see `FORMATS.md`) — don't guess the formats.
- **Improvability × tractability** → orchestrate (wrap it), reimplement (the concept,
  natively, maybe better), or leave alone. Heavy trained-model/DB tools
  (speciateIT C++, VIRGO gene catalog) lean *orchestrate*; tractable method-y pieces
  are where "better than current" could live.

---

## The workflow

### 1. Raw input
_16S amplicon reads; metagenomic reads._
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

### 2. Preprocessing & denoising
_QC, trimming; 16S → ASVs (DADA2); metagenomics → QC/host-removal._
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

### 3. Taxonomic assignment
_Who is present. speciateIT (16S, species-level)._
- What happens here: assign each ASV/sequence to a taxon, then roll up to taxon×sample.
- Current tool: speciateIT — C++, 7th-order Markov-chain models (`MC_order7`). Classifies
  **sequences, not samples**; output is ASV→(classification, posterior, #decisions).
- Where it's rigid: (a) C++ build is install pain; (b) trained reference models — not
  something you'd casually swap or retrain; (c) **sample identity isn't in the output** —
  user must join ASV→taxon against a separate ASV count table to get taxon×sample. The
  integration burden lands on the user.
- Orchestrate / reimplement / leave: **orchestrate** — trained C++ classifier, not worth
  reinventing. Value to add = wrap the install + own the ASV→sample join so the user doesn't.

### 4. Community state typing (CST)  ← richest method-development target
_What kind of community. VALENCIA, centroid-based. Scott's named example._
- What happens here: assign each sample to a vaginal community state type from its taxonomic
  composition.
- Current tool: VALENCIA (`Valencia.py`) — nearest-centroid to **13 fixed reference centroids**
  derived from >13,000 N. American women. subCST = argmax of 13 `_sim` cols; score = max sim.
- Where it's rigid: (a) **taxa names must exactly match VALENCIA's convention**
  (`Lactobacillus_crispatus`, `g_Genus`, `f_Family`) — Prevotella "often causes problems";
  brittle string-matching is real usability pain; (b) **fixed reference centroids** — a sample
  is forced to its nearest of 13 types, no native novelty/out-of-distribution detection, can't
  discover new states; (c) nearest-centroid is a *simple* method — one similarity score, no
  calibrated uncertainty, no soft/mixed-membership assignment; (d) Python script wanting an exact
  input CSV shape.
- Orchestrate / reimplement / leave: **reimplement the concept, openly.** This is where
  "beyond centroid" lives — robust taxa-name reconciliation (kills the #1 pain), then explore
  model-based clustering / classifiers / soft assignment / novelty detection. Validate any
  alternative reproduces VALENCIA on the published data first (we have it: 99.9% match), then
  show where it does better.

### 5. Functional profiling
_What genes / capacity. VIRGO (metagenomics)._
- What happens here: map metagenomic reads to a non-redundant vaginal gene catalog → per-gene
  counts per sample.
- Current tool: VIRGO — large reference gene catalog + mapping. One `<sample>.out` per sample,
  3 cols (geneID, read_count, gene_length), no header; gene→taxon/length maps in the catalog.
- Where it's rigid: (a) huge reference DB + mapping step = heavy install/run; (b) output is
  per-sample files — **user stacks them into a wide matrix themselves**; (c) interpretation
  needs the catalog annotation joins.
- Orchestrate / reimplement / leave: **orchestrate** — catalog-backed, not worth reinventing.
  Value = wrap the run + own the stack-into-matrix + annotation joins.

### 6. Normalization & compositional transforms
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

### 7. Diversity (alpha / beta)
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

### 8. Differential abundance
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

### 9. Cross-modality integration (composition × function)
_Per-sample join of the tracks into one object — the scientifically interesting questions live here._
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

### 10. Visualization
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

### 11. Clinical / sample metadata association
- What happens here:
- Current tool:
- Where it's rigid:
- Orchestrate / reimplement / leave:

---

## Synthesis (fill after mapping)

- **Where the real value is** (which steps are worth reimplementing vs. wrapping):
- **What "one-stop / easy to install" demands** (packaging, distribution, who runs it):
- **Language/structure that falls out**:
- **First spike** (one method, reimplement, validate it reproduces the real tool on real data, try one alternative):
