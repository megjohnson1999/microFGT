# microFGT — Design Vision

*Working vision doc. Captures the high-level design decisions — what microFGT is,
who it's for, and how it's shaped. Complements `fgt_workflow_map.md` (the per-step
"orchestrate vs. reimplement" analysis) and `implementation_plan.md` (the build/language
decisions). Last updated 2026-08-13.*

Refer back here whenever the thread slips. This is the map.

---

## In one sentence

A **general-purpose platform** for FGT microbiome analysis: a standard **integrated data
object** plus a **friendly, flexible analysis surface** — so researchers stop hand-rolling
scripts for every study.

## The two pillars (the intent behind the tool)

1. **The integrated object** — one standard container holding multiple data modalities
   keyed to the same samples.
2. **Streamlined, user-friendly analysis** — a tool that makes working with that object
   easy, so people don't reinvent the analysis each time.

Everything else serves these two.

## Purpose: a platform, and a general one

- **A platform, not a study.** Success is not "we found result X." Success is: a researcher
  can load their data and run the standard analyses, across modalities, without reinventing
  them. It's infrastructure.
- **General, not lab-specific.** The tool is meant to be useful to the FGT field at large.
  Our lab's questions and data are the **first concrete instance** we build and test against —
  they keep the work real, but they are **not the spec.**
- **The discipline that keeps "general" from becoming boundless:** *design for the field;
  build against our lab's data as the first instance.* When it starts feeling too wide,
  return to the concrete instance.

## The user

- **Primary user:** a researcher who wants to load their data and *explore* it easily. For
  them the integrated object is **invisible plumbing** behind a friendly surface.
- **Secondary user:** a bioinformatician who can drop down to the object **directly** if they
  want. This escape hatch is essentially *free* — if the object is clean and well-structured,
  you don't build a separate thing for them, you just don't hide the object.

## Architecture: three layers

```
Layer 3   Friendly surface (dashboard / app)   ← primary user; "easy"
Layer 2   Analysis layer (blessed methods)     ← called by BOTH surface and power users
Layer 1   Integrated object (multi-modal)      ← bioinformatician can touch directly
```

- The dashboard is a **thin surface over the same functions** a bioinformatician would call —
  you design the analysis layer **once**; the surface just presents it. Never build
  "dashboard logic" and "analysis logic" as two separate things. *Now concrete:* every analysis
  verb returns a uniform `AnalysisResult` (tidy table + headline stats + a **declarative plot
  spec**), so the surface renders results without re-deriving them — the dashboard becomes a
  spec-builder + a renderer.
- **Build order:** object → analysis → surface. (This is why the object work came first and
  wasn't premature — it's the foundation both upper layers stand on.)

## Guiding principles

1. **Design for the field; build against our lab's data first.** (The general/concrete
   discipline, above.)
2. **Consume, don't rebuild — the *general* processing.** Don't re-implement 16S denoising;
   many pipelines already do it (QIIME2, mothur, dada2, nf-core/ampliseq, or the lab's
   `16s_dada2_valencia`). microFGT is **pipeline-agnostic** — `16s_dada2_valencia` is just
   *an* example and our first test input, **not a field standard we assume everyone runs.**
   The object **begins where general processing ends.**
   - **But inputs arrive at very different levels of *FGT-completeness*.** `16s_dada2_valencia`
     is unusually complete — it runs speciateIT + VALENCIA, so its output already carries
     species taxonomy *and* CSTs. A generic QIIME2/SILVA run gives genus-level composition and
     **no CST.** So the *general* denoising stays upstream, but the **FGT-specific
     characterization — species reclassification, CST, name reconciliation — is microFGT's
     own**, and it must add it when the input lacks it. (See *What makes it FGT-specific*.)
3. **Accept standard formats, not one lab's pipeline.** Labs use different setups (QIIME2,
   mothur, different DBs), but they nearly all export the same containers — phyloseq, BIOM,
   feature+taxonomy tables. Accept *those*. Reading the standard phyloseq container (rather
   than any one pipeline's flat files) is the general move, and `import_phyloseq` is the entry
   point for it.
   - **But "reads *any* phyloseq" is aspirational, not today's code.** As built, the importer
     assumes the output *conventions* of the `16s_dada2_valencia` pipeline (the pipeline the
     FRESH data happened to be processed through — **not** anything FRESH-specific). The real
     axis of variation is the **toolchain/convention, not the dataset.** Two known limits:
     - *Sequences must be the taxa_names.* It does not read a `refseq` slot — a phyloseq named
       `ASV1…ASVn` with sequences in `refseq()` would keep the ids and lose the real sequences.
     - *Taxonomy ranks must be named `Genus_Species` / `Species` / `Genus`.* A QIIME2/mothur
       object with `Rank1…Rank7` imports structurally but lands all-`Unclassified`, collapsing
       everything to one taxon downstream.
   - This is principle #1 made concrete — *built against the first instance (FRESH), not yet
     general.* Expected, not a defect. These two are the first items on the generalization list.
4. **Standard methods, user-chosen inputs.** The computation is fixed and reproducible (one
   right way to run a diversity calc or a differential-abundance test); the *user* freely
   chooses **what to look at** — which variable to group/compare by, which subset, which taxa.
   A playground on top of a rigorous engine. *(Realized:* each verb takes `(object + which obs
   variables play which role + a subset)` — the spec IS the user's choice; the engine is fixed.)
5. **Free-form metadata; roles assigned at analysis time.** Metadata columns differ per study;
   subject/timepoint/group are not canonical fields. Nothing is hardcoded — which is exactly
   what makes flexible "look at it by *this* variable" exploration possible.
6. **Honest reconciliation.** Samples measured on different modalities overlap only partially
   (e.g. 16S on all, shotgun on a subset). Take the union; never silently drop.

## What goes in the object (data landscape — in progress)

- **Composition (16S)** — *who is there.* The **first feeder**, wired today via
  `import_phyloseq`. ✅
- **Function (shotgun metagenomics)** — *what they are doing* (genes/pathways; VIRGO/VMGC).
  Designed slot exists; **feeder still open** — no local pipeline or data yet (can source
  external data). This is the modality that makes "integrated" meaningful.
- **Candidates beyond:** metabolomics (SCFAs, amines), host/immune (cytokines, inflammation),
  clinical (Nugent, pH, diagnoses). Not yet scoped in.
- **The make-or-break design concern:** *how modalities align* — the sample/subject join key,
  and differing resolution between modalities. **Integration lives or dies here.**

## What makes it FGT-specific (and what's borrowed)

microFGT is **general analytical machinery + FGT biology encoded in the analysis layer.**
Being honest about which is which matters — it's also the honest version of "the gap."

- **Borrowed / general (not FGT-specific):** the container itself (MuData; general multi-omics
  containers like MultiAssayExperiment already exist), diversity math, differential-abundance
  methods, ordination, plotting. We do **not** claim to invent an integrated container.
- **FGT-specific (the domain knowledge that makes it an *FGT* tool):**
  - Species/subspecies taxonomy — *L. crispatus* vs. *iners* vs. *jensenii* matter clinically;
    the Gardnerella → *Bifidobacterium vaginale* mess (speciateIT / VMGC references).
  - The **inverted health axis** — in the vagina, *low* diversity + Lactobacillus dominance =
    healthy (opposite of the gut). General tools assume "diversity good."
  - FGT characterizations & clinical scores — CSTs, molecular-BV, Nugent/Amsel, Lacto depletion.
  - FGT modality semantics — D- vs. L-lactate, biogenic amines / sialidase in BV, cervicovaginal
    cytokines, pH.
  - Low-biomass & region handling — spike-ins for absolute abundance; V3V4 vs. V4.
- **Where the FGT-specificity concentrates:** because microFGT *consumes* upstream pipelines,
  most FGT-specific *processing* (species taxonomy, region handling, spike-ins) happens
  **upstream.** So microFGT's own FGT-specificity lives mainly in the **analysis +
  characterization layer being FGT-aware** — not in the container.

**The honest gap (positioning).** General containers *and* general microbiome-analysis tools
both exist. What's missing is an **FGT-specific integrated tool** — CST-aware, bundling and
reconciling the FGT toolchain, joining composition + function for FGT questions, and *easy.*
Don't claim "nobody has an integrated object"; claim "no FGT-tuned, easy, integrated analysis
tool — people assemble it from general parts plus custom scripts."

*Verified against the tool landscape (Aug 2026):* FGT tools are single-purpose steps chained by
hand — speciateIT, VALENCIA, VIRGO/VMGC, mgCST (mostly one lab). General platforms
(MicrobiomeAnalyst, QIIME2, mia) exist but are gut-oriented (no CSTs; "diversity = healthy" is
backwards for the FGT). The **one real challenger is VISTA** (Ravel lab, mBio 2026) — an
FGT-specific, function-aware framework with an interactive Streamlit app — **but it's
shotgun-only typing**, not 16S+function co-keyed, not a general analysis layer (no DA/diversity),
and not a reusable multi-modal object. **Name VISTA proactively;** microFGT's differentiation is
exactly the three things it lacks. Worth reading as prior art / design reference.

## Where it stands today (concrete)

- `import_phyloseq` built and validated on **real FRESH data** (5,772 ASVs × 893 samples):
  composition + existing CST + descriptors, biologically sensible. **Handles phyloseqs that
  follow the `16s_dada2_valencia` output conventions only so far** — sequences-as-taxa_names,
  GTDB-style rank labels; this is a per-toolchain limit, not FRESH-specific (see principle #3).
- CST computed one **blessed way** (VALENCIA-faithful, 99.9% reproduction), **augmented** with
  descriptors: dominant taxon, dominance %, # taxa > 10%.
- **Analysis layer: the statistical toolkit + the contract Layer 3 binds to.** Four
  hypothesis-test verbs over a uniform `AnalysisResult`: `compare_alpha` (group test /
  covariate-adjusted OLS / mixed model for repeated measures), `compare_beta` (PERMANOVA +
  dispersion), `associate` (CST ↔ clinical — chi-square / Fisher / Kruskal / Spearman by
  dtype), `compare_abundance` (covariate-adjusted, FDR-corrected differential abundance).
  Each takes *the object + which obs variables play which role + a subset* (principles #4/#5
  made real). Validated on FRESH: Shannon differs by CST (Kruskal p≈1e-94), CST explains
  R²=0.40 of Bray–Curtis variation (p=0.005), and adjusted DA recovers the textbook signal
  (Lactobacillus depleted, *Fannyhessea/Prevotella/Dialister* enriched in CST IV-B).
- **Engine decision (hybrid), resolved by what's real:** the portable core is pure Python
  (scipy / statsmodels / scikit-bio — including **ANCOM-BC** and a Dirichlet-multinomial mixed
  model for longitudinal DA, both brand-name adjusted methods that need no R). Orchestrating R
  (MaAsLin2 / DESeq2 / vegan's `adonis2`) stays a latent extension point behind the same verbs
  — added when a project needs those exact tools, **not** shipped as untested glue.
- Function layer: object slot exists, VIRGO orchestration wrapper exists, **no feeder run yet.**

## Open questions / what's left to decide

- **Data landscape** — which modalities the object holds, and how they align.
- **The gap** — which cross-data questions *can't* be answered today (currently done by hand in
  per-study scripts).
- **What integration actually buys** — co-storage vs. real cross-modal analysis. The value
  question. (Preview: the flexible exploration in principle #4, applied *across* modalities —
  "show function grouped by CST," "color composition by a metabolite" — is the payoff.) *The
  single-modality analysis layer now exists and is validated; the untested half is running
  those same verbs* across *modalities — which needs a second real feeder. `associate`'s
  continuous×continuous path (taxon ↔ metabolite) is the seam where that first shows up.*
- **Final scope** — what microFGT builds vs. consumes.

## What microFGT is NOT

- **Not a 16S pipeline** — that exists; we consume its output.
- **Not tied to one lab's** questions or pipeline.
- **Not a fixed-report generator** — it's for interactive, user-driven exploration.
