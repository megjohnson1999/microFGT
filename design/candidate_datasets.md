# Candidate External Datasets — for the shotgun / multi-omics phase

*Parked reference (researched 2026-08-14). We are finishing/refining the **16S** path first;
this is here to pick up when we return to adding a second modality (shotgun metagenomics) and
proving cross-modal integration. Goal at that point: a public vaginal/cervicovaginal dataset
with **16S + shotgun on the same samples** plus clinical metadata.*

---

## The key insight (why this matters, and why it's not scary)

**Partial overlap is the *normal* shape of real FGT multi-omics — including our own cohort.**
16S is cheap, so it's run on everything; shotgun is expensive, so it's run on a curated
subset. So the realistic input is "16S on all N samples, shotgun on a subset, both on the
intersection." Clean 1:1 pairing is the *unusual* case.

microFGT already anticipated this. Principle #6 (honest reconciliation) is implemented:
`build_mudata` keeps the **union** of sample IDs across modalities, never silently drops, and
emits a `Reconciliation` report (per-modality counts + how many are shared across all assays).
So a 893-16S / 122-shotgun object holds all 893 and tells you the shared subset is 122. The
object layer is built for this; it just hasn't had a second real modality to exercise it.

**Design requirement this surfaces for later:** cross-modal verbs (e.g. "function ~ CST",
"taxon ↔ metabolite") must run on the **shared subset** and report `n_used` /
`n_dropped-for-missing-modality` — the same honest bookkeeping the single-modality verbs
already do, extended across assays. The dashboard should show per-modality **coverage** up
front ("893 have 16S, 122 have shotgun, 122 have both").

---

## Shortlist at a glance

| Dataset | Same-sample 16S + shotgun? | Access | Metadata | Scope |
|---|---|---|---|---|
| **PIN cohort** (Carter 2023) | **Yes — same DNA extraction, by design** | **Fully open** (SRA) | birth outcome, gestational age, race | 72, vaginal, cross-sectional |
| **MOMS-PI** (iHMP) | Partial — paired only on a subset (~122 women) | 16S open; **shotgun dbGaP-gated** | Excellent (CSTs, cytokines, metabolomics, preterm; Black-majority) | ~597 subj, vaginal+, longitudinal |
| **VIRGO2 / VMGC** | **No** — reference catalogs, not a cohort | Open | n/a | aggregated public metagenomes |

---

## 1 — PIN cohort (top pick for the *clean* proof)

- **Citation:** Carter KA, et al. "Vaginal Microbiome Metagenome Inference Accuracy:
  Differential Measurement Error according to Community Composition." *mSystems* (2023).
  Pregnancy, Infection, and Nutrition (PIN) cohort, North Carolina.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10134888/>
- **Same-sample pairing:** **Yes, high confidence.** The *same extracted DNA* from each
  vaginal swab was used for both 16S amplification and whole-metagenome library prep — pairing
  guaranteed at the aliquot level (the study's whole point was comparing 16S-inferred vs true
  metagenome function).
- **Scope:** vaginal; **72 pregnant individuals** (35 early-preterm <32 wk + 37 term controls);
  effectively cross-sectional (one sample/participant at 24–29 wk).
- **Metadata:** birth outcome (early preterm <32 wk vs term 37–41 wk), gestational age at
  collection, race (~63% Black / ~38% White). No cytokines/metabolomics.
- **Where:** NCBI SRA, BioProject **PRJNA876771**; processed tables via the paper's GitHub.
  <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA876771>
- **Access:** **Fully open** — public SRA, no dbGaP application, immediate download.
- **Formats:** raw FASTQ for both 16S and shotgun; supplementary tables on GitHub. No
  pre-built phyloseq (we'd generate ASV/feature tables — fits our ASV-grain pipeline, more
  upfront processing than MOMS-PI's ready phyloseq).
- **Biggest risk / to verify before committing:** n=72 is small; **not yet confirmed in SRA
  that every subject has *both* a 16S run and a shotgun run** — check the PRJNA876771 run table
  (BioSample→run mapping) first. Shotgun depth may be modest (generated for inference-accuracy
  comparison). Metadata beyond birth-outcome/race may need author contact.

## 2 — MOMS-PI (best richness; the *realistic messy-overlap* rehearsal)

- **Citation:** Fettweis JM, Serrano MG, Buck GA, et al. "The vaginal microbiome and preterm
  birth." *Nature Medicine* 25:1012–1021 (2019). VCU Vaginal Microbiome Consortium; part of the
  NIH integrative Human Microbiome Project (iHMP / HMP2).
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6750801/>
- **Same-sample pairing:** **Partial, moderate confidence.** Modalities were generated on
  *overlapping but distinct* subsets, not uniformly. In the preterm case-control set: 16S on
  ~6,452 samples, shotgun on ~496 vaginal, metatranscriptomics ~243, cytokines ~1,223,
  lipidomics 63. "Paired MGS and MTS data were available for 41 women who delivered preterm and
  81 term controls" (~122). A usable same-sample 16S+shotgun subset exists, but you must
  **intersect sample IDs** — it is not the whole cohort. *(This partial overlap is a feature
  for us: it mirrors our own cohort's structure.)*
- **Scope:** vaginal (+ buccal, skin, rectum, infant sites); ~597 consented / ~1,500 women;
  ~12,000 samples; **longitudinal** (~7 visits/participant). Preterm sub-study: 45 sPTB vs 90
  matched term.
- **Metadata (richest of any candidate):** preterm/term, gestational age, race/ethnicity
  (~78% African ancestry — a Black-majority cohort, rare and valuable), CSTs/"vagitypes" (13),
  cytokines/immunoproteomics, metabolomics/lipidomics.
- **Where:** dbGaP **phs001523**; BioProject **PRJNA430482**; HMP DACC portal
  (<https://portal.hmpdacc.org/>); VMC project page (<http://vmc.vcu.edu/momspi>).
- **Access (split):** open-access **16S** (raw/trimmed FASTQ + BIOM OTU tables) from HMP DACC /
  SRA, no application. **Shotgun / metatranscriptomics is controlled-access** via dbGaP
  phs001523 — eRA Commons login, Data Use Certification, IRB/institutional sign-off
  (days–weeks). This is the main friction.
- **Formats (importer-friendly on the 16S side):** 16S raw FASTQ, trimmed FASTQ, and **QIIME
  OTU tables in BIOM**. Bioconductor **`HMP2Data`** builds a ready **phyloseq** object
  (`momspi16S`: 7,665 taxa × 9,107 samples) — directly readable by our phyloseq importer.
  <https://rdrr.io/github/dozmorovlab/HMP2Data/man/momspi16S.html> Shotgun as FASTQ (controlled).

## 3 — VIRGO2 / VMGC (reference resource, *not* a validation cohort)

- **Citation:** VIRGO — Ma B, France M, Ravel J, et al., *Nat Commun* 11:940 (2020); VIRGO2
  update 2025. VMGC — "A multi-kingdom collection of 33,804 reference genomes for the human
  vaginal microbiome," *Nat Microbiology* (2024).
- **Same-sample pairing:** **No.** These are reference catalogs aggregated from thousands of
  public vaginal metagenomes, not a same-sample 16S+shotgun cohort — cannot serve as an
  integration validation set.
- **Use to us:** the **gene/genome reference our shotgun arm maps against** (aligns with the
  existing `import_virgo` / function slot), and a way to trace back constituent cohorts.
- **Where:** VMGC genomes on **Zenodo** (<https://zenodo.org/records/10457006>) + GitHub
  (RChGO/VMGC); VIRGO/VIRGO2 from the Ravel lab. Formats: reference genomes, gene catalogs,
  gene-abundance/annotation tables.

---

## Recommendation & suggested play

1. **PIN** — *best fit for the goal.* The only candidate with **guaranteed same-sample 16S +
   shotgun** AND **fully open access**. Use it for a clean, defensible cross-modal proof.
2. **MOMS-PI** — *best for richness + credibility, and the realistic overlap rehearsal.*
   Unmatched metadata and a ready phyloseq for the 16S arm, but same-sample only on a subset
   and the shotgun modality is dbGaP-gated. Use to demonstrate scale + deep integration (and to
   stress-test reconciliation the way our own cohort will need) once dbGaP clears.
3. **VIRGO2 / VMGC** — reference for the shotgun arm, not a validation cohort.

**Play:** develop and demo integration on **PIN** (open, truly paired, fast) → scale/validate
on **MOMS-PI** for richness and the messy-overlap case once dbGaP access clears. In parallel,
our own cohort is the real target and has MOMS-PI's partial-overlap shape.

**First concrete step when we return:** pull the **PRJNA876771** SRA run table and confirm
every subject has both a 16S and a shotgun run (and check shotgun depth) before wiring PIN as
the second feeder.

---

## Sources

- dbGaP phs001523 — <https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001523.v1.p1>
- Fettweis 2019 (Nat Med, PMC6750801) — <https://pmc.ncbi.nlm.nih.gov/articles/PMC6750801/>
- HMP DACC portal — <https://portal.hmpdacc.org/>
- HMP2Data `momspi16S` phyloseq — <https://rdrr.io/github/dozmorovlab/HMP2Data/man/momspi16S.html>
- Carter 2023 PIN cohort (PMC10134888) — <https://pmc.ncbi.nlm.nih.gov/articles/PMC10134888/>
- BioProject PRJNA876771 — <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA876771>
- VIRGO (Nat Commun) — <https://www.nature.com/articles/s41467-020-14677-3>
- VMGC (Nat Microbiol) + Zenodo — <https://www.nature.com/articles/s41564-024-01751-5> · <https://zenodo.org/records/10457006>
