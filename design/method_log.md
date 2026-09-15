# Method interrogation log (both arms)

**Purpose.** Running log of *scientific/method* questions about the tools microFGT wraps —
"is this the right method, and are we using it correctly?" — as opposed to install/plumbing
snags (those live in the install runbooks, e.g. `install_runbook_16s.md`).

**Philosophy (agreed).** *Faithful by default, improve by exception, with receipts.* We run the
field-standard tools as-is unless a concrete, demonstrated problem justifies deviating — and
when it does, we document the evidence and the change here. Findings are triggered by evidence
(usually a real run), not by ambition.

**Status legend:** 🔴 open · 🟡 leaning · 🟢 decided/disposed

**Entry template:**
```
### M<n> — <title>
Arm / stage · Status
Issue:          <what, and why it matters for a field tool>
In code:        <file:line — current behavior>
Options:        <the menu>
Evidence needed: <what would let us decide>
Disposition:    <decision, or parked-pending-X>
```

---

### M1 — DADA2 error model on binned quality scores (NovaSeq/NextSeq)  🟡
**Arm/stage:** 16S · denoise

**Issue.** DADA2's `learnErrors` fits an error model *as a function of* base quality score.
Older Illumina instruments (MiSeq/HiSeq) emit ~40 distinct Q-scores; NovaSeq and newer 2-color
NextSeq chemistry emit only **~4 binned values** (e.g. 2/12/23/37). With so few bins the fitted
error rates become non-monotonic and poorly behaved, degrading ASV inference. DADA2 is not
"broken" for newer data — but its **default error-learning step needs adjustment** for binned
quality. This matters because the field is trending toward NovaSeq for amplicon data, so it's a
forward-looking concern, not a corner case.

**In code.** `microfgt/scripts/dada2_run.R:57-58` calls `learnErrors(filtFs, multithread=TRUE)`
with all defaults — the classic error model, no binned-quality handling, and no knob exposed.
As wired today, microFGT's denoise front-end would hit this on NovaSeq input.

**Mitigating factor.** The multi-entry design is an escape hatch: a user who denoises elsewhere
(or wants binned handling) enters at the ASV table and skips microFGT's denoise stage entirely.
The risk only bites users who use microFGT's *own* front-end from raw reads.

**Options.**
1. Expose `errorEstimationFunction` as a config knob (default classic; monotone variant opt-in).
2. Auto-detect binned quality (few distinct Q-values in the quality profile we already emit) and
   switch to a monotone-enforcing error function automatically.
3. Document the limitation and point NovaSeq users at the ASV-table entry point.

**Evidence needed.** A real NovaSeq/NextSeq FGT amplicon dataset run through the denoise stage,
compared against a monotone-error-function run, to confirm the effect and validate the fix.

**Disposition.** 🟡 Leaning toward option 2 (auto-detect via the quality profile) + option 1 as
the manual override — "improve by exception, with receipts." Parked pending a real binned-quality
dataset to test against. Not fixing blind.

---

### M2 — `_genus_of` mis-parses Candidatus / abbreviated names  🟢 (confirmed, fix ready)
**Arm/stage:** 16S · import (speciateIT → composition)

**Issue.** `microfgt/io/speciateit.py:_genus_of` derives genus as the first `_`-split token of a
speciateIT Classification. For `Candidatus` taxa the label is `Ca_Lachnocurva_vaginae`, so the
"genus" comes out as **`Ca`** — the Candidatus prefix, not the genus (`Lachnocurva`). Same class
of error for any abbreviated leading token. Genus is used in the taxon roll-up's `var['genus']`,
so this propagates a wrong genus into `composition_taxon`.

**Confirmed on real data (2026-09-02).** Feeding *genuine* speciateIT output (`test.fasta` →
`MC_order7_results.txt` on the real V3V4 v6 models) through `import_speciateit` produced
`ASV4 Ca_Lachnocurva_vaginae → genus 'Ca'`. The code's own docstring predicted exactly this and
deferred the fix to "first real speciateIT output" — which we now have (locally, no HTCF needed).

**Fix.** When the first token is a known abbreviation prefix (`Ca` = Candidatus; also handle
`Ca.`/`Candidatus`), take the next token as genus (`Lachnocurva`). Add `Ca_Lachnocurva_vaginae`
to the importer's unit fixtures so it regresses. Low-risk, well-scoped.

**Disposition.** 🟢 Real-output confirmed; straightforward fix + regression fixture. Do it as part
of the speciateIT importer cleanup (also fixes the region-name bug noted in the install runbook).

---

### M3 — detect quality-less / degenerate-quality input before DADA2  🟡
**Arm/stage:** 16S · denoise

**Issue.** DADA2's error model is a function of base quality, so it fails on inputs whose quality
information is degraded — the extreme case being reads stored with a **single constant quality
value**. When that happens DADA2 dies cryptically ("Error rates could not be estimated … Error
matrix is NULL"), and a user can't tell it's a *data* problem, not a tool bug.

**Confirmed on real data (2026-09-02).** PIN cohort SRA reads (PRJNA876771) have exactly **1
distinct quality character** (all Q30) — SRA stored them without real qualities. DADA2 failed on
them. This is common enough in SRA amplicon deposits to matter for a raw-reads tool.

**Fix (improve-by-exception).** A cheap preflight in the denoise stage: count distinct quality
values in the input; if ~1 (or very few), fail fast with a clear message — e.g. "these reads carry
no usable quality scores; DADA2 requires them. If from SRA, the qualities were likely stripped —
supply reads with real qualities, or enter at an ASV table instead." Related to [[M1]]: M1 =
*handle* few quality bins in the model; M3 = *detect the pathological case and message clearly*.

**Disposition.** 🟡 Clear win; do alongside the M1 work. (The multi-entry design already offers the
escape hatch — such a user can enter at an ASV table and skip denoise.)

---

### M4 — CST centroids (2020) are naming-incompatible with speciateIT v6 → silent CST degradation  🟢 (FIXED 2026-09-02)
**Arm/stage:** 16S · CST (VALENCIA)

**Issue.** VALENCIA classifies by **matching taxon names** between a sample and 13 reference
centroids (Yue–Clayton theta over the taxon union; unmatched taxa are zero-filled). So the sample's
taxonomy naming MUST agree with the centroid file's naming. It currently does not: we bundle the
**2020** centroids (`cst_centroids_012920.csv`, 199 taxa, old names) but assign taxonomy with
**speciateIT v6 (2024)**, whose names have moved on.

**Confirmed on real data (2026-09-02).** Of the 8 real vaginal taxa speciateIT v6 emitted, **5 are
absent from the 2020 centroids but present in the 2024 VALENCIA2 centroids**:
`Ca_Lachnocurva_vaginae` (2020: `BVAB1`), `Fannyhessea_vaginae` (2020: `Atopobium_vaginae`),
`Lactobacillus_mulieris`, `Megasphaera_lornae`, `Sneathia_sanguinegens`. These are **core BV
organisms** — so for CST IV (dysbiotic) communities, much of the sample's abundance is invisible to
the centroid comparison and the CST call is silently wrong. Lactobacillus-dominated samples (whose
taxa match) are unaffected. Damage is concentrated in exactly the clinically important communities.

**Fix.** Adopt the **VALENCIA2 2024 centroids** (ship in the speciateIT repo:
`VALENCIA2_CST_centroids_19Aug2024.csv`, 356 taxa, 13 subCSTs — same nearest-centroid structure).
`load_reference_centroids(reference=...)` already supports swapping the file, so the code change is
tiny: replace the bundled centroid CSV.

**Catch (must not skip).** The headline validation ("99.94% subCST agreement, 13,231 samples") was
against **2020** VALENCIA output. Swapping centroids changes the CST *version*, so that number does
NOT transfer. VALENCIA2 needs its own validation: confirm it's the same Yue–Clayton nearest-centroid
method (structure suggests yes) and re-run the gate against genuine VALENCIA2 output. Related:
[[M2]]/[[M3]] and the general taxon-namespace problem (16S vs shotgun naming won't line up either).

**Disposition.** 🔴 Open, high priority — it degrades the flagship CST output.

**Compatibility CONFIRMED (2026-09-02): it's a true drop-in.**
- The 2024 centroids have the **exact same 13 subCST labels, order, and relative-abundance format**
  as the 2020 file — `classify_cst` (and `CST_ORDER`/`_COLLAPSE`) run unchanged; only the bundled
  CSV changes.
- There is **no separate "VALENCIA2" algorithm.** VALENCIA is one nearest-centroid method (France
  2020); the official repo still ships the 2020 centroids. The `VALENCIA2_CST_centroids_2024` files
  (in the speciateIT repo) are just updated *centroids* with modern taxonomy — same method.

**Sequence to fix:**
1. Swap the bundled centroid CSV → `VALENCIA2_CST_centroids_19Aug2024.csv` (behind
   `load_reference_centroids`; keep 2020 available for reproducing the original paper).
2. Re-validate: VALENCIA's script takes a centroids file, so generate genuine VALENCIA output
   *using the 2024 centroids* and re-run the CST gate against that — preserving the
   "reproduces the real tool exactly" guarantee for the new pairing.
3. Note in docs: centroid set is coupled to the assigner's naming (speciateIT v6 ↔ 2024 centroids);
   a different assigner (SILVA/GTDB) would need its own mapping. Matched-pair principle.

**DONE (2026-09-02).** Steps 1 done, 2 (validation) done at the port-fidelity level:
- `microfgt/data/VALENCIA2_CST_centroids_19Aug2024.csv` added; `load_reference_centroids` now takes
  a set name (`"2024"` default / `"2020"`) or a path. Default is **2024**.
- The two validation tests are pinned to `reference="2020"` (their gold-standard data is 2020-named)
  and still pass — so the port's fidelity to VALENCIA's *math* is unchanged. Added a regression
  test asserting the default is 2024 and includes the BV taxa the 2020 set misses.
- Full suite: **161 passed.** Before/after on a synthetic BV sample: 2020 centroids "saw" only 50%
  of the community (score 0.43); 2024 saw 100% (score 0.68).
- **NOT committed to git yet — includes changes that supersede the uncommitted dada2 rekey edit;
  Megan to review the diff.**

**Remaining follow-ups (validation-phase, not blockers):**
- Belt-and-suspenders: run the official `Valencia.py` *with* the 2024 centroids on a few samples and
  confirm our port reproduces it exactly (our math is centroid-agnostic + the file format is
  identical, so this is confirmation, not a risk).
- Real gold-standard: a 2024-named vaginal cohort with genuine VALENCIA-2024 CST labels, to re-earn
  the 99.9% gate for the 2024 pairing.
- Docs: update the README's "validated 99.94%" line to note it is the 2020-pairing port validation,
  and that the default centroids are now 2024.
