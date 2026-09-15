# Task brief — Tier 1 fixes (16S correctness)

Two **silent-corruption** bugs from the Run 1 shakedown (Balle PRJEB30774, 6 samples). Both
produce plausible-looking wrong numbers with no error. Fix both, then **re-run the 6 Balle samples
as the acceptance check** (the `.h5mu` / workdir from Run 1 is the test set).

---

## Fix 1 — Genus: look it up from `model.tree`, do NOT parse the string

**Problem.** `_genus_of()` (in `microfgt/io/speciateit.py`, the `collapse_to_taxon` path) uses
`split("_")[0]`. On Balle, 14/76 ASVs got a nonsense genus. Root cause: speciateIT's flat label
overloads `_` three different ways, so **no positional / prefix rule can extract genus reliably.**

**Evidence** (from `~/Projects/microfgt-refdata/speciateIT/vSpeciateDB_models/vSpeciateIT_<R>/model.tree`):
- Rank-tag internal nodes present, standard ranks only: `g_` `f_` `o_` `c_` `p_` `d_`.
- GTDB polyphyly genus suffixes exist: `g_Dorea_A`, `g_Dorea_C`, `g_Eubacterium_F/G/I/J/Q`,
  `g_Clostridium_Q`, `g_Arcanobacterium_A`. So `Dorea_A_longicatena`'s real genus is `Dorea_A`;
  `split("_")[0]` silently merges `Dorea_A` + `Dorea_C` into a fake "Dorea".
- Candidatus: `Ca_Lachnocurva` is a node.

**Fix — a hybrid: tree lookup where it's authoritative, binomial fallback everywhere else.**

> **Corrected 2026-09-03 after inspecting the actual `model.tree` (V4V4).** The original spec
> ("nearest `g_` ancestor, else NA; Candidatus resolves through the tree") does not survive
> contact with the tree. Two measured facts change the rule:
> - **Only 955/1165 leaves have a `g_` ancestor. 210 (18%) do NOT** — not every genus is wrapped
>   in a `g_` clade. That NA set **includes `Gardnerella_vaginalis`** (the CST IV-B driver) and
>   `Ca_Lachnocurva_vaginae`. So "no `g_` ancestor → NA" would regress the single most important
>   vaginal genus — worse than the bug being fixed.
> - **The tree is still necessary**, though: a pure string rule (`rsplit('_',1)[0]`) disagrees
>   with the tree on **47 leaves**, and those are genuinely unresolvable from the string —
>   `Dorea_A_longicatena` = *(Dorea_A)(longicatena)* → `Dorea_A`, but `Aerococcus_urinae_A` =
>   *(Aerococcus)(urinae_A)* → `Aerococcus`. Same shape, different genus boundary; only the tree
>   (`g_Dorea_A` is a node, `g_Aerococcus_urinae` is not) can tell them apart.

Resolution order for a classification label, using the tree from the installed region
(`composition.speciateit.db`/`model.tree`):
1. **Rank-tagged backoff call** (label matches `^[dpcofgs]_`): `g_<X>` → genus `<X>`; a higher
   rank (`d_/p_/c_/o_/f_/…`) → **genus = NA** (honestly unclassified — never a fake genus like "d").
2. **Species leaf with a `g_` ancestor in the tree** → genus = strip `g_` from the nearest such
   ancestor (`Prevotella_bivia`→`Prevotella`; `Dorea_A_longicatena`→`Dorea_A`; `Aerococcus_urinae_A`
   →`Aerococcus`). This is the authoritative case the tree is *for*.
3. **Species leaf with no `g_` ancestor** (the 210, incl. `Gardnerella_vaginalis`,
   `Ca_Lachnocurva_vaginae`) → fall back to the **binomial genus = `label.rsplit('_',1)[0]`**
   (strip the species epithet, keeping any polyphyly suffix / Candidatus prefix): `Gardnerella`,
   `Ca_Lachnocurva`. **Not NA.**

Minimal, near-zero hardcoding: the only fixed list is the seven Linnaean rank codes
`{d,p,c,o,f,g,s}` (stable, not DB-specific). The genus assignments themselves come from the
installed model, so they track the DB version.

**Implementation traps (both hit while validating):**
- **skbio's Newick reader converts `_`→space by default.** You MUST pass
  `TreeNode.read(path, convert_underscores=False)` or every label silently corrupts
  (`g_Prevotella`→`g Prevotella`, so nothing matches).
- **`import_speciateit` does not currently receive the db path.** The tree lives at
  `composition.speciateit.db`; thread it through `_run_import_composition` so the importer can
  build the map. When no db/tree is available (e.g. an ASV-table entry with no model), skip the
  tree and use steps 1+3 only (rank rule + binomial fallback) — still far better than the old bug.

**Verify while building:** a label whose species leaf isn't in the tree at all still resolves via
steps 1/3 (never crashes); if you want drift-detection, log a debug note when a species label has
no tree leaf, but do **not** NA it — the binomial is a valid genus.

**Acceptance (re-run the 6 Balle samples):** the 14 previously-mangled ASVs now carry a correct
genus or NA — specifically `g_Prevotella`→`Prevotella`, `d_Bacteria`→NA, `o_Acetivibrionales`→NA,
`Ca_Lachnocurva_vaginae`→`Ca_Lachnocurva`, and `Gardnerella_vaginalis` stays `Gardnerella` (not
regressed to NA). No genus is a single letter.

---

## Fix 2 — trim_left double-trim trap

**Problem.** The V4V4 region default `trim_left:[19,20]` assumes primers are still on the reads —
but the pipeline runs cutadapt (the `primer_trim` stage) first, which already removes them. So a
user who doesn't override the `dada2` block silently loses 19/20 bp of **real sequence**. Run 1
only dodged it because the mock config set `trim_left:[0,0]` explicitly.

**Fix — make the default honor the actual pipeline contract:**
- Contract: cutadapt removes primers **before** DADA2 *when primers are configured*. So the
  `trim_left` default should be **0 whenever cutadapt actually removed the primers**.
- **Corrected trigger (2026-09-03):** key on **"primers were configured"**
  (`reads.primers.fwd`/`rev` set), NOT merely "the `primer_trim` stage is in the plan." A cutadapt
  run with no primers set is a no-op, and *then* a non-zero `trim_left` is the legitimate way to
  strip primers in DADA2 — so zeroing on stage-presence alone would break that path.
- So: default `trim_left = [0,0]` when `reads.primers.fwd`/`rev` are set; otherwise fall back to
  the region default. An explicit `dada2.trim_left` in the config always wins.
- (Note: in today's DAG `denoise`'s only input, `trimmed_reads`, is produced solely by
  `primer_trim`, so cutadapt always runs before DADA2 — the nonzero region-default `trim_left` is
  effectively always-wrong right now. This fix makes "primers configured" the deciding signal.)

**Acceptance:** with no explicit `dada2.trim_left`, a Balle run yields the same ASV lengths as the
`[0,0]` run — no lost 19/20 bp.

---

## Out of scope (don't scope-creep)
Tier 2 (`--discard-untrimmed` + plumb `extra_args`) and Tier 3 (space-in-path, download retry,
env-on-PATH). Separate tasks.
