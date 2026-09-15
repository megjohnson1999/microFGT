# 16S install runbook (manual, local Mac)

**Purpose.** Perform the 16S-arm install *by hand*, once, on an Apple-Silicon Mac, capturing
every command, URL, version, checksum, and gotcha. This document is the **spec for a future
`microfgt setup`** — nothing gets automated until it has been verified here first.

**Scope.** 16S front-end + taxonomy only: conda env (cutadapt + DADA2) → speciateIT binary →
vSpeciateDB models → `microfgt check` passing for real → a tiny end-to-end run producing ASVs.
Metagenomics is a separate runbook, later.

**Status legend:** 🔴 not started · 🟡 in progress · 🟢 done & verified · ⚠️ gotcha/finding

## Workspace layout
| Thing | Location | In git? |
|---|---|---|
| Code | `~/Projects/microFGT` | yes |
| Conda env | miniforge env `microfgt` (from `environment.yml`) | n/a |
| Reference data (speciateIT clone + models) | `~/Projects/microfgt-refdata/` | **no** |
| This runbook | `design/install_runbook_16s.md` | yes |
| Test working dir (config + tiny input + outputs) | _TBD (scratch / gitignored)_ | no |

---

## Step 1 — conda env  🟢
**Goal:** an env that builds on arm64 where `cutadapt` and R+`dada2` actually run.
**Result:** ✅ lean env at `/Users/Shared/microfgt-16s` — verified native arm64 (no Rosetta),
`cutadapt 5.2`, `dada2 1.38.0` on `R 4.5.3`, `python 3.13.15`. Tools run by name.

- ⚠️ **FINDING (2026-09-02): the bundled `environment.yml` does NOT solve on Apple Silicon.**
  `mamba env create -f environment.yml` fails: `r-randomforestsrc` (a VISTA/metagenomics R
  dep) has **no `osx-arm64` build** — conda-forge ships it for `osx-64` only (confirmed via
  `mamba search`; latest is 2.12.0, an R-4.0 build). One shotgun package thus blocks the whole
  single env, **including a 16S-only user who doesn't need VISTA at all.**
- **Decision:** for the 16S arm, use a lean env — `python>=3.10 pip cutadapt bioconductor-dada2`
  — named `microfgt-16s`. This is all the 16S front-end needs.
  - Command: `mamba create -n microfgt-16s -c conda-forge -c bioconda python">=3.10" pip cutadapt bioconductor-dada2`
  - Then: `conda activate microfgt-16s && pip install -e ".[dev]"` (from the repo)
  - Verify: `cutadapt --version`; `Rscript -e 'library(dada2); packageVersion("dada2")'`
- ⚠️ **FINDING (2026-09-02): spaces in the install path break conda tools.** The home dir is
  `/Users/Megan Johnson/` (space). An env created under it (`~/miniforge3/envs/microfgt-16s`)
  builds fine (native arm64, no Rosetta) but tools fail at runtime: `cutadapt` →
  `/Users/Megan: bad interpreter: No such file or directory` — the script's `#!` shebang points
  at `<env>/bin/python`, and Unix shebang lines split on the space. This is a known conda
  limitation (it does not support spaces in the env prefix). Since microFGT invokes `cutadapt`
  by name as a subprocess, this WOULD break a real run on this machine.
  - **Workaround / decision:** site the env at a space-free prefix →
    `mamba create -p /Users/Shared/microfgt-16s ...`. Refdata dirs may keep the space (data
    files are passed as argv, not exec'd via shebang); only the *env prefix* must be space-free.
  - For users/HTCF this usually won't bite (no space in path), but it belongs in install docs
    as a known gotcha.
- ⚠️ **FINDING: two miniforge installs present** (`/opt/homebrew/Caskroom/miniforge/base` and
  `~/miniforge3`) — `mamba create -n NAME` landed the env under `~/miniforge3/envs` while the
  sourced hook was the Homebrew one, so `conda run -n NAME` from the Homebrew hook couldn't find
  it. Using an explicit `-p <prefix>` sidesteps the ambiguity. Machine-config quirk, not
  microFGT's, but it confused env resolution.
- versions resolved (cutadapt / dada2 / r-base): _fill after space-free build verifies_

## Step 2 — speciateIT binary  🟢
**Goal:** the `classify` binary on PATH, runnable.

- Repo: `git clone --depth 1 https://github.com/ravel-lab/speciateIT` into
  `~/Projects/microfgt-refdata/` (done).
- Binary ships precompiled under `bin/macosx/classify` (single 540 KB executable).
- **Result:** ✅ `file` reports **Mach-O 64-bit arm64** — *native, no Rosetta*. Runs directly,
  prints usage (exit 0). No build needed. Good news for the Mac install story.
- Invocation (from the repo README): `classify -d <models dir> -i <fasta> -o <outdir>`; always
  writes `MC_order7_results.txt`. `--skip-err-thld` forces species-level. Count table via
  `bin/count_table.py -s <MC_order7_results.txt> -c <count table>`.

_Findings:_
- ⚠️ RAM/storage (from repo README): ~1 GB RAM to classify; each region model set is ~2.6 GB
  uncompressed (~12.5 GB for all four).

## Step 3 — vSpeciateDB models  🟢
**Goal:** region-specific model dir downloaded, placed, and proven to classify.
**Result:** ✅ V3V4 downloaded + verified working.
- figshare API (article 25254229): **open, no auth, scriptable** via `ndownloader.figshare.com`.
- **License: CC BY 4.0** → redistribution/mirroring permitted *with attribution* (resolves the
  re-hosting IOU favorably).
- **Latest version is 6** (not the v3 in the speciateIT repo README). Download sizes are the
  *zips*, far smaller than the "2.6 GB" (that's uncompressed): V1V3 142 MB · V1V9 321 MB ·
  V3V4 123 MB · V4V4 89 MB.
- Downloaded `vSpeciateIT_V3V4.zip` → **sha256 `0f6c0e198f0a6fdb5358a6c6c001367a82de10c9d3364bd427eaa303c355225a`** → unzipped to a 2.4 GB model dir
  (`MC0..7.log10cProb`, `model.tree`, `error_thlds.txt`, …).
- **classify runs: ✅** `classify -d vSpeciateDB_models/vSpeciateIT_V3V4 -i test.fasta -o test_run`
  finished in ~4 s and produced sensible vaginal-taxa calls (L. iners/crispatus, Gardnerella…).
- ⚠️ **Reproducibility finding (reinforces pin-the-version):** vs the README's documented
  expected output, **9/10 species match** but pp/nDecisions drift and **ASV9 flips**
  (`Leptotrichia_shahii` → `Sneathia_sanguinegens`). Cause: README example was made on an older
  DB version; we have v6. A dependency version change *did* change a result. → pin the exact
  figshare version + checksum; treat the README values as version-specific, not frozen truth.
- ⚠️ **Format note for `import_speciateit`:** the real `MC_order7_results.txt` has **NO header
  row** — 4 tab-separated cols: `SeqID  Classification  pp  nDecisions`. (README shows a header
  in its illustrative snippet, but the actual file omits it.)

- 🟢 **URL CONFIRMED** (resolves microFGT README's "unverified" TODO): the models live on
  **figshare, DOI `10.6084/m9.figshare.25254229` (v3)** — Holm 2024, "vSpeciateDB Models."
  Source of truth is `speciateIT/vSpeciateDB_models/README.md` in the cloned repo. Gating: TBD
  (figshare is normally open + scriptable via `ndownloader`; confirm on download).
- 🟢 **FIXED (2026-09-02) — region names in the code.** Actual figshare dirs are
  `vSpeciateIT_{V1V3,V1V9,V3V4,V4V4}`, but `REGION_DEFAULTS` used `V1V3 / V3V4 / **V4**` and
  lacked `V1V9`, so a config using the correct model-region name (`V4V4`) lost its DADA2 defaults.
  Renamed `V4`→`V4V4`, added `V1V9`, and kept `V4` as a colloquial alias (`region_defaults()`
  normalizes it). *(The region↔DB `check` itself was fine — its substring match handles the
  canonical names; earlier note overstated it.)*
- Plan: download the **V3V4** set first — it's the region of the repo's documented ground-truth
  test (see below), so it lets us validate end-to-end. Place under
  `~/Projects/microfgt-refdata/speciateIT/vSpeciateDB_models/vSpeciateIT_V3V4/`.
- 🟢 **Ground-truth test available (discharges the speciateIT real-output IOU).** The repo ships
  `test.fasta` (10 vaginal ASVs) with the *expected* `MC_order7_results.txt` printed in its
  README (ASV1→`Lactobacillus_iners` pp 0.970448; ASV2→`L_crispatus`; …). Running
  `classify -d vSpeciateDB_models/vSpeciateIT_V3V4 -i test.fasta -o test` and matching that
  output gives us a genuine `MC_order7_results.txt` to (a) confirm the binary and (b) validate
  `import_speciateit` against **real** classifier output — the exact IOU #1 in the README.

_Findings:_ URL: figshare DOI 10.6084/m9.figshare.25254229 v3 · region chosen: V3V4 (matches the
ground-truth test) · size: ~2.6 GB · checksum: _fill on download_

## Step 4 — `microfgt check` passes for real  🔴
**Goal:** the preflight doctor validates the real install, not stubs.

- Write a minimal `config.yaml` pointing at the real binary + model dir:
  ```yaml
  composition:
    speciateit:
      classify: classify                  # or explicit path
      db: ~/Projects/microfgt-refdata/vSpeciateIT_<REGION>
  ```
- Run: `microfgt check -c config.yaml`
- Verify: binary resolves, `db` path exists, region↔DB match passes. Record any check that
  fires wrongly or misses something — that's a `check` bug to fix (improve-by-exception).

_Findings:_ _fill_

## Step 5 — real-output validation & end-to-end run  🟡

### 5a — `import_speciateit` on GENUINE output  🟢 **(IOU #1 DISCHARGED, locally — no HTCF needed)**
Fed the real `test_run/MC_order7_results.txt` (headerless) + `test_count_table.csv` +
`test.fasta` through `import_speciateit`:
- ✅ auto-detected headerless format; joined 10 classifications, kept the other 1504 count-table
  ASVs as unclassified features; attached real sequences (ASV1 = 429 bp); `collapse_to_taxon` →
  clean 9-taxon roll-up. Object shape (169 samples × 1514 ASVs) as expected.
- 🔴 **Bug confirmed on real data:** `ASV4 Ca_Lachnocurva_vaginae → genus 'Ca'` (Candidatus
  prefix mis-read as genus). Logged as **method_log M2** with a ready fix. The importer docstring
  had predicted exactly this and deferred it to "first real output" — now in hand.

### 5b — full raw-reads front-end via `microfgt run`  🟢 **DONE — plumbing proven end-to-end**
Ran the whole ladder `cutadapt → DADA2 → speciateIT → CST → integrate` from raw FASTQ and wrote a
valid `.h5mu`. **Step 4 done too:** `microfgt check` passes on the real config (all 6 prereqs OK).

- **Input:** mock community `SRR3163904` (V4, 515F/806R, MiSeq, **real quality scores**) +
  `vSpeciateIT_V4V4` models. Config: `~/Projects/microfgt-refdata/mock_run/config.yaml`;
  space-free workdir `/Users/Shared/mock_run_work` (see the speciateIT-space bug below).
- **A PLUMBING test only** — biology is meaningless (non-vaginal mock through vaginal models). The
  point was to prove the pipes, not the science.
- **Result:** object has `composition` (1×54 ASVs) + `composition_taxon` (1×20), a CST call
  (subCST IV-C4 — meaningless for a mock), augment descriptors, reconciliation report, and
  provenance (`tool_runs`). DADA2 produced clean 253 bp V4 ASVs.

**Why not PIN (real vaginal):** PIN's SRA reads carry **no quality scores** (single constant Q30),
so DADA2 can't learn an error model — see method_log **M3**. Real vaginal-with-qualities data is a
validation-phase task, not a plumbing blocker.

**Bugs found & handled here:**
- 🟢 **FIXED — DADA2 single-sample crash.** `dada2_run.R` crashed on n=1 (mergePairs returns a bare
  data.frame → row-name-less table → CSV write fails). Fixed by keying the mergers list by sample id
  *before* `makeSequenceTable`. **This supersedes the earlier uncommitted row-rekey edit — Megan to
  review the diff.**
- ⚠️ **speciateIT `classify` breaks on spaces in paths** (machine artifact). It shells out to
  `mkdir` unquoted, so a workdir under `/Users/Megan Johnson/` splits at the space
  (`mkdir: /Users/Megan: Permission denied`). Workaround: space-free workdir. Robustness idea:
  microFGT could warn when the workdir path contains a space (a bundled 3rd-party binary chokes).

---

## Open IOUs surfaced during this install
_(park empirical/method questions here instead of chasing them mid-install)_
- **VISTA env on Apple Silicon** — `r-randomforestsrc` is `osx-64`-only. When the metagenomics
  runbook happens, VISTA needs one of: a `CONDA_SUBDIR=osx-64` (Rosetta) env, install from CRAN
  via R instead of conda, or a documented "VISTA not supported natively on arm64." Park until
  the metagenomics arm.

## Decisions made
_(anything we settled — e.g. env-splitting, region choice, tool version pins)_
- **Split the conda env.** The single bundled `environment.yml` is not viable on arm64 (see
  Step 1). Move toward per-arm env files (16S core vs VISTA), which `environment.yml` itself
  already anticipated. Concretely for now: a lean `microfgt-16s` env; revisit a proper
  `environment-16s.yml` / `environment-vista.yml` split as a code change later.
- **Dependency durability strategy (deferred to pre-publication).** External deps (figshare
  models, tool repos) can move/version. Standing habit adopted now, cheaply: **pin versioned
  DOIs + record sha256 for every dep as fetched** (done: vSpeciateDB v6, sha256 above). Full
  machinery — a self-owned mirror (vSpeciateDB is CC BY 4.0, so mirroring is permitted) + a
  `microfgt setup` fallback chain (upstream → our mirror → actionable manual message) — is a
  pre-publication task, not built now. Demonstrated need: the v6-vs-older DB drift on ASV9.
