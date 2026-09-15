# microFGT prototype — a thin MAE-based spine (exploratory spike)

This directory is a throwaway spike (branch `importers`), built to answer one
question: **if microFGT is built on Bioconductor's `MultiAssayExperiment` (MAE)
instead of a custom container, how thin does the FGT-specific code actually
need to be?** It is exploratory — building to learn, not shipping. Nothing in
the package's `R/` core is deleted; this lives alongside it.

## The idea

The three FGT tools share only the **sample** (a taxon and a gene don't join
directly), and they cover overlapping-but-different sample sets (amplicon broad,
metagenomics a subset). That is exactly what MAE is for. So each tool gets a
small importer to a standard object, and MAE does the cross-assay-by-sample work:

- **speciateIT** → taxonomy (`TreeSummarizedExperiment`, taxa × samples)
- **VIRGO** → gene/function abundance (`SummarizedExperiment`, genes × samples)
- **VALENCIA** → CST label per sample (`colData`)

## Files

| File | What it is |
|---|---|
| `importers.R` | `import_speciateit()` / `import_virgo()` / `import_valencia()` — each maps one tool's real output to a standard Bioc object |
| `demo.R` | end-to-end: generate real-shaped mock → import → assemble MAE → `intersectColumns()` → compare CST-IV* vs CST-I |
| `test_importers.R` | contract test: the **same** importer must parse both the real fixture and the fixed mock output |
| `real_fixtures/` | genuine tool outputs + `FORMATS.md` (the authoritative shapes) + `Valencia.py` |
| `mock/` | generated mock files — gitignored, regenerate locally |

## Run it

```r
# from the repo root
Rscript prototype/test_importers.R   # contract checks (real + mock); exits non-zero on failure
Rscript prototype/demo.R             # full import -> MAE -> payoff walkthrough
```

Requires `MultiAssayExperiment`, `TreeSummarizedExperiment`, `SummarizedExperiment`
(Bioconductor). `mia` is optional — it would add the ecology verbs (alpha/beta
diversity, `agglomerateByRank`) on top of the container, but is not needed here.

## What's validated, and how genuinely

These checks confirm only that each **importer correctly parses the tool's
output format** — *not* that the integrated pipeline has run on real co-assayed
data. It has not (see *Open / not done* below).

| Importer | Validation | Importer vs real output |
|---|---|---|
| **VIRGO** | genuine per-sample `sub*.out` from `ravel-lab/VIRGO/_test_run` | ✅ parsed real output |
| **VALENCIA** | ran the actual `Valencia.py` on 13,231 real samples; resulting `CST` matches the paper's own `Val_CST` for **99.9%** | ✅ parsed real output |
| **speciateIT** | columns confirmed vs the documented spec; ASV→sample join modeled correctly — but **no genuine run** (compiled C++ + trained model DB) | ⚠️ shape-real only |

See `real_fixtures/FORMATS.md` for the authoritative format of each tool and
every way the original mock diverged from reality.

## Key findings

- **The spine is thin.** The three importers are ~40 lines of actual logic; MAE
  gives the container, the `sampleMap`, and `intersectColumns()` for free.
- **The real FGT glue is narrow but real:** recovering speciateIT's sample
  identity from the ASV count table (speciateIT classifies ASVs, *not* samples —
  the sample is not in its output), stacking VIRGO's per-sample files, and
  knowing CST belongs in `colData`.
- **The old mock encoded an *assumed* format with no importer to check it.** This
  spike fixed the mock's writers to emit real shapes and added the round-trip
  test that forces mock and importer to agree — the contract the package lacked.

## Open / not done here

- **speciateIT genuine output** — not run (needs the C++ tool + `vSpeciateDB`).
  Its format risk is low (columns confirmed, join modeled), but it is the one
  importer not checked against real tool output.
- **A real co-assayed dataset** — 16S (→ speciateIT + VALENCIA) and metagenomics
  (→ VIRGO) on the *same* samples — is the true end-to-end test and the real
  prize. None of the per-tool fixtures can exercise the cross-tool integration.
- **The package itself is untouched** — the salvage-vs-rewrite decision and any
  deletion of the old `FGTExperiment` core are deliberately out of scope here.
