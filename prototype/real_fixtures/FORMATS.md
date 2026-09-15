# Real tool output formats — ground truth for importers + mock fix

Fixtures pulled from the tools' own repos (ravel-lab) on 2026-06-08. These are
the authoritative shapes the importers must parse and the mock must emit.
Large catalog/study files were trimmed to format samples (provenance below);
`virgo_sub1.out` / `virgo_sub2.out` are full, real per-sample outputs.

## speciateIT — github.com/ravel-lab/speciateIT
- **Classifies ASVs/sequences, NOT samples.** Input is a FASTA (`test.fasta`,
  headers `>ASV1`, `>ASV2`, …). Output `MC_order7_results.txt` is TSV, one row
  per sequence, keyed by the FASTA header:
  `Sequence ID \t Classification \t posterior probability \t number of Decisions`
- **Sample identity is NOT in this file.** It lives in the ASV count table
  (`test_count_table.csv`: rows = `sampleID`, cols = `ASV1..ASVn`, 169 samples ×
  1514 ASVs here). taxon×sample = join(ASV→Classification, ASV×sample counts).

### Divergence from current mock
| mock today | real |
|---|---|
| sample baked into seq ID `Sample_001_Seq_000072` | seq ID is the ASV id (`ASV1`); no sample in it |
| one combined `MC_order7_results.txt`, sample parsed from prefix | one output per FASTA; sample comes from the **separate ASV count table** |
| column format (4 cols) | ✅ already correct |

## VIRGO — github.com/ravel-lab/VIRGO
- **One output file per sample**, named `<sample>.out` (sample = filename).
  TSV, **no header, 3 columns**: `geneID \t read_count \t gene_length`
  (`virgo_sub1.out`: `V1593031  1417  3663`, 3294 genes).
- Catalog annotation (repo `1_VIRGO/`): `0.geneLength.txt` = `geneID\tlength`;
  `1.taxon.tbl.txt` = `Cluster\tgeneID\ttaxon\tlength`. Gene→taxon/length map.

### Divergence from current mock
| mock today | real |
|---|---|
| one `*_test.out` for all samples, **4 cols incl. a `Sample` column**, with header | **one file per sample**, **3 cols, no header, no Sample column** |
| `_counts.tsv` wide matrix as a deliverable | wide matrix is something *you* build by stacking per-sample `.out` files |

## VIRGO2 — github.com/Ravel-Laboratory/VIRGO2  (audit: ENA/PRJEB34536, 2026-08)
> Fixtures are slices from **public ENA cohort PRJEB34536** (6 samples). FRESH is
> controlled human-subjects data → no FRESH files in the repo.

- **`VIRGO2.py compile` writes exactly ONE wide matrix**, `VIRGO2_Compiled.summary.NR.txt`.
  TSV, header `Gene \t <sample> \t <sample> …`, one row per gene, **float** counts
  (fixture `virgo2_compiled.summary.NR.slice.txt`, 30 genes × 6 ERR samples). This is the
  opposite orientation to v1's per-sample files, and unlike v1 it carries **no annotation
  columns**.
- **All annotations are separate `AnnotationTables/` files, joined on `Gene` (not `Cluster`).**
  - taxon: `1.VIRGO2.taxon.txt` = `Cluster \t Gene \t Taxa \t Cat`
    (fixture `virgo2_taxon_annotation.slice.txt`). Join **on `Gene`**; `Taxa` is the label.
  - functional (KEGG shown): `Gene \t KEGG \t KEGG_Pathway \t KEGG_Module \t …`
    (fixture `virgo2_kegg_annotation.slice.txt`). Also EC / PFAM / eggNog / CAZy / AMR / phage /
    geneProduct / geneLength / VOGkey. **Partial coverage is normal** — the KEGG slice annotates
    20/30 genes → the join must tolerate missing (→ NaN / `Unannotated`).
- **Shotgun taxon composition is DERIVED, not emitted.** There is no `VIRGO2.py taxonomy` output
  file; microFGT builds taxon×sample itself by joining the gene matrix → taxon table and summing
  per `Taxa` (`import_virgo2` → `collapse_virgo2_to_taxon`), mirroring the 16S
  `import_speciateit` → `collapse_to_taxon` split. This decouples taxonomy from running VISTA.

## VISTA / mgCST — the shotgun community-type call  (audit: ENA/PRJEB34536, 2026-08)
- **`run_VISTA.R` writes six files**: `mgCSTs_*.csv`, `norm_counts_mgSs_mgCST_*.csv`,
  `norm_counts_taxa_*.csv`, `norm_counts_genes_*.csv`, `relabund_w_mgCSTs_*.csv`,
  `mgCST_heatmap_*.pdf`.
- **The authoritative per-sample call is `mgCSTs_*.csv`** (fixture `vista_mgCSTs.csv`):
  CSV, first (unnamed) column = sample id, then `mgCST` (label, e.g. `"mgCST 11"`) and
  `max_YC_theta` (YC-θ of the **best-matching** mgCST). `import_mgcst()` parses this →
  `mgCST` + `mgCST_score`.
- **No per-centroid similarities** — VISTA emits only `max_YC_theta` (best match), not θ against
  all 25 centroids. So (unlike CST's `<subCST>_sim`) there is **no `mgcst_sim` block** to route
  to `.obsm`; a low θ means "matches no reference type well," surfaced on the sample.
- **No scalar subtype in the call file.** The finer mgSs level lives in
  `norm_counts_mgSs_mgCST_*.csv` (fixture `vista_norm_counts_mgSs_mgCST.csv`) as a **feature
  matrix** (mgSs × sample), not a per-sample label — an mgSs modality is deferred.

## VALENCIA — github.com/ravel-lab/VALENCIA  (authoritative: `Valencia.py`)
> **Validated against GENUINE tool output.** Ran `Valencia.py` on the repo's real
> published composition data (13,231 samples × 212 taxa); resulting `CST` matches
> the paper's own `Val_CST` for **99.9%** of samples. `import_valencia()` parses
> that genuine output. Fixture: `valencia_genuine_output_head.csv` (head of the run).

- Input CSV: col1 `sampleID`, col2 `read_count`, then **one column per taxon**
  (taxon name = header, cells = counts).
- Output CSV = **the input, plus appended columns** (one wide file):
  - 13 per-subCST similarity columns `<subCST>_sim`
    (`I-A_sim`, `I-B_sim`, `II_sim`, `III-A_sim`, `III-B_sim`, `IV-A_sim`,
     `IV-B_sim`, `IV-C0_sim`…`IV-C4_sim`, `V_sim`)  [Valencia.py:125–132]
  - `subCST` = argmax of the 13 `_sim` columns (`_sim` stripped)
  - `score`  = max of the 13 `_sim` values
  - `CST`    = `subCST` collapsed (I-A/I-B→I, III-A/III-B→III, IV-C0..4→IV-C) [:135]

### Divergence from current mock
| mock today | real |
|---|---|
| **three files** `_cst.csv` / `_scores.csv` / `_abundance.csv` | **one wide CSV** (input + appended cols) |
| sample column `Sample` | `sampleID` |
| 7 score cols `CST-I`…`CST-V` | 13 cols `<subCST>_sim` |
| CST only (I, II, III, IV-A/B/C, V) | finer `subCST` (I-A/B, III-A/B, IV-C0–4) that collapses to `CST` |

## Cross-cutting consequence for the importers
- `import_speciateit()` must return **ASV→taxonomy** and take an **ASV count
  table** (the dada2/feature table) to produce taxon×sample. The mock's
  "sample-in-the-ID" shortcut must go.
- `import_virgo()` must read a **directory of per-sample `<sample>.out`** files
  (3 cols, no header) and stack them, not one combined file with a Sample col.
- `import_valencia()` must read **one wide CSV** and pull `CST`/`subCST`/`score`
  from the trailing columns, keyed by `sampleID`.

## Provenance (full files, if needed)
- speciateIT: `test.fasta`, `test_count_table.csv` (repo root)
- VIRGO: `_test_run/temp_mapping/sub1.out`, `sub2.out`; catalog `1_VIRGO/0.geneLength.txt`, `1_VIRGO/1.taxon.tbl.txt`
- VIRGO2 (ENA PRJEB34536): `VIRGO2_Compiled.summary.NR.txt`; `AnnotationTables/1.VIRGO2.taxon.txt`, `.../3.VIRGO2.kegg.txt`
- VISTA (ENA PRJEB34536): `mgCSTs_*.csv`, `norm_counts_mgSs_mgCST_*.csv` (from `run_VISTA.R`)
- VALENCIA: `Valencia.py`, `README.md`, `Publication_materials/Data_and_metadata/all_samples_taxonomic_composition_data.csv`
