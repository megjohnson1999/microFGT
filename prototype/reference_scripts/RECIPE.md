# Shotgun pipeline recipe (portable extraction)

The `*.sh` files in this dir are the **reference implementation** run on HTCF
(the FRESH cohort). They are **not** what microFGT ships — they carry
cluster-specific paths, conda envs, and SLURM headers a user won't have.

This file extracts the **portable recipe**: the tool + flags + ordering + gotchas,
with every HTCF path replaced by a `{config.*}` slot and every `conda activate`
replaced by "the tool, found on PATH." The orchestration stages should be built
from *this*, not from the raw scripts.

---

## What becomes user config (not hardcoded paths)

The reference scripts hardcode these; microFGT turns each into a config slot the
`microfgt check` doctor verifies (the "unopinionated about data location" half):

| Reference script path | microFGT config slot | Kind | Doctor check |
|---|---|---|---|
| `/ref/sahlab/data/GRCh38.fna.gz` | `config.host_ref` | host-removal DB | file exists & readable |
| `.../resources/VIRGO2` (has `VIRGO2.py` + `Index/VIRGO2.*.bt2`) | `config.virgo2_dir` | gene-catalog DB + index | `Index/VIRGO2.1.bt2` present |
| `.../resources/VISTA` (contains `VISTA_data/`) | `config.vista_repo` | classifier model bundle | `VISTA_data/volume/` present |

## What becomes owned tools (not conda envs)

The reference scripts `conda activate` lab envs; microFGT owns these via its own
conda/bioconda recipe (the "opinionated about software" half). The env *paths*
are irrelevant — only the tool + version matters:

- `fastp`, `minimap2`, `samtools`  (QC + host removal)
- `bowtie2`  (used internally by `VIRGO2.py map`)
- `python3` + `pandas` + `numpy`  (for `VIRGO2.py`)
- `Rscript` + `randomForestSRC`, `pheatmap`, `dplyr`, `data.table`, `R.utils`  (for VISTA)

Everything below refers to tools **on PATH** — no env path appears.

---

## Stage recipe (paths abstracted)

`N` = threads. `R1`/`R2` = input pair for the sample.

### sg_qc  (fastp — default params: auto PE adapter + quality trim)
```
fastp --in1 R1 --in2 R2 --out1 TRIM_R1 --out2 TRIM_R2 \
      --thread N --json qc.json --html qc.html
```
Output to keep for the object: `qc.json` (per-sample read/quality metrics).

### sg_host_removal  (minimap2 vs host ref -> keep read pairs where BOTH mates are unmapped)
```
minimap2 -ax sr -t N {config.host_ref} TRIM_R1 TRIM_R2 \
  | samtools view  -@ N -b -f 12 -F 256 - \   # -f 12: read+mate both unmapped; -F 256: drop secondary
  | samtools sort  -n -@ N -m 2G -T sorttmp - \   # name-sort before fastq
  | samtools fastq -@ N -1 NONHOST_R1 -2 NONHOST_R2 -0 /dev/null -s /dev/null -
```
The `-f 12 -F 256` filter *is* the host-removal logic. Vaginal samples are very
host-heavy (FRESH: ~86% of reads removed), so most depth is spent here.

### sg_virgo2_map  (VIRGO2 is SINGLE-END only -> concatenate R1+R2 first)
```
cat NONHOST_R1 NONHOST_R2 > combined.fq.gz            # <-- required: no -1/-2 path in VIRGO2
python3 {config.virgo2_dir}/VIRGO2.py map -r combined.fq.gz -o OUTDIR/<id> -p N
# -> OUTDIR/<id>.out   (per-sample gene counts)
```

### sg_virgo2_compile  (all per-sample .out -> one gene x sample matrix)
```
python3 {config.virgo2_dir}/VIRGO2.py compile -i OUTDIR -o OUTDIR/VIRGO2_Compiled
# -> OUTDIR/VIRGO2_Compiled.summary.NR.txt
```
`compile` globs `*.out` in the dir — it is sample-list-agnostic (this is why the
9-of-11 recovery "just worked" without editing anything).

### classify_mgcst  (VISTA)
```
cd OUTDIR_VISTA                                        # VISTA writes outputs to CWD (timestamped)
Rscript {config.vista_repo}/run_VISTA.R  VIRGO2_Compiled.summary.NR.txt  {config.vista_repo}
```
Outputs (six files): `mgCSTs_*.csv` (call + max_YC_theta), `norm_counts_mgSs_mgCST_*.csv`
(subtype), `norm_counts_taxa_*.csv`, `norm_counts_genes_*.csv`, `relabund_w_mgCSTs_*.csv`,
`mgCST_heatmap_*.pdf`.

---

## Gotchas the orchestration must encode (these are the "pain points")

1. **VIRGO2 is single-end only** — `map` takes one `-r` file, no `-1/-2`. Must
   `cat R1 R2` first. (Ravel-lab convention; mate info is unused by a gene catalog.)
2. **VISTA arg2 = the dir CONTAINING `VISTA_data/`**, not `VISTA_data` itself.
   `run_VISTA.R` internally prepends `/VISTA_data/...`. The README's `./VISTA_data`
   example fails.
3. **VISTA writes to the current working directory** (timestamped filenames) — the
   stage must `cd` into a dedicated output dir, not rely on an `-o` flag.
4. **Host removal keeps only pairs where BOTH mates are unmapped** (`-f 12`), and
   must name-sort (`samtools sort -n`) before `samtools fastq`.
5. **Taxonomy is not a VIRGO2 output file** — it's the gene->taxon annotation
   (`AnnotationTables/1.VIRGO2.taxon.txt`) aggregated onto the compiled gene matrix.
   (See `prototype/real_fixtures/` and the fixture correction note.)

## Resource notes (for the doctor / stage req hints, not user config)

- Host removal: minimap2 peak ~13 GB on the deepest FRESH sample.
- VIRGO2 map memory scales with clean depth (per-read gene-assignment dict). The
  two deepest FRESH samples (~175M raw pairs) OOM'd at 32 GB; 64 GB cleared them.
- These are HPC-scaling notes; on a single machine they translate to "map memory
  grows with sample depth," worth surfacing rather than hardcoding.
