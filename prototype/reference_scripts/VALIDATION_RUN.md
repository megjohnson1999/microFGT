# Metagenomics arm — HTCF validation run (transcription checklist)

The shotgun arm is stub-tested (plumbing proven, biology not). The first real HTCF run
discharges that IOU **and** doubles as validating the distributed conda env. Structure it to be
*transcribable*: everything recorded here is what `fetch-references`, the docs, and a future
container must later encode faithfully — so record it as you go, don't reconstruct it after.

This is a checklist, not a script. The runnable reference lives in the `*_fresh.sh` files here;
the portable command contracts (config-slotted) are in `RECIPE.md`.

## 1. Environment (validates `environment.yml`)
- [ ] `conda env create -f environment.yml && conda activate microfgt` — record whether the
      single-env solve is clean, or whether `r-base` + the 5 VISTA packages had to split into a
      separate `environment-vista.yml` (the one open decision in `environment.yml`).
- [ ] `pip install -e ".[dev]"`.
- [ ] Record exact versions: `fastp --version`, `minimap2 --version`, `samtools --version`,
      `bowtie2 --version`, `R --version`, and each VISTA R package version.

## 2. Reference databases (the "configure + verify, don't provision" half)
Record every DOI, the download URL, and a **sha256** of each downloaded artifact.
- [ ] **GRCh38** host reference + its minimap2/bowtie2 index. Source + sha256 → `metagenomics.host_ref`.
- [ ] **VIRGO2** from **Zenodo DOI 10.5281/zenodo.18703182** (NOT git-lfs — it breaks on HPC into
      133-byte pointer stubs). Record the **sha256 of `VIRGO2.py`** → this is the value that goes in
      `metagenomics.virgo2_sha256` so `microfgt check` catches the GitHub-vs-Zenodo divergence.
- [ ] **Build the VIRGO2 bowtie2 index once** (`Index.tar` on Zenodo is an empty placeholder):
      `VIRGO2.py install` (~2.5 GB). Record the command + the resulting `Index/VIRGO2.*.bt2` sizes.
- [ ] **VISTA** repo + `VISTA_data.tar.gz` from **figshare DOI 10.6084/m9.figshare.28684934**.
      Record the tarball sha256 and the extracted `VISTA_data/` layout (the trimmed bundle the
      container will ship). → `metagenomics.vista_repo`.

## 3. Preflight (validates the doctor)
- [ ] Fill a config from `example_metagenomics_config.yaml` with the real paths + the recorded
      `virgo2_sha256`. Run `microfgt check -c <cfg>` and confirm **every** line is OK. Note any
      slot whose real filename/layout differs from what the `req_fn`s assume (e.g. the exact
      `Index/VIRGO2.1.bt2` name, `VISTA_data/volume`) — those are the corrections to fold back.

## 4. End-to-end on 1–2 samples (discharges the stub IOU)
- [ ] `microfgt run -c <cfg> --workdir <wd>` on a couple of samples. Confirm the `.h5mu` has
      `function` + `composition_taxon_shotgun` modalities, `mgCST` + `mgCST_score` in `.obs`, and
      `shotgun_` descriptors.
- [ ] Diff each stage's real command (in `<wd>/provenance/*.json`) against `RECIPE.md`; record any
      flag/contract mismatch between the shipped orchestrators and the real tools — that's the
      other thing this run is for.

## 5. Resource notes (feed the later Snakemake/Slurm scaling — do NOT build it now)
- [ ] Record peak mem + wallclock per stage vs sample depth. Known so far (RECIPE.md): host
      removal minimap2 peak ~13 GB; VIRGO2 map OOM'd at 32 GB / cleared at 64 GB on the two
      deepest FRESH samples (~175M raw pairs). These become the array-job req hints.
