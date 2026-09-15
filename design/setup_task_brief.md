# Task brief — build `microfgt setup` (16S)

**For an implementation session.** Self-contained. The authority behind every fact here is
`design/install_runbook_16s.md` (the manual install, already verified by hand). Do not re-derive
what the runbook already nailed down; if something here is ambiguous, the runbook wins.

## Goal (one sentence)
Add a `microfgt setup --region <R> --dest <dir>` CLI command that installs the two 16S
prerequisites conda can't provide — the **speciateIT `classify` binary** and the **vSpeciateDB
model for `<R>`** — then writes the config pointing at them and runs `microfgt check`.

## Done when
`microfgt setup --region V3V4 --dest <a fresh empty dir>` completes and, on a machine that has
the `microfgt-16s` conda env active, a `microfgt check` on the config it wrote reports **all
prereqs OK** (binary resolves, `db` path exists, region↔DB match passes). Verify against the
repo's ground-truth test: `classify -d <model dir> -i speciateIT/test.fasta -o <out>` reproduces
the README's expected calls (ASV1 → `Lactobacillus_iners`, etc.; 9/10 species — ASV9 legitimately
differs on v6, that's expected and documented).

## What it must do, in order
1. **Resolve `--region`** to a canonical model-region name: one of `V1V3 | V1V9 | V3V4 | V4V4`.
   Accept `V4` as a colloquial alias for `V4V4` (the code already normalizes this in
   `region_defaults()`).
2. **Fetch the speciateIT binary.** `git clone --depth 1 https://github.com/ravel-lab/speciateIT`
   into `<dest>`. The precompiled binary ships in the repo — pick the build for the current OS:
   `bin/macosx/classify` on macOS (verified Mach-O arm64, native, no Rosetta, ~540 KB). (Linux/HTCF
   build selection is a follow-on — see "Out of scope / flag to Megan".)
3. **Download the vSpeciateDB model zip for `<R>`** from figshare — open, no auth, scriptable:
   - Article DOI `10.6084/m9.figshare.25254229`, **pin version 6**. Resolve files at runtime via
     the figshare API (`https://api.figshare.com/v2/articles/25254229/versions/6`), match the file
     named `vSpeciateIT_<R>.zip`, download its `download_url` (an `ndownloader.figshare.com` link).
     Resolve by version+filename — do **not** hardcode figshare file IDs (they change across
     versions).
   - Zip sizes (for a pre-download size/disk message): V1V3 142 MB · V1V9 321 MB · V3V4 123 MB ·
     V4V4 89 MB. Each unpacks to ~2.4–2.6 GB. Classifying needs ~1 GB RAM.
4. **Verify integrity, then unzip.** sha256-check the downloaded zip against the pinned manifest
   below before unzipping. Unzip to `<dest>/vSpeciateDB_models/vSpeciateIT_<R>/` (contains
   `MC0..7.log10cProb`, `model.tree`, `error_thlds.txt`, …).
5. **Write (and print) the config lines**, pointing at what it just installed:
   ```yaml
   composition:
     speciateit:
       classify: <resolved path to the classify binary>   # or "classify" if it's on PATH
       db: <dest>/vSpeciateDB_models/vSpeciateIT_<R>
   ```
   Support `--dest` being anywhere (a laptop dir, or a shared `/ref/...` dir on a cluster).
6. **Run `microfgt check`** on that config and surface its result as the final word.

## Pinned integrity manifest (sha256 of the zips)
Only V3V4 has been fetched + recorded by hand so far. The other three are **not yet pinned** —
capture each one's sha256 on first successful download and write it back here.
```
V3V4: 0f6c0e198f0a6fdb5358a6c6c001367a82de10c9d3364bd427eaa303c355225a
V1V3: db414c7327e72ca11974a8ca49169896a8a03a49c6dd82a5ec5c38ea1e7172b0   # pinned 2026-09-03 (figshare v6, md5-cross-checked)
V1V9: <unpinned — record on first download>   # figshare file is dated: vSpeciateIT_V1V9_14May2026.zip
V4V4: e38171df1197560fdee0f95c026397c5e955ddcde01fe72bbec27221358c1101   # pinned 2026-09-03 (figshare v6, downloaded + md5-cross-checked)
```
Because only V3V4 is pinned: build + validate the whole command against **V3V4 first**. For an
unpinned region, don't silently skip verification — record the computed sha256 and warn that it's
being trusted-on-first-use, or refuse until it's pinned (your call — see flag to Megan).

## Gotchas the runbook already found (handle, don't rediscover)
- **speciateIT `classify` breaks on spaces in a path** — it shells out to `mkdir` unquoted, so a
  path containing a space splits. On this machine `$HOME` is `/Users/Megan Johnson/` (space).
  `setup` should at least **warn** if `--dest` contains a space; a space-free `--dest` avoids it.
- v6 model output legitimately drifts from the speciateIT README's older-version example (ASV9
  flips). That's why we pin the version + checksum — treat README values as version-specific.

## Decisions already made (don't reopen)
- Env is the user's job; `setup` owns tools + DBs. The lean env is `environment-16s.yml` (done).
- License is CC BY 4.0 → mirroring is permitted, but a self-owned mirror + fallback chain is a
  **pre-publication task, not this one**. This command talks to upstream figshare only.
- One region per `setup` run.

## Out of scope / flag back to Megan (this session)
- **Non-macOS binary selection** (Linux for HTCF) — ✅ wired. `_OS_BINARY` maps `Darwin →
  bin/macosx/classify` and `Linux → bin/linux/classify` (the linux build ships in the repo too).
  Darwin is in `_OS_VALIDATED`; Linux is selected but prints a "not yet validated on this OS"
  note — **cluster validation on HTCF is still a separate follow-on.**
- **Unpinned-region policy** — ✅ **decided: trust-on-first-use + record (default)**, with
  `--require-pinned` to refuse-until-pinned. Rationale: refusing would make setup useless for any
  region not hand-pinned. TOFU computes the sha256, prints it for pasting into the manifest, and
  proceeds (figshare version is pinned to v6, so content is fixed). Reopen if you'd rather it hard-fail.
- **Mirror/fallback chain** — explicitly deferred (above).

### New findings this session (2026-09-03)
- ⚠️ **V1V9's figshare file is `vSpeciateIT_V1V9_14May2026.zip`, not `vSpeciateIT_V1V9.zip`.** An
  exact-name match would silently miss it. `resolve_figshare_file` matches on the region token
  bounded by `.`/`_`, so the dated name resolves and V1V3 never swallows V1V9.
- **V1V3 + V4V4 now pinned** (see manifest). V4V4 was downloaded fresh (89 MB, byte-identical to
  figshare's md5); V1V3's sha was taken from the local v6 zip and md5-cross-checked. V1V9 remains
  unpinned (not fetched — it's the 321 MB dated file).
