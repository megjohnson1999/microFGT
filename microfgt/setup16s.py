"""``microfgt setup`` — install the two 16S prerequisites conda can't provide.

Conda gives you cutadapt + DADA2 (the ``microfgt-16s`` env is the user's job). What it *can't*
give you is the speciateIT ``classify`` binary and the region-specific vSpeciateDB model. This
module fetches both, writes a config that points at them, and runs ``microfgt check`` as the
final word. The authority behind every constant/URL/checksum here is
``design/install_runbook_16s.md`` — the manual install, verified by hand first.

Design notes / decisions (see the runbook + the task brief):
- **One region per run.** ``setup(region, dest)`` installs exactly one model.
- **Pin the figshare version + checksum.** The models live at DOI 10.6084/m9.figshare.25254229;
  we pin **version 6** and resolve files via the figshare API at runtime (never hardcode file
  ids — they change across versions). The zip's sha256 is checked against a pinned manifest.
- **Unpinned regions: trust-on-first-use + record (with a loud warning).** Only V3V4 was pinned
  by hand. Refusing every other region would make ``setup`` useless for V1V3/V1V9/V4V4 until
  someone hand-pins them, so the default is TOFU: verify the region+version resolution, compute
  the sha256, print it prominently to be pasted into the manifest, and proceed. ``require_pinned``
  flips this to refuse-until-pinned. (Flagged to Megan in the brief — this is the pick.)
- **speciateIT ``classify`` splits on spaces in a path** (it shells out to ``mkdir`` unquoted).
  ``setup`` warns when ``--dest`` contains a space; a space-free ``--dest`` avoids it.
"""

from __future__ import annotations

import hashlib
import json
import platform
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

from microfgt.stages.registry import canonical_region

SPECIATEIT_REPO = "https://github.com/ravel-lab/speciateIT"

# figshare: pin the article + version; resolve files (and their download_urls) at runtime.
FIGSHARE_ARTICLE = "25254229"
FIGSHARE_VERSION = 6
FIGSHARE_API = "https://api.figshare.com/v2/articles/{article}/versions/{version}"

# OS -> the precompiled binary path inside the cloned repo. macOS is verified (Mach-O arm64,
# native, no Rosetta). Linux ships too (bin/linux/classify) so it's a one-line entry, but it is
# NOT yet validated on a cluster (HTCF) — that's a separate follow-on (flagged to Megan).
_OS_BINARY = {"Darwin": "bin/macosx/classify", "Linux": "bin/linux/classify"}
_OS_VALIDATED = {"Darwin"}

# Pinned sha256 of the region zips (fill in each on first successful download; see the brief).
PINNED_SHA256 = {
    "V3V4": "0f6c0e198f0a6fdb5358a6c6c001367a82de10c9d3364bd427eaa303c355225a",
    # V1V3 + V4V4 pinned 2026-09-03: fetched from figshare v6, sha256 cross-checked against the
    # article's computed_md5 for that version (byte-identical download).
    "V1V3": "db414c7327e72ca11974a8ca49169896a8a03a49c6dd82a5ec5c38ea1e7172b0",
    "V4V4": "e38171df1197560fdee0f95c026397c5e955ddcde01fe72bbec27221358c1101",
    # V1V9: unpinned — its figshare file is the dated vSpeciateIT_V1V9_14May2026.zip (321 MB),
    # not yet fetched. First download records it trust-on-first-use (see PINNED_SHA256 note above).
}

# Zip sizes (MB) for a pre-download size/disk message. Each unpacks to ~2.4-2.6 GB; classifying
# needs ~1 GB RAM.
ZIP_SIZE_MB = {"V1V3": 142, "V1V9": 321, "V3V4": 123, "V4V4": 89}

# Ground-truth: the repo ships test.fasta (10 vaginal ASVs). These are the species the speciateIT
# README documents for the V3V4 model. On the pinned v6 DB, ASV9 legitimately flips
# (Leptotrichia_shahii -> Sneathia_sanguinegens) — a documented version drift, not a regression —
# so ASV9 is an allowed mismatch. (Only meaningful for V3V4; test.fasta is a V3V4 amplicon set.)
GROUND_TRUTH_V3V4 = {
    "ASV1": "Lactobacillus_iners",
    "ASV2": "Lactobacillus_crispatus",
    "ASV3": "Lactobacillus_mulieris",
    "ASV4": "Ca_Lachnocurva_vaginae",
    "ASV5": "Gardnerella_vaginalis",
    "ASV6": "Lactobacillus_crispatus",
    "ASV7": "Lactobacillus_iners",
    "ASV8": "Fannyhessea_vaginae",
    "ASV9": "Leptotrichia_shahii",
    "ASV10": "Megasphaera_lornae",
}
GROUND_TRUTH_ALLOWED_DRIFT = {"ASV9"}


class SetupError(Exception):
    """A setup step failed in a way the user must act on (message is shown, no traceback)."""


def _echo(msg: str = "") -> None:
    print(msg, flush=True)


# --- individual steps (module-level so tests can monkeypatch the network ones) --------------
def select_binary(clone_dir: Path, system: str | None = None) -> Path:
    """Path to the precompiled ``classify`` for this OS inside the cloned repo."""
    system = system or platform.system()
    rel = _OS_BINARY.get(system)
    if rel is None:
        raise SetupError(
            f"no precompiled speciateIT binary bundled for {system!r}. "
            f"Supported: {', '.join(_OS_BINARY)}. Build from source (see the speciateIT repo) "
            "and point composition.speciateit.classify at it."
        )
    binary = clone_dir / rel
    if not binary.exists():
        raise SetupError(f"expected the speciateIT binary at {binary}, but it isn't there.")
    return binary


def clone_speciateit(dest: Path) -> Path:
    """``git clone --depth 1`` the speciateIT repo into ``dest`` (idempotent)."""
    clone_dir = dest / "speciateIT"
    if clone_dir.exists():
        _echo(f"  speciateIT already cloned at {clone_dir} — skipping clone.")
        return clone_dir
    if shutil.which("git") is None:
        raise SetupError("git is not on PATH; install it (setup clones the speciateIT repo).")
    _echo(f"  cloning {SPECIATEIT_REPO} -> {clone_dir}")
    proc = subprocess.run(
        ["git", "clone", "--depth", "1", SPECIATEIT_REPO, str(clone_dir)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise SetupError(f"git clone failed:\n{proc.stderr.strip()}")
    return clone_dir


def resolve_figshare_file(region: str) -> dict:
    """Resolve the vSpeciateDB zip for ``region`` via the pinned figshare version.

    Returns the figshare file record (``name``, ``size``, ``download_url``, ``computed_md5``).
    Matches by version + filename prefix, NOT a hardcoded file id: file ids change across
    versions. The prefix match (rather than an exact ``vSpeciateIT_<R>.zip``) is deliberate —
    the V1V9 file is published as ``vSpeciateIT_V1V9_14May2026.zip``, so an exact-name match
    would silently miss it.
    """
    api = FIGSHARE_API.format(article=FIGSHARE_ARTICLE, version=FIGSHARE_VERSION)
    with urllib.request.urlopen(api, timeout=60) as resp:  # noqa: S310 (pinned https figshare URL)
        meta = json.load(resp)
    prefix = f"vSpeciateIT_{region}"
    # Match on the region token bounded by "." or "_" so V1V3 never swallows the V1V9 file and
    # the dated V1V9 name (vSpeciateIT_V1V9_14May2026.zip) still matches.
    matches = [
        f for f in meta.get("files", [])
        if f["name"].endswith(".zip") and (
            f["name"] == f"{prefix}.zip" or f["name"].startswith(f"{prefix}_")
        )
    ]
    if not matches:
        available = ", ".join(f["name"] for f in meta.get("files", [])) or "(none)"
        raise SetupError(
            f"no vSpeciateDB zip for region {region!r} in figshare article "
            f"{FIGSHARE_ARTICLE} v{FIGSHARE_VERSION}. Files present: {available}."
        )
    if len(matches) > 1:
        names = ", ".join(f["name"] for f in matches)
        raise SetupError(f"ambiguous match for region {region!r}: {names}.")
    return matches[0]


def download(url: str, out_path: Path, *, retries: int = 3) -> Path:
    """Download ``url`` -> ``out_path`` (streamed; follows figshare's ndownloader redirect).

    Retries a few times on a network error or a short read — a transient blip that silently
    truncates the file happens in practice (hit one during Run 1). When the server sends a
    ``Content-Length``, the written size is checked against it so a truncated download is caught
    and retried rather than passed on to the (pinned-only) sha256 check.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:  # noqa: S310 (pinned https URL)
                expected = resp.headers.get("Content-Length")
                expected = int(expected) if expected and expected.isdigit() else None
                with open(out_path, "wb") as fh:
                    shutil.copyfileobj(resp, fh, length=1024 * 1024)
            got = out_path.stat().st_size
            if expected is not None and got != expected:
                raise OSError(f"incomplete download: got {got} of {expected} bytes")
            return out_path
        except (urllib.error.URLError, OSError) as e:
            last_error = e
            if attempt < retries:
                _echo(f"  download attempt {attempt}/{retries} failed ({e}); retrying…")
    raise SetupError(f"download failed after {retries} attempts: {last_error}")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_zip(region: str, zip_path: Path, *, require_pinned: bool) -> str:
    """sha256-check the zip against the pinned manifest; TOFU-record if the region is unpinned.

    Returns the computed sha256. Raises on a *pinned* mismatch (corruption / a moved file) or,
    when ``require_pinned``, on any unpinned region.
    """
    computed = sha256_file(zip_path)
    pinned = PINNED_SHA256.get(region)
    if pinned:
        if computed != pinned:
            raise SetupError(
                f"sha256 mismatch for {zip_path.name}!\n"
                f"  expected (pinned): {pinned}\n  got:               {computed}\n"
                "Refusing to unzip — the download is corrupt or the pinned file moved."
            )
        _echo(f"  sha256 OK (matches pinned manifest): {computed}")
        return computed
    if require_pinned:
        raise SetupError(
            f"region {region!r} has no pinned sha256 and --require-pinned was set.\n"
            f"  computed sha256: {computed}\n"
            "Add it to PINNED_SHA256 (and the brief's manifest), then re-run."
        )
    _echo(
        "  ⚠️  UNPINNED REGION — trusting on first use.\n"
        f"      computed sha256: {computed}\n"
        f"      Pin it: add  {region!r}: {computed!r}  to PINNED_SHA256 in "
        "microfgt/setup16s.py and to the manifest in design/setup_task_brief.md."
    )
    return computed


def unzip(zip_path: Path, models_root: Path, region: str) -> Path:
    """Unzip to ``<models_root>/vSpeciateIT_<region>/`` (idempotent)."""
    target = models_root / f"vSpeciateIT_{region}"
    if target.exists():
        _echo(f"  model dir already present at {target} — skipping unzip.")
        return target
    models_root.mkdir(parents=True, exist_ok=True)
    _echo(f"  unzipping {zip_path.name} -> {models_root}/ …")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(models_root)
    if not target.exists():
        # Some zips nest the region dir one level down (or name it slightly differently); find it.
        cand = _find_model_dir(models_root, region)
        if cand is None:
            raise SetupError(
                f"unzipped {zip_path.name} but no vSpeciateIT_{region} model dir appeared under "
                f"{models_root}. Inspect the archive layout."
            )
        target = cand
    return target


def _find_model_dir(models_root: Path, region: str) -> Path | None:
    for p in models_root.rglob(f"vSpeciateIT_{region}"):
        if p.is_dir():
            return p
    return None


def run_classify(classify: Path, db: Path, fasta: Path, outdir: Path) -> Path:
    """Run ``classify -d <db> -i <fasta> -o <outdir>``; return the results file path."""
    from microfgt.orchestrate import run_speciateit

    run_speciateit(fasta=str(fasta), db=str(db), outdir=str(outdir), executable=str(classify))
    return Path(outdir) / "MC_order7_results.txt"


def parse_results(results_path: Path) -> dict[str, str]:
    """Parse a (headerless) MC_order7_results.txt -> {ASV id: species}."""
    calls: dict[str, str] = {}
    for line in Path(results_path).read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] and parts[0].lower() not in ("seq", "sequence id"):
            calls[parts[0]] = parts[1]
    return calls


def ground_truth_report(calls: dict[str, str]) -> tuple[bool, list[str]]:
    """Compare V3V4 test.fasta calls to the README ground truth (ASV9 drift allowed)."""
    lines, matched, unexpected = [], 0, []
    for asv, expected in GROUND_TRUTH_V3V4.items():
        got = calls.get(asv)
        if got == expected:
            matched += 1
        elif asv in GROUND_TRUTH_ALLOWED_DRIFT:
            lines.append(f"    {asv}: {got}  (expected {expected} on older DB — documented v6 drift, OK)")
        else:
            unexpected.append(f"    {asv}: {got}  (expected {expected})")
    ok = not unexpected and matched >= len(GROUND_TRUTH_V3V4) - len(GROUND_TRUTH_ALLOWED_DRIFT)
    header = f"  ground-truth: {matched}/{len(GROUND_TRUTH_V3V4)} species match the README"
    if GROUND_TRUTH_ALLOWED_DRIFT:
        header += f" (+{len(GROUND_TRUTH_ALLOWED_DRIFT)} allowed v6 drift)"
    report = [header] + lines + unexpected
    return ok, report


def write_config(
    dest: Path, region: str, classify: Path, db: Path, test_fasta: Path
) -> Path:
    """Write a check-clean, immediately-runnable config pointing at the fresh install.

    It wires the bundled ``test.fasta`` (+ a generated one-sample count table) as the entry
    point, so ``microfgt check`` resolves the assign stage (binary, db, region<->DB match all
    exercised) and ``microfgt run`` on it is a real end-to-end smoke test. Swap the two
    ``asv_seqs`` / ``asv_table`` lines (or add a ``reads: {fastq_dir, primers}`` block) to point
    it at your own data.
    """
    import yaml

    selftest = dest / "selftest"
    selftest.mkdir(parents=True, exist_ok=True)
    count_table = _write_selftest_count_table(test_fasta, selftest / "asv_table.csv")

    config = {
        "composition": {
            "reads": {"region": region},
            # --- bundled self-test inputs: replace with your own ASV FASTA + count table,
            #     or delete both and add reads.fastq_dir to run from raw FASTQs ---
            "asv_seqs": str(test_fasta),
            "asv_table": str(count_table),
            "speciateit": {
                "classify": str(classify),
                "db": str(db),
            },
        },
        "cst": {"method": "centroid"},
        "output": str(selftest / "microfgt-selftest.h5mu"),
    }
    config_path = dest / "microfgt-16s.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config_path


def _write_selftest_count_table(test_fasta: Path, out_csv: Path) -> Path:
    import csv

    ids = [
        line[1:].split()[0]
        for line in Path(test_fasta).read_text().splitlines()
        if line.startswith(">")
    ]
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample"] + ids)
        w.writerow(["selftest"] + [10] * len(ids))
    return out_csv


def run_check(config_path: Path) -> bool:
    """Run ``microfgt check`` on the written config; print each line, return whether all OK."""
    from microfgt.config import load_config
    from microfgt.stages import check

    results = check(load_config(str(config_path)))
    for r in results:
        _echo(f"  {r.message}")
    return all(r.ok for r in results)


# --- the orchestrator ----------------------------------------------------------------------
def run_setup(
    region: str,
    dest: str,
    *,
    require_pinned: bool = False,
    skip_download: bool = False,
) -> int:
    """Install speciateIT + the vSpeciateDB model for ``region`` under ``dest``; write a config
    and run ``microfgt check``. Returns a process exit code (0 == all prereqs OK)."""
    try:
        region = canonical_region(region)
    except ValueError as e:
        _echo(f"error: {e}")
        return 2

    dest_path = Path(dest).expanduser().resolve()
    _echo(f"microfgt setup — region {region}, dest {dest_path}")
    if " " in str(dest_path):
        _echo(
            "  ⚠️  --dest contains a space. speciateIT's `classify` shells out to `mkdir` "
            "unquoted and WILL break on a spaced path at run time. Prefer a space-free --dest."
        )
    dest_path.mkdir(parents=True, exist_ok=True)

    try:
        # 1) speciateIT binary
        _echo("\n[1/5] speciateIT binary")
        clone_dir = clone_speciateit(dest_path)
        classify = select_binary(clone_dir)
        system = platform.system()
        note = "" if system in _OS_VALIDATED else "  (NOTE: not yet validated on this OS)"
        _echo(f"  classify: {classify}{note}")

        # 2/3/4) the vSpeciateDB model
        _echo(f"\n[2/5] vSpeciateDB model for {region}")
        models_root = dest_path / "vSpeciateDB_models"
        db_dir = _find_model_dir(models_root, region)
        if db_dir is not None:
            _echo(f"  model already installed at {db_dir} — skipping download.")
        elif skip_download:
            raise SetupError(
                f"no model dir for {region} under {models_root} and skip_download was set."
            )
        else:
            record = resolve_figshare_file(region)
            size_mb = ZIP_SIZE_MB.get(region, round(record.get("size", 0) / 1e6))
            _echo(
                f"  downloading {record['name']} (~{size_mb} MB; unpacks to ~2.4-2.6 GB, "
                "needs ~1 GB RAM to classify)"
            )
            zip_path = download(record["download_url"], models_root / record["name"])
            verify_zip(region, zip_path, require_pinned=require_pinned)
            db_dir = unzip(zip_path, models_root, region)
        _echo(f"  db: {db_dir}")

        # 5) write config, ground-truth check, microfgt check
        _echo("\n[3/5] writing config")
        test_fasta = clone_dir / "test.fasta"
        config_path = write_config(dest_path, region, classify, db_dir, test_fasta)
        _echo(f"  wrote {config_path}")

        _echo("\n[4/5] ground-truth classification test")
        if region == "V3V4" and test_fasta.exists():
            out = dest_path / "selftest" / "ground_truth"
            calls = parse_results(run_classify(classify, db_dir, test_fasta, out))
            ok, report = ground_truth_report(calls)
            for line in report:
                _echo(line)
            if not ok:
                _echo("  ⚠️  ground-truth mismatch beyond the documented ASV9 drift — investigate.")
        else:
            _echo(f"  skipped (ground-truth fixture is V3V4-specific; region is {region}).")

        _echo("\n[5/5] microfgt check")
        all_ok = run_check(config_path)
    except SetupError as e:
        _echo(f"\nerror: {e}")
        return 1

    if all_ok:
        _echo("\n✅ all prerequisites satisfied. Point the config at your data and `microfgt run`.")
        return 0
    _echo("\n❌ some prerequisites are still missing — see the check output above.")
    return 1


def main(argv: list[str] | None = None) -> int:
    """Thin argv entry (the real CLI wiring is in cli.py; this eases direct invocation/tests)."""
    import argparse

    p = argparse.ArgumentParser(prog="microfgt setup")
    p.add_argument("--region", required=True)
    p.add_argument("--dest", required=True)
    p.add_argument("--require-pinned", action="store_true")
    a = p.parse_args(argv if argv is not None else sys.argv[1:])
    return run_setup(a.region, a.dest, require_pinned=a.require_pinned)
