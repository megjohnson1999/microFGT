"""Turnkey CLI — runs the whole workflow, with the Python API underneath for power users.

Config-first (constraint A). The entry point is just which inputs the config provides:
FASTQs run the full ladder (primer-trim -> DADA2 -> speciateIT -> CST -> analysis); an ASV
table enters at speciateIT; existing outputs enter at import. ``microfgt check`` verifies the
tools/paths the resolved entry point needs, up front.
"""

from __future__ import annotations

import argparse
import sys
import tempfile

from microfgt import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="microfgt", description="One-stop FGT microbiome analysis."
    )
    parser.add_argument("--version", action="version", version=f"microfgt {__version__}")
    parser.set_defaults(_run=lambda _a: parser.print_help())
    sub = parser.add_subparsers(title="commands")

    run = sub.add_parser("run", help="Run the workflow from a config file (any entry point).")
    run.add_argument("-c", "--config", required=True, help="YAML workflow config.")
    run.add_argument("-o", "--output", help="Output .h5mu (overrides config 'output').")
    run.add_argument("--workdir", help="Directory for intermediate artifacts.")
    run.add_argument("--executor", choices=["local", "snakemake"], default="local")
    run.set_defaults(_run=_cmd_run)

    chk = sub.add_parser("check", help="Preflight: verify tools/paths for the resolved entry point.")
    chk.add_argument("-c", "--config", required=True)
    chk.set_defaults(_run=_cmd_check)

    setup = sub.add_parser(
        "setup", help="Install the 16S prerequisites conda can't (speciateIT + a vSpeciateDB model)."
    )
    setup.add_argument("--region", required=True,
                       help="16S region / model: V1V3 | V1V9 | V3V4 | V4V4 (V4 aliases V4V4).")
    setup.add_argument("--dest", required=True, help="Directory to install the binary + model into.")
    setup.add_argument("--require-pinned", action="store_true",
                       help="Refuse a region whose zip sha256 isn't pinned (default: trust-on-"
                            "first-use and record the computed hash).")
    setup.set_defaults(_run=_cmd_setup)

    clf = sub.add_parser("classify", help="Classify CST on a .h5mu's composition modality.")
    clf.add_argument("-i", "--input", required=True)
    clf.add_argument("-o", "--output", required=True)
    clf.add_argument("-m", "--method", default="centroid")
    clf.set_defaults(_run=_cmd_classify)

    cmp = sub.add_parser("compare", help="Run a hypothesis-test verb on a .h5mu; print/save the result.")
    cmp.add_argument("-i", "--input", required=True)
    cmp.add_argument("--verb", required=True, choices=["alpha", "beta", "abundance", "associate"])
    cmp.add_argument("--predictors", help="Comma-separated obs columns (alpha/beta/abundance); "
                     "first is the predictor of interest.")
    cmp.add_argument("--x"); cmp.add_argument("--y")                 # associate
    cmp.add_argument("--subject", help="obs subject-id column (repeated measures).")
    cmp.add_argument("--metric", default=None, help="diversity metric (default: shannon for "
                     "alpha, braycurtis for beta).")
    cmp.add_argument("--method", help="verb method override (e.g. ancombc / dirmult_lme / fisher).")
    cmp.add_argument("--subset", help="Restrict samples: a pandas query string.")
    cmp.add_argument("-o", "--output", help="Write the result table to this CSV.")
    cmp.add_argument("--plot", help="Render the result to this image (PNG/SVG).")
    cmp.set_defaults(_run=_cmd_compare)

    ana = sub.add_parser("analyze", help="Run analysis steps on a .h5mu's composition modality.")
    ana.add_argument("-i", "--input", required=True)
    ana.add_argument("-o", "--output", required=True)
    ana.add_argument("--transform", nargs="*", default=[], choices=["relabund", "clr"])
    ana.add_argument("--alpha", nargs="*", default=[], metavar="METRIC")
    ana.add_argument("--beta", default=None, metavar="METRIC")
    ana.add_argument("--ordinate", action="store_true")
    ana.add_argument("--diffabund-group", default=None, metavar="OBS_COL")
    ana.set_defaults(_run=_cmd_analyze)

    dash = sub.add_parser("dashboard", help="Launch the Streamlit dashboard (needs the [app] extra).")
    dash.add_argument("-i", "--input", help="Preload this .h5mu object.")
    dash.set_defaults(_run=_cmd_dashboard)

    rs = sub.add_parser("_run-stage", help=argparse.SUPPRESS)  # internal (Snakemake calls it)
    rs.add_argument("stage")
    rs.add_argument("--workdir", required=True)
    rs.add_argument("--config", required=True)
    rs.add_argument("--output")
    rs.set_defaults(_run=_cmd_run_stage)
    return parser


def _load_config(path: str) -> dict:
    from microfgt.config import load_config

    return load_config(path)


def _attach_cst(mdata, cst_df) -> None:
    from microfgt.io.integrate import attach_cst_annotations

    # Labels -> .obs; <subCST>_sim vectors -> composition_taxon.obsm['cst_sim'].
    attach_cst_annotations(mdata, cst_df)


def _cmd_run(args: argparse.Namespace) -> None:
    from microfgt.stages import (
        LocalExecutor,
        SnakemakeExecutor,
        provided_artifacts,
        resolve,
        speciateit_space_warnings,
    )

    config = _load_config(args.config)
    out = args.output or config.get("output")
    if not out:
        raise SystemExit("No output path: pass -o/--output or set 'output:' in the config.")
    workdir = args.workdir or tempfile.mkdtemp(prefix="microfgt_")
    for w in speciateit_space_warnings(config, workdir=workdir):
        print(w)
    stages = resolve("mudata", provided_artifacts(config))
    plan = " -> ".join(s.id for s in stages) or "(nothing to do)"

    if args.executor == "snakemake":
        path = SnakemakeExecutor().run(stages, args.config, workdir, out)
        print(f"entry-point plan: {plan}")
        print(f"wrote Snakefile: {path}\nsubmit on the cluster with snakemake + a Slurm profile.")
        return

    print(f"entry-point plan: {plan}")
    written = LocalExecutor().run(stages, workdir, config, out)
    print(f"wrote {written}")


def _cmd_check(args: argparse.Namespace) -> None:
    from microfgt.stages import check, speciateit_space_warnings

    config = _load_config(args.config)
    results = check(config)
    for r in results:
        print(r.message)
    for w in speciateit_space_warnings(config):
        print(w)
    missing = [r for r in results if not r.ok]
    if missing:
        hint = _conda_env_hint(missing)
        if hint:
            print(hint)
        raise SystemExit(f"{len(missing)} prerequisite(s) missing — see above.")
    print("all prerequisites satisfied.")


# The 16S tools conda provides — if these are the ones missing, the env is likely just not active.
_CONDA_TOOLS = {"cutadapt", "rscript", "dada2"}


def _conda_env_hint(missing) -> str | None:
    """A pointed hint when the *missing* prereqs are conda-provided tools: the environment is
    probably just not activated (they report MISS on PATH even though they're installed)."""
    names = {r.message.split("'")[1].split("/")[-1].lower() for r in missing if "'" in r.message}
    if names & _CONDA_TOOLS:
        return ("hint: those tools are provided by the conda env — if you installed it, activate "
                "it first (e.g. `conda activate microfgt-16s`) and re-run, before installing anything.")
    return None


def _cmd_setup(args: argparse.Namespace) -> None:
    from microfgt.setup16s import run_setup

    rc = run_setup(args.region, args.dest, require_pinned=args.require_pinned)
    if rc != 0:
        raise SystemExit(rc)


def _cmd_run_stage(args: argparse.Namespace) -> None:
    from microfgt.stages import execute_stage

    execute_stage(args.stage, args.workdir, _load_config(args.config), args.output)


def _cmd_classify(args: argparse.Namespace) -> None:
    import mudata as md

    from microfgt.cst import classify_cst

    mdata = md.read(args.input)
    # CST is computed from the taxon roll-up; fall back to composition for older objects.
    mod = "composition_taxon" if "composition_taxon" in mdata.mod else "composition"
    if mod not in mdata.mod:
        raise SystemExit("Input has no composition modality to classify.")
    cst = classify_cst(mdata[mod], method=args.method)
    _attach_cst(mdata, cst)
    mdata.write(args.output)
    print(f"wrote {args.output}: classified CST ({args.method}) for {cst.shape[0]} samples")


def _cmd_dashboard(args: argparse.Namespace) -> None:
    import os
    import shutil
    import subprocess
    from importlib import resources

    if shutil.which("streamlit") is None:
        raise SystemExit("Streamlit is not installed. Install it with: pip install 'microfgt[app]'")
    app = str(resources.files("microfgt.dashboard").joinpath("app.py"))
    env = dict(os.environ)
    if args.input:
        env["MICROFGT_H5MU"] = str(args.input)
    subprocess.run(["streamlit", "run", app], env=env)


def _cmd_compare(args: argparse.Namespace) -> None:
    import mudata as md

    from microfgt import analysis

    mdata = md.read(args.input)
    preds = [p.strip() for p in args.predictors.split(",")] if args.predictors else None
    metric = args.metric or ("braycurtis" if args.verb == "beta" else "shannon")
    try:
        result = analysis.run_verb(
            mdata, args.verb, predictors=preds, x=args.x, y=args.y, subject=args.subject,
            metric=metric, method=args.method, subset=args.subset,
        )
    except ValueError as e:
        raise SystemExit(str(e))

    print(result.summary())
    print(result.table.to_string())
    if args.output:
        result.table.to_csv(args.output)
        print(f"wrote {args.output}")
    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        from microfgt.viz import render

        ax = render(result)
        ax.figure.savefig(args.plot, dpi=110, bbox_inches="tight")
        print(f"wrote {args.plot}")


def _cmd_analyze(args: argparse.Namespace) -> None:
    import mudata as md

    from microfgt.workflow import apply_analysis

    mdata = md.read(args.input)
    if "composition" not in mdata.mod:
        raise SystemExit("Input has no 'composition' modality to analyze.")
    cfg = {
        "transforms": args.transform,
        "alpha": args.alpha,
        "beta": args.beta,
        "ordinate": args.ordinate,
    }
    if args.diffabund_group:
        cfg["diffabund"] = {"group": args.diffabund_group, "method": "ancom"}
    apply_analysis(mdata["composition"], cfg)
    mdata.write(args.output)
    print(f"wrote {args.output}: analysis applied to composition modality")


def main(argv: list[str] | None = None) -> int:
    from microfgt.orchestrate._run import ToolNotFoundError, ToolRunError

    parser = build_parser()
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    try:
        args._run(args)
    except (ToolRunError, ToolNotFoundError) as e:
        # Show the external tool's own error (or the not-found hint), not a Python traceback.
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
