"""End-to-end workflow + turnkey CLI (the P4 headline: real outputs in, .h5mu out)."""

import mudata as md
import pytest

from microfgt.cli import main
from microfgt.workflow import run_workflow

pytestmark = pytest.mark.filterwarnings("ignore::FutureWarning")


def _config(real_fixtures, test_data, out_path):
    return {
        "composition": {
            "speciateit": {
                "results": str(test_data / "speciateit_MC_order7_results.synthetic.txt"),
                "count_table": str(real_fixtures / "speciateit_test_count_table.csv"),
            }
        },
        "cst": {"method": "centroid"},
        "analysis": {
            "transforms": ["relabund", "clr"],
            "alpha": ["shannon"],
            "beta": "braycurtis",
            "ordinate": True,
        },
        "output": str(out_path),
    }


def test_run_workflow_imports_classifies_and_analyzes(real_fixtures, test_data, tmp_path):
    mdata = run_workflow(_config(real_fixtures, test_data, tmp_path / "x.h5mu"))

    comp = mdata["composition"]
    assert comp.n_obs == 169
    # CST classified from composition and attached as sample annotation.
    assert "CST" in mdata.obs.columns
    assert mdata.obs["CST"].notna().any()
    # Analysis results written onto the modality.
    assert "relabund" in comp.layers and "clr" in comp.layers
    assert "alpha_shannon" in comp.obs
    assert comp.obsm["X_pcoa"].shape[0] == 169


def test_cli_run_writes_h5mu_and_reloads(real_fixtures, test_data, tmp_path, capsys):
    out = tmp_path / "out.h5mu"
    cfg = tmp_path / "cfg.yaml"
    import yaml

    cfg.write_text(yaml.safe_dump(_config(real_fixtures, test_data, out)))

    rc = main(["run", "-c", str(cfg)])
    assert rc == 0
    assert out.exists()
    assert "wrote" in capsys.readouterr().out

    # Round-trips: the written object reloads with the classified CST + analysis.
    reloaded = md.read(out)
    assert "composition" in reloaded.mod
    assert "CST" in reloaded.obs.columns
    assert "alpha_shannon" in reloaded["composition"].obs


def test_cli_classify_then_analyze_on_existing_h5mu(real_fixtures, test_data, tmp_path):
    # First produce a composition-only object (no CST/analysis), then drive the subcommands.
    base = run_workflow({
        "composition": {"speciateit": {
            "results": str(test_data / "speciateit_MC_order7_results.synthetic.txt"),
            "count_table": str(real_fixtures / "speciateit_test_count_table.csv"),
        }},
    })
    base_path = tmp_path / "base.h5mu"
    base.write(base_path)

    clf_out = tmp_path / "clf.h5mu"
    assert main(["classify", "-i", str(base_path), "-o", str(clf_out), "-m", "centroid"]) == 0
    assert "CST" in md.read(clf_out).obs.columns

    ana_out = tmp_path / "ana.h5mu"
    assert main(["analyze", "-i", str(clf_out), "-o", str(ana_out),
                 "--transform", "clr", "--alpha", "shannon", "--ordinate"]) == 0
    reloaded = md.read(ana_out)
    assert "clr" in reloaded["composition"].layers
    assert "alpha_shannon" in reloaded["composition"].obs
