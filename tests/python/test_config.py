"""Config loading: an empty YAML section must not crash the `.get(...)` sites.

`mgcst:` with no body parses to None, not {}. The key is present, so a bare
`config.get("mgcst", {})` returns that None and the next `.get(...)` raises. The loader
normalizes empty sections to absent ones so this whole class can't recur.
"""

from microfgt.config import load_config, normalize_config
from microfgt.stages.registry import provided_artifacts


def test_normalize_drops_none_valued_sections():
    assert normalize_config({"mgcst": None, "cst": {"method": "centroid"}}) == {
        "cst": {"method": "centroid"}
    }


def test_normalize_recurses_into_nested_sections():
    # composition present, but its `reads` sub-section is empty.
    assert normalize_config({"composition": {"reads": None, "asv_table": "x.csv"}}) == {
        "composition": {"asv_table": "x.csv"}
    }


def test_empty_section_does_not_crash_provided_artifacts(tmp_path):
    # The exact bite: an empty `mgcst:` section alongside a real entry point.
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "mgcst:\n"
        "composition:\n"
        "  reads:\n"
        "metagenomics:\n"
        "  compiled: VIRGO2_Compiled.summary.NR.txt\n"
    )
    config = load_config(cfg)
    provided = provided_artifacts(config)          # must not raise
    assert provided["sg_compiled"] == "VIRGO2_Compiled.summary.NR.txt"


def test_wholly_empty_config_is_a_dict(tmp_path):
    cfg = tmp_path / "empty.yaml"
    cfg.write_text("")
    assert load_config(cfg) == {}
