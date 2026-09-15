"""microFGT dashboard — a thin Streamlit presenter over the analysis verbs.

Launch it with ``microfgt dashboard -i your_object.h5mu`` (needs the ``[app]`` extra). Every
control here just builds a *spec* of choices and hands it to ``run_verb``; the result draws
itself via ``viz.render``. There is no analysis logic in this file — that is the whole point
(Layer-2 principle: the surface presents the same calls a power user makes).
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import mudata as md
import streamlit as st

from microfgt.analysis import run_verb
from microfgt.dashboard.catalog import as_frame, continuous, groupable, variable_catalog
from microfgt.viz import render

_VERB_LABELS = {
    "alpha": "Alpha diversity ~ group",
    "beta": "Beta diversity (PERMANOVA)",
    "associate": "Associate two variables",
    "abundance": "Differential abundance",
}


@st.cache_resource(show_spinner=False)
def _load(path: str):
    # cache_resource (not cache_data): a MuData is held by reference, not pickled.
    return md.read(path)


def _sidebar_spec(mdata):
    """Build the (verb, kwargs) spec from the object's variables — the user's choices."""
    catalog = variable_catalog(mdata)
    groups = groupable(catalog)
    conts = continuous(catalog)
    allvars = [v.name for v in catalog]
    catmap = {v.name: v for v in catalog}

    verb = st.sidebar.selectbox("Analysis", list(_VERB_LABELS),
                                format_func=lambda v: _VERB_LABELS[v])
    kwargs: dict = {}

    if verb == "associate":
        kwargs["x"] = st.sidebar.selectbox("Variable X", allvars)
        kwargs["y"] = st.sidebar.selectbox("Variable Y", [v for v in allvars if v != kwargs["x"]])
    else:
        options = groups or allvars
        primary = st.sidebar.selectbox("Predictor of interest", options)
        preds = [primary]
        if verb in ("alpha", "abundance"):
            preds += st.sidebar.multiselect("Adjust for (covariates)",
                                            [v for v in allvars if v != primary])
            subj = st.sidebar.selectbox("Subject (repeated measures)", ["(none)"] + groups)
            kwargs["subject"] = None if subj == "(none)" else subj
        kwargs["predictors"] = preds
        if verb == "alpha":
            kwargs["metric"] = st.sidebar.selectbox("Alpha metric", ["shannon", "simpson"])
        if verb == "abundance":
            kwargs["method"] = st.sidebar.selectbox("Method", ["ancombc", "dirmult_lme"])

    with st.sidebar.expander("Subset (optional)"):
        scol = st.selectbox("Restrict by", ["(all samples)"] + groups)
        if scol != "(all samples)":
            keep = st.multiselect("Keep values", catmap[scol].levels, default=catmap[scol].levels)
            if keep:
                kwargs["subset"] = {scol: keep}

    return verb, kwargs, catalog


def main():
    st.set_page_config(page_title="microFGT", layout="wide")
    st.title("microFGT — explore your integrated object")

    path = st.sidebar.text_input("Object (.h5mu)", value=os.environ.get("MICROFGT_H5MU", ""))
    if not path:
        st.info("Enter the path to a microFGT `.h5mu` object in the sidebar to begin.")
        return
    try:
        mdata = _load(path)
    except Exception as e:  # noqa: BLE001 - surface any read error to the user
        st.error(f"Could not read `{path}`: {e}")
        return

    st.sidebar.success(f"{mdata.n_obs} samples · modalities: {', '.join(mdata.mod)}")

    # Parameterized descriptor, materialized on demand (not stored on the object): the cutoff is
    # an exploration knob, not a baked-in constant. Recomputed each run at the slider's value.
    if "composition_taxon" in mdata.mod:
        from microfgt.characterize import taxa_over_threshold

        pct = st.sidebar.slider("Richness cutoff (% abundance)", 1, 50, 10,
                                help="# taxa above this relative abundance — an adjustable view. "
                                     "The cutoff-free evenness lives in `effective_taxa`.")
        counts = taxa_over_threshold(mdata["composition_taxon"], pct / 100)
        mdata.obs["taxa_over_threshold"] = counts.reindex(mdata.obs_names.astype(str)).to_numpy()
        st.sidebar.caption(f"`taxa_over_threshold` = # taxa above {pct}% (adjust above)")

    verb, kwargs, catalog = _sidebar_spec(mdata)
    run = st.sidebar.button("Run", type="primary")

    home_tab, analysis_tab = st.tabs(["🏠 Home", "📊 Analysis"])
    with home_tab:
        _render_home(mdata, catalog, path)
    with analysis_tab:
        if not run:
            st.info("Configure an analysis in the sidebar, then click **Run**.")
        else:
            try:
                with st.spinner("Running…"):
                    result = run_verb(mdata, verb, **kwargs)
            except Exception as e:  # noqa: BLE001 - a bad selection should message, not crash
                st.error(str(e))
            else:
                _render_result(result)


def _render_home(mdata, catalog, path):
    """Overview tab: cohort at a glance + the standard FGT plots."""
    # Resolve annotations through the same prefix-stripped, merged view the analysis layer
    # uses. Metadata that rides on a modality's obs is stored by MuData under a "<modality>:"
    # prefix, so reading raw mdata.obs["PID"]/["read_count"] misses it and falls back to "—".
    from microfgt.analysis._frame import merged_obs

    obs = merged_obs(mdata)
    n_taxa = mdata["composition_taxon"].n_vars if "composition_taxon" in mdata.mod else mdata.n_vars
    cols = st.columns(5)
    cols[0].metric("Samples", f"{mdata.n_obs:,}")
    cols[1].metric("Taxa", f"{n_taxa:,}")
    subj_col = next((c for c in ("PID", "subject", "subject_id", "Subject_number") if c in obs.columns), None)
    subj = obs[subj_col].nunique() if subj_col else None
    cols[2].metric("Subjects", f"{subj:,}" if subj else "—")
    cols[3].metric("CSTs", obs["CST"].nunique() if "CST" in obs.columns else "—")
    depth = obs["read_count"].median() if "read_count" in obs.columns else None
    cols[4].metric("Median reads", f"{int(depth):,}" if depth else "—")

    figs = _overview_figs(path)
    if "community" in figs:
        st.markdown("**Community composition** — each column is one sample")
        st.pyplot(figs["community"])            # the iconic full-width composition stack
    row = st.columns(2)
    for name, col in zip(("cst", "diversity", "dominance", "ordination"), (row[0], row[1], row[0], row[1])):
        if name in figs:
            col.pyplot(figs[name])

    with st.expander("Variables available to explore"):
        st.dataframe(as_frame(catalog), use_container_width=True, height=240)


@st.cache_resource(show_spinner="Building overview…")
def _overview_figs(path: str):
    # cache_resource keyed on the object path so the (heavier) overview computes once, not every rerun.
    from microfgt.viz.overview import overview_figures

    return overview_figures(_load(path))


def _render_result(result):
    st.subheader(result.summary())
    fig_col, tbl_col = st.columns(2)
    with fig_col:
        st.pyplot(render(result).figure)
    with tbl_col:
        st.dataframe(result.table, use_container_width=True)
    with st.expander("stats · spec · notes"):
        st.json({"stats": result.stats, "spec": result.spec, "notes": result.notes})


if __name__ == "__main__":
    main()
