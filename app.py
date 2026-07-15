"""Streamlit UI for dd_seq.

Run with `streamlit run app.py -- --report-dir data` (after `dd_seq-run`/
`dd_seq-fetch`+`dd_seq-align` have populated that directory with
`report.json` and `aligned/*_aligned.pdb`), or just `streamlit run app.py`
and enter the directory in the sidebar.
"""
import argparse
import json
import sys
from pathlib import Path

import streamlit as st
from dd_viewer import html_with_camera_events, view3d

from dd_seq import dashboard, scene, seqplot

st.set_page_config(page_title="dd_seq", layout="wide")


def _parse_cli_defaults() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dir")
    args, _ = parser.parse_known_args(sys.argv[1:])
    return args


@st.cache_data(show_spinner=False)
def _load_report(report_dir: str) -> dict:
    return json.loads((Path(report_dir) / "report.json").read_text())


def main() -> None:
    defaults = _parse_cli_defaults()
    st.title("dd_seq -- structure & sequence comparison")

    with st.sidebar:
        st.header("Report")
        report_dir = st.text_input("Report directory (dd_seq-run/-align output)", value=defaults.report_dir or "data")

    report_path = Path(report_dir) / "report.json"
    if not report_path.exists():
        st.info(f"No report.json found in {report_dir!r}. Run `dd_seq-run UNIPROT -o {report_dir}` first.")
        st.stop()

    report = _load_report(report_dir)
    structures = report["structures"]
    labels = [s["label"] for s in structures]

    st.caption(
        f"UniProt {report['uniprot_id']} -- {len(labels)} structures -- "
        f"site_mode={report['site_mode']} (source: {report['site_source']}) -- reference: {report['reference']}"
    )

    with st.sidebar:
        st.header("Overlay display")
        selected = st.multiselect("Structures to show", labels, default=labels)
        show_ligands = st.checkbox("Show ligands", value=True)
        has_site = report["site_mode"] != "none"
        show_site = st.checkbox(
            "Highlight active-site residues", value=has_site, disabled=not has_site,
            help=None if has_site else "This report was built with --site-mode none, so no active site was defined.",
        )

        if "camera_generation" not in st.session_state:
            st.session_state.camera_generation = 0
        if st.button("Reset view"):
            st.session_state.camera_generation += 1

    tab_overview, tab_coverage, tab_overlay = st.tabs(["Overview", "Sequence coverage", "Structure overlay"])

    with tab_overview:
        st.dataframe(dashboard.summary_dataframe(report), width="stretch", hide_index=True)

    with tab_coverage:
        fig = seqplot.plot_coverage(structures, site_canonical_positions=report.get("site_canonical_positions"))
        st.pyplot(fig, width="stretch")

    with tab_overlay:
        by_label = {s["label"]: s for s in structures}
        scene_structures = []
        for label in selected:
            s = by_label[label]
            if not s.get("aligned_pdb"):
                continue  # skipped during alignment (see s.get("align_error"))
            scene_structures.append(
                {
                    "label": label,
                    "pdb_path": s["aligned_pdb"],
                    "chain_id": s["chain"],
                    "site_resseqs": s["site_resseqs"] if show_site else None,
                    "show_ligand": show_ligands,
                }
            )
        if not scene_structures:
            st.info("No selected structure has a superposed coordinate file to show (all skipped during alignment?).")
        else:
            view = scene.build_overlay_view(scene_structures)
            html = html_with_camera_events(view._make_html())
            view3d(html, height=650, reset_camera_token=st.session_state.camera_generation)

        skipped = [s for s in structures if s["label"] in selected and not s.get("aligned_pdb")]
        if skipped:
            st.warning(
                "Not shown (couldn't be superposed at the requested site): "
                + ", ".join(f"{s['label']} ({s['align_error']})" for s in skipped)
            )


if __name__ == "__main__":
    main()
