"""py3Dmol multi-structure overlay scene, built from `structalign`'s
superposed coordinate files. Every structure gets its own flat cartoon
color -- not the per-chain "spectrum" rainbow `dd_viewer.scene` uses for a
single receptor, since with a dozen structures overlaid, telling
*structures* apart matters more than telling chains within one apart.
Active-site residues (if a site was used for the fit) are drawn as sticks
in one shared highlight color across every structure, so the same site
stays visually traceable across the whole overlay; each structure's own
ligand, if it has one, can be toggled on as sticks in that structure's
color.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import py3Dmol

from dd_prep.hetero import classify_hetero_groups, pick_ligand_of_interest
from dd_prep.parse import collect_hetero_groups

PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]
AFDB_COLOR = "#444444"
SITE_COLOR = "yellow"


def assign_colors(labels: Sequence[str], afdb_label: str = "AFDB") -> Dict[str, str]:
    """One color per label, cycling through `PALETTE`; `afdb_label` (if
    present) always gets the same neutral gray so the AlphaFold model
    reads as "the reference", not just another structure in the cycle."""
    colors: Dict[str, str] = {}
    i = 0
    for label in labels:
        if label == afdb_label:
            colors[label] = AFDB_COLOR
        else:
            colors[label] = PALETTE[i % len(PALETTE)]
            i += 1
    return colors


def build_overlay_view(
    structures: Sequence[dict],
    *, colors: Optional[Dict[str, str]] = None, width: Union[int, str] = "100%", height: Union[int, str] = 600,
    focus_on_site: bool = False, focus_radius: float = 8.0,
) -> py3Dmol.view:
    """`structures`: each a dict with `label`, `pdb_path` (superposed
    coordinates, from `structalign.align_structures`'s `aligned_pdb`),
    `chain_id`, and optionally `site_resseqs` (that structure's own
    numbering, for highlighting/focusing), `highlight_site` (bool, whether
    to draw `site_resseqs` in `SITE_COLOR`; defaults to True when
    `site_resseqs` is set) and `show_ligand` (bool).

    `width` defaults to `"100%"` (py3Dmol/3Dmol.js accept a CSS size
    string, not just a pixel int) rather than a fixed pixel count -- a
    fixed-width scene silently gets cropped rather than scaled down
    whenever its embedding container (e.g. a narrower Streamlit column)
    ends up smaller than that width.

    `focus_on_site`: if True, structures with a `site_resseqs` skip the
    whole-chain cartoon and only draw residues within `focus_radius`
    angstroms of the site (as sticks), so everything outside the pocket
    is hidden rather than merely un-highlighted. Structures without a
    site fall back to the normal whole-chain cartoon, since there's
    nothing defined to focus around.
    """
    view = py3Dmol.view(width=width, height=height)
    colors = colors or assign_colors([s["label"] for s in structures])

    for model_index, s in enumerate(structures):
        pdb_text = Path(s["pdb_path"]).read_text()
        view.addModel(pdb_text, "pdb")
        color = colors.get(s["label"], "gray")
        chain_sel = {"model": model_index, "chain": s["chain_id"]}
        site = s.get("site_resseqs")
        # Only the target chain is drawn -- co-crystallized partner chains
        # (e.g. cyclin B, Cks2) were never part of the site/whole-chain fit
        # (structalign.py only fits `chain_id`), so their positions don't
        # correspond across structures and would just clutter the overlay
        # with unaligned mass around the one thing that *is* superposed.
        # Every other chain is explicitly styled to nothing first, since
        # 3Dmol.js falls back to a default line/wireframe rendering for
        # any atom left unstyled rather than hiding it.
        view.setStyle({"model": model_index}, {})
        if focus_on_site and site:
            pocket_sel = {
                **chain_sel,
                "byres": True,
                "within": {"distance": focus_radius, "sel": {**chain_sel, "resi": list(site)}},
            }
            view.setStyle(pocket_sel, {"stick": {"color": color, "radius": 0.15}})
        else:
            view.setStyle(chain_sel, {"cartoon": {"color": color}})

        if site and s.get("highlight_site", True):
            view.addStyle(
                {"model": model_index, "chain": s["chain_id"], "resi": list(site)},
                {"stick": {"color": SITE_COLOR, "radius": 0.25}},
            )

        if s.get("show_ligand"):
            groups = classify_hetero_groups(collect_hetero_groups(pdb_text))
            ligand = pick_ligand_of_interest(groups)
            if ligand is not None:
                view.addStyle(
                    {"model": model_index, "chain": ligand.chain, "resi": ligand.resseq, "resn": ligand.resname},
                    {"stick": {"color": color, "radius": 0.3}},
                )

    view.zoomTo()
    return view
