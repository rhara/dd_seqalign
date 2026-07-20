"""py3Dmol scene construction, usable from both the Streamlit app and Jupyter.

`build_view` returns a plain `py3Dmol.view` object: call `.show()` on it in a
notebook, or `._make_html()` to embed it (e.g. via `st.iframe`, as `app.py`
does). Receptor and ligand pose are each optional so the same function covers
protein-only, ligand-only, and complex views.
"""
import json
import re
from typing import Optional, Sequence

import py3Dmol
from rdkit import Chem

from .interactions import Contact, RingContact, find_contact_residues
from .io import Receptor

RECEPTOR_STYLES = ("cartoon", "stick", "surface")

LIGAND_COLORSCHEME = "cyanCarbon"
REFERENCE_COLORSCHEME = "greenCarbon"
HIGHLIGHT_COLORSCHEME = "yellowCarbon"
SELECTED_COLORSCHEME = "magentaCarbon"

# cartoon: translucent rather than fully opaque, so overlaid interaction
# lines/highlighted-residue sticks/labels near the backbone stay visible
# instead of being hidden behind the ribbon.
RECEPTOR_CARTOON_OPACITY = 0.8
# stick: thinner than the ligand's own sticks (radius 0.2, see _add_mol)
# so the receptor reads as a supporting backdrop rather than competing with
# the ligand/interaction overlays for attention -- thinner than a first
# pass at this (0.08) per user feedback, still thicker than 3Dmol's
# non-radius-configurable "line" style.
RECEPTOR_STICK_RADIUS = 0.05
# surface: translucent enough that interaction dashes at the ligand-facing
# side of the surface (right where the surface bulges closest to the
# ligand) aren't fully hidden behind the mesh -- they're always drawn
# regardless of receptor_style, but a less-transparent surface can occlude
# them from most camera angles even though they're technically present.
RECEPTOR_SURFACE_OPACITY = 0.5
RECEPTOR_SURFACE_COLOR = "white"

HBOND_COLOR = "yellow"
HBOND_RADIUS = 0.06

# Kept visually subordinate to hbonds/salt bridges (thinner radius, and
# still slightly translucent below) since hydrophobic contacts are the
# weakest, least specific signal (any nearby carbon pair) and there are
# usually many more of them -- but not washed out either, per user
# feedback that the original 0.4 read as too pale.
HYDROPHOBIC_COLOR = "orange"
HYDROPHOBIC_RADIUS = 0.03
HYDROPHOBIC_OPACITY = 0.6

SALT_BRIDGE_COLOR = "red"
SALT_BRIDGE_RADIUS = 0.07

PI_STACKING_COLOR = "purple"
PI_STACKING_RADIUS = 0.07

ELECTROSTATIC_COLOR = "blue"
ELECTROSTATIC_RADIUS = 0.05
ELECTROSTATIC_OPACITY = 0.9

PI_HALOGEN_COLOR = "teal"
PI_HALOGEN_RADIUS = 0.07

SULFUR_HALOGEN_COLOR = "brown"
SULFUR_HALOGEN_RADIUS = 0.07

LABEL_FONT_COLOR = "white"
LABEL_BACKGROUND_COLOR = "black"
LABEL_BACKGROUND_OPACITY = 0.7
LABEL_FONT_SIZE = 12


def _add_receptor(
    view, receptor: Receptor, style: str, model_index: int,
    color_by_secondary_structure: bool = False, residue_keys: Optional[Sequence[tuple]] = None,
) -> None:
    """`residue_keys` (e.g. from `find_contact_residues`), when given,
    restricts the receptor's rendered style to just those (chain, resnum)
    residues instead of the whole model -- one `setStyle`/`addSurface` call
    per residue rather than a single model-wide selector, so a multi-chain
    receptor only shows style on the requested chain/resnum combination
    (a single combined `{"chain": [...], "resi": [...]}` selector would
    instead match the cross product of every listed chain with every listed
    resnum). See `build_view`'s `only_near_ligand`.
    """
    if style not in RECEPTOR_STYLES:
        raise ValueError(f"receptor_style must be one of {RECEPTOR_STYLES}, got {style!r}")

    view.addModel(receptor.pdb_text, "pdb")

    if style == "cartoon":
        cartoon_style = {"colorscheme": "ssPyMol"} if color_by_secondary_structure else {"color": "spectrum"}
        cartoon_style["opacity"] = RECEPTOR_CARTOON_OPACITY
        style_spec = {"cartoon": cartoon_style}
    elif style == "stick":
        style_spec = {"stick": {"radius": RECEPTOR_STICK_RADIUS, "colorscheme": "grayCarbon"}}
    else:
        style_spec = {}

    sels = [{"model": model_index, "chain": chain, "resi": int(resnum)} for chain, resnum in residue_keys] \
        if residue_keys is not None else [{"model": model_index}]
    for sel in sels:
        view.setStyle(sel, style_spec)
        if style == "surface":
            view.addSurface(py3Dmol.VDW, {"opacity": RECEPTOR_SURFACE_OPACITY, "color": RECEPTOR_SURFACE_COLOR}, sel)


def _add_mol(view, mol: Chem.Mol, conf_id: int, model_index: int, colorscheme: str, radius: float) -> None:
    view.addModel(Chem.MolToMolBlock(mol, confId=conf_id), "mol")
    view.setStyle({"model": model_index}, {"stick": {"colorscheme": colorscheme, "radius": radius}})


def _style_residues(view, receptor_model_index: int, residue_keys: Sequence[tuple], colorscheme: str, radius: float) -> None:
    for chain, resnum in residue_keys:
        view.addStyle(
            {"model": receptor_model_index, "chain": chain, "resi": resnum},
            {"stick": {"colorscheme": colorscheme, "radius": radius}},
        )


def _add_residue_labels(view, receptor: Receptor, residue_labels: dict) -> None:
    """One floating text label per (chain, resnum) key in `residue_labels`
    (e.g. `{("A", 235): "S235"}`), anchored at that residue's CA atom.
    Residues missing a CA (a HETATM, or a non-standard/incomplete one) are
    silently skipped rather than raising, since a label is a minor visual
    aid, not something that should abort the whole scene over one residue.
    """
    ca = receptor.atoms[receptor.atoms["name"] == "CA"]
    for (chain, resnum), text in residue_labels.items():
        match = ca[(ca["chain"] == chain) & (ca["resnum"] == resnum)]
        if match.empty:
            continue
        row = match.iloc[0]
        view.addLabel(text, {
            "position": {"x": float(row["x"]), "y": float(row["y"]), "z": float(row["z"])},
            "fontColor": LABEL_FONT_COLOR,
            "backgroundColor": LABEL_BACKGROUND_COLOR,
            "backgroundOpacity": LABEL_BACKGROUND_OPACITY,
            "fontSize": LABEL_FONT_SIZE,
            "showBackground": True,
        })


def _draw_contact_lines(
    view, pose_mol: Chem.Mol, contacts: Sequence[Contact], color: str, conf_id: int,
    radius: float = 0.06, opacity: float = 1.0,
) -> None:
    conf = pose_mol.GetConformer(conf_id)
    for contact in contacts:
        lig_pos = conf.GetAtomPosition(contact.ligand_atom_idx)
        view.addCylinder({
            "start": {"x": lig_pos.x, "y": lig_pos.y, "z": lig_pos.z},
            "end": {"x": contact.rec_x, "y": contact.rec_y, "z": contact.rec_z},
            "radius": radius,
            "color": color,
            "opacity": opacity,
            "dashed": True,
        })


def _draw_ring_lines(view, contacts: Sequence[RingContact], color: str, radius: float, opacity: float = 1.0) -> None:
    for contact in contacts:
        view.addCylinder({
            "start": {"x": contact.lig_x, "y": contact.lig_y, "z": contact.lig_z},
            "end": {"x": contact.rec_x, "y": contact.rec_y, "z": contact.rec_z},
            "radius": radius,
            "color": color,
            "opacity": opacity,
            "dashed": True,
        })


def build_view(
    receptor: Optional[Receptor] = None,
    pose_mol: Optional[Chem.Mol] = None,
    receptor_style: str = "cartoon",
    color_by_secondary_structure: bool = False,
    highlight_residues: Optional[Sequence[tuple]] = None,
    selected_residues: Optional[Sequence[tuple]] = None,
    residue_labels: Optional[dict] = None,
    hbonds: Optional[Sequence[Contact]] = None,
    hydrophobic_contacts: Optional[Sequence[Contact]] = None,
    salt_bridges: Optional[Sequence[Contact]] = None,
    pi_stacking: Optional[Sequence[RingContact]] = None,
    electrostatic: Optional[Sequence[Contact]] = None,
    pi_halogen_bonds: Optional[Sequence[Contact]] = None,
    sulfur_halogen_bonds: Optional[Sequence[Contact]] = None,
    reference_mol: Optional[Chem.Mol] = None,
    conf_id: int = 0,
    width: int = 900,
    height: int = 600,
    only_near_ligand: bool = False,
    near_ligand_cutoff: float = 5.0,
) -> py3Dmol.view:
    """Build a scene from any combination of a receptor, a docked ligand
    pose, and a reference ligand -- at least one of `receptor` / `pose_mol`
    must be given.

    `receptor_style` picks the receptor representation ("cartoon", "stick",
    or "surface" -- a translucent surface on its own, no cartoon underneath
    it), and is ignored if `receptor` is omitted. `cartoon` and `surface`
    are both drawn partly transparent (`RECEPTOR_CARTOON_OPACITY`/
    `RECEPTOR_SURFACE_OPACITY`) so overlaid interaction lines, highlighted-
    residue sticks, and labels stay visible through them rather than being
    hidden behind an opaque receptor -- interaction lines are always drawn
    regardless of `receptor_style` (see below), so a lower surface opacity
    is what actually keeps them visible in "surface" mode, not a separate
    code path. When `receptor_style` is "cartoon", `color_by_secondary_structure`
    switches from the default rainbow (N-to-C "spectrum") coloring to
    3Dmol.js's "ssPyMol" scheme (helix / sheet / coil colored by secondary
    structure); ignored otherwise. `pose_mol` is always drawn as cyan
    sticks; `reference_mol` (e.g. a crystal pose) is drawn alongside it as
    thinner green sticks so the two stay visually distinct.

    `highlight_residues` (from `interactions.find_contact_residues`) and
    `selected_residues` (e.g. rows a user picked in a residue table) are both
    lists of (chain, resnum) pairs rendered as extra sticks on top of the
    base receptor style -- yellow and magenta respectively, so an explicit
    table selection stays visible even when it overlaps the auto-detected
    contact residues. Both require `receptor`.

    `residue_labels` (e.g. `{("A", 235): "S235"}`) draws one floating text
    label per (chain, resnum) key, anchored at that residue's CA atom --
    independent of `highlight_residues`/`selected_residues` (a caller
    typically labels only the ones it wants named, e.g. the current
    `selected_residues`, to avoid cluttering the view when there are many
    auto-detected `highlight_residues`). Requires `receptor`; a key with no
    matching CA atom (HETATM, or a non-standard/incomplete residue) is
    silently skipped.

    `hbonds`, `hydrophobic_contacts`, `salt_bridges`, `electrostatic`,
    `pi_halogen_bonds`, `sulfur_halogen_bonds` (all from the matching
    `interactions.find_*` function) and `pi_stacking` (from
    `interactions.find_pi_stacking`) are drawn as dashed cylinders between
    the ligand and receptor -- atom-to-atom for all but `pi_stacking` and
    `pi_halogen_bonds`, which go ligand-atom/ring-centroid-to-ring-centroid.
    All require both `receptor` and `pose_mol`. Hydrophobic contacts are
    drawn thin and slightly translucent since they're the weakest, most
    numerous signal, and electrostatic interactions (a softer, partial-
    charge-based signal than salt bridges) are drawn at a middling
    opacity; the rest are drawn bolder and fully opaque.

    `only_near_ligand` restricts the receptor's rendered style (cartoon/
    stick/surface) to just the residues within `near_ligand_cutoff` Å of
    `pose_mol` (via `interactions.find_contact_residues`, independent of
    whatever cutoff a caller used for `highlight_residues`), instead of the
    whole receptor -- useful for decluttering a view of a large or
    multi-chain structure down to just the binding pocket. Requires both
    `receptor` and `pose_mol`; a no-op otherwise (there's no ligand to be
    "near").
    """
    if receptor is None and pose_mol is None:
        raise ValueError("build_view needs at least one of receptor or pose_mol")

    view = py3Dmol.view(width=width, height=height)

    next_index = 0
    receptor_index = None
    pose_index = None
    if receptor is not None:
        near_ligand_residues = None
        if only_near_ligand and pose_mol is not None:
            contact_df = find_contact_residues(receptor, pose_mol, cutoff=near_ligand_cutoff, conf_id=conf_id)
            near_ligand_residues = list(zip(contact_df["chain"], contact_df["resnum"]))
        _add_receptor(view, receptor, receptor_style, next_index, color_by_secondary_structure, residue_keys=near_ligand_residues)
        receptor_index = next_index
        next_index += 1
    if pose_mol is not None:
        _add_mol(view, pose_mol, conf_id, next_index, LIGAND_COLORSCHEME, radius=0.2)
        pose_index = next_index
        next_index += 1
    if reference_mol is not None:
        _add_mol(view, reference_mol, 0, next_index, REFERENCE_COLORSCHEME, radius=0.12)
        next_index += 1

    if receptor_index is not None:
        if highlight_residues:
            _style_residues(view, receptor_index, highlight_residues, HIGHLIGHT_COLORSCHEME, radius=0.15)
        if selected_residues:
            _style_residues(view, receptor_index, selected_residues, SELECTED_COLORSCHEME, radius=0.22)
        if residue_labels:
            _add_residue_labels(view, receptor, residue_labels)
    if receptor_index is not None and pose_index is not None:
        if hbonds:
            _draw_contact_lines(view, pose_mol, hbonds, HBOND_COLOR, conf_id, radius=HBOND_RADIUS)
        if hydrophobic_contacts:
            _draw_contact_lines(
                view, pose_mol, hydrophobic_contacts, HYDROPHOBIC_COLOR, conf_id,
                radius=HYDROPHOBIC_RADIUS, opacity=HYDROPHOBIC_OPACITY,
            )
        if salt_bridges:
            _draw_contact_lines(view, pose_mol, salt_bridges, SALT_BRIDGE_COLOR, conf_id, radius=SALT_BRIDGE_RADIUS)
        if pi_stacking:
            _draw_ring_lines(view, pi_stacking, PI_STACKING_COLOR, radius=PI_STACKING_RADIUS)
        if electrostatic:
            _draw_contact_lines(
                view, pose_mol, electrostatic, ELECTROSTATIC_COLOR, conf_id,
                radius=ELECTROSTATIC_RADIUS, opacity=ELECTROSTATIC_OPACITY,
            )
        if pi_halogen_bonds:
            _draw_contact_lines(view, pose_mol, pi_halogen_bonds, PI_HALOGEN_COLOR, conf_id, radius=PI_HALOGEN_RADIUS)
        if sulfur_halogen_bonds:
            _draw_contact_lines(
                view, pose_mol, sulfur_halogen_bonds, SULFUR_HALOGEN_COLOR, conf_id, radius=SULFUR_HALOGEN_RADIUS,
            )

    if pose_index is not None:
        view.zoomTo({"model": pose_index})
    else:
        view.zoomTo({"model": receptor_index})
    return view


_VIEWER_RENDER_CALL = re.compile(r"viewer_(\d+)\.render\(\);")


def get_viewer_variable(html: str) -> Optional[str]:
    """The py3Dmol-generated JS variable name (e.g. `"viewer_123456"`) that
    a `_make_html()` string's scene is loaded into, or `None` if `html`
    doesn't look like a py3Dmol scene at all. Each `_make_html()` call
    picks a fresh random suffix, so this has to be re-extracted per HTML
    string rather than assumed constant -- used both by
    `html_with_camera_events`/`html_with_initial_view` (to know which
    variable to patch calls onto) and by callers wanting to query a
    *previously* embedded scene's live state (e.g. the desktop app calling
    `page().runJavaScript(f"{var}.getView()")` on the still-loaded page
    before replacing it, to carry the camera position across a `setHtml()`
    reload).
    """
    match = _VIEWER_RENDER_CALL.search(html)
    return f"viewer_{match.group(1)}" if match else None


def html_with_initial_view(html: str, view: list) -> str:
    """Patch a py3Dmol `_make_html()` string to apply a saved camera view
    (a `getView()`-shaped array) right after the scene renders, instead of
    leaving `build_view`'s own default `zoomTo()` auto-fit as the final
    word. For callers that reload the whole scene on every change (e.g.
    the desktop app's `QWebEngineView.setHtml()`, torn down and rebuilt
    each time rather than mutated in place) and want the camera to *not*
    visibly jump on every reload -- capture the outgoing scene's view with
    `get_viewer_variable` + `page().runJavaScript(...)` before replacing
    it, then feed that view into the next `build_view(...)._make_html()`
    through this function.

    A no-op (returns `html` unchanged) if `html` doesn't look like a
    py3Dmol scene, or if `view` is falsy (`None` or `[]` -- both mean "no
    saved view to restore", e.g. the very first render of a session).
    """
    if not view:
        return html
    viewer_var = get_viewer_variable(html)
    if viewer_var is None:
        return html
    match = _VIEWER_RENDER_CALL.search(html)
    snippet = f"\n{viewer_var}.setView({json.dumps(view)});\n"
    return html.replace(match.group(0), match.group(0) + snippet, 1)


def html_with_camera_events(html: str) -> str:
    """Patch a py3Dmol `_make_html()` string to report its camera state and
    render-readiness to the parent window via postMessage, for the
    `plviewer_3d` double-buffered component (`dd_viewer.component`):

    - On every drag/zoom/touch, posts `{plviewerCameraUpdate: true, view:
      [...]}` (the 3Dmol `getView()` array). The component's own JS keeps
      the latest one in memory (its execution context survives Streamlit
      reruns, unlike this scene's) and re-applies it to each newly-loaded
      scene -- which is what makes the camera position survive a widget
      interaction instead of snapping back to the default zoomTo fit every
      time. (An earlier version of this used `sessionStorage` instead, but
      each srcdoc-loaded scene turned out to get its own isolated storage
      partition even when nested in the same component, so the saved state
      never made it to the next scene; postMessage into the component's
      always-alive JS sidesteps that entirely.)
    - Two animation frames after the initial render (once the paint has
      actually landed), posts `{plviewerReady: true}` -- the signal the
      component waits for before swapping this scene into view, so updates
      read as a cross-fade instead of a flash to blank.

    Both are harmless no-ops if nothing is listening (e.g. a plain
    `st.iframe` embed, or Jupyter).
    """
    viewer_var = get_viewer_variable(html)
    if viewer_var is None:
        return html
    match = _VIEWER_RENDER_CALL.search(html)
    snippet = f"""
try {{
  var __plvSave = function() {{
    try {{ parent.postMessage({{plviewerCameraUpdate: true, view: {viewer_var}.getView()}}, "*"); }} catch (e) {{}}
  }};
  document.addEventListener('mouseup', __plvSave);
  document.addEventListener('wheel', __plvSave, {{passive: true}});
  document.addEventListener('touchend', __plvSave);
}} catch (e) {{}}
requestAnimationFrame(function() {{
  requestAnimationFrame(function() {{
    try {{ parent.postMessage({{plviewerReady: true}}, "*"); }} catch (e) {{}}
  }});
}});
"""
    return html.replace(match.group(0), match.group(0) + snippet, 1)


_FIXED_SIZE_DIV_STYLE = re.compile(r'style="position: relative; width: \d+px; height: \d+px;"')


def html_fill_container(html: str) -> str:
    """Patch a py3Dmol `_make_html()` string so its viewer div fills 100%
    of whatever element embeds it, instead of the fixed pixel `width`/
    `height` baked in from the `width=`/`height=` passed to
    `build_view(...)` (900x600 by default). A fixed pixel size only
    matches its embedding container by coincidence -- e.g. the desktop
    app's `QWebEngineView`, which is resized by its splitter to whatever
    the window's actual size is, ends up with the 3D scene rendered in a
    fixed-size box in one corner and a plain white margin around it.
    `html`/`body` also need an explicit reset (browsers default `body` to
    `height: auto` and a small margin, so a `height: 100%` div would
    otherwise have nothing to be 100% relative *to*).

    Safe to call on the raw `_make_html()` output or on the result of
    `html_with_camera_events` (order doesn't matter -- they touch
    different parts of the string). Not used by the Streamlit app's
    `view3d` embedding, whose component sizes itself to a fixed `height`
    passed by the caller rather than to the browser window.
    """
    html = _FIXED_SIZE_DIV_STYLE.sub('style="position: relative; width: 100%; height: 100%;"', html, count=1)
    reset = "<style>html, body { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }</style>"
    return reset + html
