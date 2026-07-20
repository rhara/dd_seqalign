from .io import Receptor, Pose, load_receptor, load_poses, load_reference_ligand, split_structure
from .scoring import detect_score, rmsd_to_reference
from .interactions import (
    find_contact_residues,
    find_hydrogen_bonds,
    find_hydrophobic_contacts,
    find_salt_bridges,
    find_pi_stacking,
    find_electrostatic_interactions,
    find_pi_halogen_bonds,
    find_sulfur_halogen_bonds,
)
from .scene import build_view, get_viewer_variable, html_fill_container, html_with_camera_events, html_with_initial_view
from .dashboard import poses_dataframe

__all__ = [
    "Receptor",
    "Pose",
    "load_receptor",
    "load_poses",
    "load_reference_ligand",
    "split_structure",
    "detect_score",
    "rmsd_to_reference",
    "find_contact_residues",
    "find_hydrogen_bonds",
    "find_hydrophobic_contacts",
    "find_salt_bridges",
    "find_pi_stacking",
    "find_electrostatic_interactions",
    "find_pi_halogen_bonds",
    "find_sulfur_halogen_bonds",
    "build_view",
    "get_viewer_variable",
    "html_fill_container",
    "html_with_camera_events",
    "html_with_initial_view",
    "view3d",
    "poses_dataframe",
]


def __getattr__(name):
    # `view3d` is the only symbol that needs `.component`, which imports
    # `streamlit` and declares a Streamlit custom component at import time --
    # that's a heavy, unnecessary cost (and it trips Streamlit's "missing
    # ScriptRunContext" warning) for callers like dd_molview's desktop app,
    # which only use `build_view`/`_make_html()` and never touch the
    # Streamlit-specific embedding. Deferred here so plain `import dd_viewer`
    # never pulls in streamlit at all.
    if name == "view3d":
        from .component import view3d as _view3d
        return _view3d
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
