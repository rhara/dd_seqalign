"""Double-buffered Streamlit embedding for py3Dmol scenes.

`st.iframe(html)` replaces the iframe's srcdoc on every Streamlit rerun; the
browser tears down the old document immediately and the new one takes a
moment to load 3Dmol.js, re-parse the model, and re-render, producing a
visible blank flash on every widget interaction (a residue-table selection,
a style checkbox, ...).

`view3d` embeds scenes through a small static Streamlit component
(`frontend/index.html`) instead. Because it's a *static* component, its own
outer iframe has a fixed local URL and is loaded exactly once; Streamlit
talks to an already-running instance via postMessage rather than reloading
it, so the component's JS can stay alive across reruns. It keeps two hidden/
visible inner iframes and only swaps which one is visible once the new scene
reports (via the postMessage `scene.html_with_camera_events` injects) that
it has actually finished painting -- so updates read as a quick cross-fade
instead of a flash to blank. That same always-alive JS context is also what
lets the camera position (rotation/zoom) survive across scenes: it's kept in
a plain JS variable in the component, not in the scene's own (short-lived,
storage-partitioned) iframe.
"""
import os

import streamlit.components.v1 as components

_FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "frontend")
_component = components.declare_component("plviewer_3d", path=_FRONTEND_DIR)


def view3d(html: str, height: int = 650, key: str = "plviewer_3d_view", reset_camera_token: int = 0) -> None:
    """Embed a py3Dmol scene's HTML (typically already passed through
    `html_with_camera_events`) via the double-buffered component.

    `key` should stay constant across reruns (the default already does) --
    it's what lets Streamlit recognize this as the same live component
    instance rather than tearing down and recreating it, which would defeat
    the whole point. `reset_camera_token` forgets the saved camera position
    (falling back to the scene's own zoomTo fit) whenever it changes from
    the previous call -- wire it to e.g. a counter incremented by a "reset
    view" button.
    """
    _component(html=html, height=height, key=key, reset_token=reset_camera_token, default=None)
