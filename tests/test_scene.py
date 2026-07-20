"""Tests for the pure-string HTML post-processing helpers in dd_viewer.scene
(get_viewer_variable, html_with_initial_view, html_fill_container). Not
covered: build_view/_make_html's actual rendered output, and
html_with_camera_events's postMessage behavior -- both require visual/
browser-runtime verification, per this project's manual-only convention
for scene.py's WebGL-rendering side.
"""
from dd_viewer.scene import get_viewer_variable, html_fill_container, html_with_initial_view

_SCENE_HTML = (
    '<div id="3dmolviewer_123" style="position: relative; width: 900px; height: 600px;"></div>\n'
    "<script>\n"
    "var viewer_123456789012345 = null;\n"
    "$3Dmolpromise.then(function() {\n"
    "viewer_123456789012345 = $3Dmol.createViewer(...);\n"
    "\tviewer_123456789012345.render();\n"
    "});\n"
    "</script>"
)


class TestGetViewerVariable:
    def test_extracts_the_viewer_variable_name(self):
        assert get_viewer_variable(_SCENE_HTML) == "viewer_123456789012345"

    def test_returns_none_for_non_scene_html(self):
        assert get_viewer_variable("<p>not a scene</p>") is None


class TestHtmlWithInitialView:
    def test_injects_setview_after_render_call(self):
        patched = html_with_initial_view(_SCENE_HTML, [1.0, 2.0, 3.0])
        assert "viewer_123456789012345.render();" in patched
        assert "viewer_123456789012345.setView([1.0, 2.0, 3.0]);" in patched
        # setView must come after render(), not before
        assert patched.index("render();") < patched.index("setView(")

    def test_noop_when_view_is_none(self):
        assert html_with_initial_view(_SCENE_HTML, None) == _SCENE_HTML

    def test_noop_when_view_is_empty(self):
        assert html_with_initial_view(_SCENE_HTML, []) == _SCENE_HTML

    def test_noop_for_non_scene_html(self):
        html = "<p>not a scene</p>"
        assert html_with_initial_view(html, [1.0, 2.0]) == html


class TestHtmlFillContainer:
    def test_rewrites_fixed_pixel_size_to_100_percent(self):
        patched = html_fill_container(_SCENE_HTML)
        assert "width: 900px" not in patched
        assert 'width: 100%; height: 100%;"' in patched

    def test_prepends_html_body_reset_style(self):
        patched = html_fill_container(_SCENE_HTML)
        assert patched.startswith("<style>html, body")
