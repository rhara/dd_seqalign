"""Shared fixtures for dd_viewer's ported pytest suite (test_interactions_*,
test_io_integration, test_scene, test_scoring). dd_seqalign's own tests
(test_sequence.py) are pure synthetic-sequence unit tests and don't use
this fixture.
"""
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).resolve().parent / "dd_viewer_data"


@pytest.fixture
def data_dir() -> Path:
    """Directory holding the small bundled 6W63 sample PDB/SDF files."""
    return DATA_DIR
