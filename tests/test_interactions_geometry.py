"""Tests for the deterministic pure-function geometry helpers in
`dd_viewer.interactions`: plane-normal computation, the linearity-angle
helper used by the halogen-bond finders, and the nearest-neighbor dedup
logic used by `find_hydrophobic_contacts`.

These are pure numpy functions on plain coordinate arrays, so they're tested
directly with small synthetic inputs rather than via RDKit Mols or the
bundled PDB/SDF sample data.
"""
import numpy as np
import pytest

from dd_viewer.interactions import (
    Contact,
    _dedupe_nearest_per_ligand_atom,
    _linearity_angle,
    _plane_normal,
)


class TestPlaneNormal:
    def test_normal_of_xy_plane_points_along_z(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        normal = _plane_normal(coords)
        assert normal == pytest.approx([0.0, 0.0, 1.0], abs=1e-9)

    def test_normal_is_unit_length(self):
        coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
        normal = _plane_normal(coords)
        assert np.linalg.norm(normal) == pytest.approx(1.0)

    def test_normal_of_tilted_plane(self):
        # Ring tilted 90 degrees from the xy-plane: lies in the xz-plane instead,
        # so its normal should point along y.
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        normal = _plane_normal(coords)
        assert abs(normal[1]) == pytest.approx(1.0)
        assert normal[0] == pytest.approx(0.0, abs=1e-9)
        assert normal[2] == pytest.approx(0.0, abs=1e-9)

    def test_degenerate_collinear_points_returns_zero_vector_without_error(self):
        # Three collinear points don't define a plane; cross product is ~0.
        # The function should not raise (division guarded by norm > 0 check).
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        normal = _plane_normal(coords)
        assert np.linalg.norm(normal) == pytest.approx(0.0, abs=1e-9)


class TestLinearityAngle:
    def test_perfectly_linear_arrangement_is_zero_degrees(self):
        # C->X and X->Acceptor point the same direction: ideal sigma-hole geometry.
        cx_vector = np.array([1.0, 0.0, 0.0])
        x_to_acceptor = np.array([1.0, 0.0, 0.0])
        assert _linearity_angle(cx_vector, x_to_acceptor) == pytest.approx(0.0, abs=1e-6)

    def test_perpendicular_arrangement_is_ninety_degrees(self):
        cx_vector = np.array([1.0, 0.0, 0.0])
        x_to_acceptor = np.array([0.0, 1.0, 0.0])
        assert _linearity_angle(cx_vector, x_to_acceptor) == pytest.approx(90.0)

    def test_opposite_arrangement_is_180_degrees(self):
        cx_vector = np.array([1.0, 0.0, 0.0])
        x_to_acceptor = np.array([-1.0, 0.0, 0.0])
        assert _linearity_angle(cx_vector, x_to_acceptor) == pytest.approx(180.0)

    def test_zero_length_vector_returns_180_without_dividing_by_zero(self):
        cx_vector = np.array([0.0, 0.0, 0.0])
        x_to_acceptor = np.array([1.0, 0.0, 0.0])
        assert _linearity_angle(cx_vector, x_to_acceptor) == 180.0

    def test_scale_invariant(self):
        # Only direction should matter, not magnitude.
        cx_vector = np.array([2.0, 0.0, 0.0])
        x_to_acceptor = np.array([5.0, 5.0, 0.0])
        angle = _linearity_angle(cx_vector, x_to_acceptor)
        assert angle == pytest.approx(45.0)


def _make_contact(ligand_atom_idx: int, distance: float) -> Contact:
    return Contact(
        ligand_atom_idx=ligand_atom_idx, chain="A", resnum=1, resname="ALA",
        atom_name="CA", distance=distance, rec_x=0.0, rec_y=0.0, rec_z=0.0,
    )


class TestDedupeNearestPerLigandAtom:
    def test_keeps_only_closest_contact_per_ligand_atom(self):
        contacts = [
            _make_contact(ligand_atom_idx=0, distance=3.0),
            _make_contact(ligand_atom_idx=0, distance=1.5),
            _make_contact(ligand_atom_idx=0, distance=2.0),
        ]
        result = _dedupe_nearest_per_ligand_atom(contacts)
        assert len(result) == 1
        assert result[0].distance == pytest.approx(1.5)

    def test_keeps_one_contact_per_distinct_ligand_atom(self):
        contacts = [
            _make_contact(ligand_atom_idx=0, distance=2.0),
            _make_contact(ligand_atom_idx=1, distance=4.0),
            _make_contact(ligand_atom_idx=0, distance=1.0),
            _make_contact(ligand_atom_idx=1, distance=3.5),
        ]
        result = _dedupe_nearest_per_ligand_atom(contacts)
        distances_by_idx = {c.ligand_atom_idx: c.distance for c in result}
        assert distances_by_idx == {0: pytest.approx(1.0), 1: pytest.approx(3.5)}

    def test_empty_input_returns_empty_list(self):
        assert _dedupe_nearest_per_ligand_atom([]) == []
