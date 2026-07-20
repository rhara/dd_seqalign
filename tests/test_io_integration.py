"""Integration-style tests for `dd_viewer.io`'s loading functions, using
the small bundled 6W63 sample data in `data/` (a real receptor PDB, a raw
co-crystal PDB, and docking-pose / reference-ligand SDFs).

These complement the pure-function/synthetic-Mol tests elsewhere by
exercising the actual Bio.PDB/RDKit parsing path end-to-end.
"""
import pandas as pd
import pytest

from dd_viewer.io import (
    ATOM_COLUMNS,
    Pose,
    Receptor,
    load_poses,
    load_receptor,
    load_reference_ligand,
    split_structure,
)


class TestLoadReceptor:
    def test_returns_receptor_with_expected_columns(self, data_dir):
        receptor = load_receptor(str(data_dir / "6W63_receptor.pdb"))
        assert isinstance(receptor, Receptor)
        assert list(receptor.atoms.columns) == ATOM_COLUMNS
        assert len(receptor.atoms) > 0

    def test_pdb_text_is_preserved(self, data_dir):
        receptor = load_receptor(str(data_dir / "6W63_receptor.pdb"))
        assert "ATOM" in receptor.pdb_text

    def test_source_path_is_recorded(self, data_dir):
        path = str(data_dir / "6W63_receptor.pdb")
        receptor = load_receptor(path)
        assert receptor.source == path


class TestLoadPoses:
    def test_loads_at_least_one_pose(self, data_dir):
        poses = load_poses(str(data_dir / "6W63_redock.sdf"))
        assert len(poses) > 0
        assert all(isinstance(p, Pose) for p in poses)

    def test_pose_indices_match_file_order(self, data_dir):
        poses = load_poses(str(data_dir / "6W63_redock.sdf"))
        assert [p.index for p in poses] == list(range(len(poses)))

    def test_each_pose_has_a_mol_with_atoms(self, data_dir):
        poses = load_poses(str(data_dir / "6W63_redock.sdf"))
        for pose in poses:
            assert pose.mol.GetNumAtoms() > 0


class TestLoadReferenceLigand:
    def test_loads_a_single_mol(self, data_dir):
        mol = load_reference_ligand(str(data_dir / "6W63_ligand_ref.sdf"))
        assert mol is not None
        assert mol.GetNumAtoms() > 0


class TestSplitStructure:
    def test_splits_raw_pdb_into_receptor_and_ligand_poses(self, data_dir):
        receptor, poses = split_structure(str(data_dir / "6W63_raw.pdb"))
        assert isinstance(receptor, Receptor)
        assert len(receptor.atoms) > 0
        # Protein-only receptor: no HETATM rows should remain.
        assert not receptor.atoms["is_hetero"].any()

    def test_extracted_poses_have_auto_extracted_label(self, data_dir):
        _, poses = split_structure(str(data_dir / "6W63_raw.pdb"))
        if poses:
            assert all("auto-extracted from PDB" in p.label for p in poses)

    def test_extracted_pose_indices_are_sequential(self, data_dir):
        _, poses = split_structure(str(data_dir / "6W63_raw.pdb"))
        assert [p.index for p in poses] == list(range(len(poses)))
