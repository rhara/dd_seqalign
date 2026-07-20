"""Tests for `dd_viewer.scoring`: docking-score property detection and
pose-vs-reference RMSD.

Poses/molecules are small synthetic RDKit Mols built in-test.
"""
import pytest
from rdkit import Chem

from dd_viewer.io import Pose
from dd_viewer.scoring import detect_score, rmsd_to_reference


def _make_pose(properties: dict, index: int = 0, label: str = "pose") -> Pose:
    mol = Chem.MolFromSmiles("CCO")
    return Pose(mol=mol, index=index, label=label, properties=properties)


class TestDetectScore:
    def test_finds_minimized_affinity_smina_style(self):
        pose = _make_pose({"minimizedAffinity": "-7.2"})
        assert detect_score(pose) == pytest.approx(-7.2)

    def test_priority_order_prefers_earlier_candidate(self):
        # minimizedAffinity outranks a generic "Score" key when both are present.
        pose = _make_pose({"Score": "-1.0", "minimizedAffinity": "-9.5"})
        assert detect_score(pose) == pytest.approx(-9.5)

    def test_falls_back_to_later_candidate_if_earlier_absent(self):
        pose = _make_pose({"docking_score": "-4.4"})
        assert detect_score(pose) == pytest.approx(-4.4)

    def test_non_numeric_value_is_skipped_in_favor_of_next_candidate(self):
        pose = _make_pose({"minimizedAffinity": "not-a-number", "Score": "-3.3"})
        assert detect_score(pose) == pytest.approx(-3.3)

    def test_no_matching_property_returns_none(self):
        pose = _make_pose({"some_other_field": "1.0"})
        assert detect_score(pose) is None

    def test_empty_properties_returns_none(self):
        pose = _make_pose({})
        assert detect_score(pose) is None


class TestRmsdToReference:
    def _embedded_ethanol(self, coord_scale: float = 1.0) -> Chem.Mol:
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(mol, randomSeed=42)
        if coord_scale != 1.0:
            conf = mol.GetConformer()
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                conf.SetAtomPosition(i, (pos.x * coord_scale, pos.y * coord_scale, pos.z * coord_scale))
        return mol

    def test_identical_pose_has_zero_rmsd(self):
        mol = self._embedded_ethanol()
        rmsd = rmsd_to_reference(mol, Chem.Mol(mol))
        assert rmsd == pytest.approx(0.0, abs=1e-3)

    def test_displaced_pose_has_positive_rmsd(self):
        reference = self._embedded_ethanol()
        probe = Chem.Mol(reference)
        conf = probe.GetConformer()
        for i in range(probe.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pos.x + 2.0, pos.y, pos.z))
        rmsd = rmsd_to_reference(probe, reference)
        assert rmsd == pytest.approx(2.0, abs=1e-3)

    def test_incompatible_molecules_return_none(self):
        probe = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(probe, randomSeed=1)
        reference = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(reference, randomSeed=1)
        assert rmsd_to_reference(probe, reference) is None
