"""Tests for the higher-level contact-finding functions in
`dd_viewer.interactions`, using small synthetic receptors/ligands built
in-test rather than the bundled PDB/SDF sample data.

A synthetic `Receptor` is built directly from a `pandas.DataFrame` (its
`pdb_text` is irrelevant to these functions -- only `.atoms` is read), and
ligand poses are built as small RDKit Mols from a molblock with an embedded
3D conformer, which keeps these tests fast and independent of external files.
"""
import pandas as pd
import pytest
from rdkit import Chem

from dd_viewer.interactions import (
    find_contact_residues,
    find_hydrogen_bonds,
    find_hydrophobic_contacts,
    find_pi_stacking,
    find_salt_bridges,
    find_sulfur_halogen_bonds,
)
from dd_viewer.io import ATOM_COLUMNS, Receptor


def _receptor_from_rows(rows) -> Receptor:
    atoms = pd.DataFrame(rows, columns=ATOM_COLUMNS)
    return Receptor(pdb_text="", atoms=atoms, source="<test>")


def _mol_from_molblock(molblock: str) -> Chem.Mol:
    mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    assert mol is not None, "test molblock failed to parse"
    return mol


# A single-atom-like "ligand" oxygen at the origin, used for straightforward
# distance-based checks against a nearby receptor atom.
_LIGAND_OH_MOLBLOCK = """
  test

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9600    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"""


class TestFindContactResidues:
    def test_finds_residue_within_cutoff(self):
        receptor = _receptor_from_rows([
            ("A", 10, "ALA", "CA", "C", 2.0, 0.0, 0.0, False),
            ("A", 50, "GLY", "CA", "C", 20.0, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(_LIGAND_OH_MOLBLOCK)
        residues = find_contact_residues(receptor, ligand, cutoff=4.5)
        assert list(residues["resnum"]) == [10]

    def test_no_contacts_when_everything_is_far(self):
        receptor = _receptor_from_rows([
            ("A", 10, "ALA", "CA", "C", 100.0, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(_LIGAND_OH_MOLBLOCK)
        residues = find_contact_residues(receptor, ligand, cutoff=4.5)
        assert len(residues) == 0


class TestFindHydrogenBonds:
    def test_polar_ligand_atom_near_polar_receptor_atom_is_flagged(self):
        receptor = _receptor_from_rows([
            ("A", 5, "SER", "OG", "O", 2.8, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(_LIGAND_OH_MOLBLOCK)
        hbonds = find_hydrogen_bonds(receptor, ligand, distance_cutoff=3.5)
        assert len(hbonds) == 1
        assert hbonds[0].resnum == 5
        assert hbonds[0].distance == pytest.approx(2.8)

    def test_carbon_receptor_atom_is_not_flagged_as_hbond(self):
        receptor = _receptor_from_rows([
            ("A", 5, "ALA", "CB", "C", 2.8, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(_LIGAND_OH_MOLBLOCK)
        hbonds = find_hydrogen_bonds(receptor, ligand, distance_cutoff=3.5)
        assert hbonds == []

    def test_beyond_cutoff_is_not_flagged(self):
        receptor = _receptor_from_rows([
            ("A", 5, "SER", "OG", "O", 10.0, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(_LIGAND_OH_MOLBLOCK)
        assert find_hydrogen_bonds(receptor, ligand, distance_cutoff=3.5) == []


class TestFindHydrophobicContacts:
    _LIGAND_METHANE_LIKE = """
  test

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""

    def test_nearest_only_dedupes_to_one_contact_per_ligand_atom(self):
        receptor = _receptor_from_rows([
            ("A", 1, "ALA", "CB", "C", 1.0, 0.0, 0.0, False),
            ("A", 2, "ALA", "CB", "C", 2.0, 0.0, 0.0, False),
            ("A", 3, "ALA", "CB", "C", 3.0, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(self._LIGAND_METHANE_LIKE)
        contacts = find_hydrophobic_contacts(receptor, ligand, cutoff=4.5, nearest_only=True)
        assert len(contacts) == 1
        assert contacts[0].resnum == 1  # closest

    def test_nearest_only_false_returns_all_raw_contacts(self):
        receptor = _receptor_from_rows([
            ("A", 1, "ALA", "CB", "C", 1.0, 0.0, 0.0, False),
            ("A", 2, "ALA", "CB", "C", 2.0, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(self._LIGAND_METHANE_LIKE)
        contacts = find_hydrophobic_contacts(receptor, ligand, cutoff=4.5, nearest_only=False)
        assert len(contacts) == 2


class TestFindSaltBridges:
    _LIGAND_CARBOXYLATE = """
  test

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.5000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -1.2000    0.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  2  0
  1  4  1  0
M  END
"""

    def test_negative_ligand_atom_near_positive_receptor_atom_is_flagged(self):
        receptor = _receptor_from_rows([
            ("A", 20, "LYS", "NZ", "N", 3.5, 0.5, 0.0, False),
        ])
        ligand = _mol_from_molblock(self._LIGAND_CARBOXYLATE)
        assert any(a.GetFormalCharge() != 0 for a in ligand.GetAtoms())
        bridges = find_salt_bridges(receptor, ligand, distance_cutoff=4.0)
        assert len(bridges) == 1
        assert bridges[0].resname == "LYS"

    def test_no_charged_ligand_atoms_returns_empty(self):
        receptor = _receptor_from_rows([
            ("A", 20, "LYS", "NZ", "N", 3.5, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(_LIGAND_OH_MOLBLOCK)
        assert find_salt_bridges(receptor, ligand, distance_cutoff=4.0) == []


class TestFindPiStacking:
    # A planar 6-membered aromatic ring (benzene-like) centered at the origin
    # in the xy-plane, matching PHE's ring atom names so the receptor side
    # can be built from plain rows.
    _LIGAND_BENZENE = """
  test

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.3900    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6950    1.2037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6950    1.2037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3900    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6950   -1.2037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6950   -1.2037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  6  1  4  0
M  END
"""

    def _phe_ring_rows(self, z_offset, resnum=30):
        names = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
        coords = [
            (1.39, 0.0), (0.695, 1.2037), (-0.695, 1.2037),
            (-1.39, 0.0), (-0.695, -1.2037), (0.695, -1.2037),
        ]
        return [
            ("A", resnum, "PHE", name, "C", x, y, z_offset)
            for name, (x, y) in zip(names, coords)
        ]

    def test_parallel_stacked_ring_is_detected(self):
        # Receptor ring directly "above" (along z) the ligand ring, both flat
        # in parallel planes -> face-to-face stacking, well within cutoff.
        receptor = _receptor_from_rows([(*row, False) for row in self._phe_ring_rows(z_offset=4.0)])
        ligand = _mol_from_molblock(self._LIGAND_BENZENE)
        contacts = find_pi_stacking(receptor, ligand, distance_cutoff=5.5)
        assert len(contacts) == 1
        assert contacts[0].resname == "PHE"
        assert contacts[0].angle == pytest.approx(0.0, abs=1.0)

    def test_ring_beyond_cutoff_is_not_detected(self):
        receptor = _receptor_from_rows([(*row, False) for row in self._phe_ring_rows(z_offset=20.0)])
        ligand = _mol_from_molblock(self._LIGAND_BENZENE)
        assert find_pi_stacking(receptor, ligand, distance_cutoff=5.5) == []

    def test_incomplete_ring_atoms_are_skipped(self):
        # Only 3 of the 6 named PHE ring atoms present -> not recognized as a ring.
        rows = self._phe_ring_rows(z_offset=4.0)[:3]
        receptor = _receptor_from_rows([(*row, False) for row in rows])
        ligand = _mol_from_molblock(self._LIGAND_BENZENE)
        assert find_pi_stacking(receptor, ligand, distance_cutoff=5.5) == []


class TestFindSulfurHalogenBonds:
    # A ligand C-Cl bond pointing along +x from the origin.
    _LIGAND_CHLOROMETHANE = """
  test

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"""

    def test_linear_approach_to_sulfur_is_detected(self):
        # Receptor sulfur placed further out along the same C->Cl direction:
        # a linear (favorable) sigma-hole geometry.
        receptor = _receptor_from_rows([
            ("A", 40, "CYS", "SG", "S", 5.0, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(self._LIGAND_CHLOROMETHANE)
        contacts = find_sulfur_halogen_bonds(receptor, ligand, distance_cutoff=4.5)
        assert len(contacts) == 1
        assert contacts[0].resname == "CYS"

    def test_perpendicular_approach_is_not_detected(self):
        # Sulfur off to the side (perpendicular to C->Cl) -> not linear enough.
        receptor = _receptor_from_rows([
            ("A", 40, "CYS", "SG", "S", 1.8, 4.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(self._LIGAND_CHLOROMETHANE)
        assert find_sulfur_halogen_bonds(receptor, ligand, distance_cutoff=4.5) == []

    def test_no_sulfur_atoms_in_receptor_returns_empty(self):
        receptor = _receptor_from_rows([
            ("A", 1, "ALA", "CB", "C", 5.0, 0.0, 0.0, False),
        ])
        ligand = _mol_from_molblock(self._LIGAND_CHLOROMETHANE)
        assert find_sulfur_halogen_bonds(receptor, ligand, distance_cutoff=4.5) == []
