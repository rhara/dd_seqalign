"""Loading receptors and ligand poses from files.

Receptor keeps both the raw PDB text (so py3Dmol/3Dmol.js can parse its own
cartoon/secondary-structure info) and a flat atom table (for the distance-based
interaction analysis in `interactions.py`, which needs plain coordinates and
residue metadata, not a JS-side structure).
"""
import logging
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select
from rdkit import Chem

logger = logging.getLogger(__name__)

ATOM_COLUMNS = ["chain", "resnum", "resname", "name", "element", "x", "y", "z", "is_hetero"]

# HETATM residue names that are essentially never the ligand of interest --
# crystallization water/buffers/cryoprotectants and common ions -- so
# `split_structure` skips them when looking for a bound ligand to extract.
COMMON_NON_LIGAND_HETATMS = {
    "HOH", "WAT", "H2O",
    "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "CD", "HG",
    "SO4", "PO4", "NO3", "IOD", "BR",
    "GOL", "EDO", "DMS", "ACT", "TRS", "PEG", "PG4", "P6G", "1PE", "MPD", "MRD", "IPA", "FMT", "CIT", "UNK",
}
MIN_LIGAND_HEAVY_ATOMS = 5


@dataclass
class Receptor:
    """A loaded protein receptor.

    Keeps both the raw PDB text (so py3Dmol/3Dmol.js can parse its own
    cartoon/secondary-structure info) and a flat, per-atom `pandas.DataFrame`
    (columns: `ATOM_COLUMNS`) used by the distance-based interaction analysis
    in `interactions.py`.
    """

    pdb_text: str
    atoms: pd.DataFrame
    source: str = ""


@dataclass
class Pose:
    """A single ligand pose (one docking result, or one embedded/reference
    ligand) plus its metadata.

    `index` is the pose's position among its siblings (e.g. its record index
    within a multi-model SDF); `label` is a short human-readable identifier
    for UI display; `properties` holds any SDF tags read alongside the
    molecule (e.g. docking-tool score fields, see `scoring.py`).
    """

    mol: Chem.Mol
    index: int
    label: str
    properties: dict = field(default_factory=dict)


def _guess_element(atom_name: str) -> str:
    """PDB files without columns 77-78 populated leave Bio.PDB's `element`
    empty; fall back to stripping digits/whitespace from the atom name
    (e.g. "1HD1" -> "H", "CA" -> "C").
    """
    stripped = "".join(c for c in atom_name if c.isalpha())
    return (stripped[:1] or "X").upper()


def load_receptor(path: str) -> Receptor:
    path = str(path)
    pdb_text = Path(path).read_text()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", path)

    rows = []
    for model in structure:
        for chain in model:
            for residue in chain:
                is_hetero = residue.id[0] != " "
                for atom in residue:
                    x, y, z = atom.coord
                    element = (atom.element or "").strip() or _guess_element(atom.get_name())
                    rows.append((
                        chain.id, residue.id[1], residue.resname, atom.get_name(),
                        element, float(x), float(y), float(z), is_hetero,
                    ))
        break  # only the first model; receptors here are single-model PDBs

    atoms = pd.DataFrame(rows, columns=ATOM_COLUMNS)
    return Receptor(pdb_text=pdb_text, atoms=atoms, source=path)


class _SelectResidue(Select):
    def __init__(self, target):
        self.target = target

    def accept_residue(self, residue):
        return residue is self.target


class _SelectProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "


def _residue_to_mol(structure, residue) -> Optional[Chem.Mol]:
    """Convert a single Bio.PDB HETATM residue to an RDKit Mol.

    Bond connectivity is perceived from interatomic distances (there's no
    CONECT-independent chemical info in a bare PDB residue), which is less
    reliable than a proper SDF -- in particular it does not reliably
    recover aromaticity/bond order, so pi-stacking and halogen-bond
    detection can under-report for a Pose built this way. Hydrogen bonds,
    hydrophobic contacts, and formal-charge-based checks are unaffected
    (they only need element/connectivity, not bond order).
    """
    io_ = PDBIO()
    io_.set_structure(structure)
    buf = StringIO()
    io_.save(buf, _SelectResidue(residue))
    try:
        mol = Chem.MolFromPDBBlock(buf.getvalue(), sanitize=True, removeHs=False)
    except (Chem.rdchem.MolSanitizeException, ValueError, RuntimeError) as exc:
        logger.warning(
            "Failed to convert residue %s%s to an RDKit Mol: %s",
            residue.get_resname(), residue.id[1], exc,
        )
        return None
    return mol


def split_structure(path: str) -> tuple[Receptor, list[Pose]]:
    """Split a single co-crystal PDB (protein + bound ligand + waters all in
    one file, e.g. downloaded straight from the PDB) into a protein-only
    `Receptor` and any bound ligand(s) as `Pose`s -- so a raw PDB works
    without a separately-prepared receptor PDB / ligand SDF.

    A HETATM residue is treated as a candidate ligand unless its name is a
    common water/ion/crystallization-buffer code (`COMMON_NON_LIGAND_HETATMS`)
    or it has fewer than `MIN_LIGAND_HEAVY_ATOMS` heavy atoms. A covalently
    bound ligand is still extracted as its own residue-based Pose (its 3D
    pose and internal geometry are captured; the covalent bond itself isn't
    modeled as a bond in the returned Mol, since PDB HETATM records don't
    carry cross-residue connectivity).
    """
    path = str(path)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", path)
    model = next(iter(structure))

    protein_rows = []
    ligand_residues = []
    for chain in model:
        for residue in chain:
            resname = residue.resname.strip()
            if residue.id[0] == " ":
                for atom in residue:
                    x, y, z = atom.coord
                    element = (atom.element or "").strip() or _guess_element(atom.get_name())
                    protein_rows.append((
                        chain.id, residue.id[1], resname, atom.get_name(),
                        element, float(x), float(y), float(z), False,
                    ))
                continue
            if resname in COMMON_NON_LIGAND_HETATMS:
                continue
            heavy_atoms = sum(1 for a in residue if ((a.element or "").strip() or _guess_element(a.get_name())).upper() != "H")
            if heavy_atoms < MIN_LIGAND_HEAVY_ATOMS:
                continue
            ligand_residues.append((chain, residue))

    io_ = PDBIO()
    io_.set_structure(structure)
    buf = StringIO()
    io_.save(buf, _SelectProteinOnly())
    receptor = Receptor(pdb_text=buf.getvalue(), atoms=pd.DataFrame(protein_rows, columns=ATOM_COLUMNS), source=path)

    poses = []
    for chain, residue in ligand_residues:
        mol = _residue_to_mol(structure, residue)
        if mol is None:
            continue
        label = f"{residue.resname.strip()} {chain.id}{residue.id[1]} (auto-extracted from PDB)"
        poses.append(Pose(mol=mol, index=len(poses), label=label, properties={}))

    return receptor, poses


def load_poses(path: str) -> list[Pose]:
    """Read every record of a multi-model SDF as one docking pose each,
    in file order (the order most docking tools already rank poses by score).
    """
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    poses = []
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name") else f"pose {i}"
        poses.append(Pose(mol=mol, index=i, label=name, properties=mol.GetPropsAsDict()))
    return poses


def load_reference_ligand(path: str) -> Optional[Chem.Mol]:
    """Load a single reference/crystal ligand pose (first record of the file)."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    for mol in supplier:
        if mol is not None:
            return mol
    return None
