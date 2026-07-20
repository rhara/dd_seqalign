"""Distance-based binding-site contact detection.

No PLIP/ProDy dependency (full interaction typing needs protonation
states/geometry beyond what a receptor.pdb + docked SDF reliably give you).
Instead this uses the same kind of simplified heuristics most quick-look
docking viewers use:

- "contact residue": any residue with an atom within `cutoff` of the ligand.
- "hydrogen bond": an N/O ligand atom within `distance_cutoff` of an N/O
  receptor atom (donor/acceptor role and angle are not resolved  --  this
  flags plausible polar contacts, not confirmed H-bonds).
- "hydrophobic contact": a ligand carbon within `cutoff` of a receptor carbon,
  kept to the nearest receptor atom per ligand atom by default (there are
  usually many more raw carbon-carbon pairs within range than there are
  meaningfully distinct contacts).
- "salt bridge": a formally-charged ligand atom (RDKit `GetFormalCharge`)
  within `distance_cutoff` of an oppositely-charged receptor side-chain atom
  (ARG/LYS as positive, ASP/GLU as negative; HIS is skipped since its charge
  depends on an unknown protonation state).
- "pi-stacking": an aromatic ring on the ligand (RDKit ring perception) with
  its centroid within `distance_cutoff` of an aromatic receptor side-chain
  ring (PHE/TYR/TRP/HIS) centroid, and the two ring-plane normals roughly
  parallel (face-to-face) or roughly perpendicular (edge-to-face/T-shaped);
  in-between tilt angles are treated as ambiguous and dropped.
- "electrostatic interaction": broader than a salt bridge -- a ligand atom
  with a significant RDKit Gasteiger partial charge near a receptor atom
  with an approximate partial charge of the opposite sign (backbone
  carbonyl O / amide N, plus a small lookup table of polar/charged
  side-chain atoms; see `RECEPTOR_APPROX_CHARGES`). Unlike `find_salt_bridges`
  this doesn't require either side to be formally/fully charged, so it also
  catches e.g. a ligand carbonyl drawn charge-neutral near a lysine.
- "halogen bond": a ligand halogen (Cl/Br/I; F is excluded -- it's a poor
  halogen-bond donor) within `distance_cutoff` of a receptor aromatic ring
  centroid (`find_pi_halogen_bonds`) or sulfur atom (CYS SG / MET SD,
  `find_sulfur_halogen_bonds`), with the C-X...acceptor angle required to be
  roughly linear (<= `HALOGEN_BOND_MAX_LINEARITY_ANGLE`) since halogen bonds
  are directional (they come from the sigma-hole opposite the C-X bond).
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from .io import Receptor

POLAR_ELEMENTS = {"N", "O"}
HYDROPHOBIC_ELEMENTS = {"C"}

RECEPTOR_POSITIVE_ATOMS = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("LYS", "NZ")}
RECEPTOR_NEGATIVE_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}

AROMATIC_RESIDUE_RING_ATOMS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),  # six-membered ring only
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}

PI_STACKING_PARALLEL_MAX_ANGLE = 35.0
PI_STACKING_PERPENDICULAR_MIN_ANGLE = 55.0

# Approximate partial charges for standard-residue atoms, used by
# find_electrostatic_interactions. Not a real force field -- just enough to
# rank "clearly positive" / "clearly negative" / "not really charged" for a
# quick-look view. Backbone amide N and carbonyl O are handled separately
# below since their names ("N"/"O") are shared by every residue type.
RECEPTOR_APPROX_CHARGES = {
    ("ARG", "NH1"): 0.5, ("ARG", "NH2"): 0.5, ("ARG", "NE"): 0.3, ("LYS", "NZ"): 0.6,
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5, ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("SER", "OG"): -0.4, ("THR", "OG1"): -0.4, ("TYR", "OH"): -0.35,
    ("ASN", "OD1"): -0.4, ("ASN", "ND2"): 0.3, ("GLN", "OE1"): -0.4, ("GLN", "NE2"): 0.3,
    ("HIS", "ND1"): -0.2, ("HIS", "NE2"): -0.2,
    ("CYS", "SG"): -0.2, ("TRP", "NE1"): 0.2,
}
RECEPTOR_BACKBONE_CARBONYL_CHARGE = -0.45
RECEPTOR_BACKBONE_AMIDE_CHARGE = 0.3
ELECTROSTATIC_CHARGE_THRESHOLD = 0.2

HALOGEN_ELEMENTS = {"Cl", "Br", "I"}
RECEPTOR_SULFUR_ATOMS = {("CYS", "SG"), ("MET", "SD")}
HALOGEN_BOND_MAX_LINEARITY_ANGLE = 60.0


@dataclass
class Contact:
    """A single ligand-atom-to-receptor-atom contact (H-bond, hydrophobic,
    salt bridge, electrostatic, or halogen-bond candidate, depending on which
    `find_*` function produced it).

    `ligand_atom_idx` is the RDKit atom index into the ligand `Chem.Mol`;
    `chain`/`resnum`/`resname`/`atom_name` identify the receptor atom;
    `distance` is the Euclidean distance between the two in Angstroms;
    `rec_x`/`rec_y`/`rec_z` are the receptor atom's coordinates (handy for
    drawing the contact without re-looking-up the atom).
    """

    ligand_atom_idx: int
    chain: str
    resnum: int
    resname: str
    atom_name: str
    distance: float
    rec_x: float
    rec_y: float
    rec_z: float


@dataclass
class RingContact:
    """A single ligand-ring-to-receptor-ring pi-stacking contact (as
    produced by `find_pi_stacking`).

    `ligand_ring_atoms` are the RDKit atom indices making up the ligand's
    aromatic ring; `chain`/`resnum`/`resname` identify the receptor's
    aromatic side chain; `distance` is the centroid-to-centroid distance in
    Angstroms; `angle` is the angle between the two ring-plane normals in
    degrees (near 0 = face-to-face, near 90 = edge-to-face/T-shaped);
    `lig_x/y/z` and `rec_x/y/z` are the two ring centroids' coordinates.
    """

    ligand_ring_atoms: list
    chain: str
    resnum: int
    resname: str
    distance: float
    angle: float
    lig_x: float
    lig_y: float
    lig_z: float
    rec_x: float
    rec_y: float
    rec_z: float


def _ligand_atom_coords(pose_mol: Chem.Mol, conf_id: int = 0):
    conf = pose_mol.GetConformer(conf_id)
    elements = np.array([atom.GetSymbol() for atom in pose_mol.GetAtoms()])
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(pose_mol.GetNumAtoms())])
    return elements, coords


def _receptor_coords(receptor: Receptor) -> np.ndarray:
    return receptor.atoms[["x", "y", "z"]].to_numpy()


def find_contact_residues(receptor: Receptor, pose_mol: Chem.Mol, cutoff: float = 4.5, conf_id: int = 0) -> pd.DataFrame:
    """Residues with at least one atom within `cutoff` of any ligand atom,
    sorted by closest approach.
    """
    _, ligand_coords = _ligand_atom_coords(pose_mol, conf_id)
    receptor_coords = _receptor_coords(receptor)

    dists = np.linalg.norm(receptor_coords[:, None, :] - ligand_coords[None, :, :], axis=-1)
    min_dist = dists.min(axis=1)

    atoms = receptor.atoms.copy()
    atoms["min_dist"] = min_dist
    contacts = atoms[atoms["min_dist"] <= cutoff]
    residues = (
        contacts.groupby(["chain", "resnum", "resname"], as_index=False)["min_dist"]
        .min()
        .sort_values("min_dist")
        .reset_index(drop=True)
    )
    return residues


def _dedupe_nearest_per_ligand_atom(contacts: list) -> list:
    best = {}
    for contact in contacts:
        current = best.get(contact.ligand_atom_idx)
        if current is None or contact.distance < current.distance:
            best[contact.ligand_atom_idx] = contact
    return list(best.values())


def _receptor_atom_subset(receptor: Receptor, rec_mask: np.ndarray):
    """Return (coords, rows) for the receptor atoms selected by `rec_mask`.

    Shared by `_atom_contacts` and `find_sulfur_halogen_bonds`: both need the
    same "pull out the masked receptor atoms' coordinates and metadata rows"
    step before computing distances to ligand atoms/halogens.
    """
    atoms = receptor.atoms
    rec_coords = atoms.loc[rec_mask, ["x", "y", "z"]].to_numpy()
    rec_rows = atoms.loc[rec_mask].reset_index(drop=True)
    return rec_coords, rec_rows


def _build_receptor_atom_contacts(rec_rows: pd.DataFrame, ligand_atom_idx: int, distance: float, rec_row_idx: int) -> Contact:
    """Build a single `Contact` from a receptor atom row and a precomputed distance.

    Shared by `_atom_contacts` and `find_sulfur_halogen_bonds` so the
    `Contact(...)` field mapping (chain/resnum/resname/atom_name/rec_x/y/z)
    only needs to be kept in sync in one place.
    """
    row = rec_rows.iloc[rec_row_idx]
    return Contact(
        ligand_atom_idx=ligand_atom_idx,
        chain=row["chain"], resnum=int(row["resnum"]), resname=row["resname"],
        atom_name=row["name"], distance=float(distance),
        rec_x=float(row["x"]), rec_y=float(row["y"]), rec_z=float(row["z"]),
    )


def _atom_contacts(receptor: Receptor, rec_mask: np.ndarray, lig_coords: np.ndarray, lig_idx: np.ndarray, distance_cutoff: float) -> list:
    """All receptor atoms selected by `rec_mask` within `distance_cutoff` of
    any ligand atom in `lig_coords[lig_idx]`, as `Contact` objects (one per
    close receptor/ligand atom pair; a ligand atom close to several receptor
    atoms yields several `Contact`s, and vice versa).
    """
    if not rec_mask.any() or len(lig_idx) == 0:
        return []
    rec_coords, rec_rows = _receptor_atom_subset(receptor, rec_mask)

    dists = np.linalg.norm(rec_coords[:, None, :] - lig_coords[lig_idx][None, :, :], axis=-1)
    contacts = []
    close_i, close_j = np.nonzero(dists <= distance_cutoff)
    for i, j in zip(close_i, close_j):
        contacts.append(_build_receptor_atom_contacts(rec_rows, int(lig_idx[j]), dists[i, j], i))
    return contacts


def _polar_contacts(receptor: Receptor, pose_mol: Chem.Mol, distance_cutoff: float, elements: set, conf_id: int) -> list:
    lig_elements, lig_coords = _ligand_atom_coords(pose_mol, conf_id)
    lig_idx = np.nonzero(np.isin(lig_elements, list(elements)))[0]
    rec_mask = receptor.atoms["element"].isin(elements).to_numpy()
    return _atom_contacts(receptor, rec_mask, lig_coords, lig_idx, distance_cutoff)


def find_hydrogen_bonds(receptor: Receptor, pose_mol: Chem.Mol, distance_cutoff: float = 3.5, conf_id: int = 0) -> list:
    """Ligand N/O atoms within `distance_cutoff` of a receptor N/O atom.

    Donor/acceptor role and angle are not resolved, so this flags plausible
    polar contacts rather than confirmed hydrogen bonds (see module
    docstring).
    """
    return _polar_contacts(receptor, pose_mol, distance_cutoff, POLAR_ELEMENTS, conf_id)


def find_hydrophobic_contacts(receptor: Receptor, pose_mol: Chem.Mol, cutoff: float = 4.5, conf_id: int = 0, nearest_only: bool = True) -> list:
    """Ligand carbons within `cutoff` of a receptor carbon.

    When `nearest_only` is True (the default), only the closest receptor
    contact per ligand atom is kept (`_dedupe_nearest_per_ligand_atom`) --
    there are usually many more raw carbon-carbon pairs within range than
    there are meaningfully distinct contacts.
    """
    contacts = _polar_contacts(receptor, pose_mol, cutoff, HYDROPHOBIC_ELEMENTS, conf_id)
    return _dedupe_nearest_per_ligand_atom(contacts) if nearest_only else contacts


def find_salt_bridges(receptor: Receptor, pose_mol: Chem.Mol, distance_cutoff: float = 4.0, conf_id: int = 0) -> list:
    """Formally-charged ligand atoms within `distance_cutoff` of an
    oppositely-charged receptor side-chain atom. Ligand charge comes from
    RDKit's `GetFormalCharge` on the input structure as given -- a ligand
    drawn charge-neutral (e.g. a carboxylic acid instead of a carboxylate)
    won't be flagged even if it would be ionized at physiological pH.
    """
    _, lig_coords = _ligand_atom_coords(pose_mol, conf_id)
    lig_charges = np.array([atom.GetFormalCharge() for atom in pose_mol.GetAtoms()])
    if not (lig_charges != 0).any():
        return []

    atoms = receptor.atoms
    rec_keys = list(zip(atoms["resname"], atoms["name"]))
    rec_pos_mask = np.array([k in RECEPTOR_POSITIVE_ATOMS for k in rec_keys])
    rec_neg_mask = np.array([k in RECEPTOR_NEGATIVE_ATOMS for k in rec_keys])

    lig_pos_idx = np.nonzero(lig_charges > 0)[0]
    lig_neg_idx = np.nonzero(lig_charges < 0)[0]

    contacts = _atom_contacts(receptor, rec_neg_mask, lig_coords, lig_pos_idx, distance_cutoff)
    contacts += _atom_contacts(receptor, rec_pos_mask, lig_coords, lig_neg_idx, distance_cutoff)
    return contacts


def _plane_normal(coords: np.ndarray) -> np.ndarray:
    normal = np.cross(coords[1] - coords[0], coords[2] - coords[0])
    norm = np.linalg.norm(normal)
    return normal / norm if norm > 0 else normal


def _receptor_aromatic_rings(receptor: Receptor):
    atoms = receptor.atoms
    rings = []
    for (chain, resnum, resname), group in atoms.groupby(["chain", "resnum", "resname"]):
        ring_atom_names = AROMATIC_RESIDUE_RING_ATOMS.get(resname)
        if ring_atom_names is None:
            continue
        sub = group[group["name"].isin(ring_atom_names)]
        if len(sub) < len(ring_atom_names):
            continue
        coords = sub[["x", "y", "z"]].to_numpy()
        rings.append((chain, int(resnum), resname, coords.mean(axis=0), _plane_normal(coords)))
    return rings


def _ligand_aromatic_rings(pose_mol: Chem.Mol, conf_id: int):
    conf = pose_mol.GetConformer(conf_id)
    rings = []
    for atom_ring in pose_mol.GetRingInfo().AtomRings():
        if not all(pose_mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ring):
            continue
        coords = np.array([list(conf.GetAtomPosition(i)) for i in atom_ring])
        rings.append((list(atom_ring), coords.mean(axis=0), _plane_normal(coords)))
    return rings


def find_pi_stacking(receptor: Receptor, pose_mol: Chem.Mol, distance_cutoff: float = 5.5, conf_id: int = 0) -> list:
    """Aromatic ligand rings with a centroid within `distance_cutoff` of an
    aromatic receptor side-chain ring (PHE/TYR/TRP/HIS) centroid, restricted
    to roughly parallel (face-to-face, angle <= `PI_STACKING_PARALLEL_MAX_ANGLE`)
    or roughly perpendicular (edge-to-face/T-shaped, angle >=
    `PI_STACKING_PERPENDICULAR_MIN_ANGLE`) ring-plane orientations;
    in-between tilt angles are ambiguous and dropped.
    """
    rec_rings = _receptor_aromatic_rings(receptor)
    if not rec_rings:
        return []
    lig_rings = _ligand_aromatic_rings(pose_mol, conf_id)

    contacts = []
    for lig_atoms, lig_centroid, lig_normal in lig_rings:
        for chain, resnum, resname, rec_centroid, rec_normal in rec_rings:
            distance = float(np.linalg.norm(lig_centroid - rec_centroid))
            if distance > distance_cutoff:
                continue
            cos_angle = min(1.0, max(-1.0, float(abs(np.dot(lig_normal, rec_normal)))))
            angle = float(np.degrees(np.arccos(cos_angle)))
            if PI_STACKING_PARALLEL_MAX_ANGLE < angle < PI_STACKING_PERPENDICULAR_MIN_ANGLE:
                continue  # ambiguous tilt -- neither face-to-face nor edge-to-face
            contacts.append(RingContact(
                ligand_ring_atoms=lig_atoms, chain=chain, resnum=resnum, resname=resname,
                distance=distance, angle=angle,
                lig_x=float(lig_centroid[0]), lig_y=float(lig_centroid[1]), lig_z=float(lig_centroid[2]),
                rec_x=float(rec_centroid[0]), rec_y=float(rec_centroid[1]), rec_z=float(rec_centroid[2]),
            ))
    return contacts


def _receptor_approx_charges(atoms: pd.DataFrame) -> np.ndarray:
    charges = np.zeros(len(atoms))
    for i, (resname, name) in enumerate(zip(atoms["resname"], atoms["name"])):
        if name == "O":
            charges[i] = RECEPTOR_BACKBONE_CARBONYL_CHARGE
        elif name == "N":
            charges[i] = RECEPTOR_BACKBONE_AMIDE_CHARGE
        else:
            charges[i] = RECEPTOR_APPROX_CHARGES.get((resname, name), 0.0)
    return charges


def find_electrostatic_interactions(
    receptor: Receptor, pose_mol: Chem.Mol, distance_cutoff: float = 5.0,
    charge_threshold: float = ELECTROSTATIC_CHARGE_THRESHOLD, conf_id: int = 0,
) -> list:
    """Ligand atoms with a Gasteiger partial charge of magnitude >=
    `charge_threshold` within `distance_cutoff` of a receptor atom with an
    approximate partial charge of the opposite sign and magnitude >=
    `charge_threshold` (see `RECEPTOR_APPROX_CHARGES`). Broader than
    `find_salt_bridges`: catches polar-but-not-formally-charged contacts
    (e.g. a neutral carbonyl near a lysine) that a strict formal-charge
    check misses.
    """
    _, lig_coords = _ligand_atom_coords(pose_mol, conf_id)
    charge_mol = Chem.Mol(pose_mol)
    AllChem.ComputeGasteigerCharges(charge_mol)
    lig_charges = np.array([
        charge_mol.GetAtomWithIdx(i).GetDoubleProp("_GasteigerCharge") for i in range(charge_mol.GetNumAtoms())
    ])
    lig_charges = np.nan_to_num(lig_charges, nan=0.0, posinf=0.0, neginf=0.0)

    rec_charges = _receptor_approx_charges(receptor.atoms)
    rec_pos_mask = rec_charges >= charge_threshold
    rec_neg_mask = rec_charges <= -charge_threshold

    lig_pos_idx = np.nonzero(lig_charges >= charge_threshold)[0]
    lig_neg_idx = np.nonzero(lig_charges <= -charge_threshold)[0]

    contacts = _atom_contacts(receptor, rec_neg_mask, lig_coords, lig_pos_idx, distance_cutoff)
    contacts += _atom_contacts(receptor, rec_pos_mask, lig_coords, lig_neg_idx, distance_cutoff)
    return contacts


def _ligand_halogen_atoms(pose_mol: Chem.Mol, conf_id: int):
    conf = pose_mol.GetConformer(conf_id)
    halogens = []
    for atom in pose_mol.GetAtoms():
        if atom.GetSymbol() not in HALOGEN_ELEMENTS:
            continue
        neighbors = atom.GetNeighbors()
        if not neighbors:
            continue
        x_pos = np.array(conf.GetAtomPosition(atom.GetIdx()))
        c_pos = np.array(conf.GetAtomPosition(neighbors[0].GetIdx()))
        halogens.append((atom.GetIdx(), x_pos, x_pos - c_pos))
    return halogens


def _linearity_angle(cx_vector: np.ndarray, x_to_acceptor: np.ndarray) -> float:
    """Deviation from ideal sigma-hole geometry, in degrees.

    For a linear C-X...Acceptor arrangement (the acceptor approaching along
    the extension of the C-X bond, where the sigma hole sits), the C->X
    vector and X->Acceptor vector point the same way -- 0 degrees apart, not
    180. A right-angle approach (acceptor beside X rather than in front of
    it) gives ~90 degrees; smaller is more linear/favorable.
    """
    denom = np.linalg.norm(cx_vector) * np.linalg.norm(x_to_acceptor)
    if denom == 0:
        return 180.0
    cos_angle = min(1.0, max(-1.0, float(np.dot(cx_vector, x_to_acceptor) / denom)))
    return float(np.degrees(np.arccos(cos_angle)))


def find_pi_halogen_bonds(receptor: Receptor, pose_mol: Chem.Mol, distance_cutoff: float = 4.5, conf_id: int = 0) -> list:
    """Ligand halogens (Cl/Br/I) within `distance_cutoff` of a receptor
    aromatic ring (PHE/TYR/TRP/HIS) centroid, with the C-X...ring-centroid
    angle required to be roughly linear (`_linearity_angle` <=
    `HALOGEN_BOND_MAX_LINEARITY_ANGLE`) since halogen bonds are directional
    (they come from the sigma-hole opposite the C-X bond).
    """
    halogens = _ligand_halogen_atoms(pose_mol, conf_id)
    if not halogens:
        return []
    rec_rings = _receptor_aromatic_rings(receptor)

    contacts = []
    for lig_idx, x_pos, cx_vec in halogens:
        for chain, resnum, resname, rec_centroid, _rec_normal in rec_rings:
            distance = float(np.linalg.norm(x_pos - rec_centroid))
            if distance > distance_cutoff:
                continue
            if _linearity_angle(cx_vec, rec_centroid - x_pos) > HALOGEN_BOND_MAX_LINEARITY_ANGLE:
                continue
            contacts.append(Contact(
                ligand_atom_idx=lig_idx, chain=chain, resnum=resnum, resname=resname,
                atom_name="ring_centroid", distance=distance,
                rec_x=float(rec_centroid[0]), rec_y=float(rec_centroid[1]), rec_z=float(rec_centroid[2]),
            ))
    return contacts


def find_sulfur_halogen_bonds(receptor: Receptor, pose_mol: Chem.Mol, distance_cutoff: float = 4.5, conf_id: int = 0) -> list:
    """Ligand halogens (Cl/Br/I) within `distance_cutoff` of a receptor
    sulfur atom (CYS SG / MET SD), with the C-X...S angle required to be
    roughly linear (`_linearity_angle` <= `HALOGEN_BOND_MAX_LINEARITY_ANGLE`).
    """
    halogens = _ligand_halogen_atoms(pose_mol, conf_id)
    if not halogens:
        return []
    atoms = receptor.atoms
    rec_keys = list(zip(atoms["resname"], atoms["name"]))
    rec_mask = np.array([k in RECEPTOR_SULFUR_ATOMS for k in rec_keys])
    if not rec_mask.any():
        return []
    rec_coords, rec_rows = _receptor_atom_subset(receptor, rec_mask)

    contacts = []
    for lig_idx, x_pos, cx_vec in halogens:
        dists = np.linalg.norm(rec_coords - x_pos, axis=1)
        for i in np.nonzero(dists <= distance_cutoff)[0]:
            if _linearity_angle(cx_vec, rec_coords[i] - x_pos) > HALOGEN_BOND_MAX_LINEARITY_ANGLE:
                continue
            contacts.append(_build_receptor_atom_contacts(rec_rows, lig_idx, dists[i], i))
    return contacts
