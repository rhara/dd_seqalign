"""Score extraction and pose-vs-reference RMSD."""
from typing import Optional

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from .io import Pose

# Common SDF property names used by popular docking tools for their pose
# score, checked in priority order (lower-is-better docking scores first).
SCORE_PROPERTY_CANDIDATES = [
    "minimizedAffinity",  # smina
    "REMARK.VINARESULT",
    "vina_score",
    "VINA_RESULT",
    "Score",
    "score",
    "affinity",
    "r_i_docking_score",  # Glide
    "docking_score",
]


def detect_score(pose: Pose) -> Optional[float]:
    for key in SCORE_PROPERTY_CANDIDATES:
        if key in pose.properties:
            try:
                return float(pose.properties[key])
            except (TypeError, ValueError):
                continue
    return None


def rmsd_to_reference(pose_mol: Chem.Mol, reference_mol: Chem.Mol) -> Optional[float]:
    """Heavy-atom RMSD against `reference_mol`, without any re-superposition
    (both poses are already docked into the same receptor frame — aligning
    them first would hide genuine pose displacement). Uses RDKit's atom
    mapping with symmetry perception, since docking tools can permute atom
    order and swap chemically-equivalent atoms (e.g. a symmetric phenyl ring).
    """
    try:
        probe = Chem.RemoveHs(pose_mol)
        ref = Chem.RemoveHs(reference_mol)
        return rdMolAlign.CalcRMS(probe, ref)
    except (RuntimeError, ValueError):
        return None
