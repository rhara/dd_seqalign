"""Summary table across all poses, for the score/RMSD/interaction-count dashboard."""
from typing import Optional, Sequence

import pandas as pd
from rdkit import Chem

from .interactions import (
    find_contact_residues,
    find_hydrogen_bonds,
    find_hydrophobic_contacts,
    find_salt_bridges,
    find_pi_stacking,
    find_electrostatic_interactions,
    find_pi_halogen_bonds,
    find_sulfur_halogen_bonds,
)
from .io import Pose, Receptor
from .scoring import detect_score, rmsd_to_reference


def poses_dataframe(
    poses: Sequence[Pose],
    receptor: Optional[Receptor] = None,
    reference_mol: Optional[Chem.Mol] = None,
) -> pd.DataFrame:
    """One row per pose: index, label, docking score, RMSD to `reference_mol`
    (if given), and interaction counts (contact residues, H-bonds,
    hydrophobic, salt bridges, pi-stacking, electrostatic, pi-halogen bonds,
    sulfur-halogen bonds) against `receptor` (if given). Either extra
    argument can be omitted to skip its columns -- e.g. call with just
    `poses` for a bare score table.
    """
    rows = []
    for pose in poses:
        row = {
            "index": pose.index,
            "label": pose.label,
            "score": detect_score(pose),
        }
        if reference_mol is not None:
            row["rmsd"] = rmsd_to_reference(pose.mol, reference_mol)
        if receptor is not None:
            row["n_contact_residues"] = len(find_contact_residues(receptor, pose.mol))
            row["n_hbonds"] = len(find_hydrogen_bonds(receptor, pose.mol))
            row["n_hydrophobic"] = len(find_hydrophobic_contacts(receptor, pose.mol))
            row["n_salt_bridges"] = len(find_salt_bridges(receptor, pose.mol))
            row["n_pi_stacking"] = len(find_pi_stacking(receptor, pose.mol))
            row["n_electrostatic"] = len(find_electrostatic_interactions(receptor, pose.mol))
            row["n_pi_halogen"] = len(find_pi_halogen_bonds(receptor, pose.mol))
            row["n_sulfur_halogen"] = len(find_sulfur_halogen_bonds(receptor, pose.mol))
        rows.append(row)
    return pd.DataFrame(rows)
