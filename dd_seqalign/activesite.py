"""Active-site residue detection, in two selectable modes, plus translating
a detected site into any other structure's own numbering via the canonical-
UniProt-position coordinate system `sequence.py` establishes.

- `site_from_ligand`: distance-based, around the structure's own bound
  ligand (reuses `pdbio`'s water/additive/cofactor/unknown classification
  to find the real ligand rather than a cryoprotectant). Only usable on
  structures that actually have a ligand.
- `site_from_pocket`: fpocket-based auto-detection (reuses `pocket`),
  usable on any structure including apo ones and the AlphaFold model.

Both return residues as (chain_id, author_resseq) pairs in the *input
structure's own* numbering -- `map_site_to_structure` is what makes them
comparable across structures with different numbering schemes/chain
compositions, by round-tripping through each structure's `ChainAlignment`
(structure resseq -> canonical UniProt position -> other structure's
resseq).
"""
from __future__ import annotations

import shutil
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from Bio.PDB import NeighborSearch, PDBParser

from .pdbio import classify_hetero_groups, collect_hetero_groups, group_coords, pick_ligand_of_interest, select_protein
from .pocket import find_druggable_pocket
from .sequence import ChainAlignment

SiteResidue = Tuple[str, int]  # (chain_id, author resseq)


def site_from_ligand(
    pdb_path: Union[str, Path], *, chain_id: Optional[str] = None, cutoff: float = 5.0, min_ligand_atoms: int = 5,
) -> List[SiteResidue]:
    """Protein residues with any atom within `cutoff` Angstrom of the
    structure's auto-picked ligand of interest (see
    `pdbio.pick_ligand_of_interest`). Returns `[]` if this
    structure has no plausible ligand (apo structures, the AlphaFold
    model) -- callers should fall back to `site_from_pocket` in that case.
    `chain_id`, if given, restricts the result to that chain (the target
    protein chain, since a bound partner chain, e.g. cyclin, can also have
    ligand-proximal residues that aren't part of the site of interest).
    """
    text = Path(pdb_path).read_text()
    groups = classify_hetero_groups(collect_hetero_groups(text))
    ligand = pick_ligand_of_interest(groups, min_atoms=min_ligand_atoms)
    if ligand is None:
        return []
    ligand_coords = group_coords(ligand.lines)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", str(pdb_path))
    model = next(iter(structure))
    protein_atoms = [
        atom
        for chain in model
        for res in chain
        for atom in res
        if res.id[0] == " " and (chain_id is None or chain.id == chain_id)
    ]
    ns = NeighborSearch(protein_atoms)

    seen = set()
    site: List[SiteResidue] = []
    for coord in ligand_coords:
        for atom in ns.search(coord, cutoff):
            res = atom.get_parent()
            key = (res.get_parent().id, res.id[1])
            if key not in seen:
                seen.add(key)
                site.append(key)
    return sorted(site)


def site_from_pocket(
    pdb_path: Union[str, Path], *, chain_id: str, work_dir: Optional[Union[str, Path]] = None, pocket_rank: int = 1,
) -> List[SiteResidue]:
    """Auto-detected druggable pocket (fpocket, via `pocket`) on the
    given chain in isolation -- the input is first stripped to that
    chain's protein atoms only (a temp file) so fpocket sees a single
    kinase domain rather than e.g. a CDK1/CyclinB/Cks2 assembly, which
    would let it detect an inter-chain groove instead of the intended
    (single-chain) active site. Works on apo structures and the AlphaFold
    model, unlike `site_from_ligand`.
    """
    text = Path(pdb_path).read_text()
    protein_lines = select_protein(text, chains=[chain_id])
    if not protein_lines:
        raise ValueError(f"{pdb_path}: no protein atoms found for chain {chain_id!r}")

    own_tmp = work_dir is None
    work_dir = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="dd_seqalign_pocket_"))
    work_dir.mkdir(parents=True, exist_ok=True)
    chain_pdb = work_dir / f"{Path(pdb_path).stem}_{chain_id}_protein.pdb"
    chain_pdb.write_text("\n".join(protein_lines) + "\nEND\n")

    try:
        selection = find_druggable_pocket(chain_pdb, work_dir, pocket_rank=pocket_rank, show_progress=False)
    finally:
        if own_tmp:
            shutil.rmtree(work_dir, ignore_errors=True)

    return sorted((r.chain, r.resnum) for r in selection.residues)


def map_site_to_structure(
    site: Sequence[SiteResidue], site_chain_alignment: ChainAlignment, target_chain_alignment: ChainAlignment,
) -> List[int]:
    """Translate a site detected on one structure (`site`, in that
    structure's own numbering, restricted to the chain
    `site_chain_alignment` describes) into `target_chain_alignment`'s
    structure's own residue numbers, by round-tripping through canonical
    UniProt positions. Site residues at a canonical position the target
    structure doesn't resolve (missing density there) are silently
    dropped -- the caller ends up with however much of the site actually
    overlaps what's modeled in the target."""
    canonical_positions = (site_chain_alignment.canonical_for_resseq(resseq) for _chain, resseq in site)
    resseqs = (
        target_chain_alignment.resseq_for_canonical(pos) for pos in canonical_positions if pos is not None
    )
    return sorted({r for r in resseqs if r is not None})
