"""Minimal fixed-column PDB parsing and HETATM classification, vendored (not
imported) from `dd_prep.parse`/`dd_prep.hetero` -- dd_seqalign only needs a
handful of read-only helpers from those modules, and duplicating them here
(instead of depending on the whole dd_prep package) keeps dd_seqalign
installable on its own. No external dependency, plain text, same convention
as dd_docking/receptor_prep.py.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

Coord = Tuple[float, float, float]
GroupKey = Tuple[str, str, int]  # (resname, chain, resseq)

WATERS = {"HOH", "WAT", "DOD", "H2O"}

CRYO_ADDITIVES = {
    # polyols / cryoprotectants
    "GOL", "EDO", "PEG", "PG4", "1PE", "PGE", "P6G", "PE4", "PE8", "MPD", "DIO", "BU3", "12P",
    # sulfoxide / alcohols
    "DMS", "IPA", "MOH",
    # buffer / salt ions and small anions-cations from crystallization liquor
    "SO4", "PO4", "ACT", "CIT", "FMT", "TRS", "MES", "EPE", "CAC", "BCT", "NO3", "OXL",
    # reducing agents
    "BME", "DTT", "TCEP",
    # misc additives
    "IOD", "BR", "AZI", "IMD", "GSH",
}

COFACTOR_METALS = {"ZN", "MG", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "CD", "CA", "NA", "K"}

COFACTORS_ORGANIC = {
    "HEM", "HEC", "NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "FMNH", "PLP",
    "SAM", "SAH", "COA", "ATP", "ADP", "AMP", "GTP", "GDP", "GNP", "TPP", "BIO",
}

COFACTORS = COFACTOR_METALS | COFACTORS_ORGANIC


def altloc_ok(line: str) -> bool:
    a = line[16] if len(line) > 16 else " "
    return a in (" ", "A")


def select_protein(text: str, chains: Optional[Sequence[str]] = None) -> List[str]:
    """ATOM (+ TER) lines for the given chains (None = all chains present),
    keeping only the primary altloc."""
    keep = set(chains) if chains is not None else None
    out: List[str] = []
    for ln in text.splitlines():
        rec = ln[:6]
        if rec == "ATOM  " and altloc_ok(ln) and (keep is None or ln[21] in keep):
            out.append(ln)
        elif rec == "TER   " and (keep is None or ln[21:22] in keep):
            out.append(ln)
    return out


def collect_hetero_groups(text: str, chains: Optional[Sequence[str]] = None) -> Dict[GroupKey, List[str]]:
    """All HETATM lines (primary altloc only), grouped by (resname, chain,
    resseq). `chains=None` collects hetero groups from every chain."""
    keep = set(chains) if chains is not None else None
    groups: Dict[GroupKey, List[str]] = {}
    for ln in text.splitlines():
        if ln[:6] != "HETATM" or not altloc_ok(ln):
            continue
        chain = ln[21]
        if keep is not None and chain not in keep:
            continue
        try:
            resseq = int(ln[22:26])
        except ValueError:
            continue
        resname = ln[17:20].strip().upper()
        groups.setdefault((resname, chain, resseq), []).append(ln)
    return groups


def group_coords(lines: Sequence[str]) -> List[Coord]:
    return [(float(ln[30:38]), float(ln[38:46]), float(ln[46:54])) for ln in lines]


def classify_group(resname: str) -> str:
    """One of 'water', 'additive', 'cofactor', 'unknown'."""
    r = resname.strip().upper()
    if r in WATERS:
        return "water"
    if r in CRYO_ADDITIVES:
        return "additive"
    if r in COFACTORS:
        return "cofactor"
    return "unknown"


@dataclass
class HeteroGroup:
    resname: str
    chain: str
    resseq: int
    lines: List[str]
    category: str

    @property
    def n_atoms(self) -> int:
        return len(self.lines)

    @property
    def label(self) -> str:
        return f"{self.resname}_{self.chain}{self.resseq}"


def classify_hetero_groups(groups: Dict[GroupKey, List[str]]) -> List[HeteroGroup]:
    out = []
    for (resname, chain, resseq), lines in groups.items():
        out.append(HeteroGroup(resname, chain, resseq, lines, classify_group(resname)))
    return sorted(out, key=lambda g: (g.chain, g.resseq))


def pick_ligand_of_interest(
    hetero_groups: Sequence[HeteroGroup], resname: str = "", min_atoms: int = 5,
) -> Optional[HeteroGroup]:
    """Pick the group to treat as "the ligand". If `resname` is given, the
    first group matching it exactly (any category). Otherwise, the largest
    'unknown'-category group with at least `min_atoms` heavy atoms (ties
    broken by atom count, then first-seen) -- water/additives/cofactors are
    never picked implicitly."""
    if resname:
        for g in hetero_groups:
            if g.resname == resname.upper():
                return g
        return None
    candidates = [g for g in hetero_groups if g.category == "unknown" and g.n_atoms >= min_atoms]
    if not candidates:
        return None
    return max(candidates, key=lambda g: g.n_atoms)
