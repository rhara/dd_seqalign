"""Structural superposition of every structure onto one reference
structure, in one shared PyMOL session (`pymol2`, not the classic global
`pymol`/`cmd` singleton -- keeps this importable as a library without
launching a GUI or fighting other code over global state).

Two site modes:
- `"pocket"` / `"ligand"`: fit on an explicit, already-known residue
  correspondence (`cmd.pair_fit`) -- both structures' active-site residues
  were already mapped onto the same canonical UniProt positions by
  `activesite.map_site_to_structure`, so there is a known 1:1 pairing and
  no need for (or risk from) an automatic re-matching step.
- `"none"`: no residue correspondence assumed; PyMOL's CE algorithm
  (`cealign`) finds the best structural superposition on its own, which is
  what makes it usable across oligomeric state/fragment/numbering
  differences, at the cost of being a whole-chain (not site-focused) fit.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

SITE_MODES = ("pocket", "ligand", "none")


@dataclass
class StructureInput:
    label: str
    pdb_path: str
    chain_id: str
    site_resseqs: Optional[List[int]] = None  # this structure's own numbering; None/[] => whole-chain cealign


@dataclass
class AlignmentResult:
    label: str
    reference_label: str
    site_mode: str
    rmsd: Optional[float]
    n_atoms: int  # atom-pair count (site mode) or aligned-residue count (cealign)
    aligned_pdb: str
    error: Optional[str] = None  # set instead of raising when this one structure can't be fit (see align_structures)


def _pair_fit(cmd, mobile: StructureInput, ref: StructureInput) -> AlignmentResult:
    if len(mobile.site_resseqs) != len(ref.site_resseqs):
        raise ValueError(
            f"{mobile.label}: site residue count ({len(mobile.site_resseqs)}) doesn't match "
            f"reference {ref.label}'s ({len(ref.site_resseqs)}) -- both should come from the "
            f"same canonical-position site mapped through each structure's own alignment"
        )
    if len(mobile.site_resseqs) < 3:
        raise ValueError(
            f"{mobile.label}: only {len(mobile.site_resseqs)} corresponding site residue(s) resolved "
            f"against {ref.label} -- need >=3 for a stable fit"
        )
    pairs = []
    for m_resi, r_resi in zip(mobile.site_resseqs, ref.site_resseqs):
        pairs.append(f"{mobile.label} and chain {mobile.chain_id} and resi {m_resi} and name CA and alt ''+A")
        pairs.append(f"{ref.label} and chain {ref.chain_id} and resi {r_resi} and name CA and alt ''+A")
    rmsd = cmd.pair_fit(*pairs)
    return AlignmentResult(mobile.label, ref.label, "site", rmsd, len(mobile.site_resseqs), "")


def _cealign(cmd, mobile: StructureInput, ref: StructureInput) -> AlignmentResult:
    result = cmd.cealign(
        f"{ref.label} and chain {ref.chain_id} and alt ''+A",
        f"{mobile.label} and chain {mobile.chain_id} and alt ''+A",
    )
    return AlignmentResult(mobile.label, ref.label, "none", result["RMSD"], int(result["alignment_length"]), "")


def align_structures(
    structures: Sequence[StructureInput], reference_label: str, out_dir: Union[str, Path], *,
    site_mode: str = "none", show_progress: bool = True,
) -> List[AlignmentResult]:
    """Superpose every structure in `structures` onto `reference_label`
    and save each (including the unmoved reference, for a consistent
    output set) to `out_dir/{label}_aligned.pdb`. `site_mode` selects
    `_pair_fit` (`"pocket"`/`"ligand"` -- requires `site_resseqs` set on
    every non-reference `StructureInput`) or `_cealign` (`"none"`).

    `show_progress` prints one line per completed fit (`print(...,
    flush=True)`) -- PyMOL's own `cmd.pair_fit` already writes an
    unlabeled `ExecutiveRMSPairs: RMSD = ...` line to stdout on its own,
    which is not disabled here, but doesn't say *which* structure it was
    for when fitting several in a row, hence this project's own labeled
    line alongside it.
    """
    if site_mode not in SITE_MODES:
        raise ValueError(f"site_mode must be one of {SITE_MODES}, got {site_mode!r}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    by_label = {s.label: s for s in structures}
    if reference_label not in by_label:
        raise ValueError(f"reference_label {reference_label!r} not found in structures")
    reference = by_label[reference_label]

    import pymol2

    results: List[AlignmentResult] = []
    with pymol2.PyMOL() as session:
        cmd = session.cmd
        for s in structures:
            cmd.load(s.pdb_path, s.label)

        ref_out = out_dir / f"{reference.label}_aligned.pdb"
        cmd.save(str(ref_out), reference.label)
        results.append(AlignmentResult(reference.label, reference.label, site_mode, 0.0, 0, str(ref_out)))
        if show_progress:
            print(f"[structalign] {reference.label}: reference (unmoved) -> {ref_out.name}", flush=True)

        mobiles = [s for s in structures if s.label != reference_label]
        for i, s in enumerate(mobiles, start=1):
            try:
                result = _cealign(cmd, s, reference) if site_mode == "none" else _pair_fit(cmd, s, reference)
            except Exception as e:
                # One structure genuinely not comparable at the requested site
                # (e.g. a co-complex where the site's canonical positions
                # simply aren't resolved in this particular entry) shouldn't
                # abort the fit for every other structure in the batch.
                results.append(AlignmentResult(s.label, reference_label, site_mode, None, 0, "", error=str(e)))
                if show_progress:
                    print(f"[structalign] ({i}/{len(mobiles)}) {s.label}: SKIPPED ({e})", flush=True)
                continue
            out_pdb = out_dir / f"{s.label}_aligned.pdb"
            cmd.save(str(out_pdb), s.label)
            result.aligned_pdb = str(out_pdb)
            results.append(result)
            if show_progress:
                print(f"[structalign] ({i}/{len(mobiles)}) {s.label}: rmsd={result.rmsd:.3f} ({result.n_atoms} atoms) -> {out_pdb.name}", flush=True)

    return results
