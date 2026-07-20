"""Per-chain sequence extraction from a structure file, and a glocal
(free-end-gap) pairwise alignment of each chain against the UniProt
canonical sequence.

A full multiple-sequence-alignment tool (mafft/clustalo) is deliberately
not used here: every input is a fragment/oligomer/isoform of the *same*
protein (co-crystal fragments, AlphaFold's full-length model), not a set of
divergent homologs, so there is a single natural reference -- the UniProt
canonical sequence -- and comparing every structure to that one reference
with Biopython's pairwise aligner is both simpler and gives a more directly
useful result (a per-canonical-residue coverage/mismatch table) than a
generic multi-sequence alignment would.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

Residue = Tuple[int, str]  # (author resseq, one-letter code)

# Bio.PDB.Polypeptide.three_to_one was removed in recent Biopython; this is
# the same strict standard-20-only mapping it used to provide (deliberately
# not Bio.SeqUtils.seq1, which silently maps any unrecognized/modified
# residue, e.g. MSE, SEP, to 'X' -- we'd rather skip a residue we can't
# confidently call than let a fabricated 'X' pollute the sequence used for
# alignment).
_THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


@dataclass
class ChainSequence:
    chain_id: str
    residues: List[Residue]  # ordered by structure position, not necessarily contiguous resseq

    @property
    def sequence(self) -> str:
        return "".join(code for _, code in self.residues)


def extract_chain_sequences(pdb_path: Union[str, Path]) -> Dict[str, ChainSequence]:
    """One `ChainSequence` per protein chain in the structure, in
    structure (not necessarily author-numbering-contiguous) order. Non-
    standard residues that `Bio.PDB.Polypeptide.is_aa(..., standard=False)`
    still recognizes as an amino acid but `three_to_one` cannot map
    (unusual modified residues) are skipped rather than raising, so one
    oddly-modified residue doesn't abort extraction for the whole chain.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(Path(pdb_path).stem, str(pdb_path))
    model = next(iter(structure))  # first model (NMR/EM multi-model files: take the first)

    chains: Dict[str, ChainSequence] = {}
    for chain in model:
        residues: List[Residue] = []
        for res in chain:
            if not is_aa(res, standard=False):
                continue
            code = _THREE_TO_ONE.get(res.get_resname().strip().upper())
            if code is None:
                continue
            resseq = res.id[1]
            residues.append((resseq, code))
        if residues:
            chains[chain.id] = ChainSequence(chain_id=chain.id, residues=residues)
    return chains


@dataclass
class ResidueAlignment:
    canonical_pos: int  # 1-indexed position in the canonical sequence
    canonical_code: str
    structure_resseq: Optional[int]  # None if this canonical position is absent from the structure
    structure_code: Optional[str]
    status: str  # "match" | "mismatch" | "missing"


@dataclass
class ChainAlignment:
    chain_id: str
    residues: List[ResidueAlignment] = field(default_factory=list)

    @property
    def n_covered(self) -> int:
        return sum(1 for r in self.residues if r.status != "missing")

    @property
    def n_mismatch(self) -> int:
        return sum(1 for r in self.residues if r.status == "mismatch")

    @property
    def coverage(self) -> float:
        return self.n_covered / len(self.residues) if self.residues else 0.0

    def resseq_for_canonical(self, canonical_pos: int) -> Optional[int]:
        """Author residue number in this structure for a given canonical
        UniProt position, or None if that position isn't resolved here --
        the lookup `activesite.py` needs to translate a reference-defined
        site into each structure's own numbering."""
        idx = canonical_pos - 1
        if 0 <= idx < len(self.residues):
            return self.residues[idx].structure_resseq
        return None

    def canonical_for_resseq(self, resseq: int) -> Optional[int]:
        """Inverse of `resseq_for_canonical`: the canonical UniProt position
        for one of this structure's own author residue numbers -- used to
        translate a site detected directly on a structure (e.g. ligand-
        proximal residues, fpocket lining residues) into the common
        canonical-numbering coordinate system before comparing it across
        structures."""
        if self._reverse is None:
            self._reverse = {
                r.structure_resseq: r.canonical_pos for r in self.residues if r.structure_resseq is not None
            }
        return self._reverse.get(resseq)

    _reverse: Optional[Dict[int, int]] = field(default=None, repr=False, compare=False)


def pick_target_chain(chain_alignments: Dict[str, "ChainAlignment"]) -> str:
    """The chain that is actually the canonical protein, ranked by number
    of identical (matching) residues, not raw coverage -- e.g. picks the
    CDK1 chain out of a CDK1/CyclinB/Cks2 co-crystal so downstream active-
    site/structural-alignment steps operate on the protein of interest,
    not a bound partner chain that happens to also be a polypeptide.
    Coverage alone is not enough to tell apart the real ortholog from a
    homologous paralog: in a CAK-CDK1-cyclinB1 assembly, the CDK7 chain
    (CAK's own kinase) aligns to the CDK1 canonical sequence with even
    *higher* coverage than the true CDK1 chain, but the vast majority of
    those aligned positions are mismatches -- so this ranks by
    `n_covered - n_mismatch` (matching residues) instead."""
    return max(chain_alignments, key=lambda cid: chain_alignments[cid].n_covered - chain_alignments[cid].n_mismatch)


_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALIGNER.open_gap_score = -10
_ALIGNER.extend_gap_score = -0.5
_ALIGNER.end_insertion_score = 0.0  # canonical residues outside the structure's span: free
_ALIGNER.end_deletion_score = 0.0  # structure residues outside... (shouldn't happen, same reason)


def align_to_canonical(chain_seq: ChainSequence, canonical_seq: str) -> ChainAlignment:
    """Glocal-align one chain's observed residues against the full-length
    canonical sequence and return a per-canonical-position table: which
    structures cover which residues, and where a resolved residue differs
    from canonical (a real point mutation in the crystallized construct,
    not a modeling artifact -- missing density shows up as "missing", not
    "mismatch")."""
    alignment = _ALIGNER.align(canonical_seq, chain_seq.sequence)[0]
    result = ChainAlignment(chain_id=chain_seq.chain_id)
    covered = [False] * len(canonical_seq)
    slot: List[Optional[Residue]] = [None] * len(canonical_seq)

    t_blocks, q_blocks = alignment.aligned
    for (t_start, t_end), (q_start, q_end) in zip(t_blocks, q_blocks):
        for i in range(t_end - t_start):
            canon_idx = t_start + i
            chain_idx = q_start + i
            covered[canon_idx] = True
            slot[canon_idx] = chain_seq.residues[chain_idx]

    for canon_idx, canon_code in enumerate(canonical_seq):
        if covered[canon_idx]:
            resseq, code = slot[canon_idx]
            status = "match" if code == canon_code else "mismatch"
            result.residues.append(ResidueAlignment(canon_idx + 1, canon_code, resseq, code, status))
        else:
            result.residues.append(ResidueAlignment(canon_idx + 1, canon_code, None, None, "missing"))
    return result
