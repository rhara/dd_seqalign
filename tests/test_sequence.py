"""Fast, offline unit tests for sequence.py's alignment logic (no
network, no PyMOL/fpocket) -- the parts of dd_seq that don't need a real
structure file or external tool to exercise meaningfully."""
from dd_seq.sequence import ChainSequence, align_to_canonical, pick_target_chain

CANONICAL = "MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRHPNIVSLQDVLMQDSRLYLIFEFLSMDLKKYLDSI"


def test_align_to_canonical_full_coverage_no_mismatch():
    chain = ChainSequence(chain_id="A", residues=list(enumerate(CANONICAL, start=1)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.coverage == 1.0
    assert aln.n_mismatch == 0
    assert aln.resseq_for_canonical(1) == 1
    assert aln.canonical_for_resseq(1) == 1


def test_align_to_canonical_fragment_has_free_end_gaps():
    fragment = CANONICAL[10:40]  # a co-crystal that only resolved residues 11-40
    chain = ChainSequence(chain_id="A", residues=list(enumerate(fragment, start=11)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_covered == 30
    assert aln.n_mismatch == 0
    assert aln.resseq_for_canonical(1) is None  # outside the fragment: missing, not a bad alignment
    assert aln.resseq_for_canonical(11) == 11
    assert aln.resseq_for_canonical(40) == 40


def test_align_to_canonical_flags_point_mutation_not_missing():
    mutant = CANONICAL[:20] + "A" + CANONICAL[21:]  # one substituted residue
    chain = ChainSequence(chain_id="A", residues=list(enumerate(mutant, start=1)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_mismatch == 1
    mismatch = next(r for r in aln.residues if r.status == "mismatch")
    assert mismatch.canonical_pos == 21
    assert mismatch.structure_code == "A"


def test_pick_target_chain_prefers_identity_over_raw_coverage():
    # A homologous paralog chain can align across the same span (same raw
    # coverage) as the true target chain while being mostly mismatches (see
    # CAK-CDK1-cyclinB1, 9SKQ, in the README) -- pick_target_chain must
    # rank by matching residues, not just how much of the canonical
    # sequence got covered.
    true_chain = ChainSequence(chain_id="A", residues=list(enumerate(CANONICAL, start=1)))
    paralog_seq = "".join("X" if i % 2 == 0 else c for i, c in enumerate(CANONICAL))
    paralog_chain = ChainSequence(chain_id="B", residues=list(enumerate(paralog_seq, start=1)))

    alignments = {
        "A": align_to_canonical(true_chain, CANONICAL),
        "B": align_to_canonical(paralog_chain, CANONICAL),
    }
    assert pick_target_chain(alignments) == "A"
