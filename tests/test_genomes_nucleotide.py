"""Tests for nucleotide genomes — slice 1: the coordinate representation + split + inversion.

The genome is a signed segment list over a nucleotide coordinate axis. The strongest check is the
oracle: apply the same inversions to a plain per-nucleotide array (the ground truth) and to the
``NucleotideGenome``, and require ``trace_back`` to agree.
"""

import numpy as np
import pytest

from zombi2.genomes.nucleotide import NucleotideGenome, Segment, _split_segment


# --- seeding + the readers -----------------------------------------------------------------------

def test_seed_is_one_forward_segment_of_the_root():
    g = NucleotideGenome.seed(5)
    assert g.length == 5
    assert g.mosaic() == [(0, 0, 5, 1)]
    assert g.trace_back() == [(0, 0, 1), (0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)]


def test_seed_rejects_empty():
    with pytest.raises(ValueError):
        NucleotideGenome.seed(0)


# --- split_at ------------------------------------------------------------------------------------

def test_split_at_creates_a_boundary_without_changing_the_sequence():
    g = NucleotideGenome.seed(5)
    before = g.trace_back()
    g.split_at(2)
    assert g.mosaic() == [(0, 0, 2, 1), (0, 2, 5, 1)]     # one boundary added
    assert g.trace_back() == before                       # the nucleotides are unchanged
    assert g.length == 5


def test_split_at_is_a_noop_at_ends_and_existing_boundaries():
    g = NucleotideGenome.seed(5)
    g.split_at(0)
    g.split_at(5)
    assert g.mosaic() == [(0, 0, 5, 1)]                    # nothing split at the two ends
    g.split_at(2)
    g.split_at(2)                                          # second split at the same place is a no-op
    assert g.mosaic() == [(0, 0, 2, 1), (0, 2, 5, 1)]


def test_split_segment_is_strand_aware():
    # forward: the first o physical positions are the LOW end of the source
    assert _split_segment(Segment(0, 10, 20, 1), 3) == [Segment(0, 10, 13, 1), Segment(0, 13, 20, 1)]
    # reversed: the first o physical positions are the HIGH end of the source
    assert _split_segment(Segment(0, 10, 20, -1), 3) == [Segment(0, 17, 20, -1), Segment(0, 10, 17, -1)]


# --- inversion -----------------------------------------------------------------------------------

def test_inversion_hand_example():
    g = NucleotideGenome.seed(5)
    g.invert(1, 3)                                         # reverse source positions 1,2,3
    assert g.mosaic() == [(0, 0, 1, 1), (0, 1, 4, -1), (0, 4, 5, 1)]
    assert g.trace_back() == [(0, 0, 1), (0, 3, -1), (0, 2, -1), (0, 1, -1), (0, 4, 1)]


def _ancestry(g):
    # the multiset of ancestral (source, position) — orientation-agnostic (inversion flips strand)
    return sorted((src, pos) for (src, pos, _strand) in g.trace_back())


def test_inversion_conserves_length_and_ancestry():
    g = NucleotideGenome.seed(20)
    origin = _ancestry(g)
    g.invert(4, 9)
    g.invert(0, 20)
    g.invert(7, 5)
    assert g.length == 20                                  # length is conserved
    assert _ancestry(g) == origin                          # every nucleotide is still present, once


def test_inversion_flips_strands_in_the_arc_only():
    g = NucleotideGenome.seed(10)
    g.invert(3, 4)                                         # positions 3,4,5,6
    tb = g.trace_back()
    assert [s for (_, _, s) in tb] == [1, 1, 1, -1, -1, -1, -1, 1, 1, 1]


def test_an_inversion_is_its_own_inverse():
    g = NucleotideGenome.seed(12)
    before = g.trace_back()
    g.invert(2, 6)
    g.invert(2, 6)                                         # the same inversion again restores it
    assert g.trace_back() == before                        # identical at the nucleotide level


def test_inversion_of_a_reversed_region_uses_the_strand_aware_split():
    g = NucleotideGenome.seed(8)
    g.invert(0, 8)                                          # whole genome now reversed
    g.invert(2, 3)                                          # invert a sub-arc of the reversed genome
    # the strand-aware split must keep every ancestral nucleotide present exactly once
    assert g.length == 8
    assert _ancestry(g) == _ancestry(NucleotideGenome.seed(8))


# --- the oracle: random inversions vs a per-nucleotide array -------------------------------------

def _oracle_invert(arr, start, length):
    start = max(0, min(start, len(arr)))
    end = min(start + max(0, length), len(arr))
    if end <= start:
        return
    arr[start:end] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[start:end])]


def test_trace_back_matches_the_oracle_under_random_inversions():
    rng = np.random.default_rng(0)
    for trial in range(50):
        L = int(rng.integers(5, 40))
        g = NucleotideGenome.seed(L)
        arr = [(0, p, 1) for p in range(L)]
        for _ in range(int(rng.integers(1, 15))):
            s = int(rng.integers(0, L))
            ell = int(rng.integers(0, L))
            g.invert(s, ell)
            _oracle_invert(arr, s, ell)
            assert g.trace_back() == arr                   # exact agreement, every step
        assert g.length == L
