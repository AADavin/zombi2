"""Tests for nucleotide genomes — slice 1: the coordinate representation + split + inversion.

The genome is a signed segment list over a nucleotide coordinate axis. The strongest check is the
oracle: apply the same inversions to a plain per-nucleotide array (the ground truth) and to the
``NucleotideGenome``, and require ``trace_back`` to agree.
"""

import numpy as np
import pytest

from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_nucleotide
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


# --- topology: circular (default) vs linear ------------------------------------------------------

def test_seed_topology_default_is_circular_and_validates():
    assert NucleotideGenome.seed(5).topology == "circular"          # canonical bacterial default
    assert NucleotideGenome.seed(5, "linear").topology == "linear"
    with pytest.raises(ValueError):
        NucleotideGenome.seed(5, "loop")


def test_circular_inversion_wraps_the_origin():
    # seed(6): positions 0..5. invert(4, 4) spans 4,5,0,1 — it wraps the origin. The ring rotates so
    # the arc is contiguous (the physical origin drifts), then the arc reverses & flips.
    g = NucleotideGenome.seed(6, "circular")
    g.invert(4, 4)
    assert g.trace_back() == [(0, 1, -1), (0, 0, -1), (0, 5, -1), (0, 4, -1), (0, 2, 1), (0, 3, 1)]
    assert g.length == 6


def test_linear_inversion_clamps_and_never_wraps():
    g = NucleotideGenome.seed(6, "linear")
    g.invert(4, 10)                                          # length overruns the end -> clamp to [4, 6)
    assert g.trace_back() == [(0, 0, 1), (0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 5, -1), (0, 4, -1)]


# --- the oracle: random inversions vs a per-nucleotide array (both topologies) --------------------

def _oracle_invert(arr, start, length, topology):
    n = len(arr)
    if n == 0:
        return
    if topology == "linear":
        start = max(0, min(start, n))
        end = min(start + max(0, length), n)
        if end <= start:
            return
        arr[start:end] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[start:end])]
        return
    ell = max(1, min(length, n))
    s = start % n
    if s + ell <= n:                                         # non-wrapping: reverse in place
        arr[s:s + ell] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[s:s + ell])]
    else:                                                    # wrapping: rotate arc to front, then reverse
        arr[:] = arr[s:] + arr[:s]
        arr[:ell] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[:ell])]


@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_trace_back_matches_the_oracle_under_random_inversions(topology):
    rng = np.random.default_rng(0)
    for _trial in range(50):
        length_total = int(rng.integers(5, 40))
        g = NucleotideGenome.seed(length_total, topology)
        arr = [(0, p, 1) for p in range(length_total)]
        for _ in range(int(rng.integers(1, 15))):
            s = int(rng.integers(0, length_total))
            ell = int(rng.integers(0, length_total))
            g.invert(s, ell)
            _oracle_invert(arr, s, ell, topology)
            assert g.trace_back() == arr                     # exact agreement, every step
        assert g.length == length_total


# --- slice A: the tree wiring (inversions along the species tree) ---------------------------------

def _run(seed=1, inversion=0.02, root_length=200, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
    params = dict(inversion=inversion, root_length=root_length, inversion_length=20, seed=seed)
    params.update(kw)
    return sp, simulate_genomes_nucleotide(sp, **params)


def test_every_node_has_a_genome_of_conserved_length():
    sp, r = _run(seed=2)
    assert set(r.genomes) == set(sp.complete_tree.nodes)     # a genome at every node of the complete tree
    assert all(g.length == 200 for g in r.genomes.values())  # inversion conserves length everywhere


def test_every_node_carries_the_whole_root_sequence_permuted():
    # the strong invariant: inheritance copies the ancestry and an inversion changes none of it, so
    # EVERY node's trace_back is a signed permutation of the root's nucleotides — all present once
    sp, r = _run(seed=3)
    root_positions = sorted((0, p) for p in range(200))
    for node_id in r.genomes:
        assert sorted((src, pos) for (src, pos, _s) in r.trace_back(node_id)) == root_positions


def test_zero_inversion_leaves_every_genome_as_the_seed():
    sp, r = _run(seed=4, inversion=0.0)
    assert r.rearrangements == []
    assert all(r.mosaic(n) == [(0, 0, 200, 1)] for n in r.genomes)


def test_inversions_fire_and_are_recorded_within_their_branch():
    sp, r = _run(seed=5, inversion=0.05)
    assert r.rearrangements                                   # inversions really happened
    for inv in r.rearrangements:
        node = sp.complete_tree.nodes[inv.lineage]
        assert node.birth_time <= inv.time <= node.end_time   # each fired within its branch's interval
        assert 0 <= inv.start < 200 and inv.length >= 1


def test_deterministic_given_seed():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=6)
    kw = dict(inversion=0.03, root_length=200, inversion_length=20, seed=6)
    a = simulate_genomes_nucleotide(sp, **kw)
    b = simulate_genomes_nucleotide(sp, **kw)
    assert all(a.mosaic(n) == b.mosaic(n) for n in a.genomes)
    assert a.rearrangements == b.rearrangements


def test_linear_topology_passes_through():
    _sp, r = _run(seed=7, topology="linear")
    assert all(g.topology == "linear" for g in r.genomes.values())


def test_simulate_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=5, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, inversion=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, inversion_length=0)
