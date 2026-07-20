"""Tests for nucleotide genomes — the coordinate representation, inversion, and the multi-chromosome,
identity-bearing karyotype wired along the species tree.

A chromosome is an ordered list of blocks (maximal runs of one ancestry). The strongest checks are
the oracle (apply the same inversions to a plain per-nucleotide array and require ``trace_back`` to
agree) and the strong invariant (inversion conserves ancestry, so every node carries the whole root
sequence, permuted across its chromosomes).
"""

import numpy as np
import pytest

from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_nucleotide
from zombi2.genomes.nucleotide import Block, Chromosome, NucleotideGenome, _split_block


def _chrom(length, topology="circular", cid=0, source=0):
    return Chromosome(cid, topology, [Block(source, 0, length, 1)])


def _ancestry(chrom):
    return sorted((src, pos) for (src, pos, _s) in chrom.trace_back())


def _is_maximal(chrom):
    for a, b in zip(chrom.blocks, chrom.blocks[1:]):
        if a.source == b.source and a.strand == b.strand:
            if (a.strand == 1 and a.end == b.start) or (a.strand == -1 and a.start == b.end):
                return False
    return True


# --- blocks + split -------------------------------------------------------------------------------

def test_split_block_is_strand_aware():
    assert _split_block(Block(0, 10, 20, 1), 3) == [Block(0, 10, 13, 1), Block(0, 13, 20, 1)]
    # reversed: the first o physical positions are the HIGH end of the source
    assert _split_block(Block(0, 10, 20, -1), 3) == [Block(0, 17, 20, -1), Block(0, 10, 17, -1)]


def test_chromosome_seed_and_readers():
    ch = _chrom(5)
    assert ch.length == 5
    assert ch.mosaic() == [(0, 0, 5, 1)]
    assert ch.trace_back() == [(0, 0, 1), (0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)]


def test_topology_validation():
    with pytest.raises(ValueError):
        Chromosome(0, "loop", [Block(0, 0, 5, 1)])


# --- inversion ------------------------------------------------------------------------------------

def test_inversion_hand_example():
    ch = _chrom(5)
    ch.invert(1, 3)
    assert ch.mosaic() == [(0, 0, 1, 1), (0, 1, 4, -1), (0, 4, 5, 1)]
    assert ch.trace_back() == [(0, 0, 1), (0, 3, -1), (0, 2, -1), (0, 1, -1), (0, 4, 1)]


def test_circular_inversion_wraps_the_origin():
    ch = _chrom(6, "circular")
    ch.invert(4, 4)                                          # arc 4,5,0,1 wraps the origin
    assert ch.trace_back() == [(0, 1, -1), (0, 0, -1), (0, 5, -1), (0, 4, -1), (0, 2, 1), (0, 3, 1)]
    assert ch.length == 6


def test_linear_inversion_clamps_and_never_wraps():
    ch = _chrom(6, "linear")
    ch.invert(4, 10)                                         # overruns the end -> clamp to [4, 6)
    assert ch.trace_back() == [(0, 0, 1), (0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 5, -1), (0, 4, -1)]


def test_inversion_flips_strands_in_the_arc_only():
    ch = _chrom(10)
    ch.invert(3, 4)
    assert [s for (_x, _p, s) in ch.trace_back()] == [1, 1, 1, -1, -1, -1, -1, 1, 1, 1]


# --- blocks stay maximal (canonicalisation) -------------------------------------------------------

def test_invert_keeps_blocks_maximal():
    ch = _chrom(20)
    rng = np.random.default_rng(0)
    for _ in range(25):
        ch.invert(int(rng.integers(20)), int(rng.integers(1, 10)))
        assert _is_maximal(ch)                               # never two adjacent collinear blocks
        assert ch.length == 20


def test_inversion_then_inverse_canonicalises_to_one_block():
    ch = _chrom(12)
    ch.invert(2, 6)
    ch.invert(2, 6)                                          # undone -> the seam heals
    assert ch.mosaic() == [(0, 0, 12, 1)]                    # a single maximal block again


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
    if s + ell <= n:
        arr[s:s + ell] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[s:s + ell])]
    else:
        arr[:] = arr[s:] + arr[:s]
        arr[:ell] = [(src, pos, -strand) for (src, pos, strand) in reversed(arr[:ell])]


@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_trace_back_matches_the_oracle_under_random_inversions(topology):
    rng = np.random.default_rng(0)
    for _trial in range(50):
        n = int(rng.integers(5, 40))
        ch = _chrom(n, topology)
        arr = [(0, p, 1) for p in range(n)]
        for _ in range(int(rng.integers(1, 15))):
            s, ell = int(rng.integers(0, n)), int(rng.integers(0, n))
            ch.invert(s, ell)
            _oracle_invert(arr, s, ell, topology)
            assert ch.trace_back() == arr                    # exact, every step (canonicalise is trace-neutral)
        assert ch.length == n


# --- the karyotype (NucleotideGenome) -------------------------------------------------------------

def test_karyotype_length_and_readers():
    g = NucleotideGenome([Chromosome(0, "circular", [Block(0, 0, 100, 1)]),
                          Chromosome(1, "linear", [Block(1, 0, 40, 1)])])
    assert g.length == 140
    assert set(g.mosaic()) == {0, 1}                         # keyed by chromosome id
    assert g.ancestry() == sorted([(0, p) for p in range(100)] + [(1, p) for p in range(40)])


# --- the tree wiring: heterogeneous seeding, identity, the strong invariant -----------------------

def _run(seed=1, inversion=0.04, chromosomes=1, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=seed)
    params = dict(inversion=inversion, inversion_length=15, root_length=120, chromosomes=chromosomes,
                  seed=seed)
    params.update(kw)
    return sp, simulate_genomes_nucleotide(sp, **params)


def test_every_node_has_a_karyotype_of_conserved_length():
    sp, r = _run(seed=2, chromosomes=3, root_length=100)
    assert set(r.genomes) == set(sp.complete_tree.nodes)
    assert all(g.length == 300 for g in r.genomes.values())  # 3 x 100, inversion conserves length
    assert all(len(g.chromosomes) == 3 for g in r.genomes.values())   # no tier yet: count conserved


def test_heterogeneous_seeding_sizes_and_shapes():
    specs = [(100, "circular"), (40, "circular"), (25, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=4, seed=5)
    r = simulate_genomes_nucleotide(sp, inversion=0.0, chromosomes=specs, seed=5)
    tip = sorted(n.id for n in sp.complete_tree.extant())[0]
    chroms = r.genomes[tip].chromosomes
    assert [c.topology for c in chroms] == ["circular", "circular", "linear"]
    assert sorted(c.length for c in chroms) == [25, 40, 100]     # the three sizes, preserved
    # with no inversion, each replicon is a single block under its own source
    assert {c.blocks[0].source for c in chroms} == {0, 1, 2}


def test_every_node_carries_the_whole_root_sequence_permuted():
    # the strong invariant across chromosomes: inheritance copies ancestry, inversion changes none
    specs = [(100, "circular"), (40, "circular"), (25, "linear")]
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=3)
    r = simulate_genomes_nucleotide(sp, inversion=0.06, inversion_length=10, chromosomes=specs, seed=3)
    full = sorted((s, p) for s, (length, _t) in enumerate(specs) for p in range(length))
    for node_id in r.genomes:
        assert r.ancestry(node_id) == full


def test_chromosome_identity_network_is_well_formed():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=4)
    r = simulate_genomes_nucleotide(sp, inversion=0.02,
                                    chromosomes=[(80, "circular"), (30, "linear")], seed=4)
    ce = r.chromosome_events
    roots = [e for e in ce if e.kind == "origination"]
    spec = [e for e in ce if e.kind == "speciation"]
    assert len(roots) == 2                                    # one origination per seed replicon
    assert all(e.parents == () and len(e.children) == 1 for e in roots)
    assert spec and all(len(e.parents) == 1 and len(e.children) == 2 for e in spec)
    born = {}
    for e in ce:                                             # every chromosome id minted exactly once
        for ch in e.children:
            assert ch not in born
            born[ch] = e
    assert all(p in born for e in ce for p in e.parents)     # every speciation parent born earlier


def test_zero_inversion_leaves_each_replicon_a_single_block():
    sp, r = _run(seed=6, inversion=0.0, chromosomes=2, root_length=50)
    assert r.rearrangements == []
    for g in r.genomes.values():
        assert all(len(c.blocks) == 1 for c in g.chromosomes)


def test_inversions_recorded_within_their_branch_and_chromosome():
    sp, r = _run(seed=7, inversion=0.08, chromosomes=2, root_length=100)
    assert r.rearrangements
    node_chroms = {n: {c.id for c in r.genomes[n].chromosomes} for n in r.genomes}
    for inv in r.rearrangements:
        node = sp.complete_tree.nodes[inv.lineage]
        assert node.birth_time <= inv.time <= node.end_time
        assert inv.chromosome in node_chroms[inv.lineage]    # fired on a chromosome that node carried
        assert inv.length >= 1


def test_deterministic_given_seed():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=8)
    kw = dict(inversion=0.05, inversion_length=15, chromosomes=[(120, "circular"), (30, "linear")], seed=8)
    a = simulate_genomes_nucleotide(sp, **kw)
    b = simulate_genomes_nucleotide(sp, **kw)
    assert all(a.mosaic(n) == b.mosaic(n) for n in a.genomes)
    assert a.rearrangements == b.rearrangements
    assert a.chromosome_events == b.chromosome_events


def test_simulate_validation():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=5, seed=1)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, inversion=-1.0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, inversion_length=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, chromosomes=0)
    with pytest.raises(ValueError):
        simulate_genomes_nucleotide(sp, chromosomes=[(0, "circular")])   # a bad replicon length
