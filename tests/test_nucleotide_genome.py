"""Nucleotide-level circular genome, M1 (inversions only).

The workhorse is a dead-simple O(L) **array oracle** that applies each inversion by
literal circular slice-reverse over one cell per nucleotide. It is obviously correct and
hopelessly slow; running it against the efficient segment structure on the *same*
(s, length) stream pins down every circular-geometry bug. On top of that we assert the
strong invariants inversion gives us (content conservation, bijection, involution) and
that the trace-back / atom decomposition is self-consistent.
"""

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.events import EventType, TargetParams
from zombi2.genome import IdManager
from zombi2.nucleotide_genome import NucleotideGenome
from zombi2.nucleotide_sim import simulate_nucleotide_genomes


# --------------------------------------------------------------------------- #
# Independent brute-force oracle: one cell per nucleotide, fixed origin.
# --------------------------------------------------------------------------- #
class ArrayGenome:
    """One cell per nucleotide; same rotate-on-wrap origin convention as the genome."""

    def __init__(self, source, length):
        self.cells = [(source, i, 1) for i in range(length)]

    def invert(self, s, ell):
        L = len(self.cells)
        ell = max(1, min(ell, L))
        s %= L
        if s + ell > L:  # wrapping: rotate the origin to s, then reverse the prefix
            self.cells = self.cells[s:] + self.cells[:s]
            s = 0
        seg = [(src, p, -st) for (src, p, st) in self.cells[s:s + ell][::-1]]
        self.cells[s:s + ell] = seg


def _fresh(length, ext=0.99):
    g = NucleotideGenome(IdManager(), root_length=length, extension=ext)
    g.originate(np.random.default_rng(0), TargetParams())  # seed the root chromosome
    return g


# --------------------------------------------------------------------------- #
# Seeding
# --------------------------------------------------------------------------- #
def test_seed_is_single_forward_segment():
    g = _fresh(1000)
    assert g.size() == 1000
    assert g.n_segments() == 1
    assert g.to_cells() == [("1", i, 1) for i in range(1000)]


# --------------------------------------------------------------------------- #
# Geometry vs the oracle (the workhorse)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(25))
def test_inversion_matches_oracle_random(seed):
    rng = np.random.default_rng(seed)
    L = int(rng.integers(20, 200))
    g = _fresh(L, ext=0.9)
    o = ArrayGenome("1", L)
    for _ in range(80):
        s = int(rng.integers(L))
        ell = int(rng.integers(1, L + 1))
        g._apply_inversion(s, ell)
        o.invert(s, ell)
        assert g.to_cells() == o.cells       # identical content & layout
        assert g.size() == L                  # inversion preserves length


def test_wrapping_and_edge_cases_match_oracle():
    L = 12
    g = _fresh(L, ext=0.7)
    o = ArrayGenome("1", L)
    for s, ell in [(8, 6), (0, 12), (11, 3), (5, 12), (3, 9), (7, 5),
                   (0, 1), (11, 1), (10, 4), (6, 12)]:
        g._apply_inversion(s, ell)  # wrapping arcs, whole-genome, length-1
        o.invert(s, ell)
        assert g.to_cells() == o.cells


def test_whole_genome_inversion_reverses_all():
    L = 8
    g = _fresh(L)
    g._apply_inversion(0, L)
    assert g.to_cells() == [("1", p, -1) for p in range(L - 1, -1, -1)]


def test_length_one_inversion_flips_one_strand():
    g = _fresh(6)
    before = g.to_cells()
    g._apply_inversion(3, 1)
    after = g.to_cells()
    assert after[3] == (before[3][0], before[3][1], -before[3][2])
    assert after[:3] == before[:3] and after[4:] == before[4:]


def test_inversion_involution():
    g = _fresh(50, ext=0.9)
    before = g.to_cells()
    g._apply_inversion(12, 20)
    g._apply_inversion(12, 20)  # same arc twice == identity on content
    assert g.to_cells() == before


def test_split_preserves_content():
    g = _fresh(30)
    before = g.to_cells()
    g._split_at(7)
    g._split_at(19)
    g._split_at(7)  # idempotent
    assert g.n_segments() == 3
    assert g.to_cells() == before  # a breakpoint reorders nothing


# --------------------------------------------------------------------------- #
# Tree-level: content conservation, bijection, reproducibility
# --------------------------------------------------------------------------- #
def test_content_conserved_and_bijective_at_leaves():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    res = simulate_nucleotide_genomes(tree, inversion=0.02, root_length=200,
                                      extension=0.9, seed=7)
    assert any(r.event is EventType.INVERSION for r in res.event_log)
    root = {("1", i) for i in range(200)}
    for genome in res.leaf_genomes.values():
        cells = genome.to_cells()
        origins = [(src, p) for (src, p, _st) in cells]
        assert len(cells) == 200
        assert set(origins) == root            # content conserved
        assert len(origins) == len(set(origins))  # bijection: each origin once


def test_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=2)
    a = simulate_nucleotide_genomes(tree, inversion=0.02, root_length=200, seed=9)
    b = simulate_nucleotide_genomes(tree, inversion=0.02, root_length=200, seed=9)
    assert len(a.event_log) == len(b.event_log)
    for leaf, ga in a.leaf_genomes.items():
        assert ga.to_cells() == b.leaf_genomes[leaf].to_cells()


# --------------------------------------------------------------------------- #
# Trace-back / atoms
# --------------------------------------------------------------------------- #
def test_atoms_tile_the_root():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=3)
    res = simulate_nucleotide_genomes(tree, inversion=0.03, root_length=150, seed=4)
    atoms = sorted((a for a in res.atoms if a.source == "1"), key=lambda a: a.start)
    assert atoms[0].start == 0 and atoms[-1].end == 150
    for x, y in zip(atoms, atoms[1:]):
        assert x.end == y.start                # contiguous, no gaps/overlaps
    n_inv = sum(1 for r in res.event_log if r.event is EventType.INVERSION)
    assert len(atoms) <= 2 * n_inv + 1         # each inversion adds <= 2 breakpoints


def test_mosaic_reassembles_each_leaf():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=5)
    res = simulate_nucleotide_genomes(tree, inversion=0.03, root_length=150,
                                      extension=0.9, seed=6)
    amap = {a.atom_id: a for a in res.atoms}
    for leaf, genome in res.leaf_genomes.items():
        mosaic = res.leaf_mosaic(leaf)
        # every atom appears exactly once per leaf (inversion loses nothing)
        assert sorted(aid for aid, _ in mosaic) == sorted(a.atom_id for a in res.atoms)
        cells = []
        for aid, strand in mosaic:
            a = amap[aid]
            if strand == 1:
                cells.extend((a.source, p, 1) for p in range(a.start, a.end))
            else:
                cells.extend((a.source, p, -1) for p in range(a.end - 1, a.start - 1, -1))
        assert cells == genome.to_cells()      # mosaic reconstructs the leaf exactly


def test_atom_histories_track_inversions():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=3.0, seed=8)
    res = simulate_nucleotide_genomes(tree, inversion=0.05, root_length=120,
                                      extension=0.9, seed=8)
    histories = res.atom_histories()
    branches = {n.name for n in tree.nodes_preorder()}
    n_inv = sum(1 for r in res.event_log if r.event is EventType.INVERSION)
    if n_inv:
        assert any(histories[a.atom_id] for a in res.atoms)  # something recorded
    for events in histories.values():
        for branch, _t in events:
            assert branch in branches           # every entry is a real branch
