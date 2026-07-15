"""Orientation-flip on transposition for the OrderedGenome model.

A transposed segment can reinsert reverse-complemented (gene order reversed, every strand
flipped) with probability ``transposition_flip``. The default ``0.0`` preserves orientation
and must leave every existing ordered run byte-identical (it never even draws an ``rng``
value — the same guard as biased gene conversion).
"""

import numpy as np

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    Rates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.genomes.genome import IdManager, OrderedGene
from zombi2.genomes.events import EventType, Region, Selection, TargetParams


def _ordered(extension=0.6, flip=0.0):
    return lambda ids: OrderedGenome(ids, extension=extension, transposition_flip=flip)


def _chromosomes(genomes):
    """Leaf-name -> list of (family, orientation), sorted by leaf name (a stable fingerprint)."""
    return {
        leaf.name: [(g.family, g.orientation) for g in genome.chromosome]
        for leaf, genome in sorted(genomes.leaf_genomes.items(), key=lambda kv: kv[0].name)
    }


class _FakeRNG:
    """A minimal Generator stand-in for the transposition branch.

    ``integers(high)`` returns a fixed insertion index; ``random()`` returns a fixed draw, or
    raises if ``forbid_random`` — so a test can *prove* the flip==0 guard never draws.
    """

    def __init__(self, *, random_value=0.0, insert_at=0, forbid_random=False):
        self._random_value = random_value
        self._insert_at = insert_at
        self._forbid_random = forbid_random

    def random(self):
        if self._forbid_random:
            raise AssertionError("rng.random() drawn when transposition_flip == 0")
        return self._random_value

    def integers(self, high):
        assert 0 <= self._insert_at < high
        return self._insert_at


def _segment(genome, start, length):
    cid = next(iter(genome.chromosomes))  # the (single) chromosome's chrom_id
    genes = tuple(genome.chromosome[start:start + length])
    return Selection(genes=genes, region=Region(chromosome=cid, start=start, length=length))


# --- default (0.0) preserves current behaviour EXACTLY -----------------------------

def test_default_matches_explicit_zero_flip():
    """Passing the new default explicitly changes nothing (same seed -> identical leaves)."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)
    rates = Rates(transposition=0.3)
    a = simulate_genomes(tree, rates, initial_families=15, seed=1,
                         genome_factory=lambda ids: OrderedGenome(ids, extension=0.6))
    b = simulate_genomes(tree, rates, initial_families=15, seed=1,
                         genome_factory=_ordered(extension=0.6, flip=0.0))
    assert _chromosomes(a) == _chromosomes(b)


def test_flip_zero_never_draws_from_rng():
    """The flip==0 guard short-circuits: the transposition branch must not draw ``rng.random()``
    (this is what keeps an orientation-preserving run byte-identical to a flip-free engine)."""
    ids = IdManager()
    genome = OrderedGenome(ids, extension=0.6, transposition_flip=0.0)
    next(iter(genome.chromosomes.values())).genes = [
        OrderedGene("g1", "f1", 1), OrderedGene("g2", "f2", -1),
        OrderedGene("g3", "f3", 1), OrderedGene("g4", "f4", 1),
    ]
    before = [(g.gid, g.family, g.orientation) for g in genome.chromosome]
    sel = _segment(genome, start=0, length=2)
    rng = _FakeRNG(insert_at=0, forbid_random=True)  # .random() would raise
    genome.apply(EventType.TRANSPOSITION, sel, rng, TargetParams())  # must not raise
    # nothing was flipped; the same genes, same orientations, are still present
    assert sorted((g.gid, g.family, g.orientation) for g in genome.chromosome) == sorted(before)


# --- flip == 1.0 reinserts reverse-complemented ------------------------------------

def test_flip_one_reverses_order_and_flips_strands():
    """A deterministic single transposition with flip==1.0: the moved block is reversed and every
    strand is flipped, exactly like an inversion applied to the relocated segment."""
    ids = IdManager()
    genome = OrderedGenome(ids, extension=0.6, transposition_flip=1.0)
    a, b = OrderedGene("g1", "f1", 1), OrderedGene("g2", "f2", -1)
    c, d = OrderedGene("g3", "f3", 1), OrderedGene("g4", "f4", 1)
    next(iter(genome.chromosomes.values())).genes = [a, b, c, d]
    sel = _segment(genome, start=0, length=2)     # the block {a(+1), b(-1)}
    rng = _FakeRNG(random_value=0.0, insert_at=2)  # 0.0 < 1.0 -> flip; reinsert after [c, d]
    groups = genome.apply(EventType.TRANSPOSITION, sel, rng, TargetParams())

    # remaining genes keep their strands; the block is reversed with both strands flipped
    assert genome.chromosome == [c, d, b, a]
    assert a.orientation == -1 and b.orientation == 1
    assert c.orientation == 1 and d.orientation == 1
    # the logged "transposed" group follows the reinserted (reversed) order
    assert [op.gid for op in groups[0]] == ["g2", "g1"]


def test_flip_one_differs_from_no_flip_in_simulation():
    """Over a real ordered sim with transpositions, flip==1.0 produces a clearly different genome
    layout than flip==0.0 (reversed/strand-flipped relocated blocks)."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)
    rates = Rates(transposition=0.5)
    no_flip = simulate_genomes(tree, rates, initial_families=15, seed=1,
                               genome_factory=_ordered(flip=0.0))
    flipped = simulate_genomes(tree, rates, initial_families=15, seed=1,
                               genome_factory=_ordered(flip=1.0))
    # transpositions actually fired (otherwise the comparison is vacuous)
    assert any(r.event is EventType.TRANSPOSITION for r in no_flip.event_log)
    assert _chromosomes(no_flip) != _chromosomes(flipped)


def test_flip_inherited_through_speciation_clone():
    """clone_reminting() (speciation) must carry ``transposition_flip`` to the child genome."""
    ids = IdManager()
    parent = OrderedGenome(ids, extension=0.6, transposition_flip=0.7)
    next(iter(parent.chromosomes.values())).genes = [OrderedGene("g1", "f1", 1)]
    child, _ = parent.clone_reminting()
    assert isinstance(child, OrderedGenome)
    assert child.transposition_flip == 0.7
    assert child.extension == 0.6
