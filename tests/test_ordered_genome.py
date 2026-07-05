"""OrderedGenome (ZOMBI-1 basic model) — segment events, inversions, transpositions.

These run a genuinely different genome representation through the *unchanged* simulator,
rate model, sampler, gene-tree reconstruction and profile matrix.
"""

import numpy as np

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    SharedRates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.events import EventType


def _ordered(extension=0.5):
    return lambda ids: OrderedGenome(ids, extension=extension)


def _n_leaves(newick):
    return newick.count(",") + 1


def test_ordered_genome_runs_with_rearrangements():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)
    rates = SharedRates(duplication=0.1, transfer=0.1, loss=0.15, origination=0.4,
                         inversion=0.3, transposition=0.3)
    g = simulate_genomes(tree, rates, initial_families=15, seed=2, genome_factory=_ordered(0.5))
    kinds = {r.event for r in g.event_log}
    assert EventType.INVERSION in kinds and EventType.TRANSPOSITION in kinds
    for leaf in g.leaf_genomes.values():
        assert isinstance(leaf, OrderedGenome)
        assert leaf.size() == len(leaf.chromosome)
    assert g.profiles.matrix.shape[1] == 10


def test_ordered_gene_trees_reconstruct():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=3)
    rates = SharedRates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.5,
                         inversion=0.3, transposition=0.3)
    g = simulate_genomes(tree, rates, initial_families=12, seed=3, genome_factory=_ordered(0.6))
    fam_row = {f: i for i, f in enumerate(g.profiles.families)}
    for fam, (complete, extant) in g.gene_trees().items():
        expected = int(g.profiles.matrix[fam_row[fam]].sum()) if fam in fam_row else 0
        if extant is None:
            assert expected == 0
        else:
            assert _n_leaves(extant) == expected  # rearrangements don't corrupt genealogy


def test_inversion_conserves_content_and_flips_orientation():
    tree = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=6, age=2.0, seed=1)
    # only rearrangements: no D/T/L/O beyond the seeds
    rates = SharedRates(inversion=1.0, transposition=0.5, origination=0.0)
    g = simulate_genomes(tree, rates, initial_families=8, seed=2, genome_factory=_ordered(0.7))
    for leaf in g.leaf_genomes.values():  # content unchanged (I/P only reorder/flip)
        assert leaf.size() == 8
        assert all(leaf.copy_number(str(i)) == 1 for i in range(1, 9))
    assert any(r.event is EventType.INVERSION for r in g.event_log)
    flipped = any(gene.orientation == -1
                  for leaf in g.leaf_genomes.values() for gene in leaf.chromosome)
    assert flipped  # inversions actually reversed some strands


def test_ordered_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=5)
    rates = SharedRates(duplication=0.15, transfer=0.1, loss=0.2, origination=0.4,
                         inversion=0.2, transposition=0.2)
    a = simulate_genomes(tree, rates, initial_families=8, seed=6, genome_factory=_ordered(0.5))
    b = simulate_genomes(tree, rates, initial_families=8, seed=6, genome_factory=_ordered(0.5))
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)
    assert len(a.event_log) == len(b.event_log)


def test_unordered_genome_ignores_rearrangement_rates():
    # SharedRates emits I/P weights, but UnorderedGenome doesn't support them -> filtered
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    rates = SharedRates(duplication=0.1, loss=0.15, origination=0.4,
                         inversion=1.0, transposition=1.0)
    g = simulate_genomes(tree, rates, initial_families=8, seed=2)  # default UnorderedGenome
    assert not any(r.event in (EventType.INVERSION, EventType.TRANSPOSITION)
                   for r in g.event_log)
