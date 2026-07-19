"""OrderedGenome (ZOMBI-1 basic model) — segment events, inversions, transpositions.

These run a genuinely different genome representation through the *unchanged* simulator,
rate model, sampler, gene-tree reconstruction and profile matrix.
"""

import numpy as np

from zombi2 import (
    BirthDeath,
    OrderedGenome,
    Rates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.genomes.events import EventType


def _ordered(extension=0.5):
    return lambda ids: OrderedGenome(ids, extension=extension)


def _n_leaves(newick):
    return newick.count(",") + 1


def test_ordered_genome_runs_with_rearrangements():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)
    rates = Rates(duplication=0.1, transfer=0.1, loss=0.15, origination=0.4,
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
    rates = Rates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.5,
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
    rates = Rates(inversion=1.0, transposition=0.5, origination=0.0)
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
    rates = Rates(duplication=0.15, transfer=0.1, loss=0.2, origination=0.4,
                         inversion=0.2, transposition=0.2)
    a = simulate_genomes(tree, rates, initial_families=8, seed=6, genome_factory=_ordered(0.5))
    b = simulate_genomes(tree, rates, initial_families=8, seed=6, genome_factory=_ordered(0.5))
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)
    assert len(a.event_log) == len(b.event_log)


def test_unordered_genome_ignores_rearrangement_rates():
    # Rates emits I/P weights, but UnorderedGenome doesn't support them -> filtered
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    rates = Rates(duplication=0.1, loss=0.15, origination=0.4,
                         inversion=1.0, transposition=1.0)
    g = simulate_genomes(tree, rates, initial_families=8, seed=2)  # default UnorderedGenome
    assert not any(r.event in (EventType.INVERSION, EventType.TRANSPOSITION)
                   for r in g.event_log)


def _total_branch_length(tree):
    return sum(n.branch_length() for n in tree.nodes() if n.parent is not None)


def _count_inversions(g):
    return sum(1 for r in g.event_log if r.event is EventType.INVERSION)


def test_inversion_count_matches_poisson_mean():
    # A FIXED species tree: every replicate walks the same branches, so the theoretical
    # intensity is a single constant, not a per-run random quantity.
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=12, age=3.0, seed=17)
    total_bl = _total_branch_length(tree)

    inversion_rate = 0.4
    initial_families = 10

    # Only rearrangements fire (no D/T/L/O). Inversions/transpositions conserve genome
    # CONTENT, so the genome size stays exactly `initial_families` on every branch. The
    # per-branch inversion hazard is `inversion_rate * size`, so integrated over the whole
    # tree the number of INVERSION events is Poisson with mean:
    #     lambda = inversion_rate * initial_families * total_branch_length
    rates = Rates(inversion=inversion_rate, transposition=0.2, origination=0.0)

    expected_mean = inversion_rate * initial_families * total_bl

    reps = 400
    counts = []
    for s in range(reps):
        g = simulate_genomes(tree, rates, initial_families=initial_families,
                             seed=1000 + s, genome_factory=_ordered(0.5))
        # content really is conserved -> size constant -> the oracle's `n` is exact
        for leaf in g.leaf_genomes.values():
            assert leaf.size() == initial_families
        counts.append(_count_inversions(g))

    observed_mean = np.asarray(counts).mean()

    # The Poisson-mean estimator over `reps` runs has s.e. = sqrt(lambda / reps); assert
    # within 4 sigma (deterministic seeds make this a fixed, several-sigma-margin check).
    sigma = np.sqrt(expected_mean / reps)
    assert abs(observed_mean - expected_mean) < 4 * sigma, (
        f"observed {observed_mean:.3f} vs expected {expected_mean:.3f} "
        f"({abs(observed_mean - expected_mean) / sigma:.2f} sigma)")
