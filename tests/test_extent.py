"""Tests for the extent — how much a segmental event takes once it has started (SPEC §6).

An extent is the second axis of a segmental event: the rate says how often one *starts*, the extent
how much it takes. These cover the written form (a bare number is the **mean**), the one place it
parts company with :func:`as_distribution`, and each resolution's gate on what it can honour.
"""

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import simulate_genomes_nucleotide, simulate_genomes_ordered
from zombi2.rates.distributions import Fixed, Geometric, as_distribution, as_extent


@pytest.fixture
def tree():
    return species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)


# --- the written form: a bare number is the mean -------------------------------------------------

def test_a_bare_number_is_the_mean_not_a_fixed_size():
    """SPEC §6: ``500`` means runs of *about* 500, not exactly 500. This is the whole reason
    ``as_extent`` exists rather than reusing ``as_distribution``, which reads a bare number as
    ``Fixed`` — the right call for a sampled rate, the wrong one for a segment size."""
    d = as_extent(5)
    assert isinstance(d, Geometric) and d.mean == 5.0
    assert isinstance(as_distribution(5), Fixed)          # the sibling reading, deliberately unchanged

    rng = np.random.default_rng(0)
    draws = [d.sample(rng) for _ in range(4000)]
    assert len(set(draws)) > 1                            # it varies — not degenerate
    assert 4.5 < float(np.mean(draws)) < 5.5              # and it varies *around the number given*


def test_none_is_a_single_unit():
    d = as_extent(None)
    assert isinstance(d, Geometric) and d.mean == 1.0
    rng = np.random.default_rng(0)
    assert {d.sample(rng) for _ in range(200)} == {1.0}    # Geometric(mean=1) is degenerate at 1


def test_a_distribution_passes_through_untouched():
    """An explicit distribution is the escape hatch for an exact size: ``Fixed(3)`` still means
    exactly three, so the bare-number reading takes nothing away."""
    f = Fixed(3)
    assert as_extent(f) is f
    g = Geometric(mean=7)
    assert as_extent(g) is g


def test_a_callable_still_works():
    d = as_extent(lambda rng: 4.0)
    assert d.sample(np.random.default_rng(0)) == 4.0


# --- ordered: extents are in genes ----------------------------------------------------------------

def test_ordered_bare_number_gives_varying_segments(tree):
    """The behaviour the bare-number reading buys: segments spread around the mean instead of all
    being the same size."""
    g = simulate_genomes_ordered(tree, duplication=0.4, loss=0.2, origination=0.3, inversion=0.8,
                                 inversion_extent=4, initial_families=10, seed=3)
    lengths = {i.length for i in g.rearrangements}
    assert len(lengths) > 1, "a mean of 4 should produce a spread of run lengths, not one size"


def test_ordered_fixed_gives_exactly_that_many(tree):
    g = simulate_genomes_ordered(tree, duplication=0.3, origination=0.3, inversion=0.8,
                                 inversion_extent=Fixed(2), initial_families=12, seed=3)
    # a run is clamped by what the chromosome can carry, so 2 is the ceiling, never exceeded
    assert g.rearrangements and {i.length for i in g.rearrangements} <= {1, 2}


def test_ordered_default_extent_is_one_gene(tree):
    g = simulate_genomes_ordered(tree, duplication=0.3, origination=0.3, inversion=0.8,
                                 initial_families=10, seed=3)
    assert {i.length for i in g.rearrangements} == {1}


# --- nucleotide: extents are in base pairs, and only a geometric shape is wired -------------------

def test_nucleotide_refuses_a_shape_it_cannot_honour(tree):
    """It draws each arc's far end straight from the legal breakpoints, so a non-geometric shape
    would have to be re-weighted over that set rather than drawn. SPEC's rule: raise, never
    silently approximate."""
    with pytest.raises(ValueError, match="geometric extent only"):
        simulate_genomes_nucleotide(tree, root_length=2000, inversion=1.0,
                                    inversion_extent=Fixed(300), seed=1)


def test_nucleotide_rejects_a_sub_base_pair_extent(tree):
    with pytest.raises(ValueError, match=r"loss_extent must be >= 1 bp"):
        simulate_genomes_nucleotide(tree, root_length=2000, loss=1.0, loss_extent=0, seed=1)


def test_nucleotide_accepts_an_explicit_geometric_identically(tree):
    """A number and the geometric it stands for are the same run, byte for byte — the bare form is
    shorthand, not a different path."""
    kw = dict(root_length=3000, genes=2, gene_length=200, inversion=1.5, seed=11)
    a = simulate_genomes_nucleotide(tree, inversion_extent=400, **kw)
    b = simulate_genomes_nucleotide(tree, inversion_extent=Geometric(mean=400), **kw)
    assert [(r.time, r.lineage, r.start, r.length) for r in a.rearrangements] == \
           [(r.time, r.lineage, r.start, r.length) for r in b.rearrangements]
