"""Tests for the extent — how much a segmental event takes once it has started (SPEC §6).

An extent is the second axis of a segmental event: the rate says how often one *starts*, the extent
how much it takes. These cover the written form (a bare number is the **mean**), the one place it
parts company with :func:`as_distribution`, and each resolution's gate on what it can honour.
"""

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import simulate_genomes_nucleotide, simulate_genomes_ordered
from zombi2.rates.distributions import Fixed, Geometric, as_distribution
from zombi2.rates.extent import Extent, as_extent


@pytest.fixture
def tree():
    return species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)


# --- the written form: a bare number is the mean -------------------------------------------------

def test_a_bare_number_is_the_mean_not_a_fixed_size():
    """SPEC §6: ``500`` means runs of *about* 500, not exactly 500. This is the whole reason
    ``as_extent`` exists rather than reusing ``as_distribution``, which reads a bare number as
    ``Fixed`` — the right call for a sampled rate, the wrong one for a segment size."""
    d = as_extent(5)
    assert isinstance(d, Extent)
    assert isinstance(d.base, Geometric) and d.base.mean() == 5.0
    assert isinstance(as_distribution(5), Fixed)          # the sibling reading, deliberately unchanged

    rng = np.random.default_rng(0)
    draws = [d.sample(rng) for _ in range(4000)]
    assert len(set(draws)) > 1                            # it varies — not degenerate
    assert 4.5 < float(np.mean(draws)) < 5.5              # and it varies *around the number given*


def test_none_is_a_single_unit():
    d = as_extent(None)
    assert isinstance(d.base, Geometric) and d.base.mean() == 1.0
    rng = np.random.default_rng(0)
    assert {d.sample(rng) for _ in range(200)} == {1.0}    # Geometric(mean=1) is degenerate at 1


def test_a_distribution_passes_through_untouched():
    """An explicit distribution is the escape hatch for an exact size: ``Fixed(3)`` still means
    exactly three, so the bare-number reading takes nothing away."""
    f = Fixed(3)
    assert as_extent(f).base is f
    g = Geometric(mean=7)
    assert as_extent(g).base is g


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


# --- extent × modifiers (SPEC §6): how much, not how often ---------------------------------------

def test_a_modifier_scales_the_size():
    """``500 * OnTime(...)`` is an extent, not a rate: the factor scales the size drawn."""
    from zombi2.rates import modifiers as mod
    e = as_extent(400 * mod.OnTime({0: 1.0, 5.0: 3.0}))
    assert isinstance(e, Extent) and e.has_modifiers
    assert e.mean(time=0.0) == pytest.approx(400.0)
    assert e.mean(time=6.0) == pytest.approx(1200.0)     # past the step, three times as long

    rng = np.random.default_rng(0)
    early = [e.sample(rng, time=0.0) for _ in range(2000)]
    late = [e.sample(rng, time=6.0) for _ in range(2000)]
    assert float(np.mean(late)) > 2.5 * float(np.mean(early))


def test_an_extent_takes_no_scope():
    """An extent is already an absolute size, so there is no 'per what?' to answer (SPEC §6)."""
    from zombi2.rates import modifiers as mod
    from zombi2.rates.scope import PerLineage
    with pytest.raises(ValueError, match="an extent takes no scope"):
        as_extent(PerLineage(500) * mod.OnTime({0: 1.0}))


def test_a_plain_extent_is_not_driven():
    """The common case reads no context at all, which is what lets an engine skip building one."""
    assert as_extent(500).has_modifiers is False
    assert as_extent(None).has_modifiers is False


def test_mean_needs_a_geometric_base():
    e = as_extent(Fixed(50))
    with pytest.raises(ValueError, match="no mean to scale"):
        e.mean()


# --- the engines ----------------------------------------------------------------------------------

def test_nucleotide_extent_can_be_driven_by_a_trait(tree):
    """The genome-reduction statement people actually mean: host-restricted lineages delete in
    **bigger chunks**, which is a different model from deleting more often. The rate is left
    undriven here, so only the sizes may differ."""
    from zombi2 import traits
    from zombi2.rates import modifiers as mod
    from zombi2.rates.driver import driver_from_result

    # the trait and the genome take *different* seeds, as two levels of one study should: each level
    # draws its own stream, so a shared number no longer makes them the same draws (see zombi2.rng)
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=1)
    traj = driver_from_result(habitat)
    res = simulate_genomes_nucleotide(
        tree, root_length=20000, genes=4, gene_length=300, loss=0.6,
        loss_extent=150 * mod.DrivenBy(habitat, {"host": 6.0, "free": 1.0}), seed=5)

    sizes = {"host": [], "free": []}
    for e in res.events:
        if type(e).__name__ == "Loss":
            sizes[traj.value(e.lineage, e.time)].append(sum(hi - lo for (_s, _c, lo, hi) in e.lost))
    assert sizes["host"] and sizes["free"], "both states should have deleted something"
    assert np.mean(sizes["host"]) > np.mean(sizes["free"])


def test_ordered_extent_takes_a_skyline(tree):
    """``OnTime`` on an extent works at the ordered resolution, where the unit is genes."""
    from zombi2.rates import modifiers as mod
    g = simulate_genomes_ordered(
        tree, duplication=0.5, loss=0.2, origination=0.4, inversion=1.0,
        inversion_extent=1 * mod.OnTime({0: 1.0, 0.5: 6.0}),
        initial_families=15, chromosomes=1, seed=5)
    early = [r.length for r in g.rearrangements if r.time < 0.5]
    late = [r.length for r in g.rearrangements if r.time >= 0.5]
    assert early and late
    assert max(late) > max(early), "runs should get longer after the step"


def test_ordered_extent_can_be_driven_by_a_trait(tree):
    """The same statement as the nucleotide test above, in genes rather than base pairs: a
    host-restricted lineage duplicates in **longer runs**, not more often. The duplication rate is
    left undriven, so only the sizes may differ."""
    from zombi2 import traits
    from zombi2.rates import modifiers as mod
    from zombi2.rates.driver import driver_from_result

    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=2)
    traj = driver_from_result(habitat)
    g = simulate_genomes_ordered(
        tree, duplication=0.5, chromosomes=1, initial_families=30,
        duplication_extent=2 * mod.DrivenBy(habitat, {"host": 6.0, "free": 1.0}), seed=2)

    sizes = {"host": [], "free": []}
    for p in g.event_positions:
        if p.kind == "duplication":
            sizes[traj.value(p.lineage, p.time)].append(p.length)
    assert sizes["host"] and sizes["free"], "both states should have duplicated something"
    assert np.mean(sizes["host"]) > np.mean(sizes["free"])


def test_nucleotide_refuses_an_unwired_extent_modifier(tree):
    from zombi2.rates import modifiers as mod
    with pytest.raises(ValueError, match="an extent takes the same modifiers a rate does"):
        simulate_genomes_nucleotide(tree, root_length=2000, loss=0.5,
                                    loss_extent=200 * mod.Drawn(per='family', spread=0.5), seed=1)


def test_an_extent_is_read_in_the_same_context_as_a_rate():
    """One gate admits a modifier onto a rate and onto an extent, so both have to be handed the
    context `Modifier.implemented_for` promises. The extent sites passed `{"time": t}` alone, so a
    modifier written the documented way (`**_`, with defaults) silently read zeros there while
    reading real counts on a rate, and one with a required keyword died mid-run."""
    from zombi2 import genomes, species
    from zombi2.rates.modifiers import Modifier

    seen: dict[str, list[tuple]] = {"rate": [], "extent": []}

    class Records(Modifier):
        implemented_for = ("genomes.ordered", "genomes.nucleotide")

        def __init__(self, where):
            self.where = where

        def factor(self, *, time=0.0, copies=0, lineages=0, chromosomes=0, **_):
            seen[self.where].append((copies, lineages, chromosomes, time))
            return 1.0

    tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1).complete_tree
    genomes.simulate_genomes_ordered(tree, initial_families=20, loss=0.4 * Records("rate"),
                                     loss_extent=3 * Records("extent"), seed=2)

    assert seen["rate"] and seen["extent"]
    # the counts an extent is given are the run's, not zeros
    assert any(c or ln or ch for c, ln, ch, _ in seen["extent"])


def test_an_extent_is_read_at_the_instant_the_event_fires():
    """An extent is sampled when the event fires, not when the Gillespie stretch it fell in began.
    An extent's own breakpoints are deliberately kept out of the horizon, so a schedule's breakpoint
    routinely falls *inside* a stretch, and reading the stretch's start time would size the event on
    the wrong side of it. Threading the rate loop's whole context to the extent is what made that
    possible: the context is snapshotted before `t` advances to the firing instant."""
    from zombi2 import genomes, species
    from zombi2.rates import modifiers as mod

    tree = species.simulate_species_tree(birth=0.4, death=0.05, total_time=12.0,
                                         seed=3).complete_tree
    times: list[float] = []

    class RecordsTime(mod.Modifier):
        implemented_for = ("genomes.ordered",)

        def factor(self, *, time=0.0, **_):
            times.append(time)
            return 1.0

    run = genomes.simulate_genomes_ordered(
        tree, inversion=0.3, inversion_extent=8 * RecordsTime(),
        initial_families=40, chromosomes=1, seed=11)

    fired = {round(e.time, 9) for e in run.rearrangements}
    assert times, "the extent's modifier was never read"
    assert all(round(t, 9) in fired for t in times), \
        "an extent was read at a time no event fired at — the stretch's start, not the event's"


def test_a_scheduled_extent_changes_size_on_the_right_side_of_its_breakpoint():
    """The consequence of the above, on a shipped modifier. A schedule that collapses the extent
    after time 2 must give small events after time 2 and full-sized ones before it — which is only
    true if the schedule is read at the firing instant."""
    from zombi2 import genomes, species
    from zombi2.rates import modifiers as mod

    tree = species.simulate_species_tree(birth=0.4, death=0.05, total_time=12.0,
                                         seed=3).complete_tree
    run = genomes.simulate_genomes_ordered(
        tree, inversion=0.3, inversion_extent=20 * mod.OnTime({0: 1.0, 2.0: 0.05}),
        initial_families=40, chromosomes=1, seed=11)

    early = [e.length for e in run.rearrangements if e.time < 2.0]
    late = [e.length for e in run.rearrangements if e.time > 2.0]
    assert early and late
    assert sum(early) / len(early) > 3 * (sum(late) / len(late))


def test_the_nucleotide_engine_reads_an_extent_the_same_way():
    """The same fix landed at two sites, and only the ordered one was covered. A future edit could
    put the thin context back at the nucleotide site and nothing would say so."""
    from zombi2 import genomes, species
    from zombi2.rates.modifiers import Modifier

    seen: list[tuple] = []

    class Records(Modifier):
        implemented_for = ("genomes.nucleotide",)

        def factor(self, *, time=0.0, lineages=0, chromosomes=0, **_):
            seen.append((lineages, chromosomes, time))
            return 1.0

    tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1).complete_tree
    genomes.simulate_genomes_nucleotide(
        tree, genes=4, gene_length=200, root_length=3000,
        loss=0.6, loss_extent=100 * Records(), seed=2)

    assert seen, "the extent's modifier was never read"
    # the counts are the run's, not zeros — `copies` is excluded, being 0 by design here
    assert any(ln or ch for ln, ch, _ in seen)
    assert all(t > 0.0 for _, _, t in seen), "an extent was read at time 0, i.e. the stretch's start"
