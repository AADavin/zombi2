"""Tests for the unordered D/L/O gene-family core (zombi2.genomes_unordered)."""

import pytest

from zombi2 import modifiers as mod
from zombi2.species_tree import simulate_species_tree
from zombi2.genomes_unordered import GeneCopy, simulate_unordered


def _tree(seed=1, n_extant=12, death=0.3):
    return simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)


# --- the walk covers the whole complete tree -------------------------------

def test_genomes_on_every_node_including_extinct():
    sp = _tree(seed=3, death=0.5)
    g = simulate_unordered(sp, duplication=0.2, loss=0.2, origination=0.5, initial_families=4, seed=1)
    assert set(g.genomes) == set(sp.complete_tree.nodes)          # every node has a genome
    extinct = {n.id for n in sp.complete_tree.extinct()}
    assert extinct and extinct <= set(g.genomes)                 # extinct lineages included


def test_extant_genomes_are_exactly_the_tips():
    sp = _tree(seed=5)
    g = simulate_unordered(sp, duplication=0.1, loss=0.1, origination=0.4, initial_families=3, seed=1)
    assert set(g.extant_genomes) == {n.id for n in sp.complete_tree.extant()}
    assert len(g.extant_genomes) == sp.n_extant


def test_accepts_a_result_or_a_bare_tree():
    sp = _tree(seed=7)
    a = simulate_unordered(sp, origination=0.5, initial_families=3, seed=1)
    b = simulate_unordered(sp.complete_tree, origination=0.5, initial_families=3, seed=1)
    assert [(e.time, e.kind, e.copy) for e in a.events] == [(e.time, e.kind, e.copy) for e in b.events]


# --- determinism -----------------------------------------------------------

def test_deterministic_given_seed():
    sp = _tree(seed=2)
    kw = dict(duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=9)
    a, b = simulate_unordered(sp, **kw), simulate_unordered(sp, **kw)
    assert [(e.time, e.kind, e.lineage, e.copy, e.parent) for e in a.events] == \
           [(e.time, e.kind, e.lineage, e.copy, e.parent) for e in b.events]
    assert a.genomes == b.genomes


def test_different_seeds_differ():
    sp = _tree(seed=2)
    a = simulate_unordered(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=1)
    b = simulate_unordered(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=2)
    assert len(a.events) != len(b.events) or \
        [e.copy for e in a.events] != [e.copy for e in b.events]


# --- the three events behave --------------------------------------------------

def test_initial_families_seed_originations_at_the_crown():
    sp = _tree(seed=1)
    root = sp.complete_tree.root
    t0 = sp.complete_tree.nodes[root].birth_time
    g = simulate_unordered(sp, initial_families=6, seed=1)          # no D/L/O rates
    crown = [e for e in g.events if e.time == t0]
    assert len(crown) == 6 and all(e.kind == "origination" for e in crown)
    assert len(g.genomes[root]) == 6                               # the root carries all 6


def test_no_rates_means_pure_inheritance():
    # with every rate 0, no event fires after the crown seeding and every node's genome is the
    # root's, unchanged (the same GeneCopy objects threaded down) — inheritance in isolation
    sp = _tree(seed=4)
    g = simulate_unordered(sp, initial_families=3, seed=1)
    root_genome = g.genomes[sp.complete_tree.root]
    assert all(g.genomes[i] == root_genome for i in g.genomes)     # every node equals the root
    assert all(e.time == 0.0 for e in g.events)                    # only the crown originations


def test_origination_only_families_never_exceed_one_copy():
    sp = _tree(seed=6)
    g = simulate_unordered(sp, origination=0.8, initial_families=2, seed=1)  # no duplication
    for node_id in g.genomes:
        assert all(count == 1 for count in g.family_counts(node_id).values())


def test_duplication_grows_a_family():
    sp = _tree(seed=6)
    g = simulate_unordered(sp, duplication=0.8, initial_families=3, seed=1)  # no loss, no origination
    biggest = max((max(g.family_counts(i).values(), default=0)) for i in g.genomes)
    assert biggest > 1                                            # some family reached >1 copy
    # duplication never introduces a new family (only origination does)
    assert {e.kind for e in g.events} <= {"origination", "duplication"}


def test_loss_can_shrink_and_empty_a_genome():
    # high loss, no origination/duplication except the seeded families -> some lineage loses all
    sp = _tree(seed=6)
    g = simulate_unordered(sp, loss=2.0, initial_families=4, seed=1)
    sizes = [len(g.genomes[i]) for i in g.genomes]
    assert min(sizes) < 4                                          # at least one node shrank
    assert any(e.kind == "loss" for e in g.events)


def test_duplication_parent_survives_and_is_a_real_copy():
    sp = _tree(seed=8)
    g = simulate_unordered(sp, duplication=0.6, initial_families=3, seed=1)
    born = {e.copy for e in g.events if e.kind in ("origination", "duplication")}
    for e in g.events:
        if e.kind == "duplication":
            assert e.parent in born                               # parent is a real, earlier copy
            assert e.family == next(x.family for x in _copies_born(g) if x.id == e.parent)  # same family


def _copies_born(result):
    return [GeneCopy(e.copy, e.family) for e in result.events if e.kind in ("origination", "duplication")]


def test_every_born_copy_id_is_unique():
    sp = _tree(seed=8, death=0.5)
    g = simulate_unordered(sp, duplication=0.4, loss=0.3, origination=0.6, initial_families=5, seed=2)
    born = [e.copy for e in g.events if e.kind in ("origination", "duplication")]
    assert len(born) == len(set(born))


def test_family_counts_matches_the_genome():
    sp = _tree(seed=9)
    g = simulate_unordered(sp, duplication=0.3, loss=0.2, origination=0.5, initial_families=4, seed=1)
    for node_id in g.genomes:
        assert sum(g.family_counts(node_id).values()) == len(g.genomes[node_id])


# --- empty / validation ----------------------------------------------------

def test_empty_run_has_no_events_or_content():
    sp = _tree(seed=1)
    g = simulate_unordered(sp, seed=1)                            # no families, no rates
    assert g.events == []
    assert all(genome == () for genome in g.genomes.values())


def test_validation():
    sp = _tree(seed=1)
    with pytest.raises(ValueError):
        simulate_unordered(sp, initial_families=-1, seed=1)
    with pytest.raises(ValueError):
        simulate_unordered(sp, initial_families=True, seed=1)     # bool is not a valid count
    with pytest.raises(ValueError):
        simulate_unordered(sp, origination=-1.0, initial_families=1, seed=1)   # negative rate (via scope)


# --- modifiers: Time (skyline) is wired; the rest are rejected, not silently dropped ---

def test_time_skyline_modifier_is_supported():
    # Time reads only `time`, which the walk supplies, so a skyline origination works: the rate
    # drops to 0 at t=1.5, so no family originates after it
    sp = simulate_species_tree(birth=1.0, death=0.2, total_time=4.0, seed=3)
    r = simulate_unordered(sp, origination=1.0 * mod.Time({0: 1.0, 1.5: 0.0}), seed=1)
    orig_times = [e.time for e in r.events if e.kind == "origination"]
    assert orig_times and max(orig_times) < 1.5


def test_unsupported_modifiers_are_rejected_not_silently_dropped():
    sp = _tree(seed=1)
    # clade drift would need per-lineage threading the walk doesn't do → reject, don't no-op
    with pytest.raises(ValueError, match="does not support"):
        simulate_unordered(sp, duplication=0.5 * mod.Inherited(spread=0.8), initial_families=3, seed=1)
    # Diversity reads a `diversity` context the genome walk doesn't supply → reject, don't crash raw
    with pytest.raises(ValueError, match="does not support"):
        simulate_unordered(sp, loss=0.25 * mod.Diversity(cap=100), initial_families=3, seed=1)
