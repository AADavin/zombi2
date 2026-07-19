"""Tests for the forward birth-death engine (species)."""

import pytest

from zombi2 import modifiers as mod
from zombi2 import scope
from zombi2.species_tree import Event, simulate_species_tree


def test_yule_reaches_n_extant_with_no_extinction():
    r = simulate_species_tree(birth=1.0, death=0.0, n_extant=50, seed=1)
    assert r.n_extant == 50
    assert all(e.kind == "speciation" for e in r.events)  # Yule → no deaths
    assert all(n.fate == "extant" for n in r.complete.leaves())


def test_birth_death_has_extinctions_and_survivors():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=60, seed=7)
    assert r.n_extant == 60
    assert "extinction" in {e.kind for e in r.events}
    assert len(r.complete.extinct()) > 0


def test_age_stopping_makes_extant_lineages_end_at_age():
    r = simulate_species_tree(birth=1.0, death=0.2, age=4.0, seed=3)
    for n in r.complete.extant():
        assert n.end_time == pytest.approx(4.0)


def test_deterministic_given_seed():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=42)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=42)
    assert len(a.complete.nodes) == len(b.complete.nodes)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]


def test_different_seeds_differ():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=1)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=2)
    assert len(a.complete.nodes) != len(b.complete.nodes) or \
        [e.time for e in a.events] != [e.time for e in b.events]


def test_tree_structure_invariants():
    r = simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=5)
    for node in r.complete.nodes.values():
        if node.children is None:
            assert node.fate in ("extant", "extinct")
        else:
            c1, c2 = node.children
            assert c1 != c2
            assert r.complete.nodes[c1].parent == node.id
            assert r.complete.nodes[c2].parent == node.id
            assert node.fate == "speciation"


def test_events_are_time_ordered():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=50, seed=9)
    times = [e.time for e in r.events]
    assert times == sorted(times)


def test_n_extant_conditions_on_survival():
    # seed 7's first attempt dies on the first event; conditioning restarts until 60 living
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=60, seed=7)
    assert r.n_extant == 60


def test_raises_when_n_extant_unreachable():
    with pytest.raises(RuntimeError):
        simulate_species_tree(birth=0.1, death=3.0, n_extant=200, seed=1)


# --- the scope × modifiers grammar actually drives the engine -------------

def test_global_grows_slower_than_per_lineage():
    # per-lineage birth compounds (exponential); Global is a constant tree-wide budget (linear)
    per = simulate_species_tree(birth=1.0, death=0.0, age=4.0, seed=1)
    glob = simulate_species_tree(birth=scope.Global(2.0), death=0.0, age=4.0, seed=1)
    assert per.n_extant > glob.n_extant


def test_diversity_caps_growth():
    # Diversity(cap=20): the birth factor falls to 0 at 20 lineages, so the tree saturates
    r = simulate_species_tree(birth=1.0 * mod.Diversity(cap=20), death=0.0, age=100.0, seed=1)
    assert r.n_extant <= 20   # never exceeds the cap
    assert r.n_extant >= 15   # but it grew toward it


def test_no_scope_default_is_per_lineage():
    # a bare number and an explicit PerLineage wrapper give the same tree
    a = simulate_species_tree(birth=1.0, death=0.2, n_extant=40, seed=11)
    b = simulate_species_tree(birth=scope.PerLineage(1.0), death=scope.PerLineage(0.2), n_extant=40, seed=11)
    assert [(e.time, e.kind) for e in a.events] == [(e.time, e.kind) for e in b.events]


# --- validation -----------------------------------------------------------

def test_validation():
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0)                          # neither n_extant nor age
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, n_extant=10, age=5.0)    # both
    with pytest.raises(ValueError):
        simulate_species_tree(birth=-1.0, n_extant=10)            # negative rate (via scope)
    with pytest.raises(TypeError):
        simulate_species_tree(birth="fast", n_extant=10)          # non-numeric rate
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, n_extant=0)              # non-positive n_extant
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, age=-2.0)               # non-positive age


def test_event_is_frozen_record():
    e = Event(1.0, "extinction", 3)
    with pytest.raises(Exception):
        e.time = 2.0  # type: ignore[misc]
