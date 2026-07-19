"""Tests for the forward birth-death engine (species slice 1)."""

import pytest

from zombi2.species_tree import Event, simulate_species_tree


def test_yule_reaches_n_tips_with_no_extinction():
    r = simulate_species_tree(birth=1.0, death=0.0, n_tips=50, seed=1)
    assert r.n_extant == 50
    assert all(e.kind == "speciation" for e in r.events)  # Yule → no deaths
    assert all(n.fate == "extant" for n in r.complete.tips())


def test_birth_death_has_extinctions_and_survivors():
    r = simulate_species_tree(birth=1.0, death=0.4, n_tips=60, seed=7)
    assert r.n_extant == 60
    kinds = {e.kind for e in r.events}
    assert "extinction" in kinds  # death > 0 → some lineages die
    assert len(r.complete.extinct()) > 0


def test_age_stopping_makes_extant_tips_end_at_age():
    r = simulate_species_tree(birth=1.0, death=0.2, age=4.0, seed=3)
    for n in r.complete.extant():
        assert n.end_time == pytest.approx(4.0)


def test_deterministic_given_seed():
    a = simulate_species_tree(birth=1.0, death=0.3, n_tips=40, seed=42)
    b = simulate_species_tree(birth=1.0, death=0.3, n_tips=40, seed=42)
    assert len(a.complete.nodes) == len(b.complete.nodes)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]


def test_different_seeds_differ():
    a = simulate_species_tree(birth=1.0, death=0.3, n_tips=40, seed=1)
    b = simulate_species_tree(birth=1.0, death=0.3, n_tips=40, seed=2)
    assert len(a.complete.nodes) != len(b.complete.nodes) or \
        [e.time for e in a.events] != [e.time for e in b.events]


def test_tree_structure_invariants():
    r = simulate_species_tree(birth=1.0, death=0.2, n_tips=30, seed=5)
    for node in r.complete.nodes.values():
        if node.children is None:
            assert node.fate in ("extant", "extinct")
        else:
            c1, c2 = node.children
            assert c1 != c2
            assert r.complete.nodes[c1].parent == node.id
            assert r.complete.nodes[c2].parent == node.id
            assert node.fate == "speciation"
    # every speciation event lands as a 2-child internal node
    for e in r.events:
        if e.kind == "speciation":
            assert r.complete.nodes[e.node].children is not None


def test_events_are_time_ordered():
    r = simulate_species_tree(birth=1.0, death=0.3, n_tips=50, seed=9)
    times = [e.time for e in r.events]
    assert times == sorted(times)


def test_validation():
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0)                       # neither n_tips nor age
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, n_tips=10, age=5.0)   # both
    with pytest.raises(ValueError):
        simulate_species_tree(birth=-1.0, n_tips=10)           # negative rate
    with pytest.raises(TypeError):
        simulate_species_tree(birth="fast", n_tips=10)         # non-numeric rate
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, n_tips=0)             # non-positive n_tips
    with pytest.raises(ValueError):
        simulate_species_tree(birth=1.0, age=-2.0)             # non-positive age


def test_event_is_frozen_record():
    e = Event(1.0, "extinction", 3)
    with pytest.raises(Exception):
        e.time = 2.0  # type: ignore[misc]


def test_n_tips_conditions_on_survival():
    # seed 7's first attempt dies on the first event; conditioning restarts until 60 tips
    r = simulate_species_tree(birth=1.0, death=0.4, n_tips=60, seed=7)
    assert r.n_extant == 60


def test_raises_when_n_tips_unreachable():
    # death swamps birth → 200 tips essentially never happens → give up, don't hang
    with pytest.raises(RuntimeError):
        simulate_species_tree(birth=0.1, death=3.0, n_tips=200, seed=1)
