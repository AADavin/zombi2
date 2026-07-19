"""Tests for the forward birth-death engine (species)."""

import pytest

from zombi2 import modifiers as mod
from zombi2 import scope
from zombi2.species_tree import Event, simulate_species_tree


def test_yule_reaches_n_extant_with_no_extinction():
    r = simulate_species_tree(birth=1.0, death=0.0, n_extant=50, seed=1)
    assert r.n_extant == 50
    assert all(e.kind == "speciation" for e in r.events)  # Yule → no deaths
    assert all(n.fate == "extant" for n in r.complete_tree.leaves())


def test_birth_death_has_extinctions_and_survivors():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=60, seed=7)
    assert r.n_extant == 60
    assert "extinction" in {e.kind for e in r.events}
    assert len(r.complete_tree.extinct()) > 0


def test_age_stopping_makes_extant_lineages_end_at_age():
    r = simulate_species_tree(birth=1.0, death=0.2, age=4.0, seed=3)
    for n in r.complete_tree.extant():
        assert n.end_time == pytest.approx(4.0)


def test_deterministic_given_seed():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=42)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=42)
    assert len(a.complete_tree.nodes) == len(b.complete_tree.nodes)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]


def test_different_seeds_differ():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=1)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=40, seed=2)
    assert len(a.complete_tree.nodes) != len(b.complete_tree.nodes) or \
        [e.time for e in a.events] != [e.time for e in b.events]


def test_tree_structure_invariants():
    r = simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=5)
    for node in r.complete_tree.nodes.values():
        if node.children is None:
            assert node.fate in ("extant", "extinct")
        else:
            c1, c2 = node.children
            assert c1 != c2
            assert r.complete_tree.nodes[c1].parent == node.id
            assert r.complete_tree.nodes[c2].parent == node.id
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


def test_skyline_stops_births_after_a_zero_breakpoint():
    # birth 1.0 on [0, 2), then 0 → the interval-aware sampler must forbid births at/after t=2
    r = simulate_species_tree(birth=1.0 * mod.Time({0: 1.0, 2.0: 0.0}), death=0.0, age=10.0, seed=1)
    spec_times = [e.time for e in r.events if e.kind == "speciation"]
    assert spec_times                 # growth happened before the breakpoint
    assert max(spec_times) < 2.0      # and nothing after the rate dropped to 0


def test_skyline_is_deterministic():
    kw = dict(birth=1.0 * mod.Time({0: 2.0, 3: 0.2}), death=0.1, age=6.0, seed=4)
    a = simulate_species_tree(**kw)
    b = simulate_species_tree(**kw)
    assert [(e.time, e.kind) for e in a.events] == [(e.time, e.kind) for e in b.events]


# --- the extant tree + Newick output --------------------------------------

def test_extant_tree_prunes_to_survivors():
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, seed=3)
    ext = r.extant_tree
    assert len(ext.leaves()) == r.n_extant                       # exactly the survivors
    assert all(n.fate == "extant" for n in ext.leaves())         # no extinct tips
    assert all(n.fate != "extinct" for n in ext.nodes.values())  # no extinct nodes remain
    for n in ext.nodes.values():                                 # bifurcating
        assert n.children is None or len(n.children) == 2
    # every branch is now strictly positive — the present sits after the last split (n_extant fix)
    for n in ext.nodes.values():
        assert n.end_time - n.birth_time > 0


def test_yule_extant_equals_complete_leaves():
    r = simulate_species_tree(birth=1.0, death=0.0, n_extant=30, seed=1)
    assert len(r.extant_tree.leaves()) == len(r.complete_tree.leaves()) == 30


def test_newick_is_wellformed():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=2)
    nwk = r.extant_tree.to_newick()
    assert nwk.endswith(";")
    assert nwk.count("(") == nwk.count(")")   # balanced parens
    assert nwk.count(",") == 25 - 1           # a bifurcating tree of 25 tips has 24 joins


def test_write_produces_newick_files(tmp_path):
    r = simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=5)
    r.write(tmp_path)
    assert (tmp_path / "complete.nwk").read_text().strip().endswith(";")
    assert (tmp_path / "extant.nwk").read_text().strip().endswith(";")


def test_extant_tree_is_deterministic():
    a = simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=8)
    b = simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=8)
    assert a.extant_tree.to_newick() == b.extant_tree.to_newick()


def test_dead_tree_has_no_extant_tree():
    r = simulate_species_tree(birth=0.1, death=10.0, age=5.0, seed=1)
    assert r.n_extant == 0
    assert r.extant_tree is None


# --- Inherited (ClaDS): rates drift down the tree, picking is rate-weighted ---

def test_weighted_index_respects_weights():
    import numpy as np

    from zombi2.species_tree import _weighted_index
    rng = np.random.default_rng(0)
    weights = [1.0, 1.0, 8.0]           # index 2 carries 80% of the total rate
    counts = [0, 0, 0]
    for _ in range(20000):
        counts[_weighted_index(rng, weights, sum(weights))] += 1
    assert 0.77 < counts[2] / 20000 < 0.83   # ≈ 0.8
    assert counts[0] > 0 and counts[1] > 0    # the light lineages still get picked sometimes


def test_clads_is_deterministic_given_seed():
    kw = dict(birth=1.0 * mod.Inherited(spread=0.5), death=0.1, n_extant=40, seed=3)
    a = simulate_species_tree(**kw)
    b = simulate_species_tree(**kw)
    assert [(e.time, e.kind, e.node) for e in a.events] == [(e.time, e.kind, e.node) for e in b.events]


def test_inherited_zero_spread_reaches_target():
    # spread 0 → every step is ×1, so no drift; still a valid birth-death tree
    r = simulate_species_tree(birth=1.0 * mod.Inherited(spread=0.0), death=0.2, n_extant=40, seed=5)
    assert r.n_extant == 40


def test_death_can_drift_independently():
    # drift lives on death, not birth; birth and death are bent independently
    r = simulate_species_tree(birth=1.0, death=0.4 * mod.Inherited(spread=0.5), n_extant=50, seed=4)
    assert r.n_extant == 50
    assert len(r.complete_tree.extinct()) > 0


def test_clads_composes_with_diversity_cap():
    # ClaDS drift × diversity-dependence: the cap still bounds the tree
    r = simulate_species_tree(
        birth=1.0 * mod.Inherited(spread=0.4) * mod.Diversity(cap=25), death=0.0, age=100.0, seed=1)
    assert r.n_extant <= 25          # the cap is a hard ceiling even with drift
    assert r.n_extant >= 12          # and the tree grew toward it


def test_inherited_requires_per_lineage_scope():
    # per-lineage drift on a Global (tree-wide) budget is contradictory — reject it clearly
    with pytest.raises(ValueError, match="per lineage"):
        simulate_species_tree(birth=scope.Global(1.0) * mod.Inherited(spread=0.2), n_extant=10, seed=1)


def test_drifting_birth_with_non_drifting_global_death_is_allowed():
    # only the drifting rate must be per lineage; a Global death budget alongside it is fine
    r = simulate_species_tree(
        birth=1.0 * mod.Inherited(spread=0.3), death=scope.Global(0.2), n_extant=30, seed=2)
    assert r.n_extant == 30


def _colless(result):
    """Colless imbalance of the extant tree: Σ over internal nodes of |left tips − right tips|.
    Higher = more lopsided. A pure Yule tree is comparatively balanced; heritable rate drift
    concentrates tips in the fast clades, so it climbs."""
    tree = result.extant_tree
    size = {}
    for i in sorted(tree.nodes, reverse=True):        # children (higher ids) before parents
        nd = tree.nodes[i]
        size[i] = 1 if nd.children is None else sum(size[c] for c in nd.children)
    return sum(abs(size[nd.children[0]] - size[nd.children[1]])
               for nd in tree.nodes.values() if nd.children is not None)


def test_clads_is_more_imbalanced_than_yule():
    # the signature of heritable rate drift, at a fixed tip count so it is shape not size: fast
    # clades are inherited, so they hoard the tips and the tree is far more lopsided than Yule
    import statistics
    seeds = range(40)
    yule = [_colless(simulate_species_tree(birth=1.0, death=0.0, n_extant=64, seed=s)) for s in seeds]
    clads = [_colless(simulate_species_tree(birth=1.0 * mod.Inherited(spread=0.9), death=0.0, n_extant=64, seed=s))
             for s in seeds]
    assert statistics.mean(clads) > 1.5 * statistics.mean(yule)   # observed ≈ 2.7× (margin to spare)
