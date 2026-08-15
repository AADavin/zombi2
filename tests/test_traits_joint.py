"""Two traits reading each other on one tree — the trait level joined to itself (design note §6).

The pair is a Markov chain on the product of the two state spaces, so this is **exact**: no
approximation and no thinning, just the branch walk a single trait already takes, over more states.
"""

import collections

import pytest

from zombi2 import species, traits
from zombi2.params import PerLineage
from zombi2.traits import TraitsResult


def _tree(n=40, seed=4):
    return species.simulate_species_tree(birth=1.0, n_extant=n, seed=seed).complete_tree


def _pair(cave_makes_large=6.0, large_makes_cave=8.0, seed=1, tree=None):
    """Big animals go underground, and underground makes you big."""
    return traits.simulate_traits(
        tree if tree is not None else _tree(),
        [traits.discrete(name="habitat", states=["surface", "cave"], start="surface",
                         switch={"surface->cave": PerLineage(0.05).scaled_by(
                                     "traits:size", {"small": 1.0, "large": large_makes_cave}),
                                 "cave->surface": 0.1}),
         traits.discrete(name="size", states=["small", "large"], start="small",
                         switch={"small->large": PerLineage(0.05).scaled_by(
                                     "traits:habitat", {"surface": 1.0, "cave": cave_makes_large}),
                                 "large->small": 0.1})],
        joint=True, seed=seed)


# --- what comes back ------------------------------------------------------------------------------

def test_one_result_per_trait():
    tree = _tree()
    r = _pair(tree=tree)
    assert sorted(r) == ["habitat", "size"]
    for name, res in r.items():
        assert isinstance(res, TraitsResult) and res.kind == "discrete"
        assert set(res.node_values) == set(tree.nodes)


def test_each_result_is_an_ordinary_trait_result():
    """Each one is exactly what `simulate_discrete` returns, so everything downstream reads it —
    the stochastic map covers each branch, and the log opens with the `initial` row a driver needs."""
    tree = _tree()
    r = _pair(tree=tree)
    for res in r.values():
        first = res.events[0]
        assert (first.kind, first.lineage, first.from_state) == ("initial", tree.root, None)
        for i, node in tree.nodes.items():
            assert sum(d for _s, d in res.history[i]) == pytest.approx(
                node.end_time - node.birth_time)


def test_it_is_deterministic():
    a, b = _pair(seed=7), _pair(seed=7)
    for name in a:
        assert [(c.time, c.kind, c.to_state) for c in a[name].events] == \
               [(c.time, c.kind, c.to_state) for c in b[name].events]


# --- the loop actually couples them ---------------------------------------------------------------

def test_the_loop_associates_the_two_traits():
    """The claim, measured as a lift rather than as raw agreement.

    Both traits drift the same way here — cave and large are each favoured once the other arrives —
    so counting tips where the two "agree" mostly counts that drift. What the loop actually claims is
    conditional: a lineage in the cave is likelier to be large **than one on the surface is**. That is
    what this measures, over ten seeds, against the same run with both factors flattened to 1.0.
    """
    tree = _tree()
    tips = [tree.labels()[i] for i in tree.extant_leaves()]

    def lift(run):
        pairs = collections.Counter((run["habitat"].values[n], run["size"].values[n]) for n in tips)
        cave = pairs[("cave", "large")] + pairs[("cave", "small")]
        surface = pairs[("surface", "large")] + pairs[("surface", "small")]
        if not cave or not surface:
            return None                       # a run that reached one habitat says nothing here
        return pairs[("cave", "large")] / cave - pairs[("surface", "large")] / surface

    def mean_lift(**kw):
        vals = [lift(_pair(tree=tree, seed=s, **kw)) for s in range(1, 11)]
        got = [v for v in vals if v is not None]
        assert len(got) >= 6, "too few runs reached both habitats to measure anything"
        return sum(got) / len(got)

    coupled = mean_lift()
    flat = mean_lift(cave_makes_large=1.0, large_makes_cave=1.0)
    assert coupled > 0.2, f"the loop left no association: {coupled:+.3f}"
    assert coupled > flat + 0.2, f"coupled {coupled:+.3f} against flat {flat:+.3f}"


def test_at_speciation_hops_the_trait_that_carries_it():
    """A hop at the split belongs to one trait, not to the pair: the trait that was given
    `at_speciation` gets `on_speciation` rows, and the other one gets none."""
    tree = _tree(30)
    r = traits.simulate_traits(
        tree,
        [traits.discrete(name="habitat", states=["surface", "cave"], start="surface",
                         at_speciation=0.6,
                         switch={"surface->cave": PerLineage(0.05).scaled_by(
                                     "traits:size", {"small": 1.0, "large": 8.0}),
                                 "cave->surface": 0.1}),
         traits.discrete(name="size", states=["small", "large"], start="small",
                         switch=PerLineage(0.05).scaled_by(
                             "traits:habitat", {"surface": 1.0, "cave": 6.0}))],
        joint=True, seed=3)

    hops = [c for c in r["habitat"].events if c.kind == "on_speciation"]
    assert hops, "at_speciation=0.6 over a 30-tip tree left no hop at all"
    assert all(c.from_state != c.to_state for c in hops)
    assert [c for c in r["size"].events if c.kind == "on_speciation"] == []
    assert all(c.time == tree.nodes[c.lineage].birth_time for c in hops)


# --- what is refused ------------------------------------------------------------------------------

def test_one_trait_is_not_this_function():
    with pytest.raises(ValueError, match="simulate_discrete"):
        traits.simulate_traits(_tree(10), [traits.discrete(name="a", states=["x", "y"], switch=0.1)],
                               seed=1)


def test_every_trait_needs_a_name():
    with pytest.raises(ValueError, match="needs a name"):
        traits.simulate_traits(_tree(10), [traits.discrete(states=["x", "y"], switch=0.1),
                                           traits.discrete(name="b", states=["x", "y"], switch=0.1)],
                               seed=1)


def test_names_must_be_unique():
    with pytest.raises(ValueError, match="unique"):
        traits.simulate_traits(_tree(10), [traits.discrete(name="a", states=["x", "y"], switch=0.1),
                                           traits.discrete(name="a", states=["x", "y"], switch=0.1)],
                               seed=1)


def test_reading_another_trait_needs_joint_true():
    with pytest.raises(ValueError, match="joint=True"):
        traits.simulate_traits(
            _tree(10),
            [traits.discrete(name="a", states=["x", "y"],
                             switch=PerLineage(0.1).scaled_by("traits:b", {"x": 2.0, "y": 1.0})),
             traits.discrete(name="b", states=["x", "y"], switch=0.1)], seed=1)


def test_joint_true_needs_something_to_read():
    with pytest.raises(ValueError, match="none reads another"):
        traits.simulate_traits(_tree(10), [traits.discrete(name="a", states=["x", "y"], switch=0.1),
                                           traits.discrete(name="b", states=["x", "y"], switch=0.1)],
                               joint=True, seed=1)


def test_a_driver_naming_no_trait_in_the_run_is_refused():
    with pytest.raises(ValueError, match="not a trait in this run"):
        traits.simulate_traits(
            _tree(10),
            [traits.discrete(name="a", states=["x", "y"],
                             switch=PerLineage(0.1).scaled_by("traits:nope", {"x": 2.0})),
             traits.discrete(name="b", states=["x", "y"], switch=0.1)], joint=True, seed=1)


def test_a_continuous_trait_in_a_cycle_says_it_is_not_built():
    with pytest.raises(TypeError, match="continuous trait in a cycle"):
        traits.simulate_traits(_tree(10), [traits.discrete(name="a", states=["x", "y"], switch=0.1),
                                           "not a spec"], joint=True, seed=1)
