"""Per-lineage weights kept between events (issue #487).

A rate that differs between lineages (a driver, a per-lineage draw) is summed over the living
lineages, and the affected one is drawn with the same weights. The family engines used to rebuild
every lineage's weight at every event, so a driven run slowed in proportion to its tree. They keep
the weights now and rebuild one lineage when its genome or its driver changes. What has to hold:

- a kept weight is the float a fresh one would be, at every event (`_CHECK_KEPT`);
- with the weights added left to right, a driven run is byte for byte the run 0.57.1 made, on both
  engines, so keeping them changed nothing but the cost;
- a run that reads no driver is unchanged as shipped, because it adds left to right (`SummedWeights`);
- a run that reads a driver adds in pairs (`WeightedIndex`): the same process, and an event that
  costs no more on a larger tree.
"""

import hashlib
import math
import sys

import numpy as np
import pytest

from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Curve, Gamma, LogNormal, PerCopy, PerLineage
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete


def _up(y):        # module-level rather than lambdas, so the per-family engine could ship them
    return math.exp(y)


def _down(y):
    return math.exp(-y)


@pytest.fixture(scope="module")
def world():
    tree = simulate_species_tree(birth=1.0, death=0.4, n_extant=25, seed=4).complete_tree
    return (tree, simulate_discrete(tree, states=["cave", "surface"], switch=1.5, seed=1),
            simulate_continuous(tree, rate=0.3, seed=2))


def _driven(hab, bm):
    """Runs that read a driver, one per way a rate can end up read lineage by lineage."""
    tab, back = {"cave": 3.0, "surface": 0.5}, {"cave": 0.5, "surface": 2.0}
    draw = Gamma(shape=2.0, scale=0.5)
    return {
        "a trait on duplication and loss": dict(
            duplication=PerCopy(0.3).scaled_by(hab, tab), loss=PerCopy(0.3).scaled_by(hab, back),
            transfer=0.2),
        "a continuous trait, stepped": dict(
            duplication=PerCopy(0.2).scaled_by(bm, Curve(_up), step=0.05),
            loss=PerCopy(0.3).scaled_by(bm, Curve(_down), step=0.05), transfer=0.2),
        "transfer and origination": dict(
            duplication=0.2, loss=0.3, transfer=PerCopy(0.4).scaled_by(hab, tab),
            origination=PerLineage(0.5).scaled_by(hab, tab)),
        "a family draw beside it": dict(
            duplication=PerCopy(0.3).varying_among("families", draw),
            loss=PerCopy(0.3).scaled_by(hab, tab), transfer=0.2),
        "a family draw on it": dict(
            duplication=PerCopy(0.3).scaled_by(hab, tab).varying_among("families", draw),
            loss=0.3, transfer=PerCopy(0.2).scaled_by(hab, tab)),
        "a schedule on it": dict(
            duplication=PerCopy(0.3).scaled_by(hab, tab).changing_at({0: 1.0, 1.0: 0.3, 2.0: 2.0}),
            loss=PerCopy(0.3).changing_at({0: 1.0, 1.5: 2.0}), transfer=0.2),
        "a per-lineage scope": dict(
            duplication=PerLineage(0.8).scaled_by(hab, tab),
            loss=PerLineage(0.5).scaled_by(hab, tab), transfer=0.2),
        "replacing transfers to oneself": dict(
            duplication=0.2, loss=0.2, transfer=PerCopy(0.6).scaled_by(hab, tab),
            replacement=0.5, self_transfer=True),
        "a per-lineage draw on it": dict(
            duplication=PerCopy(0.3).scaled_by(hab, tab).varying_among("lineages", LogNormal(0.0, 0.5)),
            loss=PerCopy(0.3).varying_among("lineages", LogNormal(0.0, 0.5)), transfer=0.2),
    }


def _driven_serial_only(hab):
    """A family writing its own rate runs on the serial engine alone."""
    tab = {"cave": 3.0, "surface": 0.5}
    own = PerCopy(0.5).scaled_by(hab, {"cave": 0.0, "surface": 1.0})
    return {
        "a family's own number": dict(
            duplication=PerCopy(0.3).scaled_by(hab, tab), loss=0.3, transfer=0.2,
            families=[family("core", loss=0.01), family("mobile", transfer=1.5)]),
        "a family's own driven rate": dict(
            duplication=0.3, loss=0.3, transfer=0.2, families=[family("x", loss=own)]),
        "both driven": dict(
            duplication=0.3, loss=PerCopy(0.3).scaled_by(hab, tab), transfer=0.2,
            families=[family("x", loss=own),
                      family("y", loss=PerCopy(0.4).changing_at({0: 0.2, 1.0: 2.0}))]),
    }


def _undriven():
    """Runs that read no driver: these stay what they were, as shipped."""
    lin, draw = LogNormal(0.0, 0.5), Gamma(shape=2.0, scale=0.5)
    return {
        "per-lineage draws": dict(
            duplication=PerCopy(0.3).varying_among("lineages", lin),
            loss=PerCopy(0.3).varying_among("lineages", lin),
            transfer=PerCopy(0.2).varying_among("lineages", lin),
            origination=PerLineage(0.3).varying_among("lineages", lin)),
        "a per-lineage and a per-family draw": dict(
            duplication=PerCopy(0.3).varying_among("lineages", lin).varying_among("families", draw),
            loss=0.3, transfer=0.2),
        "plain numbers": dict(duplication=0.3, loss=0.3, transfer=0.2, origination=0.2),
        "a per-family draw": dict(
            duplication=PerCopy(0.3).varying_among("families", draw), loss=0.3, transfer=0.2),
        "a schedule": dict(
            duplication=PerCopy(0.3).changing_at({0: 1.0, 1.0: 0.3}), loss=0.3, transfer=0.2),
    }


def _digest(run) -> str:
    return hashlib.sha256("\n".join(repr(e) for e in run.events).encode()).hexdigest()[:12]


def _run(tree, kw, **engine):
    return simulate_genomes_family(tree, initial_families=15, seed=7, **kw, **engine)


def _modules():
    """The two engine modules. `zombi2.genomes.family` the attribute is the `family()` function,
    so the module is fetched by its name."""
    import zombi2.genomes._perfamily  # noqa: F401
    return sys.modules["zombi2.genomes.family"], sys.modules["zombi2.genomes._perfamily"]


# --- the two small structures -------------------------------------------------------------------

def test_summed_weights_adds_and_draws_as_a_list_would():
    from zombi2._runtime.draw import SummedWeights, weighted_index
    rng = np.random.default_rng(1)
    values = [float(v) for v in rng.random(37)]
    w = SummedWeights()
    for v in values:
        w.append(v)
    w.set(5, 0.25)
    values[5] = 0.25
    w.remove(3)                                   # the swap-remove `_live.retire` makes
    values[3] = values[-1]
    values.pop()
    assert list(w) == values and len(w) == 36 and w[3] == values[3]
    assert w.total == sum(values)                 # the same float, not an approximation of it
    a, b = np.random.default_rng(9), np.random.default_rng(9)
    assert [w.pick(a) for _ in range(200)] == [weighted_index(b, values, sum(values)) for _ in range(200)]


def test_driver_switches_names_the_next_switch_and_who_switched(world):
    """Against asking every living lineage, which is what the loop used to do — with lineages
    leaving and coming back, as one family's footprint does."""
    from zombi2.genomes._live import DriverSwitches
    from zombi2.params.conditioned import resolve_driver
    tree, hab, _bm = world
    trajs = {"h": resolve_driver(hab, tree)}
    rng = np.random.default_rng(3)
    nodes = sorted(tree.nodes)
    t = tree.nodes[tree.root].birth_time
    pos: dict[int, int] = {}
    switches = DriverSwitches(trajs)
    seen = 0
    for _ in range(400):
        node = nodes[int(rng.integers(len(nodes)))]
        if node in pos:
            del pos[node]                         # it leaves; its entry is dropped when it surfaces
        else:
            pos[node] = len(pos)
            switches.entered(node, t)
        expected = min((trajs["h"].next_change(n, t) for n in pos), default=math.inf)
        assert switches.next(pos) == expected
        if expected == math.inf or rng.random() < 0.5:
            continue
        who = sorted(n for n in pos if trajs["h"].next_change(n, t) == expected)
        t = expected
        assert sorted(switches.due(t, pos)) == who
        seen += len(who)
    assert seen >= 10                             # the walk did cross switches, not only churn


# --- a kept weight is the float a fresh one would be ----------------------------------------------

def test_every_kept_weight_matches_one_built_afresh(world, monkeypatch):
    tree, hab, bm = world
    serial, _ = _modules()
    monkeypatch.setattr(serial, "_CHECK_KEPT", True)
    for kw in ({**_driven(hab, bm), **_driven_serial_only(hab), **_undriven()}).values():
        _run(tree, kw)                            # raises on the first kept value that differs


# --- and so the run is the run it was -------------------------------------------------------------

#: Taken on an export of main at 0.57.1 (a66b2dff), serial engine then per-family engine.
_BEFORE = {
    "a trait on duplication and loss": ("e2517460a321", "b5d56ffce2ee"),
    "a continuous trait, stepped": ("d8bf0831620d", "2e302820edd2"),
    "transfer and origination": ("f0be955f1627", "b2595fc91b67"),
    "a family draw beside it": ("2fa196e14efa", "07a7d069dc9b"),
    "a family draw on it": ("3f8d9db9382e", "4e79c6cc6538"),
    "a schedule on it": ("00feabcd3f9e", "93821f415508"),
    "a per-lineage scope": ("6ad4c4b2dd9c", "96d727220cbd"),
    "replacing transfers to oneself": ("4d3e6c7ab8b9", "c14da29d0645"),
    "a per-lineage draw on it": ("c8203c50d32d", "85547f40a927"),
    "per-lineage draws": ("676f26ab97bf", "a4b69b8a6610"),
    "a per-lineage and a per-family draw": ("b5cb06d8af3b", "3b218906e3e1"),
    "plain numbers": ("153ae959d10d", "ffa8d8dd72aa"),
    "a per-family draw": ("930fb71bee3e", "315c4c508cca"),
    "a schedule": ("12d6ded0bea5", "f8f202276e58"),
    "a family's own number": ("c49e5db9e468", None),
    "a family's own driven rate": ("f4f0483a13ee", None),
    "both driven": ("d15df3f3172e", None),
}


def test_with_sums_left_to_right_a_driven_run_is_the_run_it_was(world, monkeypatch):
    """Keeping the weights changes the cost and nothing else. Shown by taking the one other change
    away: hold the kept weights in `SummedWeights`, which adds as the old loop did, and every driven
    run is the 0.57.1 run to the last digit, on both engines."""
    from zombi2._runtime.draw import SummedWeights
    tree, hab, bm = world
    for module in _modules():
        monkeypatch.setattr(module, "WeightedIndex", SummedWeights)
    for label, kw in _driven(hab, bm).items():
        assert (_digest(_run(tree, kw)), _digest(_run(tree, kw, parallel=1))) == _BEFORE[label], label
    for label, kw in _driven_serial_only(hab).items():
        assert _digest(_run(tree, kw)) == _BEFORE[label][0], label


def test_a_run_reading_no_driver_is_unchanged(world):
    tree, _hab, _bm = world
    for label, kw in _undriven().items():
        assert (_digest(_run(tree, kw)), _digest(_run(tree, kw, parallel=1))) == _BEFORE[label], label


def test_a_driven_run_on_the_tree_is_the_same_process(world):
    """Added in pairs, a total can differ from the left-to-right one in its last digit, and an event
    time with it. It is the same process: on these runs the same events, in the same order, at
    times that agree to a part in a million million."""
    from zombi2._runtime.draw import SummedWeights
    tree, hab, bm = world
    kw = _driven(hab, bm)["a continuous trait, stepped"]
    on_tree = _run(tree, kw)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(_modules()[0], "WeightedIndex", SummedWeights)
        left_to_right = _run(tree, kw)
    assert len(on_tree.events) == len(left_to_right.events)
    for a, b in zip(on_tree.events, left_to_right.events):
        assert (a.kind, a.family, a.parents, a.children) == (b.kind, b.family, b.parents, b.children)
        assert a.time == pytest.approx(b.time, rel=1e-12)


# --- and an event costs no more on a larger tree ----------------------------------------------------

@pytest.mark.parametrize("engine", [{}, {"parallel": 1}], ids=["serial", "per-family"])
def test_a_driven_event_costs_no_more_on_a_larger_tree(monkeypatch, engine):
    """The regression #487 reported, counted rather than timed: rate evaluations per row of the
    event log. Rebuilding every weight made that number grow with the living lineages (57 per event
    at 100 tips); kept, it is a few, whatever the tree."""
    from zombi2.params.parameter import Rate
    calls = 0
    effective = Rate.effective

    def counted(self, **context):
        nonlocal calls
        calls += 1
        return effective(self, **context)

    monkeypatch.setattr(Rate, "effective", counted)
    per_row = []
    for tips in (30, 120):
        tree = simulate_species_tree(birth=1.0, n_extant=tips, seed=1).complete_tree
        hab = simulate_discrete(tree, states=["a", "b"], switch=0.5, seed=2)
        calls = 0
        run = simulate_genomes_family(
            tree, duplication=PerCopy(0.2).scaled_by(hab, {"a": 1.5, "b": 0.7}),
            loss=PerCopy(0.3).scaled_by(hab, {"a": 0.7, "b": 1.4}), transfer=0.2,
            initial_families=60, max_family_size=20, seed=3, **engine)
        per_row.append(calls / len(run.events))
    assert per_row[1] < 10                         # it was past 60 at this size
    assert per_row[1] < 1.5 * per_row[0]           # and four times the tips is not four times the cost
