"""Coupling slice 1 — conditioned trait→gene loss (SPEC §2, coupling-api.md).

The one mechanism (``mod.DrivenBy``) and its conditioned use: a discrete trait grown first, written
to a driver file, then read by a genome run whose loss rate is driven by it. Covers the mapping
shapes, the DrivenBy modifier, the driver trajectory + file round-trip, the traits driver writer,
and the end-to-end coupling with a seed-independent correctness invariant.
"""

import math

import pytest

from zombi2 import genomes, traits
from zombi2.rates.driver import DriverTrajectory, load_driver
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Curve, Scalar, Table, as_mapping
from zombi2.species import simulate_species_tree


# --- the mapping shapes (Table / Curve / Scalar) --------------------------------------------------

def test_table_lookup_and_default():
    m = Table({"aquatic": 3.0, "terrestrial": 1.0})
    assert m.multiplier("aquatic") == 3.0
    assert m.multiplier("terrestrial") == 1.0
    assert m.multiplier("unlisted") == 1.0            # default
    assert Table({"a": 2.0}, default=0.5).multiplier("b") == 0.5


def test_table_rejects_bad_factors():
    with pytest.raises(ValueError):
        Table({"a": -1.0})                            # negative factor
    with pytest.raises(ValueError):
        Table({})                                     # empty
    with pytest.raises(ValueError):
        Table({"a": 1.0}, default=-2.0)


def test_curve_and_bound():
    m = Curve(lambda x: math.exp(0.5 * x))
    assert m.multiplier(0.0) == pytest.approx(1.0)
    assert m.multiplier(2.0) == pytest.approx(math.exp(1.0))
    assert Curve(lambda x: x, bound=5.0).multiplier(100.0) == 5.0
    with pytest.raises(ValueError):
        Curve(lambda x: -1.0).multiplier(3.0)         # a rate multiplier cannot be negative


def test_scalar_log_link():
    assert Scalar(0.0).multiplier(1.0) == pytest.approx(1.0)       # strength 0 → null
    assert Scalar(0.7).multiplier(1.0) == pytest.approx(math.exp(0.7))
    assert Scalar(1.0).multiplier(1000.0) == pytest.approx(math.exp(40.0))  # clamped, no overflow


def test_table_matches_states_by_string_form():
    # an int-labelled trait keyed with int factors must NOT silently miss (the driver file is text)
    m = Table({0: 3.0, 1: 1.0})
    assert m.multiplier(0) == 3.0          # native int value
    assert m.multiplier("0") == 3.0        # same value read back from a file as a string
    assert m.multiplier(1) == 1.0
    assert Table({"0": 3.0}).multiplier(0) == 3.0   # string key, int value → still matches
    with pytest.raises(ValueError, match="collide"):
        Table({0: 1.0, "0": 2.0})          # two keys collide as strings


def test_continuous_mapping_on_discrete_driver_errors_clearly():
    with pytest.raises(ValueError, match="continuous-driver response"):
        Scalar(0.7).multiplier("cave")
    with pytest.raises(ValueError, match="continuous-driver response"):
        Curve(lambda x: x).multiplier("cave")


def test_as_mapping_coercion():
    assert isinstance(as_mapping({"a": 2.0}), Table)
    assert isinstance(as_mapping(lambda x: x), Curve)
    assert isinstance(as_mapping(0.5), Scalar)
    t = Table({"a": 2.0})
    assert as_mapping(t) is t                           # already a mapping → unchanged
    with pytest.raises(TypeError):
        as_mapping(True)
    with pytest.raises(TypeError):
        as_mapping("nope")


# --- the DrivenBy modifier ------------------------------------------------------------------------

def test_drivenby_factor_reads_threaded_value():
    d = mod.DrivenBy("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})
    assert d.factor(drivers={"habitat.tsv": "aquatic"}) == 3.0
    assert d.factor(drivers={"habitat.tsv": "terrestrial"}) == 1.0


def test_drivenby_inert_without_driver():
    d = mod.DrivenBy("trait", {"a": 5.0})
    assert d.factor() == 1.0                            # no drivers threaded → inert
    assert d.factor(drivers={"other": "a"}) == 1.0      # this source absent → inert
    assert d.next_change(0.0) == math.inf               # the engine owns per-lineage switching


def test_drivenby_builds_a_rate():
    los = 0.25 * mod.DrivenBy("f.tsv", {"hi": 4.0})
    from zombi2.rates.scope import PerCopy
    from zombi2.rates.rate import as_rate
    r = as_rate(los, default_scope=PerCopy)
    # base × copies × mapped factor
    assert r.effective(copies=10, drivers={"f.tsv": "hi"}) == pytest.approx(0.25 * 10 * 4.0)
    assert r.effective(copies=10, drivers={"f.tsv": "lo"}) == pytest.approx(0.25 * 10 * 1.0)  # default


def test_drivenby_validates_source():
    with pytest.raises(ValueError):
        mod.DrivenBy("", {"a": 1.0})
    with pytest.raises(ValueError):
        mod.DrivenBy("   ", {"a": 1.0})


# --- the driver trajectory + file round-trip ------------------------------------------------------

def test_driver_trajectory_lookup():
    # lineage 0 switches lo→hi at t=1.5; lineage 1 is constant "lo"
    traj = DriverTrajectory({0: [(0.0, "lo"), (1.5, "hi")], 1: [(0.0, "lo")]})
    assert traj.value(0, 0.0) == "lo"
    assert traj.value(0, 1.49) == "lo"
    assert traj.value(0, 1.5) == "hi"                   # right-continuous at the switch
    assert traj.value(0, 3.0) == "hi"
    assert traj.next_change(0, 0.0) == pytest.approx(1.5)
    assert traj.next_change(0, 1.5) == math.inf         # no further switch
    assert traj.value(1, 9.0) == "lo"
    assert traj.next_change(1, 0.0) == math.inf
    with pytest.raises(KeyError):
        traj.value(99, 0.0)                             # a lineage not in the file


def test_driver_file_round_trip(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.2, seed=7).complete_tree
    hab = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=1.5, seed=1)
    hab.write(tmp_path, outputs=("driver",))
    traj = load_driver(tmp_path / "trait_driver.tsv")
    # the reconstructed trajectory agrees with the trait's own node values at each node's end time
    for i, node in tree.nodes.items():
        assert traj.value(i, node.end_time - 1e-9) == hab.node_values[i]


def test_driver_write_requires_discrete(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=2).complete_tree
    cont = traits.simulate_continuous(tree, rate=1.0, seed=1)
    with pytest.raises(ValueError):
        cont.write(tmp_path, outputs=("driver",))       # a diffusion has no stochastic map


# --- end-to-end conditioned coupling: a trait drives gene loss ------------------------------------

def _write_driver(path, tree, state_of):
    """Write a one-segment-per-branch driver file assigning ``state_of[node]`` to each lineage."""
    rows = ["node\tstart\tend\tstate"]
    for i in sorted(tree.nodes):
        node = tree.nodes[i]
        rows.append(f"n{i}\t{node.birth_time:.6g}\t{node.end_time:.6g}\t{state_of[i]}")
    path.write_text("\n".join(rows) + "\n")


def test_zero_factor_lineages_never_lose(tmp_path):
    """The sharp invariant: a lineage whose loss factor is 0 cannot lose a gene, at any seed. We
    assign states by hand so half the tree has loss factor 0 (``lo``) and half a large factor
    (``hi``); every loss event must land on a ``hi`` lineage and every ``lo`` lineage keeps all its
    inherited families."""
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    state_of[tree.root] = "lo"                          # keep the crown families intact to inherit
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, state_of)

    res = genomes.simulate_genomes_unordered(
        tree,
        loss=0.25 * mod.DrivenBy(str(driver), {"lo": 0.0, "hi": 40.0}),
        initial_families=6, seed=3,
    )
    # every loss event lands on a hi lineage (a lo lineage's loss rate is exactly 0)
    losses = [e for e in res.events if e.kind == "loss"]
    assert losses, "expected some loss on the hi lineages"
    assert all(state_of[e.lineage] == "hi" for e in losses)
    # a lo lineage changes nothing on its own branch (loss factor 0, no dup/transfer/origination
    # here): its family set equals what it inherited from its parent, whatever that was.
    for i, node in tree.nodes.items():
        if node.parent is not None and state_of[i] == "lo":
            fams = {c.family for c in res.genomes[i]}
            parent_fams = {c.family for c in res.genomes[node.parent]}
            assert fams == parent_fams, f"lo lineage n{i} changed vs parent: {fams} != {parent_fams}"


def test_driven_loss_is_deterministic(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=5).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, state_of)
    kw = dict(loss=0.3 * mod.DrivenBy(str(driver), {"lo": 1.0, "hi": 5.0}),
              initial_families=4, seed=9)
    a = genomes.simulate_genomes_unordered(tree, **kw)
    b = genomes.simulate_genomes_unordered(tree, **kw)
    assert [(e.time, e.kind, e.lineage, e.family) for e in a.events] == \
           [(e.time, e.kind, e.lineage, e.family) for e in b.events]


def test_end_to_end_trait_drives_loss(tmp_path):
    """The full coupling-api.md workflow: grow a habitat trait, write it, drive gene loss by it.
    Across the tree, lineages in the high-loss state should carry fewer copies on average than
    lineages in the low-loss state."""
    tree = simulate_species_tree(birth=1.1, total_time=3.0, seed=3).complete_tree
    hab = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.5, seed=1)
    hab.write(tmp_path, outputs=("driver",))
    res = genomes.simulate_genomes_unordered(
        tree,
        loss=0.15 * mod.DrivenBy(str(tmp_path / "trait_driver.tsv"),
                                 {"cave": 6.0, "surface": 1.0}),
        origination=0.2, initial_families=5, seed=2,
    )
    # compare mean copy count of extant tips by their (end-of-branch) habitat
    cave = [len(res.genomes[n.id]) for n in tree.extant() if hab.node_values[n.id] == "cave"]
    surface = [len(res.genomes[n.id]) for n in tree.extant() if hab.node_values[n.id] == "surface"]
    assert cave and surface, "need both habitats represented among the tips"
    assert sum(cave) / len(cave) < sum(surface) / len(surface)


def test_int_state_trait_drives_loss_end_to_end(tmp_path):
    """Finding-#1 regression, full pipeline: an int-labelled binary trait (states=[0, 1]) must drive
    loss through the file round-trip, not silently no-op. State 1 loses fast, state 0 never loses."""
    tree = simulate_species_tree(birth=1.1, total_time=2.5, seed=6).complete_tree
    trait = traits.simulate_discrete(tree, states=[0, 1], switch=0.5, seed=1)
    trait.write(tmp_path, outputs=("driver",))
    res = genomes.simulate_genomes_unordered(
        tree,
        loss=0.3 * mod.DrivenBy(str(tmp_path / "trait_driver.tsv"), {0: 0.0, 1: 30.0}),
        initial_families=5, seed=2,
    )
    losses = [e for e in res.events if e.kind == "loss"]
    assert losses, "the int-keyed mapping must actually bite (not silently default to 1.0)"
    # every loss is on a state-1 branch (state 0 has loss factor exactly 0)
    driver = load_driver(tmp_path / "trait_driver.tsv")
    assert all(driver.value(e.lineage, e.time) == "1" for e in losses)


# --- the guard: what is not wired this slice ------------------------------------------------------

def test_driven_transfer_rejected():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="transfer"):
        genomes.simulate_genomes_unordered(
            tree, transfer=0.1 * mod.DrivenBy("f.tsv", {"a": 2.0}), initial_families=1, seed=1)
