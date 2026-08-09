"""Conditioning — a trait drives a genome or a sequence rate (SPEC §2, §3, §5).

The one mechanism (``mod.Driven``) and its conditioned uses: a discrete trait grown first, written
to a driver file, then read by a later run. Covers the mapping shapes, the Driven modifier, the
driver trajectory + file round-trip, the traits driver writer, the end-to-end trait→loss drive
with a seed-independent correctness invariant, both halves of trait-driven transfer — the donor
**rate** (how much HGT) and the ``transfer_to`` recipient **weight** (where it lands) — and the
trait→substitution drive at the sequences level, where the branch length in substitutions is the
driver integrated across the branch.
"""

import hashlib
import math
import re

import pytest

from zombi2 import genomes, sequences, traits
from zombi2.rates.driver import DriverTrajectory, load_driver
from zombi2.rates import LogNormal, ScaledBy, Weights, modifiers as mod
from zombi2.rates.mapping import Between, Curve, Scalar, Table, as_mapping, check_kernel_fires
from zombi2.sequences.substitution_models import hky85, jc69
from zombi2.species import simulate_species_tree
from zombi2.tree import Node, Tree, read_newick


def _branch_lengths(newick: str) -> dict:
    """``{species node id: branch length}`` from a species phylogram — read straight out of the text
    rather than through a tree, so the numbers are the ones the file carries."""
    return {int(i): float(bl) for _, i, bl in re.findall(r"([ne])(\d+):([0-9.eE+-]+)", newick)}


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
    with pytest.raises(ValueError, match="continuous-driver mapping"):
        Scalar(0.7).multiplier("cave")
    with pytest.raises(ValueError, match="continuous-driver mapping"):
        Curve(lambda x: x).multiplier("cave")


def test_as_mapping_coercion():
    assert isinstance(as_mapping({"a": 2.0}), Table)
    assert isinstance(as_mapping(lambda x: x), Curve)
    assert isinstance(as_mapping(0.5), Scalar)
    t = Table({"a": 2.0})
    assert as_mapping(t) is t                           # already a mapping → unchanged
    b = Between({("a", "b"): 2.0})
    assert as_mapping(b) is b                           # a kernel passes through (for ScaledBy(., Between))
    with pytest.raises(TypeError):
        as_mapping(True)
    with pytest.raises(TypeError):
        as_mapping("nope")


# --- the Between kernel: the 2-D, donor-conditioned choice-slot weight -----------------------------

def test_between_weight_pair_lookup_and_default():
    k = Between({("A", "B"): 3.0, ("B", "A"): 2.0}, default=0.5)
    assert k.weight("A", "B") == 3.0
    assert k.weight("B", "A") == 2.0
    assert k.weight("A", "A") == 0.5                    # unlisted pair → default
    assert k.weight("A", "rest") == 0.5
    assert Between({("A", "B"): 1.0}).weight("A", "A") == 1.0  # default default is 1.0 (baseline)
    assert k.groups() == {"A", "B"}


def test_between_matches_groups_by_string_form():
    k = Between({(0, 1): 3.0})
    assert k.weight(0, 1) == 3.0 and k.weight("0", "1") == 3.0   # int and str keys agree, like Table


def test_between_rejects_bad_input():
    with pytest.raises(ValueError, match="pairs"):
        Between({"A": 1.0})                             # a bare key, not a (from, to) pair
    with pytest.raises(ValueError, match="non-empty"):
        Between({})
    with pytest.raises(ValueError, match="weight"):
        Between({("A", "B"): -1.0})                     # negative weight
    with pytest.raises(ValueError, match="default"):
        Between({("A", "B"): 1.0}, default=-2.0)
    with pytest.raises(ValueError, match="collide"):
        Between({(0, 1): 1.0, ("0", "1"): 2.0})         # two keys collide as strings


def test_between_repr_round_trips_through_the_parser():
    from zombi2.rates.parse import parse_rate
    d = ScaledBy("f.tsv", Between({("A", "B"): 3.0, ("B", "A"): 3.0}, default=0.0))
    assert parse_rate(repr(d)) == d                     # the log line pastes back into a flag


def test_check_kernel_fires_needs_a_named_pair_to_occur():
    check_kernel_fires(Between({("A", "B"): 1.0}), {"A", "B"}, driver_label="x")  # both present → ok
    check_kernel_fires(Between({("A", "B"): 1.0, ("Q", "Z"): 1.0}), {"A", "B"}, driver_label="x")  # partial ok
    with pytest.raises(ValueError, match="silently do nothing"):
        check_kernel_fires(Between({("Q", "Z"): 1.0}), {"A", "B"}, driver_label="x")  # none present


# --- the Driven modifier ------------------------------------------------------------------------

def test_drivenby_factor_reads_threaded_value():
    d = ScaledBy("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})
    assert d.factor(drivers={"habitat.tsv": "aquatic"}) == 3.0
    assert d.factor(drivers={"habitat.tsv": "terrestrial"}) == 1.0


def test_drivenby_inert_without_driver():
    d = ScaledBy("trait", {"a": 5.0})
    assert d.factor() == 1.0                            # no drivers threaded → inert
    assert d.factor(drivers={"other": "a"}) == 1.0      # this source absent → inert
    assert d.next_change(0.0) == math.inf               # the engine owns per-lineage switching


def test_drivenby_builds_a_rate():
    los = 0.25 * ScaledBy("f.tsv", {"hi": 4.0})
    from zombi2.rates.scope import PerCopy
    from zombi2.rates.rate import as_rate
    r = as_rate(los, default_scope=PerCopy)
    # base × copies × mapped factor
    assert r.effective(copies=10, drivers={"f.tsv": "hi"}) == pytest.approx(0.25 * 10 * 4.0)
    assert r.effective(copies=10, drivers={"f.tsv": "lo"}) == pytest.approx(0.25 * 10 * 1.0)  # default


def test_drivenby_validates_driver():
    with pytest.raises(ValueError):
        ScaledBy("", {"a": 1.0})
    with pytest.raises(ValueError):
        ScaledBy("   ", {"a": 1.0})


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
    hab.write(tmp_path, outputs=("events",))
    traj = load_driver(tmp_path / "trait_events.tsv", tree)   # the log, replayed against the tree
    # the reconstructed trajectory agrees with the trait's own node values at each node's end time
    for i, node in tree.nodes.items():
        assert traj.value(i, node.end_time - 1e-9) == hab.node_values[i]


def test_a_continuous_log_is_only_its_origin(tmp_path):
    # a diffusion cannot be reconstructed from events, so its log carries only the t=0 marker — there
    # is no discrete map to drive a rate with (driving on a continuous trait is a later slice anyway)
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=2).complete_tree
    cont = traits.simulate_continuous(tree, rate=1.0, seed=1)
    cont.write(tmp_path, outputs=("events",))
    lines = (tmp_path / "trait_events.tsv").read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2 and lines[1].split("\t")[1] == "initial"


# --- end-to-end conditioning: a trait drives gene loss --------------------------------------------

def _write_driver(path, tree, state_of):
    """Write a trait **event log** assigning ``state_of[node]`` to each lineage for its whole branch:
    an ``initial`` row for the crown, then one ``on_speciation`` row per other node fixing its start
    state (no ``on_branch`` switches, so every branch is constant). Replayed against ``tree`` this
    rebuilds exactly ``state_of`` — the format a conditioned run reads now."""
    root = tree.root
    rows = ["time\tkind\tlineage\tfrom\tto",
            f"{tree.nodes[root].birth_time!r}\tinitial\tn{root}\t\t{state_of[root]}"]
    for i in sorted(tree.nodes):
        if i != root:
            rows.append(f"{tree.nodes[i].birth_time!r}\ton_speciation\tn{i}\t\t{state_of[i]}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


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

    res = genomes.simulate_genomes_family(
        tree,
        loss=0.25 * ScaledBy(str(driver), {"lo": 0.0, "hi": 40.0}),
        initial_families=6, seed=3,
    )
    # every loss event lands on a hi lineage (a lo lineage's loss rate is exactly 0)
    losses = [e for e in res.edges if e.kind == "loss"]
    assert losses, "expected some loss on the hi lineages"
    assert all(state_of[e.lineage] == "hi" for e in losses)
    # a lo lineage changes nothing on its own branch (loss factor 0, no dup/transfer/origination
    # here): its family set equals what it inherited from its parent, whatever that was.
    for i, node in tree.nodes.items():
        if node.parent is not None and state_of[i] == "lo":
            fams = {c.family for c in res.genomes[i]}
            parent_fams = {c.family for c in res.genomes[node.parent]}
            assert fams == parent_fams, f"lo lineage n{i} changed vs parent: {fams} != {parent_fams}"


def test_mapping_matching_no_driver_state_is_refused(tmp_path):
    # a mapping whose keys occur nowhere in the driver would leave every lineage at the default factor
    # — a silently undriven run — so it must be refused, not run as if it were driven
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, {i: ("hi" if i % 2 else "lo") for i in tree.nodes})
    with pytest.raises(ValueError, match="match none of the driver's states"):
        genomes.simulate_genomes_family(
            tree, loss=0.25 * ScaledBy(str(driver), {"cave": 4.0}),  # 'cave' is never a driver state
            initial_families=6, seed=3)


def test_partial_mapping_with_one_matching_state_still_runs(tmp_path):
    # ≥1 overlap is enough: a mapping may name a state this realisation never reached, as long as at
    # least one of its states does occur — that is a legitimate partial mapping, not a mistake
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, {i: "lo" for i in tree.nodes})   # only 'lo' ever occurs
    res = genomes.simulate_genomes_family(
        tree, loss=0.25 * ScaledBy(str(driver), {"lo": 2.0, "hi": 9.0}),  # 'hi' listed but absent
        initial_families=6, seed=3)
    assert res.events is not None                                # it ran; the absent 'hi' key is fine


def test_driven_loss_is_deterministic(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=5).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, state_of)
    kw = dict(loss=0.3 * ScaledBy(str(driver), {"lo": 1.0, "hi": 5.0}),
              initial_families=4, seed=9)
    a = genomes.simulate_genomes_family(tree, **kw)
    b = genomes.simulate_genomes_family(tree, **kw)
    assert [(e.time, e.kind, e.lineage, e.family) for e in a.edges] == \
           [(e.time, e.kind, e.lineage, e.family) for e in b.edges]


def test_end_to_end_trait_drives_loss(tmp_path):
    """The full conditioning workflow: grow a habitat trait, write it, drive gene loss by it.
    Across the tree, lineages in the high-loss state should carry fewer copies on average than
    lineages in the low-loss state."""
    tree = simulate_species_tree(birth=1.1, total_time=3.0, seed=3).complete_tree
    hab = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.5, seed=1)
    hab.write(tmp_path, outputs=("events",))
    res = genomes.simulate_genomes_family(
        tree,
        loss=0.15 * ScaledBy(str(tmp_path / "trait_events.tsv"),
                                 {"cave": 6.0, "surface": 1.0}),
        origination=0.2, initial_families=5, seed=2,
    )
    # compare mean copy count of extant tips by their (end-of-branch) habitat
    cave = [len(res.genomes[n.id]) for n in tree.extant_leaves() if hab.node_values[n.id] == "cave"]
    surface = [len(res.genomes[n.id]) for n in tree.extant_leaves() if hab.node_values[n.id] == "surface"]
    assert cave and surface, "need both habitats represented among the tips"
    assert sum(cave) / len(cave) < sum(surface) / len(surface)


def test_int_state_trait_drives_loss_end_to_end(tmp_path):
    """Finding-#1 regression, full pipeline: an int-labelled binary trait (states=[0, 1]) must drive
    loss through the file round-trip, not silently no-op. State 1 loses fast, state 0 never loses."""
    tree = simulate_species_tree(birth=1.1, total_time=2.5, seed=6).complete_tree
    trait = traits.simulate_discrete(tree, states=[0, 1], switch=0.5, seed=1)
    trait.write(tmp_path, outputs=("events",))
    res = genomes.simulate_genomes_family(
        tree,
        loss=0.3 * ScaledBy(str(tmp_path / "trait_events.tsv"), {0: 0.0, 1: 30.0}),
        initial_families=5, seed=2,
    )
    losses = [e for e in res.edges if e.kind == "loss"]
    assert losses, "the int-keyed mapping must actually bite (not silently default to 1.0)"
    # every loss is on a state-1 branch (state 0 has loss factor exactly 0)
    driver = load_driver(tmp_path / "trait_events.tsv", tree)
    assert all(driver.value(e.lineage, e.time) == "1" for e in losses)


# --- conditioning in-memory: Driven accepts the trait result object (no file step) --------------

def test_drivenby_accepts_traits_result_object(tmp_path):
    """Passing the discrete TraitsResult directly is the same conditioning as writing a file, and
    gives an IDENTICAL run (the file round-trip is lossless)."""
    tree = simulate_species_tree(birth=1.1, total_time=2.5, seed=3).complete_tree
    habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.5, seed=1)
    kw = dict(loss=0.5 * ScaledBy(habitat, {"aquatic": 3.0, "terrestrial": 1.0}),
              origination=0.2, initial_families=8, seed=2)
    by_object = genomes.simulate_genomes_family(tree, **kw)

    habitat.write(tmp_path, outputs=("events",))
    by_file = genomes.simulate_genomes_family(
        tree,
        loss=0.5 * ScaledBy(str(tmp_path / "trait_events.tsv"), {"aquatic": 3.0, "terrestrial": 1.0}),
        origination=0.2, initial_families=8, seed=2)
    key = lambda r: [(e.time, e.kind, e.lineage, e.copy) for e in r.edges]
    assert key(by_object) == key(by_file)


def test_drivenby_object_must_be_a_trait_result(tmp_path):
    # a conditioned driver object must be a grown trait result (discrete or continuous), carrying its
    # own complete tree — anything else (here a bare object) is refused
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=2).complete_tree
    with pytest.raises(ValueError, match="grown trait result"):
        genomes.simulate_genomes_family(
            tree, loss=0.5 * ScaledBy(object(), {"a": 2.0}), initial_families=3, seed=1)


# --- a trait drives transfer, side 1: the DONOR rate ----------------------------------------------
# Driving `transfer` says how often a lineage DONATES, so it changes the total amount of HGT.

_UNDRIVEN_TRANSFER_DIGESTS = {
    # captured from the engine BEFORE driven transfer was wired: the whole event log of a seeded
    # run, under each transfer_to rule. An undriven transfer must stay byte-identical — same rng
    # draw order, same results — however the driven path is built around it.
    #
    # Re-pinned when the log went to one row per event: a speciation's two edges are now recorded
    # one gene at a time (both daughters of a gene together) instead of one daughter at a time, which
    # moves this digest because it hashes the list in order. Nothing the run *produced* moved — the
    # ids are minted in the same order as before, and the genomes, the gene trees and the events as a
    # set are byte-identical against the previous engine.
    "uniform": "2f068999dc990e8cdfe17997c6c418e71372b4e35c4cbf3666513cfd43935ced",
    "distance": "a52f01064f17c8c1706065214d409326896fd98727dead3359fb8ccc739a4042",
}


def _event_digest(result) -> str:
    key = repr([(round(e.time, 12), e.kind, e.lineage, e.family, e.copy, e.parent, e.recipient)
                for e in result.edges])
    return hashlib.sha256(key.encode()).hexdigest()


@pytest.mark.parametrize("rule", ["uniform", "distance"])
def test_undriven_transfer_is_unchanged(rule):
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=2.5, seed=17).complete_tree
    res = genomes.simulate_genomes_family(
        tree, duplication=0.2, transfer=0.4, loss=0.15, origination=0.3,
        transfer_to=rule, initial_families=5, seed=23)
    assert _event_digest(res) == _UNDRIVEN_TRANSFER_DIGESTS[rule], (
        "an undriven transfer changed: the rng draw order of the undriven path must not move")


def test_undriven_transfer_is_unchanged_under_replacement_and_self_transfer():
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=2.5, seed=17).complete_tree
    res = genomes.simulate_genomes_family(
        tree, duplication=0.2, transfer=0.4, loss=0.15, origination=0.3, replacement=True,
        self_transfer=True, initial_families=5, seed=23)
    assert _event_digest(res) == "e14c05d3480058943566532526710b8ac482d7816739c3749080a51cfdfabe21"


def test_driven_transfer_picks_the_donor(tmp_path):
    """The sharp invariant on the donor side: a lineage whose transfer factor is 0 never donates."""
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, state_of)

    res = genomes.simulate_genomes_family(
        tree, transfer=0.2 * ScaledBy(str(driver), {"lo": 0.0, "hi": 20.0}),
        initial_families=6, seed=3)
    donations = [e for e in res.edges if e.kind == "transfer" and e.recipient is None]
    assert donations, "expected some donation from the hi lineages"
    assert all(state_of[e.lineage] == "hi" for e in donations)
    # the recipients are still drawn uniformly, so lo lineages do receive — the drive is on the
    # donor side only
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert any(state_of[e.lineage] == "lo" for e in arrivals)


def test_driven_transfer_changes_how_much_transfer_happens(tmp_path):
    """A driven transfer rate scales the amount of HGT: a flat factor of 3 gives about 3× the
    transfers. ``replacement`` holds the copy pool fixed so the count is linear in the rate. Pooled
    over seeds rather than trusting one: a driver switch is a Gillespie horizon, and how many an
    individual run hits is an rng-path detail that shifts the single-seed count without touching the
    rate — the mean is what the factor governs."""
    n_plain = n_driven = 0
    for seed in range(20):
        tree = simulate_species_tree(birth=1.1, total_time=2.0, seed=seed).complete_tree
        driver = tmp_path / f"flat{seed}.tsv"
        _write_driver(driver, tree, {i: "any" for i in tree.nodes})
        kw = dict(replacement=True, initial_families=8, seed=seed)
        plain = genomes.simulate_genomes_family(tree, transfer=0.2, **kw)
        driven = genomes.simulate_genomes_family(
            tree, transfer=0.2 * ScaledBy(str(driver), {"any": 3.0}), **kw)
        n_plain += sum(1 for e in plain.edges if e.kind == "transfer" and e.recipient is not None)
        n_driven += sum(1 for e in driven.edges if e.kind == "transfer" and e.recipient is not None)
    assert 2.7 < n_driven / n_plain < 3.3


def test_driven_transfer_is_deterministic(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=5).complete_tree
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, {i: ("hi" if i % 2 else "lo") for i in tree.nodes})
    kw = dict(transfer=0.3 * ScaledBy(str(driver), {"lo": 1.0, "hi": 5.0}),
              initial_families=4, seed=9)
    a = genomes.simulate_genomes_family(tree, **kw)
    b = genomes.simulate_genomes_family(tree, **kw)
    assert _event_digest(a) == _event_digest(b)


# --- a trait drives transfer, side 2: the RECIPIENT weight ----------------------------------------
# `transfer_to = Weights(...)` is the choice slot (SPEC §5): the mapping's numbers are weights
# over the contemporaneous candidates, so the same transfers are redistributed, not multiplied.

def _flat_tree_and_driver(tmp_path, competent):
    """Eight lineages alive together for essentially the whole run (a balanced tree whose internal
    branches are 1e-6 long), with ``competent`` of the eight tips in state ``competent``. Holding the
    candidate set fixed at all eight makes the expected recipient share exact rather than an average
    over changing alive sets."""
    e, tiny = 1e-6, 1e-6
    length = 1.0 - 2 * e
    nwk = (f"(((A:{length!r},B:{length!r}):{tiny!r},(C:{length!r},D:{length!r}):{tiny!r}):{tiny!r},"
           f"((E:{length!r},F:{length!r}):{tiny!r},(G:{length!r},H:{length!r}):{tiny!r}):{tiny!r});")
    tree, _ = read_newick(nwk)
    tips = [i for i, n in sorted(tree.nodes.items()) if n.children is None]
    hot = set(tips[:competent])
    driver = tmp_path / "competence.tsv"
    _write_driver(driver, tree, {i: ("competent" if i in hot else "normal") for i in tree.nodes})
    return tree, tips, hot, driver


def test_recipient_weight_splits_transfers_two_to_one(tmp_path):
    """Four candidates at weight 2 and four at weight 1 send 2/3 of transfers to the weight-2 group.
    ``self_transfer`` keeps the donor in the candidate set, so the normaliser is always 4·2 + 4·1."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    # no cap: a recipient already full of a family turns that arrival away, which at transfer=4.0
    # thins most of them (429 arrivals instead of 3000+) and biases what is left toward whoever is
    # not yet full. The cap is a confound for a measurement about the weights, so lift it — the same
    # move anyone estimating rates from a run has to make.
    res = genomes.simulate_genomes_family(
        tree, transfer=4.0, initial_families=6, self_transfer=True, max_family_size=None,
        transfer_to=Weights(str(driver), {"competent": 2.0, "normal": 1.0}), seed=5)
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert len(arrivals) > 1500                        # enough events for a 0.03 tolerance
    assert all(e.lineage in tips for e in arrivals)     # the internal branches are 1e-6 long
    share = sum(1 for e in arrivals if e.lineage in hot) / len(arrivals)
    assert share == pytest.approx(2 / 3, abs=0.03)


def test_recipient_weight_zero_cannot_receive(tmp_path):
    """Weight 0 means "cannot receive": every transfer lands on a competent lineage."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_family(
        tree, transfer=1.0, initial_families=6,
        transfer_to=Weights(str(driver), {"competent": 1.0, "normal": 0.0}), seed=5)
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert arrivals
    assert all(e.lineage in hot for e in arrivals)


def test_no_eligible_recipient_means_no_transfer_at_all(tmp_path):
    """When every candidate weighs 0 the transfer cannot happen, so the event is dropped whole —
    no donor continuation, no arrival, no copy minted. The same run under 'uniform' transfers
    freely, so the difference is the weighting and not the setup."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=0)   # nobody is competent
    kw = dict(transfer=1.0, initial_families=6, seed=5)
    blocked = genomes.simulate_genomes_family(
        tree, transfer_to=Weights(str(driver), {"competent": 1.0, "normal": 0.0}), **kw)
    free = genomes.simulate_genomes_family(tree, transfer_to="uniform", **kw)
    assert not [e for e in blocked.edges if e.kind == "transfer"]
    assert [e for e in free.edges if e.kind == "transfer"]
    # a dropped event leaves the genomes untouched: the six crown families are simply inherited
    assert all(len(g) == 6 for g in blocked.genomes.values())


def test_both_drivers_compose(tmp_path):
    """The donor rate and the recipient weight are independent models and may be used together."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_family(
        tree, transfer=0.5 * ScaledBy(str(driver), {"competent": 5.0, "normal": 0.0}),
        transfer_to=Weights(str(driver), {"competent": 0.0, "normal": 1.0}),
        initial_families=6, seed=5)
    donations = [e for e in res.edges if e.kind == "transfer" and e.recipient is None]
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert donations and arrivals
    assert all(e.lineage in hot for e in donations)          # only competent lineages donate
    assert all(e.lineage not in hot for e in arrivals)       # only non-competent lineages receive


# --- a Between kernel makes transfer_to donor-conditioned (assortative by a trait) -----------------

def test_recipient_kernel_keeps_transfer_within_the_donor_state(tmp_path):
    """ScaledBy(trait, Between(...)) reads the driver on the DONOR too, so a same-state kernel keeps
    every transfer within one habitat — the thing a 1-D recipient weight cannot express."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_family(
        tree, transfer=4.0, initial_families=6, self_transfer=True,
        transfer_to=Weights(str(driver),
                                 Between({("competent", "competent"): 1.0,
                                          ("normal", "normal"): 1.0}, default=0.0)), seed=5)
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert arrivals
    assert all((e.donor in hot) == (e.lineage in hot) for e in arrivals)  # donor and recipient agree


def test_recipient_kernel_fires_check_catches_absent_groups(tmp_path):
    """A kernel naming groups the driver never takes would weight every candidate at the default —
    secretly uniform — so it is refused, like a Table that names no occurring state."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    with pytest.raises(ValueError, match="silently do nothing"):
        genomes.simulate_genomes_family(
            tree, transfer=1.0, initial_families=6,
            transfer_to=Weights(str(driver), Between({("x", "y"): 1.0})), seed=5)


def test_between_is_rejected_in_a_rate_slot():
    """A rate has no donor to condition on, so a Between kernel there is a category error."""
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="donor-conditioned"):
        genomes.simulate_genomes_family(
            tree, loss=0.5 * ScaledBy("f.tsv", Between({("a", "b"): 1.0})),
            initial_families=1, seed=1)


def test_every_rate_and_extent_slot_rejects_a_between_kernel():
    """The same category error, refused at every slot that takes a driven rate or extent.

    Two of these used to reach the engine instead. `Between` deliberately implements no
    ``multiplier`` — it answers for a (donor, recipient) *pair* — so a slot without the guard did not
    quietly do the wrong thing, it died part-way through the run with
    ``AttributeError: 'Between' object has no attribute 'multiplier'``: a traceback from inside the
    engine naming neither the rate nor the mistake, after however long the run had taken to get
    there. The nucleotide engine had no guard on either its rates or its extents, and the ordered
    engine had one on its rates but not its extents."""
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1).complete_tree
    kernel = Between({("a", "a"): 3.0, ("b", "b"): 3.0})
    driven = ScaledBy("f.tsv", kernel)

    slots = [
        ("family rate", genomes.simulate_genomes_family, {"initial_families": 4, "loss": 0.2 * driven}),
        ("ordered rate", genomes.simulate_genomes_ordered, {"initial_families": 6, "loss": 0.2 * driven}),
        ("ordered extent", genomes.simulate_genomes_ordered,
         {"initial_families": 6, "inversion": 0.4, "inversion_extent": 3 * driven}),
        ("nucleotide rate", genomes.simulate_genomes_nucleotide,
         {"root_length": 2000, "genes": 4, "loss": 0.2 * driven}),
        ("nucleotide extent", genomes.simulate_genomes_nucleotide,
         {"root_length": 2000, "genes": 4, "inversion": 0.3, "inversion_extent": 200 * driven}),
    ]
    for label, fn, kwargs in slots:
        with pytest.raises(ValueError, match="donor-conditioned"):
            fn(tree, seed=1, **kwargs)
        # the message names the slot, so a user with several driven rates knows which one to change
        try:
            fn(tree, seed=1, **kwargs)
        except ValueError as e:
            assert "transfer_to" in str(e), label


# --- the guard: what is not wired ------------------------------------------------------------------

def test_transfer_to_rejects_a_rate():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="on its own, not a rate"):
        genomes.simulate_genomes_family(
            tree, transfer=0.1, transfer_to=1.0 * Weights("f.tsv", {"a": 2.0}),
            initial_families=1, seed=1)


def test_transfer_to_rejects_combining_distance_with_a_driven_weight():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="one recipient rule"):
        genomes.simulate_genomes_family(
            tree, transfer=0.1,
            transfer_to=(genomes.Distance(decay=1.0), Weights("f.tsv", {"a": 2.0})),
            initial_families=1, seed=1)


def test_transfer_to_rejects_an_unknown_rule():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="transfer_to must be"):
        genomes.simulate_genomes_family(tree, transfer=0.1, transfer_to="closest",
                                        initial_families=1, seed=1)


# --- the same choice slot at the ordered and nucleotide resolutions -------------------------------
# `transfer_to` is one slot with one kernel, and a block of genes or an arc of DNA is chosen a
# recipient exactly as a single copy is. These are the tests that say so at the other two resolutions.

def test_ordered_engine_takes_a_driven_transfer_to(tmp_path):
    """Weight 0 means "cannot receive" at the ordered resolution too: every arriving block lands on a
    competent lineage. What moves is a run of genes rather than one copy, which is a statement about
    the extent and not about who receives."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_ordered(
        tree, transfer=1.0, initial_families=6,
        transfer_to=Weights(str(driver), {"competent": 1.0, "normal": 0.0}), seed=5)
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert arrivals
    assert all(e.lineage in hot for e in arrivals)


def test_nucleotide_engine_takes_a_driven_transfer_to(tmp_path):
    """The same at the nucleotide resolution, read off the `Transfer` records rather than the gene
    genealogy. A nucleotide transfer is additive, so steering says only where the arc lands."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_nucleotide(
        tree, transfer=4.0, root_length=2000, genes=3, gene_length=100,
        transfer_to=Weights(str(driver), {"competent": 1.0, "normal": 0.0}), seed=5)
    arrivals = [e for e in res.events if type(e).__name__ == "Transfer"]
    assert arrivals
    assert all(e.recipient in hot for e in arrivals)


def test_recipient_weight_share_is_two_to_one_at_ordered(tmp_path):
    """The choice slot is a normalised WEIGHT, not a rate multiplier — the invariant that separates
    the two driven transfer slots (SPEC §5) — checked where it was never checked before. Four
    candidates at weight 2 and four at weight 1 take 2/3 of the arrivals. As at the family resolution,
    ``self_transfer`` keeps the donor in the candidate set so the normaliser is always 4·2 + 4·1, and
    the family cap is lifted because a recipient already full of a family turns that arrival away and
    biases what is left toward whoever is not yet full."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_ordered(
        tree, transfer=4.0, initial_families=6, self_transfer=True, max_family_size=None,
        transfer_to=Weights(str(driver), {"competent": 2.0, "normal": 1.0}), seed=5)
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert len(arrivals) > 1500                        # enough events for a 0.03 tolerance
    assert all(e.lineage in tips for e in arrivals)    # the internal branches are 1e-6 long
    share = sum(1 for e in arrivals if e.lineage in hot) / len(arrivals)
    assert share == pytest.approx(2 / 3, abs=0.03)


def test_driven_transfer_to_leaves_how_much_transfer_alone_at_ordered(tmp_path):
    """The other half of the same invariant: a flat driven weight redistributes nothing and must
    therefore change nothing at all about how many transfers happen. This is the test a threading bug
    would fail — putting the transfer_to trajectory into ``trajs`` would add a Gillespie breakpoint at
    every driver switch, moving the run while every assertion about *who* received still passed.
    Pooled over seeds, because a single run's count is an rng-path detail."""
    n_plain = n_driven = 0
    for seed in range(20):
        tree = simulate_species_tree(birth=1.1, total_time=2.0, seed=seed).complete_tree
        driver = tmp_path / f"flat{seed}.tsv"
        _write_driver(driver, tree, {i: "any" for i in tree.nodes})
        kw = dict(transfer=0.4, initial_families=8, max_family_size=None, seed=seed)
        plain = genomes.simulate_genomes_ordered(tree, **kw)
        driven = genomes.simulate_genomes_ordered(
            tree, transfer_to=Weights(str(driver), {"any": 3.0}), **kw)
        n_plain += sum(1 for e in plain.edges if e.kind == "transfer" and e.recipient is not None)
        n_driven += sum(1 for e in driven.edges if e.kind == "transfer" and e.recipient is not None)
    assert n_plain > 200
    assert 0.9 < n_driven / n_plain < 1.1


def test_recipient_kernel_keeps_transfer_within_the_donor_state_at_nucleotide(tmp_path):
    """A ``Between`` mapping reads the driver on the DONOR too, so a same-state kernel keeps every arc
    inside one guild. Proves the donor-conditioned kernel — not just the per-recipient weight —
    reaches the nucleotide engine."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_nucleotide(
        tree, transfer=6.0, root_length=2000, genes=3, gene_length=100, self_transfer=True,
        transfer_to=Weights(str(driver),
                                 Between({("competent", "competent"): 1.0,
                                          ("normal", "normal"): 1.0}, default=0.0)), seed=5)
    arrivals = [e for e in res.events if type(e).__name__ == "Transfer"]
    assert arrivals
    assert all((e.lineage in hot) == (e.recipient in hot) for e in arrivals)


def test_a_driven_transfer_to_composes_with_a_driven_transfer_rate_at_ordered(tmp_path):
    """The two slots are independent models — how often a lineage donates, and who receives — and the
    ordered engine now has both. The rate driver sets the Gillespie horizon; the weight driver does
    not, and they share one loaded trajectory."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_ordered(
        tree, transfer=0.5 * ScaledBy(str(driver), {"competent": 5.0, "normal": 0.0}),
        transfer_to=Weights(str(driver), {"competent": 0.0, "normal": 1.0}),
        initial_families=6, seed=5)
    donations = [e for e in res.edges if e.kind == "transfer" and e.recipient is None]
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert donations and arrivals
    assert all(e.lineage in hot for e in donations)          # only competent lineages donate
    assert all(e.lineage not in hot for e in arrivals)       # only non-competent lineages receive


def test_ordered_engine_still_refuses_a_between_kernel_on_a_rate():
    """The change reaches the choice slot and stops there. A rate has no donor to condition on, so a
    kernel in a rate slot stays a category error at the ordered resolution."""
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="donor-conditioned"):
        genomes.simulate_genomes_ordered(
            tree, loss=0.5 * ScaledBy("f.tsv", Between({("a", "b"): 1.0})),
            initial_families=1, seed=1)


def test_the_trait_and_joint_engines_refuse_a_between_kernel_on_a_rate():
    """The kernel guard, at the three slots it did not reach.

    A ``Between`` weights a recipient by the (donor, recipient) pair, so it belongs in
    ``transfer_to``. The genome engines refuse it in a rate; the two trait engines and the joint
    engine did not check, and ``Between`` implements no ``multiplier`` — so instead of a modelling
    error the run died part-way through with ``AttributeError: 'Between' object has no attribute
    'multiplier'``, a traceback from inside the engine naming neither the rate nor the mistake."""
    from zombi2 import joint, traits
    tree = simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1).complete_tree
    hab = traits.simulate_discrete(tree, states=["a", "b"], switch=0.3, seed=2)
    kernel = Between({("a", "b"): 2.0})

    with pytest.raises(ValueError, match="donor-conditioned"):
        traits.simulate_continuous(tree, start=0.0, rate=1.0 * ScaledBy(hab, kernel), seed=3)
    with pytest.raises(ValueError, match="donor-conditioned"):
        traits.simulate_discrete(tree, states=["x", "y"],
                                 switch=0.2 * ScaledBy(hab, kernel), seed=3)
    with pytest.raises(ValueError, match="donor-conditioned"):
        joint.simulate_joint(birth=1.0 * ScaledBy("trait", kernel), death=0.3,
                             trait=traits.discrete(states=["a", "b"], switch=0.3),
                             n_extant=20, seed=3)


def test_a_family_draw_on_one_rate_beside_a_driven_rate_is_refused(tmp_path):
    """Regression. The guard has to see every per-family draw in the run, not the one on the driven
    rate: a draw it missed was accepted, and then the loss total was summed WITHOUT the family
    multipliers while the copy was still drawn WITH them — a total saying one thing and a pick doing
    another."""
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, {i: ("hi" if i % 2 else "lo") for i in tree.nodes})
    with pytest.raises(ValueError, match="per-family draw and a driver on the same run"):
        genomes.simulate_genomes_family(
            tree, duplication=0.3 * mod.Drawn(per='family', dist=LogNormal(0.0, 0.5)),
            loss=0.2 * ScaledBy(str(driver), {"lo": 0.0, "hi": 5.0}),
            initial_families=6, seed=3)


def test_missing_driver_file_is_named_as_a_drivenby_driver(tmp_path):
    # a conditioned rate pointing at a file that is not there must say so as a Driven driver, not
    # leak a bare errno — so the user knows which kind of input is missing and how to fix it
    with pytest.raises(FileNotFoundError, match="driver file not found"):
        load_driver(str(tmp_path / "not_here.tsv"), None)


# --- the nucleotide resolution: a trait drives DNA loss ------------------------------------------
#
# This is the Chapter 1 genome-reduction story at the resolution where it is actually biology: a
# lifestyle trait drives how much DNA a lineage sheds, not how many family tokens it drops.

def _nuc_tree():
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=10, seed=4).complete_tree


def test_nucleotide_loss_is_driven_by_a_trait():
    """The seed-independent invariant: with the free-living factor at zero, **every** recorded loss
    must fall on a lineage that was host-restricted at the moment it happened. Anything else means
    the driver is being read for the total but not for the pick, or not at all."""
    tree = _nuc_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=2)
    res = genomes.simulate_genomes_nucleotide(
        tree, root_length=6000, genes=4, gene_length=300,
        loss=0.9 * ScaledBy(habitat, {"host": 25.0, "free": 0.0}),
        loss_extent=200, seed=2)

    from zombi2.rates.driver import driver_from_result
    traj = driver_from_result(habitat)
    losses = [e for e in res.events if type(e).__name__ == "Loss"]
    assert losses, "the run should produce losses at all"
    assert all(traj.value(e.lineage, e.time) == "host" for e in losses), \
        "a loss fell on a free-living lineage, whose factor is zero"


def test_nucleotide_driven_rate_matches_between_an_object_and_a_file(tmp_path):
    """A driver read from an in-memory trait result and the same driver read back from its own
    ``trait_events.tsv`` must give the identical run — the file is a serialisation, not a model."""
    tree = _nuc_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.6, seed=5)
    kw = dict(root_length=4000, genes=3, gene_length=250, loss_extent=180, seed=5)

    by_object = genomes.simulate_genomes_nucleotide(
        tree, loss=0.6 * ScaledBy(habitat, {"host": 8.0, "free": 1.0}), **kw)
    habitat.write(str(tmp_path), outputs=("events",))
    by_file = genomes.simulate_genomes_nucleotide(
        tree, loss=0.6 * ScaledBy(str(tmp_path / "trait_events.tsv"),
                                      {"host": 8.0, "free": 1.0}), **kw)
    assert [(e.time, type(e).__name__, e.lineage) for e in by_object.events] == \
           [(e.time, type(e).__name__, e.lineage) for e in by_file.events]


def test_nucleotide_driving_changes_the_run():
    """A driven rate must actually bite: the same seed with the driver switched off is a different
    run. Guards against a driver that resolves but never reaches the rate."""
    tree = _nuc_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=6)
    kw = dict(root_length=5000, genes=3, gene_length=250, loss_extent=200, seed=6)
    driven = genomes.simulate_genomes_nucleotide(
        tree, loss=0.6 * ScaledBy(habitat, {"host": 20.0, "free": 0.05}), **kw)
    flat = genomes.simulate_genomes_nucleotide(tree, loss=0.6, **kw)
    assert len(driven.events) != len(flat.events)


def test_nucleotide_undriven_run_is_unchanged():
    """No driver ⇒ the pooled path, byte for byte. The driven machinery must cost an undriven run
    nothing."""
    tree = _nuc_tree()
    kw = dict(root_length=4000, genes=3, gene_length=250, inversion=1.5, loss=0.4, seed=9)
    a = genomes.simulate_genomes_nucleotide(tree, **kw)
    b = genomes.simulate_genomes_nucleotide(tree, **kw)
    assert [(e.time, type(e).__name__, e.lineage) for e in a.events] == \
           [(e.time, type(e).__name__, e.lineage) for e in b.events]


def test_nucleotide_refuses_a_mapping_that_never_fires():
    """A mapping naming states the driver never takes leaves every lineage on the default factor, so
    the run would secretly be the undriven model."""
    tree = _nuc_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=7)
    with pytest.raises(ValueError):
        genomes.simulate_genomes_nucleotide(
            tree, root_length=3000, loss=0.5 * ScaledBy(habitat, {"aquatic": 3.0}), seed=7)


# --- the ordered resolution: a trait drives gene order -------------------------------------------
#
# The gene-level rates here are PER COPY, so a driven rate's per-lineage weight carries that
# lineage's own gene count and the pick is two-stage — a lineage by its weight, then a gene inside
# it. That is the family resolution's shape, not the nucleotide one (whose gene rates are per
# lineage and whose driven pick is a single "which lineage" draw). The count test below is what
# separates the two: a weight that ignored genome size would give ~N times fewer events.


def _ord_tree():
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=10, seed=4).complete_tree


def _ord_driver(tmp_path, tree, name="habitat.tsv"):
    """Half the lineages ``host``, half ``free``, constant along every branch — the same by-hand
    driver the family tests use, so a factor of 0 is an exact statement about which lineages can act."""
    state_of = {i: ("host" if i % 2 else "free") for i in tree.nodes}
    path = tmp_path / name
    _write_driver(path, tree, state_of)
    return state_of, str(path)


def test_ordered_loss_is_driven_by_a_trait(tmp_path):
    """The seed-independent invariant: with the free-living factor at zero, every recorded loss falls
    on a lineage that was host-restricted at that instant. Anything else means the driver reaches the
    total but not the pick, or neither."""
    tree = _ord_tree()
    state_of, driver = _ord_driver(tmp_path, tree)
    res = genomes.simulate_genomes_ordered(
        tree, loss=0.5 * ScaledBy(driver, {"host": 25.0, "free": 0.0}),
        loss_extent=2, initial_families=20, seed=2)
    losses = [e for e in res.edges if e.kind == "loss"]
    assert losses, "the run should produce losses at all"
    assert all(state_of[e.lineage] == "host" for e in losses)


@pytest.mark.parametrize("rate", ["inversion", "transposition", "translocation"])
def test_ordered_rearrangements_are_driven_by_a_trait(tmp_path, rate):
    """The headline case this feature exists for: a trait drives how often a lineage reshuffles its
    gene order. With the free-living factor at zero every inversion, transposition and translocation
    must land on a host lineage — and the rearrangement log is the only place they appear, since they
    begin and end no gene lineage."""
    tree = _ord_tree()
    state_of, driver = _ord_driver(tmp_path, tree)
    res = genomes.simulate_genomes_ordered(
        tree, chromosomes=3, initial_families=24,
        **{rate: 0.6 * ScaledBy(driver, {"host": 25.0, "free": 0.0}),
           f"{rate}_extent": 3}, seed=2)
    assert res.rearrangements, f"the run should produce {rate}s at all"
    assert all(state_of[x.lineage] == "host" for x in res.rearrangements)


def test_ordered_chromosome_tier_is_driven(tmp_path):
    """The per-chromosome scope path: a driven ``fission`` weights each lineage by its own chromosome
    count times its driver factor, so the pick is the weighted lineage and then a chromosome inside
    it. With the free factor at zero every fission is on a host lineage."""
    tree = _ord_tree()
    state_of, driver = _ord_driver(tmp_path, tree)
    res = genomes.simulate_genomes_ordered(
        tree, chromosomes=4, initial_families=24,
        fission=0.5 * ScaledBy(driver, {"host": 20.0, "free": 0.0}), seed=3)
    fissions = [c for c in res.chromosome_events if c.kind == "fission"]
    assert fissions, "the run should produce fissions at all"
    assert all(state_of[c.lineage] == "host" for c in fissions)


@pytest.mark.parametrize("rate", ["origination", "chromosome_origination"])
def test_ordered_per_lineage_rates_are_driven(tmp_path, rate):
    """The per-lineage scope path. Origination makes something from nothing, so it is counted per
    lineage and the driven pick is a single weighted-lineage draw — no gene or chromosome to choose
    inside it. A zero factor still means "never"."""
    tree = _ord_tree()
    state_of, driver = _ord_driver(tmp_path, tree)
    res = genomes.simulate_genomes_ordered(
        tree, initial_families=4,
        **{rate: 0.8 * ScaledBy(driver, {"host": 20.0, "free": 0.0})}, seed=4)
    if rate == "origination":
        t0 = tree.nodes[tree.root].birth_time     # the initial genome is logged at the origin
        acted = [e for e in res.edges if e.kind == "origination" and e.time > t0]
    else:
        acted = [c for c in res.chromosome_events if c.kind == "origination"]
    assert acted, f"the run should produce {rate}s at all"
    assert all(state_of[x.lineage] == "host" for x in acted)


def test_ordered_driven_transfer_picks_the_donor(tmp_path):
    """A driven ``transfer`` says how often a lineage **donates**. Free-living lineages never donate,
    but they still receive: the recipient is drawn by ``transfer_to``, which is a separate slot and
    takes no driver here."""
    tree = _ord_tree()
    state_of, driver = _ord_driver(tmp_path, tree)
    res = genomes.simulate_genomes_ordered(
        tree, transfer=0.4 * ScaledBy(driver, {"host": 20.0, "free": 0.0}),
        transfer_extent=2, initial_families=20, seed=5)
    donations = [e for e in res.edges if e.kind == "transfer" and e.recipient is None]
    arrivals = [e for e in res.edges if e.kind == "transfer" and e.recipient is not None]
    assert donations, "expected some donation from the host lineages"
    assert all(state_of[e.lineage] == "host" for e in donations)
    assert any(state_of[e.lineage] == "free" for e in arrivals)


def test_ordered_driven_rate_matches_the_per_copy_theory(tmp_path):
    """The quantitative check, against theory derived from SPEC §5 rather than from the code.

    ``inversion`` is per copy, so its total rate at any instant is ``base × (genes alive) × factor``.
    Inversions create and destroy no gene, and with no duplication, loss or origination every living
    lineage holds exactly ``initial_families`` genes — so the total is ``base × N × (lineages alive)``
    and the expected number of inversions over the whole run is

        E[#inversions] = factor × base × N × L,        L = the complete tree's total branch length.

    Both the flat driver (factor 1) and the tripled one are checked against that number. This is what
    pins the two-stage pick: a driven weight that forgot the lineage's gene count — the per-lineage
    shape the nucleotide engine uses — would give N times fewer events and sail past a ratio test."""
    base, families, factor = 0.4, 12, 3.0
    seen = {1.0: 0, factor: 0}
    expected = {1.0: 0.0, factor: 0.0}
    for seed in range(25):
        tree = simulate_species_tree(birth=1.1, total_time=2.0, seed=seed).complete_tree
        driver = tmp_path / f"flat{seed}.tsv"
        _write_driver(driver, tree, {i: "any" for i in tree.nodes})
        length = sum(n.end_time - n.birth_time for n in tree.nodes.values())
        for f in seen:
            res = genomes.simulate_genomes_ordered(
                tree, inversion=base * ScaledBy(str(driver), {"any": f}),
                inversion_extent=2, initial_families=families, seed=seed)
            seen[f] += len(res.rearrangements)
            expected[f] += f * base * families * length
    for f in seen:
        assert seen[f] == pytest.approx(expected[f], rel=0.06), (
            f"factor {f}: {seen[f]} inversions against a theoretical {expected[f]:.0f}")


def test_ordered_extent_is_driven_by_a_trait(tmp_path):
    """The other axis (SPEC §6). An extent is ``base × modifiers``, and a driven extent's modifier
    scales the size drawn, so the mean run length on a host lineage is the factor times the mean on a
    free one. The inversion **rate** is flat here, so the two axes are shown to move independently."""
    tree = _ord_tree()
    state_of, driver = _ord_driver(tmp_path, tree)
    res = genomes.simulate_genomes_ordered(
        tree, inversion=1.0, chromosomes=1, initial_families=60,
        inversion_extent=2 * ScaledBy(driver, {"host": 6.0, "free": 1.0}), seed=6)
    host = [x.length for x in res.rearrangements if state_of[x.lineage] == "host"]
    free = [x.length for x in res.rearrangements if state_of[x.lineage] == "free"]
    assert len(host) > 50 and len(free) > 50
    assert sum(host) / len(host) > 3.0 * (sum(free) / len(free))
    # ...and the *rate* is untouched. Nothing here creates or destroys a gene, so every lineage holds
    # the same 60 genes throughout and its share of the inversions is its share of the branch length.
    # A driver that leaked from the extent into the rate would shift that share; it does not.
    span = {"host": 0.0, "free": 0.0}
    for i, node in tree.nodes.items():
        span[state_of[i]] += node.end_time - node.birth_time
    assert len(host) / len(free) == pytest.approx(span["host"] / span["free"], rel=0.15)


def test_ordered_driven_rate_matches_between_an_object_and_a_file(tmp_path):
    """A driver read from an in-memory trait result and the same driver read back from its own
    ``trait_events.tsv`` must give the identical run — the file is a serialisation, not a model."""
    tree = _ord_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.6, seed=5)
    kw = dict(inversion_extent=3, loss_extent=2, chromosomes=2, initial_families=20, seed=5)

    by_object = genomes.simulate_genomes_ordered(
        tree, loss=0.4 * ScaledBy(habitat, {"host": 8.0, "free": 1.0}),
        inversion=0.5 * ScaledBy(habitat, {"host": 4.0, "free": 1.0}), **kw)
    habitat.write(str(tmp_path), outputs=("events",))
    path = str(tmp_path / "trait_events.tsv")
    by_file = genomes.simulate_genomes_ordered(
        tree, loss=0.4 * ScaledBy(path, {"host": 8.0, "free": 1.0}),
        inversion=0.5 * ScaledBy(path, {"host": 4.0, "free": 1.0}), **kw)
    assert [(e.time, e.kind, e.lineage, e.copy) for e in by_object.edges] == \
           [(e.time, e.kind, e.lineage, e.copy) for e in by_file.edges]
    assert by_object.rearrangements == by_file.rearrangements


def test_ordered_driving_changes_the_run():
    """A driven rate must actually bite: the same seed with the driver switched off is a different
    run. Guards against a driver that resolves but never reaches the rate."""
    tree = _ord_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=6)
    kw = dict(inversion_extent=3, initial_families=20, seed=6)
    driven = genomes.simulate_genomes_ordered(
        tree, inversion=0.6 * ScaledBy(habitat, {"host": 20.0, "free": 0.05}), **kw)
    flat = genomes.simulate_genomes_ordered(tree, inversion=0.6, **kw)
    assert len(driven.rearrangements) != len(flat.rearrangements)


def test_ordered_refuses_a_mapping_that_never_fires():
    """A mapping naming states the driver never takes leaves every lineage on the default factor, so
    the run would secretly be the undriven model."""
    tree = _ord_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=7)
    with pytest.raises(ValueError, match="match none of the driver's states"):
        genomes.simulate_genomes_ordered(
            tree, inversion=0.5 * ScaledBy(habitat, {"aquatic": 3.0}),
            initial_families=6, seed=7)


def test_ordered_refuses_a_mapping_that_never_fires_on_an_extent():
    """The same check on the extent axis: an extent's driver is resolved and validated with the
    rates', so a stale mapping there cannot slip through either."""
    tree = _ord_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=7)
    with pytest.raises(ValueError, match="match none of the driver's states"):
        genomes.simulate_genomes_ordered(
            tree, inversion=0.5,
            inversion_extent=3 * ScaledBy(habitat, {"aquatic": 3.0}),
            initial_families=6, seed=7)


def test_ordered_refuses_byfamily_and_a_driver_together(tmp_path):
    """The one combination the ordered engine refuses, and it says why: one weights lineages by a
    driver, the other weights the segment by what it covers."""
    tree = _ord_tree()
    _state_of, driver = _ord_driver(tmp_path, tree)
    with pytest.raises(ValueError, match="per-family draw and a driver on the same run"):
        genomes.simulate_genomes_ordered(
            tree, duplication=0.2 * mod.Drawn(per='family', dist=LogNormal(0.0, 0.5)),
            loss=0.2 * ScaledBy(driver, {"host": 3.0, "free": 1.0}),
            initial_families=6, seed=1)


def test_ordered_refuses_a_family_draw_and_a_driver_together(tmp_path):
    """The same refusal at the ordered resolution, where a per-family draw weights the segment by
    what it covers rather than the lineage."""
    tree = _ord_tree()
    _state_of, driver = _ord_driver(tmp_path, tree)
    with pytest.raises(ValueError, match="per-family draw and a driver on the same run"):
        genomes.simulate_genomes_ordered(
            tree, duplication=0.3 * mod.Drawn(per='family', dist=LogNormal(0.0, 0.5)),
            loss=0.2 * ScaledBy(driver, {"host": 3.0, "free": 1.0}),
            initial_families=6, seed=1)


def test_ordered_refuses_between_in_a_rate_slot():
    """A kernel weights a recipient by the (donor, recipient) pair, and a rate has no donor to
    condition on — so it is a category error here as much as at the family resolution."""
    tree = _ord_tree()
    with pytest.raises(ValueError, match="donor-conditioned"):
        genomes.simulate_genomes_ordered(
            tree, inversion=0.5 * ScaledBy("f.tsv", Between({("a", "b"): 1.0})),
            initial_families=1, seed=1)


def test_ordered_refuses_byfamily_on_an_extent():
    """The 'meaningless' half of SPEC §5's two rejections: an extent is drawn before the run's genes
    are known, so there is no one family to draw a factor for — and the message must name the slot
    ByFamily does belong in rather than merely refuse."""
    tree = _ord_tree()
    with pytest.raises(ValueError, match="Put it on inversion"):
        genomes.simulate_genomes_ordered(
            tree, inversion=0.5, inversion_extent=3 * mod.Drawn(per='family', dist=LogNormal(0.0, 0.5)),
            initial_families=6, seed=1)


# --- a CONTINUOUS trait as the driver (approximate: each branch cut into constant sub-steps) --------

def test_continuous_driver_trajectory_interpolates():
    """`driver_from_continuous_result` cuts each branch into stretches of at most ``step`` time units
    whose value is the trait linearly interpolated from the parent's value to the node's, sampled at
    each stretch's midpoint — so `value()` returns exactly that, and `next_change()` steps within a
    branch."""
    import math

    from zombi2.rates.driver import driver_from_continuous_result
    ct = simulate_species_tree(birth=1.0, n_extant=6, seed=1).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=2)

    node = next(n for n in ct.nodes.values() if n.parent is not None and n.end_time > n.birth_time)
    dt = node.end_time - node.birth_time
    step = dt / 10                                       # ten stretches on *this* branch
    traj = driver_from_continuous_result(met, step=step)

    start_v, end_v = met.node_values[node.parent], met.node_values[node.id]
    n = max(1, math.ceil(dt / step))
    for k in (0, 3, n - 1):                              # the value in the k-th stretch = its midpoint
        expected = start_v + (end_v - start_v) * (k + 0.5) / n
        assert traj.value(node.id, node.birth_time + k * dt / n) == pytest.approx(expected)
    nxt = traj.next_change(node.id, node.birth_time)     # a within-branch breakpoint, not inf
    assert node.birth_time < nxt <= node.end_time


def test_the_continuous_driver_step_is_a_duration_not_a_count_of_pieces():
    """The resolution is per unit of **time**, so a stretch means the same thing on every branch.

    Cutting each branch into a fixed number of pieces makes the approximation as coarse as the branch
    is long — the error is then worst exactly where the driver has had most time to move, and refining
    it on the branches that need it also refines every branch that does not. A duration gives every
    stretch the same length wherever it sits: a branch twice as long gets twice as many, and halving
    the step doubles them everywhere."""
    from zombi2.rates.driver import default_step, driver_from_continuous_result, tree_height
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=25, seed=3).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=4)

    step = tree_height(ct) / 40
    traj = driver_from_continuous_result(met, step=step)
    for node in ct.nodes.values():                       # every stretch is at most `step` long
        starts = traj._starts[node.id]
        bounds = [*starts[1:], node.end_time]
        assert all(b - s <= step + 1e-12 for s, b in zip(starts, bounds)), node.id
        assert len(starts) == max(1, math.ceil((node.end_time - node.birth_time) / step))

    # halving the step doubles the work; a per-branch scheme would not move at all
    fine = driver_from_continuous_result(met, step=step / 2)
    n_coarse = sum(len(v) for v in traj._starts.values())
    n_fine = sum(len(v) for v in fine._starts.values())
    assert 1.9 < n_fine / n_coarse < 2.1, (n_coarse, n_fine)

    assert default_step(ct) == pytest.approx(tree_height(ct) * 0.01)
    for bad in (0.0, -1.0, float("inf")):
        with pytest.raises(ValueError, match="finite and positive"):
            driver_from_continuous_result(met, step=bad)


def test_a_continuous_trait_is_conditioned_on_from_its_values_file(tmp_path):
    """Conditioning from disk on a diffusion, which used to be a silent constant.

    A continuous trait has no switches, so its event log holds only the ``initial`` row. Replaying
    that log gave a driver frozen at the root value on every lineage — accepted without complaint, so
    a run that looked conditioned was the undriven model with one constant factor. The event log now
    refuses and names the value table, and the value table reproduces the in-memory driver exactly."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=5).complete_tree
    bm = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=6)
    bm.write(tmp_path, outputs=("values", "events"))
    curve = (lambda x: 1.0 + max(0.0, x))

    with pytest.raises(ValueError, match="CONTINUOUS trait's event log"):
        genomes.simulate_genomes_family(
            ct, initial_families=30, seed=7,
            loss=0.1 * ScaledBy(str(tmp_path / "trait_events.tsv"), curve))

    from_file = genomes.simulate_genomes_family(
        ct, initial_families=30, seed=7,
        loss=0.1 * ScaledBy(str(tmp_path / "trait_values.tsv"), curve))
    in_memory = genomes.simulate_genomes_family(
        ct, initial_families=30, seed=7, loss=0.1 * ScaledBy(bm, curve))
    assert [e.kind for e in from_file.events] == [e.kind for e in in_memory.events]

    # and the same trait read at two resolutions stays two drivers rather than being shared
    coarse = ScaledBy(bm, curve, step=0.5)
    fine = ScaledBy(bm, curve, step=0.01)
    assert coarse.key != fine.key


def test_continuous_driver_drives_a_rate_and_is_deterministic():
    """A continuous trait drives the duplication rate through a Curve: the driven run differs from an
    otherwise-identical undriven one (the driver is threaded and the mapping applied), and repeats
    exactly under the same seed."""
    ct = simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.2, seed=3)

    def run(fn):
        return genomes.simulate_genomes_family(
            ct, initial_families=15, loss=0.04,
            duplication=0.06 * ScaledBy(met, Curve(fn)), seed=9)

    driven = run(lambda v: 3.0 ** v)
    control = run(lambda v: 1.0)                          # a flat Curve == the undriven model
    tips = list(ct.extant_leaves())
    assert any(len(driven.genomes[n.id]) != len(control.genomes[n.id]) for n in tips)
    again = run(lambda v: 3.0 ** v)
    assert all(len(driven.genomes[n.id]) == len(again.genomes[n.id]) for n in tips)


def test_continuous_driver_takes_a_scalar_link():
    ct = simulate_species_tree(birth=1.0, n_extant=12, seed=1).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=2)
    genomes.simulate_genomes_family(                      # Scalar (log-link) is a valid continuous mapping
        ct, initial_families=10, duplication=0.05 * ScaledBy(met, Scalar(0.5)), seed=1)


def test_discrete_table_on_continuous_driver_is_refused():
    """A ``{state: factor}`` table names discrete states a continuous value never equals — refuse it
    with a message that points at Curve / Scalar, rather than silently leaving the rate undriven."""
    ct = simulate_species_tree(birth=1.0, n_extant=12, seed=1).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=2)
    with pytest.raises(ValueError, match="CONTINUOUS"):
        genomes.simulate_genomes_family(
            ct, initial_families=10, duplication=0.05 * ScaledBy(met, {"hi": 2.0}), seed=1)


def test_multitrait_continuous_driver_is_refused():
    from zombi2.rates.driver import driver_from_continuous_result
    ct = simulate_species_tree(birth=1.0, n_extant=8, seed=1).complete_tree
    two = traits.simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0}, seed=2)
    with pytest.raises(ValueError, match="SINGLE-trait"):
        driver_from_continuous_result(two)


# --- a trait drives the SUBSTITUTION rate (Traits → Sequences) ------------------------------------
#
# SPEC §3 allows Traits–Sequences to be conditioned and never joined: a trait can be grown first and
# held fixed, and a sequence never feeds back into it. So the driver is an ordinary modifier on
# `substitution`, and the branch length in substitutions/site is the driver integrated across the
# branch — a discrete trait switches mid-branch, and a single sample per branch would be a different
# model.

def _one_branch_run(total_time: float = 2.0):
    """A species tree of one lineage running from 0 to ``total_time``, carrying one gene family whose
    gene tree is that single branch. The smallest run in which a branch length can be read off by
    hand: the phylogram is one number."""
    tree = Tree({0: Node(0, None, 0.0, total_time, None, "extant")}, 0)
    return tree, genomes.simulate_genomes_family(tree, initial_families=1, seed=1)


def _switch_driver(path, at: float, before: str, after: str, total_time: float):
    """A driver on the one-lineage tree that is ``before`` on ``[0, at)`` and ``after`` from ``at`` on
    — one switch, strictly inside the branch."""
    path.write_text(
        "time\tkind\tlineage\tfrom\tto\n"
        f"0.0\tinitial\tn0\t\t{before}\n"
        f"{at!r}\ton_branch\tn0\t{before}\t{after}\n", encoding="utf-8")
    assert 0.0 < at < total_time                      # the switch has to be mid-branch to test anything


def _gene_sequences(result, gene_tree, labels, family=0):
    """``{id(gene-tree node): its sequence}`` for one family, from the two halves the level splits them
    into — the extant tips are in ``alignments`` and everything else in ``ancestral``. Keyed by
    identity because a gene-tree node is not hashable, which is also how the engine keys them."""
    from zombi2.genomes.events import gene_label

    out = {}
    stack = [gene_tree.complete]
    while stack:
        n = stack.pop()
        key = f"{labels[n.species]}_{gene_label(n.copy)}"
        out[id(n)] = result.alignments[family].get(key) or result.ancestral[family][key]
        stack.extend(n.children)
    return out


def test_a_trait_drives_the_substitution_rate(tmp_path):
    """The sharp, seed-independent invariant, the sequence twin of `test_zero_factor_lineages_never_lose`:
    a species branch whose driver factor is 0 has a branch length of 0 in substitutions, so every gene
    passing through it arrives at its far end unchanged. States are assigned by hand so half the tree
    is ``lo`` (factor 0) and half ``hi`` (factor 20)."""
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, state_of)
    run = genomes.simulate_genomes_family(tree, duplication=0.3, loss=0.2, initial_families=6, seed=3)

    res = sequences.simulate_sequences(
        run, model=jc69(), length=400,
        substitution=0.5 * ScaledBy(str(driver), {"lo": 0.0, "hi": 20.0}), seed=4)

    labels = tree.labels()
    changed_on_hi = False
    for family, gt in run.gene_trees.items():
        seq = _gene_sequences(res, gt, labels, family)
        pairs = [(res.founding[family], gt.complete)]
        stack = [gt.complete]
        while stack:
            n = stack.pop()
            for c in n.children:
                pairs.append((seq[id(n)], c))
                stack.append(c)
        for parent_seq, node in pairs:
            if state_of[node.species] == "lo":
                assert seq[id(node)] == parent_seq, \
                    f"a gene changed on lineage n{node.species}, whose substitution rate is 0"
            elif seq[id(node)] != parent_seq:
                changed_on_hi = True
    assert changed_on_hi, "nothing changed anywhere — the run would pass vacuously"


def test_a_mid_branch_switch_is_integrated_not_sampled(tmp_path):
    """A driver that switches halfway along a branch: the branch length in substitutions must be the
    driver **integrated** over the branch, not one sample of it.

    One lineage of length 2, factor 0 on ``[0, 1)`` and 2 on ``[1, 2]``, base 0.3. The integral is
    ``0.3 × (0×1 + 2×1) = 0.6``. Sampling the driver once per branch would give ``0.3 × 0 × 2 = 0``
    or ``0.3 × 2 × 2 = 1.2`` depending on which state was read — both are asserted against, because
    either would be a different model rather than a rounding difference."""
    tree, run = _one_branch_run(2.0)
    driver = tmp_path / "d.tsv"
    _switch_driver(driver, at=1.0, before="lo", after="hi", total_time=2.0)
    res = sequences.simulate_sequences(
        run, model=jc69(), length=10,
        substitution=0.3 * ScaledBy(str(driver), {"lo": 0.0, "hi": 2.0}), seed=2)

    length = float(res.phylograms[0]["complete"].split(":")[1].rstrip(";"))
    assert length == pytest.approx(0.6, rel=1e-6)
    assert length != pytest.approx(0.0, abs=1e-6)     # the branch read as 'lo' throughout
    assert length != pytest.approx(1.2, rel=1e-6)     # the branch read as 'hi' throughout
    # the species phylogram is the same branch, so it must carry the same number
    assert float(res.species_phylogram["complete"].split(":")[1].rstrip(";")) == \
        pytest.approx(0.6, rel=1e-6)


def test_the_driven_branch_length_is_the_exact_integral(tmp_path):
    """The integral itself, at full precision and by hand — the phylogram above is written to seven
    significant figures, which is not enough to tell an exact integral from a good approximation.

    Three stretches on one branch (factors 0.5, 4.0, 1.0 over [0,1), [1,1.25), [1.25,2]) give
    ``0.5×1 + 4×0.25 + 1×0.75 = 2.25``; a sub-stretch is checked too, so the running total is not
    only right at the breakpoints."""
    from zombi2.rates.driver import resolve_driver
    from zombi2.sequences.clock import resolve_clock

    tree = Tree({0: Node(0, None, 0.0, 2.0, None, "extant")}, 0)
    driver = tmp_path / "d.tsv"
    driver.write_text("time\tkind\tlineage\tfrom\tto\n"
                      "0.0\tinitial\tn0\t\ta\n"
                      "1.0\ton_branch\tn0\ta\tb\n"
                      "1.25\ton_branch\tn0\tb\tc\n", encoding="utf-8")
    m = ScaledBy(str(driver), {"a": 0.5, "b": 4.0, "c": 1.0})
    clock = resolve_clock(None, [(m, resolve_driver(m.driver, tree))], tree, {}, None)

    assert clock.branch_length(1.0, 0, 0.0, 2.0) == pytest.approx(2.25, abs=1e-12)
    # a sub-stretch ending inside the middle piece: 0.5×1 + 4×0.1
    assert clock.branch_length(1.0, 0, 0.0, 1.1) == pytest.approx(0.9, abs=1e-12)
    # and one starting inside it: 4×0.15 + 1×0.75
    assert clock.branch_length(2.0, 0, 1.1, 2.0) == pytest.approx(2.0 * 1.35, abs=1e-12)


def test_driven_divergence_matches_jukes_cantor(tmp_path):
    """The model check against theory. Under JC69 two sequences separated by ``d`` substitutions per
    site differ at a proportion ``0.75 · (1 − e^(−4d/3))`` of sites — derived from the model, not read
    off the implementation. With the driver integrating to ``d = 0.6`` that is 0.4130.

    The point of running it on a mid-branch switch is that the two per-branch-constant readings give
    ``d = 0`` (0.0) and ``d = 1.2`` (0.5628), both far outside the sampling window of 20 000 sites."""
    tree, run = _one_branch_run(2.0)
    driver = tmp_path / "d.tsv"
    _switch_driver(driver, at=1.0, before="lo", after="hi", total_time=2.0)
    res = sequences.simulate_sequences(
        run, model=jc69(), length=20000,
        substitution=0.3 * ScaledBy(str(driver), {"lo": 0.0, "hi": 2.0}), seed=17)

    start, end = res.founding[0], res.alignments[0]["n0_g0"]
    observed = sum(a != b for a, b in zip(start, end)) / len(start)
    expected = 0.75 * (1.0 - math.exp(-4.0 * 0.6 / 3.0))
    assert expected == pytest.approx(0.4130, abs=1e-4)          # the number this pins, spelled out
    assert observed == pytest.approx(expected, abs=0.014)       # ~4 standard errors at 20 000 sites


def test_a_driver_composes_with_a_lineage_clock():
    """SPEC §5: modifiers multiply. A `a per-lineage draw` clock and a ``Driven`` driver on one rate give a
    branch ``base × clock × ∫driver``, and each factor is recovered here independently — the clock
    from an otherwise-identical undriven run at the same seed (the draw comes first and consumes the
    same randomness either way), the integral by walking the trait's own trajectory."""
    from zombi2.rates.driver import driver_from_result

    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=21).complete_tree
    habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=2.0, seed=22)
    run = genomes.simulate_genomes_family(tree, duplication=0.2, loss=0.2, initial_families=3, seed=23)
    table = {"cave": 0.25, "surface": 3.0}
    base, kw = 0.4, dict(model=jc69(), length=20, seed=24)

    clocked = sequences.simulate_sequences(
        run, substitution=base * mod.Drawn(per='lineage', dist=LogNormal(0.0, 0.5)), **kw)
    both = sequences.simulate_sequences(
        run, substitution=base * mod.Drawn(per='lineage', dist=LogNormal(0.0, 0.5)) * ScaledBy(habitat, table), **kw)

    by_clock = _branch_lengths(clocked.species_phylogram["complete"])
    by_both = _branch_lengths(both.species_phylogram["complete"])
    traj = driver_from_result(habitat)
    moved = 0
    for i, node in tree.nodes.items():
        span = node.end_time - node.birth_time
        clock_factor = by_clock[i] / (base * span) if span else 1.0
        integral, t = 0.0, node.birth_time
        while t < node.end_time:
            nxt = min(traj.next_change(i, t), node.end_time)
            integral += table[traj.value(i, t)] * (nxt - t)
            t = nxt
        assert by_both[i] == pytest.approx(base * clock_factor * integral, rel=1e-6, abs=1e-9)
        if by_both[i] != pytest.approx(by_clock[i], rel=1e-6, abs=1e-12):
            moved += 1
    assert moved, "the driver changed no branch — the composition test would pass vacuously"


def test_the_species_phylogram_shows_the_driver():
    """``clock_species_tree_*.nwk`` is the molecular clock made visible, so it has to show the driven
    part too — a run whose tree disagreed with its own alignments would be worse than no tree."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=31).complete_tree
    habitat = traits.simulate_discrete(tree, states=["fast", "slow"], switch=1.0, seed=32)
    run = genomes.simulate_genomes_family(tree, duplication=0.2, initial_families=2, seed=33)
    kw = dict(model=jc69(), length=20, seed=34)

    flat = sequences.simulate_sequences(run, substitution=0.2, **kw)
    driven = sequences.simulate_sequences(
        run, substitution=0.2 * ScaledBy(habitat, {"fast": 8.0, "slow": 1.0}), **kw)
    plain = _branch_lengths(flat.species_phylogram["complete"])
    shown = _branch_lengths(driven.species_phylogram["complete"])
    assert any(shown[i] > plain[i] * 1.5 for i in plain), "no branch was stretched by the driver"
    assert all(shown[i] >= plain[i] - 1e-12 for i in plain)   # every factor here is ≥ 1


def test_driven_substitution_matches_between_an_object_and_a_file(tmp_path):
    """A driver handed over in memory and the same driver read back from its own ``trait_events.tsv``
    must give the identical run — the file is a serialisation, not a model (SPEC §2)."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=41).complete_tree
    habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=1.5, seed=42)
    run = genomes.simulate_genomes_family(tree, duplication=0.2, loss=0.2, initial_families=4, seed=43)
    table = {"cave": 0.3, "surface": 2.0}
    kw = dict(model=hky85(2.0), length=150, seed=44)

    by_object = sequences.simulate_sequences(run, substitution=0.3 * ScaledBy(habitat, table), **kw)
    habitat.write(str(tmp_path), outputs=("events",))
    by_file = sequences.simulate_sequences(
        run, substitution=0.3 * ScaledBy(str(tmp_path / "trait_events.tsv"), table), **kw)

    assert by_object.alignments == by_file.alignments
    assert by_object.ancestral == by_file.ancestral
    assert by_object.phylograms == by_file.phylograms
    assert by_object.species_phylogram == by_file.species_phylogram


def test_driven_substitution_is_deterministic():
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=51).complete_tree
    habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=1.5, seed=52)
    run = genomes.simulate_genomes_family(tree, duplication=0.2, initial_families=4, seed=53)
    spec = 0.3 * mod.Drawn(per='lineage', dist=LogNormal(0.0, 0.3)) * ScaledBy(habitat, {"cave": 0.5, "surface": 2.0})
    kw = dict(model=jc69(), length=120, substitution=spec, seed=54)
    a = sequences.simulate_sequences(run, **kw)
    b = sequences.simulate_sequences(run, **kw)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral
    assert a.phylograms == b.phylograms and a.species_phylogram == b.species_phylogram


def test_driven_sequences_refuse_a_mapping_that_never_fires():
    """The fires-check is wired at this level too, so a secretly-undriven run cannot pass as driven."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=61).complete_tree
    habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=1.5, seed=62)
    run = genomes.simulate_genomes_family(tree, duplication=0.2, initial_families=2, seed=63)
    with pytest.raises(ValueError, match="match none of the driver's states"):
        sequences.simulate_sequences(run, model=jc69(), length=20, seed=64,
                                     substitution=0.3 * ScaledBy(habitat, {"aquatic": 4.0}))


def test_a_continuous_trait_drives_the_substitution_rate():
    """A continuous driver needs a continuous mapping — a ``Curve`` or a ``Scalar``. Its trajectory is
    the piecewise-constant approximation `driver_from_continuous_result` builds, so the branch length
    is still an integral over stretches; a discrete ``{state: factor}`` table names states a float
    never equals and is refused."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=10, seed=71).complete_tree
    metabolism = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=72)
    run = genomes.simulate_genomes_family(tree, duplication=0.2, initial_families=3, seed=73)
    kw = dict(model=jc69(), length=60, seed=74)

    flat = sequences.simulate_sequences(run, substitution=0.3, **kw)
    curved = sequences.simulate_sequences(
        run, substitution=0.3 * ScaledBy(metabolism, Curve(lambda x: math.exp(0.3 * x))), **kw)
    linked = sequences.simulate_sequences(
        run, substitution=0.3 * ScaledBy(metabolism, Scalar(0.7)), **kw)
    assert curved.phylograms != flat.phylograms
    assert linked.phylograms != flat.phylograms
    with pytest.raises(ValueError, match="CONTINUOUS"):
        sequences.simulate_sequences(
            run, substitution=0.3 * ScaledBy(metabolism, {"hi": 2.0}), **kw)


def test_a_driven_run_is_the_same_in_parallel():
    """The parallel engine ships the resolved clock to its workers, and a ``Curve`` mapping is usually
    a lambda — which does not pickle. Precomputing the driver's integral to plain numbers is what lets
    it cross; this is the regression that catches storing the mapping on the clock instead."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=81).complete_tree
    metabolism = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=82)
    run = genomes.simulate_genomes_family(tree, duplication=0.3, initial_families=5, seed=83)
    spec = 0.3 * ScaledBy(metabolism, Curve(lambda x: math.exp(0.4 * x)))
    kw = dict(model=jc69(), length=80, substitution=spec, seed=84)
    one = sequences.simulate_sequences(run, parallel=1, **kw)
    two = sequences.simulate_sequences(run, parallel=2, **kw)
    assert one.alignments == two.alignments and one.phylograms == two.phylograms


def test_a_mapping_key_the_driver_never_takes_warns(tmp_path):
    """A typo in ONE state of a mapping used to pass in complete silence.

    The existing guard refuses a mapping where *nothing* matches. But a mapping with one good key and
    one typo fires (the good key drives the rate), so it sailed through — and the factor the user
    actually cared about was applied to nobody, while the run reported itself as driven. That is the
    shape of failure that gets a wrong result published: it does not look like a failure at all.

    It warns rather than raises because a mapping may legitimately name a state that this particular
    realisation never reached — the observed states are all we know when replaying a written driver.
    (A joint run, whose alphabet is declared up front, already raises: `exhaustive=True`.)"""
    import warnings as _w

    tree = simulate_species_tree(birth=1.0, n_extant=8, seed=1)
    habitat = traits.simulate_discrete(tree, states=("cave", "surface"), switch=0.6, seed=1)

    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        genomes.simulate_genomes_family(
            tree.complete_tree, initial_families=5,
            loss=0.3 * ScaledBy(habitat, {"caves": 5.0, "surface": 1.0}), seed=1)
    messages = [str(w.message) for w in caught]
    assert any("caves" in m and "never takes" in m for m in messages), messages
    # the state that IS right must not be reported as stray
    assert not any("'surface'" in m.split("never takes")[0] for m in messages if "never takes" in m)


def test_a_mapping_whose_keys_all_match_is_silent(tmp_path):
    # the guard must not cry wolf on the ordinary case
    import warnings as _w

    tree = simulate_species_tree(birth=1.0, n_extant=8, seed=1)
    habitat = traits.simulate_discrete(tree, states=("cave", "surface"), switch=0.6, seed=1)
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        genomes.simulate_genomes_family(
            tree.complete_tree, initial_families=5,
            loss=0.3 * ScaledBy(habitat, {"cave": 5.0, "surface": 1.0}), seed=1)
    assert [str(w.message) for w in caught if "never takes" in str(w.message)] == []
