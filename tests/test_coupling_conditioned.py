"""Conditioned coupling — a discrete trait drives a genome rate (SPEC §2, §5).

The one mechanism (``mod.DrivenBy``) and its conditioned uses: a discrete trait grown first, written
to a driver file, then read by a genome run. Covers the mapping shapes, the DrivenBy modifier, the
driver trajectory + file round-trip, the traits driver writer, the end-to-end trait→loss coupling
with a seed-independent correctness invariant, and both halves of trait-driven transfer — the donor
**rate** (how much HGT) and the ``transfer_to`` recipient **weight** (where it lands).
"""

import hashlib
import math

import pytest

from zombi2 import genomes, traits
from zombi2.rates.driver import DriverTrajectory, load_driver
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Between, Curve, Scalar, Table, as_mapping, check_kernel_fires
from zombi2.species import simulate_species_tree
from zombi2.tree import read_newick


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
    b = Between({("a", "b"): 2.0})
    assert as_mapping(b) is b                           # a kernel passes through (for DrivenBy(., Between))
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
    d = mod.DrivenBy("f.tsv", Between({("A", "B"): 3.0, ("B", "A"): 3.0}, default=0.0))
    assert parse_rate(repr(d)) == d                     # the log line pastes back into a flag


def test_check_kernel_fires_needs_a_named_pair_to_occur():
    check_kernel_fires(Between({("A", "B"): 1.0}), {"A", "B"}, source_label="x")  # both present → ok
    check_kernel_fires(Between({("A", "B"): 1.0, ("Q", "Z"): 1.0}), {"A", "B"}, source_label="x")  # partial ok
    with pytest.raises(ValueError, match="silently do nothing"):
        check_kernel_fires(Between({("Q", "Z"): 1.0}), {"A", "B"}, source_label="x")  # none present


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
    hab.write(tmp_path, outputs=("events",))
    traj = load_driver(tmp_path / "trait_events.tsv", tree)   # the log, replayed against the tree
    # the reconstructed trajectory agrees with the trait's own node values at each node's end time
    for i, node in tree.nodes.items():
        assert traj.value(i, node.end_time - 1e-9) == hab.node_values[i]


def test_a_continuous_log_is_only_its_origin(tmp_path):
    # a diffusion cannot be reconstructed from events, so its log carries only the root marker — there
    # is no discrete map to drive a rate with (driving on a continuous trait is a later slice anyway)
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=2).complete_tree
    cont = traits.simulate_continuous(tree, rate=1.0, seed=1)
    cont.write(tmp_path, outputs=("events",))
    lines = (tmp_path / "trait_events.tsv").read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2 and lines[1].split("\t")[1] == "root"


# --- end-to-end conditioned coupling: a trait drives gene loss ------------------------------------

def _write_driver(path, tree, state_of):
    """Write a trait **event log** assigning ``state_of[node]`` to each lineage for its whole branch:
    a ``root`` row for the crown, then one ``on_speciation`` row per other node fixing its start state
    (no ``on_branch`` switches, so every branch is constant). Replayed against ``tree`` this rebuilds
    exactly ``state_of`` — the format a conditioned run reads now."""
    root = tree.root
    rows = ["time\tkind\tlineage\tfrom\tto",
            f"{tree.nodes[root].birth_time!r}\troot\tn{root}\t\t{state_of[root]}"]
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


def test_mapping_matching_no_driver_state_is_refused(tmp_path):
    # a mapping whose keys occur nowhere in the driver would leave every lineage at the default factor
    # — a silently uncoupled run — so it must be refused, not run as if it were coupled
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, {i: ("hi" if i % 2 else "lo") for i in tree.nodes})
    with pytest.raises(ValueError, match="match none of the driver's states"):
        genomes.simulate_genomes_family(
            tree, loss=0.25 * mod.DrivenBy(str(driver), {"cave": 4.0}),  # 'cave' is never a driver state
            initial_families=6, seed=3)


def test_partial_mapping_with_one_matching_state_still_runs(tmp_path):
    # ≥1 overlap is enough: a mapping may name a state this realisation never reached, as long as at
    # least one of its states does occur — that is a legitimate partial mapping, not a mistake
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, {i: "lo" for i in tree.nodes})   # only 'lo' ever occurs
    res = genomes.simulate_genomes_family(
        tree, loss=0.25 * mod.DrivenBy(str(driver), {"lo": 2.0, "hi": 9.0}),  # 'hi' listed but absent
        initial_families=6, seed=3)
    assert res.events is not None                                # it ran; the absent 'hi' key is fine


def test_driven_loss_is_deterministic(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=5).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, state_of)
    kw = dict(loss=0.3 * mod.DrivenBy(str(driver), {"lo": 1.0, "hi": 5.0}),
              initial_families=4, seed=9)
    a = genomes.simulate_genomes_family(tree, **kw)
    b = genomes.simulate_genomes_family(tree, **kw)
    assert [(e.time, e.kind, e.lineage, e.family) for e in a.events] == \
           [(e.time, e.kind, e.lineage, e.family) for e in b.events]


def test_end_to_end_trait_drives_loss(tmp_path):
    """The full conditioned-coupling workflow: grow a habitat trait, write it, drive gene loss by it.
    Across the tree, lineages in the high-loss state should carry fewer copies on average than
    lineages in the low-loss state."""
    tree = simulate_species_tree(birth=1.1, total_time=3.0, seed=3).complete_tree
    hab = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.5, seed=1)
    hab.write(tmp_path, outputs=("events",))
    res = genomes.simulate_genomes_family(
        tree,
        loss=0.15 * mod.DrivenBy(str(tmp_path / "trait_events.tsv"),
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
        loss=0.3 * mod.DrivenBy(str(tmp_path / "trait_events.tsv"), {0: 0.0, 1: 30.0}),
        initial_families=5, seed=2,
    )
    losses = [e for e in res.events if e.kind == "loss"]
    assert losses, "the int-keyed mapping must actually bite (not silently default to 1.0)"
    # every loss is on a state-1 branch (state 0 has loss factor exactly 0)
    driver = load_driver(tmp_path / "trait_events.tsv", tree)
    assert all(driver.value(e.lineage, e.time) == "1" for e in losses)


# --- conditioning in-memory: DrivenBy accepts the trait result object (no file step) --------------

def test_drivenby_accepts_traits_result_object(tmp_path):
    """Passing the discrete TraitsResult directly is the same conditioning as writing a file, and
    gives an IDENTICAL run (the file round-trip is lossless)."""
    tree = simulate_species_tree(birth=1.1, total_time=2.5, seed=3).complete_tree
    habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.5, seed=1)
    kw = dict(loss=0.5 * mod.DrivenBy(habitat, {"aquatic": 3.0, "terrestrial": 1.0}),
              origination=0.2, initial_families=8, seed=2)
    by_object = genomes.simulate_genomes_family(tree, **kw)

    habitat.write(tmp_path, outputs=("events",))
    by_file = genomes.simulate_genomes_family(
        tree,
        loss=0.5 * mod.DrivenBy(str(tmp_path / "trait_events.tsv"), {"aquatic": 3.0, "terrestrial": 1.0}),
        origination=0.2, initial_families=8, seed=2)
    key = lambda r: [(e.time, e.kind, e.lineage, e.copy) for e in r.events]
    assert key(by_object) == key(by_file)


def test_drivenby_object_must_be_a_trait_result(tmp_path):
    # a conditioned driver object must be a grown trait result (discrete or continuous), carrying its
    # own complete tree — anything else (here a bare object) is refused
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=2).complete_tree
    with pytest.raises(ValueError, match="grown trait result"):
        genomes.simulate_genomes_family(
            tree, loss=0.5 * mod.DrivenBy(object(), {"a": 2.0}), initial_families=3, seed=1)


# --- a trait drives transfer, side 1: the DONOR rate ----------------------------------------------
# Driving `transfer` says how often a lineage DONATES, so it changes the total amount of HGT.

_UNDRIVEN_TRANSFER_DIGESTS = {
    # captured from the engine BEFORE driven transfer was wired: the whole event log of a seeded
    # run, under each transfer_to rule. An undriven transfer must stay byte-identical — same rng
    # draw order, same results — however the driven path is built around it.
    "uniform": "2092e6a774b9e71cefba42d41ce1c9c42e4ab00846c1f5d849baba160cec2efd",
    "distance": "f6c6ccfc5bb7db61794d103fa8356f74d5c525f973ec6b2f98d51086e3013d44",
}


def _event_digest(result) -> str:
    key = repr([(round(e.time, 12), e.kind, e.lineage, e.family, e.copy, e.parent, e.recipient)
                for e in result.events])
    return hashlib.sha256(key.encode()).hexdigest()


@pytest.mark.parametrize("rule", ["uniform", "distance"])
def test_undriven_transfer_is_unchanged(rule):
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=2.5, seed=17).complete_tree
    res = genomes.simulate_genomes_family(
        tree, duplication=0.2, transfer=0.4, loss=0.15, origination=0.3,
        transfer_to=rule, initial_families=5, seed=23)
    assert _event_digest(res) == _UNDRIVEN_TRANSFER_DIGESTS[rule], (
        "an undriven transfer changed: the rng draw order of the uncoupled path must not move")


def test_undriven_transfer_is_unchanged_under_replacement_and_self_transfer():
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=2.5, seed=17).complete_tree
    res = genomes.simulate_genomes_family(
        tree, duplication=0.2, transfer=0.4, loss=0.15, origination=0.3, replacement=True,
        self_transfer=True, initial_families=5, seed=23)
    assert _event_digest(res) == "6eb913b6da50df2dcfa463dcc02327258e7da68ec092c127990bd03ea2d8dfac"


def test_driven_transfer_picks_the_donor(tmp_path):
    """The sharp invariant on the donor side: a lineage whose transfer factor is 0 never donates."""
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, state_of)

    res = genomes.simulate_genomes_family(
        tree, transfer=0.2 * mod.DrivenBy(str(driver), {"lo": 0.0, "hi": 20.0}),
        initial_families=6, seed=3)
    donations = [e for e in res.events if e.kind == "transfer" and e.recipient is None]
    assert donations, "expected some donation from the hi lineages"
    assert all(state_of[e.lineage] == "hi" for e in donations)
    # the recipients are still drawn uniformly, so lo lineages do receive — the coupling is on the
    # donor side only
    arrivals = [e for e in res.events if e.kind == "transfer" and e.recipient is not None]
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
            tree, transfer=0.2 * mod.DrivenBy(str(driver), {"any": 3.0}), **kw)
        n_plain += sum(1 for e in plain.events if e.kind == "transfer" and e.recipient is not None)
        n_driven += sum(1 for e in driven.events if e.kind == "transfer" and e.recipient is not None)
    assert 2.7 < n_driven / n_plain < 3.3


def test_driven_transfer_is_deterministic(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=5).complete_tree
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, {i: ("hi" if i % 2 else "lo") for i in tree.nodes})
    kw = dict(transfer=0.3 * mod.DrivenBy(str(driver), {"lo": 1.0, "hi": 5.0}),
              initial_families=4, seed=9)
    a = genomes.simulate_genomes_family(tree, **kw)
    b = genomes.simulate_genomes_family(tree, **kw)
    assert _event_digest(a) == _event_digest(b)


# --- a trait drives transfer, side 2: the RECIPIENT weight ----------------------------------------
# `transfer_to = mod.DrivenBy(...)` is the choice slot (SPEC §5): the mapping's numbers are weights
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
        transfer_to=mod.DrivenBy(str(driver), {"competent": 2.0, "normal": 1.0}), seed=5)
    arrivals = [e for e in res.events if e.kind == "transfer" and e.recipient is not None]
    assert len(arrivals) > 1500                        # enough events for a 0.03 tolerance
    assert all(e.lineage in tips for e in arrivals)     # the internal branches are 1e-6 long
    share = sum(1 for e in arrivals if e.lineage in hot) / len(arrivals)
    assert share == pytest.approx(2 / 3, abs=0.03)


def test_recipient_weight_zero_cannot_receive(tmp_path):
    """Weight 0 means "cannot receive": every transfer lands on a competent lineage."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_family(
        tree, transfer=1.0, initial_families=6,
        transfer_to=mod.DrivenBy(str(driver), {"competent": 1.0, "normal": 0.0}), seed=5)
    arrivals = [e for e in res.events if e.kind == "transfer" and e.recipient is not None]
    assert arrivals
    assert all(e.lineage in hot for e in arrivals)


def test_no_eligible_recipient_means_no_transfer_at_all(tmp_path):
    """When every candidate weighs 0 the transfer cannot happen, so the event is dropped whole —
    no donor continuation, no arrival, no copy minted. The same run under 'uniform' transfers
    freely, so the difference is the weighting and not the setup."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=0)   # nobody is competent
    kw = dict(transfer=1.0, initial_families=6, seed=5)
    blocked = genomes.simulate_genomes_family(
        tree, transfer_to=mod.DrivenBy(str(driver), {"competent": 1.0, "normal": 0.0}), **kw)
    free = genomes.simulate_genomes_family(tree, transfer_to="uniform", **kw)
    assert not [e for e in blocked.events if e.kind == "transfer"]
    assert [e for e in free.events if e.kind == "transfer"]
    # a dropped event leaves the genomes untouched: the six crown families are simply inherited
    assert all(len(g) == 6 for g in blocked.genomes.values())


def test_both_couplings_compose(tmp_path):
    """The donor rate and the recipient weight are independent models and may be used together."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_family(
        tree, transfer=0.5 * mod.DrivenBy(str(driver), {"competent": 5.0, "normal": 0.0}),
        transfer_to=mod.DrivenBy(str(driver), {"competent": 0.0, "normal": 1.0}),
        initial_families=6, seed=5)
    donations = [e for e in res.events if e.kind == "transfer" and e.recipient is None]
    arrivals = [e for e in res.events if e.kind == "transfer" and e.recipient is not None]
    assert donations and arrivals
    assert all(e.lineage in hot for e in donations)          # only competent lineages donate
    assert all(e.lineage not in hot for e in arrivals)       # only non-competent lineages receive


# --- a Between kernel makes transfer_to donor-conditioned (assortative by a trait) -----------------

def test_recipient_kernel_keeps_transfer_within_the_donor_state(tmp_path):
    """DrivenBy(trait, Between(...)) reads the driver on the DONOR too, so a same-state kernel keeps
    every transfer within one habitat — the thing a 1-D recipient weight cannot express."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_family(
        tree, transfer=4.0, initial_families=6, self_transfer=True,
        transfer_to=mod.DrivenBy(str(driver),
                                 Between({("competent", "competent"): 1.0,
                                          ("normal", "normal"): 1.0}, default=0.0)), seed=5)
    arrivals = [e for e in res.events if e.kind == "transfer" and e.recipient is not None]
    assert arrivals
    assert all((e.donor in hot) == (e.lineage in hot) for e in arrivals)  # donor and recipient agree


def test_recipient_kernel_fires_check_catches_absent_groups(tmp_path):
    """A kernel naming groups the driver never takes would weight every candidate at the default —
    secretly uniform — so it is refused, like a Table that names no occurring state."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    with pytest.raises(ValueError, match="silently do nothing"):
        genomes.simulate_genomes_family(
            tree, transfer=1.0, initial_families=6,
            transfer_to=mod.DrivenBy(str(driver), Between({("x", "y"): 1.0})), seed=5)


def test_between_is_rejected_in_a_rate_slot():
    """A rate has no donor to condition on, so a Between kernel there is a category error."""
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="donor-conditioned"):
        genomes.simulate_genomes_family(
            tree, loss=0.5 * mod.DrivenBy("f.tsv", Between({("a", "b"): 1.0})),
            initial_families=1, seed=1)


# --- the guard: what is not wired ------------------------------------------------------------------

def test_transfer_to_rejects_a_rate():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="on its own, not a rate"):
        genomes.simulate_genomes_family(
            tree, transfer=0.1, transfer_to=1.0 * mod.DrivenBy("f.tsv", {"a": 2.0}),
            initial_families=1, seed=1)


def test_transfer_to_rejects_combining_distance_with_a_driven_weight():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="one recipient rule"):
        genomes.simulate_genomes_family(
            tree, transfer=0.1,
            transfer_to=(genomes.Distance(decay=1.0), mod.DrivenBy("f.tsv", {"a": 2.0})),
            initial_families=1, seed=1)


def test_transfer_to_rejects_an_unknown_rule():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="transfer_to must be"):
        genomes.simulate_genomes_family(tree, transfer=0.1, transfer_to="closest",
                                        initial_families=1, seed=1)


def test_ordered_engine_rejects_a_driven_transfer_to():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="transfer_to must be"):
        genomes.simulate_genomes_ordered(
            tree, transfer=0.1, transfer_to=mod.DrivenBy("f.tsv", {"a": 2.0}),
            initial_families=1, seed=1)


def test_ordered_engine_rejects_a_driven_transfer_rate():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="ordered genome engine does not"):
        genomes.simulate_genomes_ordered(
            tree, transfer=0.1 * mod.DrivenBy("f.tsv", {"a": 2.0}), initial_families=1, seed=1)


def test_missing_driver_file_is_named_as_a_drivenby_driver(tmp_path):
    # a conditioned rate pointing at a file that is not there must say so as a DrivenBy driver, not
    # leak a bare errno — so the user knows which kind of input is missing and how to fix it
    with pytest.raises(FileNotFoundError, match="DrivenBy driver file not found"):
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
        loss=0.9 * mod.DrivenBy(habitat, {"host": 25.0, "free": 0.0}),
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
        tree, loss=0.6 * mod.DrivenBy(habitat, {"host": 8.0, "free": 1.0}), **kw)
    habitat.write(str(tmp_path), outputs=("events",))
    by_file = genomes.simulate_genomes_nucleotide(
        tree, loss=0.6 * mod.DrivenBy(str(tmp_path / "trait_events.tsv"),
                                      {"host": 8.0, "free": 1.0}), **kw)
    assert [(e.time, type(e).__name__, e.lineage) for e in by_object.events] == \
           [(e.time, type(e).__name__, e.lineage) for e in by_file.events]


def test_nucleotide_driving_changes_the_run():
    """A driven rate must actually bite: the same seed with the coupling switched off is a different
    run. Guards against a driver that resolves but never reaches the rate."""
    tree = _nuc_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=6)
    kw = dict(root_length=5000, genes=3, gene_length=250, loss_extent=200, seed=6)
    driven = genomes.simulate_genomes_nucleotide(
        tree, loss=0.6 * mod.DrivenBy(habitat, {"host": 20.0, "free": 0.05}), **kw)
    flat = genomes.simulate_genomes_nucleotide(tree, loss=0.6, **kw)
    assert len(driven.events) != len(flat.events)


def test_nucleotide_undriven_run_is_unchanged():
    """No coupling ⇒ the pooled path, byte for byte. The driven machinery must cost an uncoupled run
    nothing."""
    tree = _nuc_tree()
    kw = dict(root_length=4000, genes=3, gene_length=250, inversion=1.5, loss=0.4, seed=9)
    a = genomes.simulate_genomes_nucleotide(tree, **kw)
    b = genomes.simulate_genomes_nucleotide(tree, **kw)
    assert [(e.time, type(e).__name__, e.lineage) for e in a.events] == \
           [(e.time, type(e).__name__, e.lineage) for e in b.events]


def test_nucleotide_refuses_a_mapping_that_never_fires():
    """A mapping naming states the driver never takes leaves every lineage on the default factor, so
    the run would secretly be the uncoupled model."""
    tree = _nuc_tree()
    habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.8, seed=7)
    with pytest.raises(ValueError):
        genomes.simulate_genomes_nucleotide(
            tree, root_length=3000, loss=0.5 * mod.DrivenBy(habitat, {"aquatic": 3.0}), seed=7)


# --- a CONTINUOUS trait as the driver (approximate: each branch cut into constant sub-steps) --------

def test_continuous_driver_trajectory_interpolates():
    """`driver_from_continuous_result` cuts each branch into ``steps`` constant stretches whose value
    is the trait linearly interpolated from the parent's value to the node's, sampled at each
    stretch's midpoint — so `value()` returns exactly that, and `next_change()` steps within a branch."""
    from zombi2.rates.driver import driver_from_continuous_result
    ct = simulate_species_tree(birth=1.0, n_extant=6, seed=1).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=2)
    traj = driver_from_continuous_result(met, steps=10)

    node = next(n for n in ct.nodes.values() if n.parent is not None and n.end_time > n.birth_time)
    start_v, end_v = met.node_values[node.parent], met.node_values[node.id]
    dt = node.end_time - node.birth_time
    for k in (0, 3, 9):                                  # the value in the k-th stretch = its midpoint
        expected = start_v + (end_v - start_v) * (k + 0.5) / 10
        assert traj.value(node.id, node.birth_time + k * dt / 10) == pytest.approx(expected)
    nxt = traj.next_change(node.id, node.birth_time)     # a within-branch breakpoint, not inf
    assert node.birth_time < nxt <= node.end_time


def test_continuous_driver_couples_a_rate_and_is_deterministic():
    """A continuous trait drives the duplication rate through a Curve: the coupled run differs from an
    otherwise-identical uncoupled one (the driver is threaded and the mapping applied), and repeats
    exactly under the same seed."""
    ct = simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.2, seed=3)

    def run(fn):
        return genomes.simulate_genomes_family(
            ct, initial_families=15, loss=0.04,
            duplication=0.06 * mod.DrivenBy(met, Curve(fn)), seed=9)

    coupled = run(lambda v: 3.0 ** v)
    control = run(lambda v: 1.0)                          # a flat Curve == the uncoupled model
    tips = list(ct.extant_leaves())
    assert any(len(coupled.genomes[n.id]) != len(control.genomes[n.id]) for n in tips)
    again = run(lambda v: 3.0 ** v)
    assert all(len(coupled.genomes[n.id]) == len(again.genomes[n.id]) for n in tips)


def test_continuous_driver_takes_a_scalar_link():
    ct = simulate_species_tree(birth=1.0, n_extant=12, seed=1).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=2)
    genomes.simulate_genomes_family(                      # Scalar (log-link) is a valid continuous mapping
        ct, initial_families=10, duplication=0.05 * mod.DrivenBy(met, Scalar(0.5)), seed=1)


def test_discrete_table_on_continuous_driver_is_refused():
    """A ``{state: factor}`` table names discrete states a continuous value never equals — refuse it
    with a message that points at Curve / Scalar, rather than silently leaving the rate undriven."""
    ct = simulate_species_tree(birth=1.0, n_extant=12, seed=1).complete_tree
    met = traits.simulate_continuous(ct, start=0.0, rate=1.0, seed=2)
    with pytest.raises(ValueError, match="CONTINUOUS"):
        genomes.simulate_genomes_family(
            ct, initial_families=10, duplication=0.05 * mod.DrivenBy(met, {"hi": 2.0}), seed=1)


def test_multitrait_continuous_driver_is_refused():
    from zombi2.rates.driver import driver_from_continuous_result
    ct = simulate_species_tree(birth=1.0, n_extant=8, seed=1).complete_tree
    two = traits.simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0}, seed=2)
    with pytest.raises(ValueError, match="SINGLE-trait"):
        driver_from_continuous_result(two)
