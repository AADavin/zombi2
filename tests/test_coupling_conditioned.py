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
from zombi2.rates.mapping import Curve, Scalar, Table, as_mapping
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
    lines = (tmp_path / "trait_events.tsv").read_text().splitlines()
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


def test_mapping_matching_no_driver_state_is_refused(tmp_path):
    # a mapping whose keys occur nowhere in the driver would leave every lineage at the default factor
    # — a silently uncoupled run — so it must be refused, not run as if it were coupled
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, {i: ("hi" if i % 2 else "lo") for i in tree.nodes})
    with pytest.raises(ValueError, match="match none of the driver's states"):
        genomes.simulate_genomes_unordered(
            tree, loss=0.25 * mod.DrivenBy(str(driver), {"cave": 4.0}),  # 'cave' is never a driver state
            initial_families=6, seed=3)


def test_partial_mapping_with_one_matching_state_still_runs(tmp_path):
    # ≥1 overlap is enough: a mapping may name a state this realisation never reached, as long as at
    # least one of its states does occur — that is a legitimate partial mapping, not a mistake
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    driver = tmp_path / "habitat.tsv"
    _write_driver(driver, tree, {i: "lo" for i in tree.nodes})   # only 'lo' ever occurs
    res = genomes.simulate_genomes_unordered(
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
    a = genomes.simulate_genomes_unordered(tree, **kw)
    b = genomes.simulate_genomes_unordered(tree, **kw)
    assert [(e.time, e.kind, e.lineage, e.family) for e in a.events] == \
           [(e.time, e.kind, e.lineage, e.family) for e in b.events]


def test_end_to_end_trait_drives_loss(tmp_path):
    """The full conditioned-coupling workflow: grow a habitat trait, write it, drive gene loss by it.
    Across the tree, lineages in the high-loss state should carry fewer copies on average than
    lineages in the low-loss state."""
    tree = simulate_species_tree(birth=1.1, total_time=3.0, seed=3).complete_tree
    hab = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.5, seed=1)
    hab.write(tmp_path, outputs=("events",))
    res = genomes.simulate_genomes_unordered(
        tree,
        loss=0.15 * mod.DrivenBy(str(tmp_path / "trait_events.tsv"),
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
    trait.write(tmp_path, outputs=("events",))
    res = genomes.simulate_genomes_unordered(
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
    by_object = genomes.simulate_genomes_unordered(tree, **kw)

    habitat.write(tmp_path, outputs=("events",))
    by_file = genomes.simulate_genomes_unordered(
        tree,
        loss=0.5 * mod.DrivenBy(str(tmp_path / "trait_events.tsv"), {"aquatic": 3.0, "terrestrial": 1.0}),
        origination=0.2, initial_families=8, seed=2)
    key = lambda r: [(e.time, e.kind, e.lineage, e.copy) for e in r.events]
    assert key(by_object) == key(by_file)


def test_drivenby_object_must_be_discrete(tmp_path):
    tree = simulate_species_tree(birth=1.0, total_time=1.5, seed=2).complete_tree
    cont = traits.simulate_continuous(tree, rate=1.0, seed=1)   # a diffusion has no stochastic map
    with pytest.raises(ValueError, match="DISCRETE"):
        genomes.simulate_genomes_unordered(
            tree, loss=0.5 * mod.DrivenBy(cont, {"a": 2.0}), initial_families=3, seed=1)


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
    res = genomes.simulate_genomes_unordered(
        tree, duplication=0.2, transfer=0.4, loss=0.15, origination=0.3,
        transfer_to=rule, initial_families=5, seed=23)
    assert _event_digest(res) == _UNDRIVEN_TRANSFER_DIGESTS[rule], (
        "an undriven transfer changed: the rng draw order of the uncoupled path must not move")


def test_undriven_transfer_is_unchanged_under_replacement_and_self_transfer():
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=2.5, seed=17).complete_tree
    res = genomes.simulate_genomes_unordered(
        tree, duplication=0.2, transfer=0.4, loss=0.15, origination=0.3, replacement=True,
        self_transfer=True, initial_families=5, seed=23)
    assert _event_digest(res) == "6eb913b6da50df2dcfa463dcc02327258e7da68ec092c127990bd03ea2d8dfac"


def test_driven_transfer_picks_the_donor(tmp_path):
    """The sharp invariant on the donor side: a lineage whose transfer factor is 0 never donates."""
    tree = simulate_species_tree(birth=1.2, death=0.2, total_time=1.5, seed=11).complete_tree
    state_of = {i: ("hi" if i % 2 else "lo") for i in tree.nodes}
    driver = tmp_path / "d.tsv"
    _write_driver(driver, tree, state_of)

    res = genomes.simulate_genomes_unordered(
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
        plain = genomes.simulate_genomes_unordered(tree, transfer=0.2, **kw)
        driven = genomes.simulate_genomes_unordered(
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
    a = genomes.simulate_genomes_unordered(tree, **kw)
    b = genomes.simulate_genomes_unordered(tree, **kw)
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
    res = genomes.simulate_genomes_unordered(
        tree, transfer=4.0, initial_families=6, self_transfer=True,
        transfer_to=mod.DrivenBy(str(driver), {"competent": 2.0, "normal": 1.0}), seed=5)
    arrivals = [e for e in res.events if e.kind == "transfer" and e.recipient is not None]
    assert len(arrivals) > 1500                        # enough events for a 0.03 tolerance
    assert all(e.lineage in tips for e in arrivals)     # the internal branches are 1e-6 long
    share = sum(1 for e in arrivals if e.lineage in hot) / len(arrivals)
    assert share == pytest.approx(2 / 3, abs=0.03)


def test_recipient_weight_zero_cannot_receive(tmp_path):
    """Weight 0 means "cannot receive": every transfer lands on a competent lineage."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_unordered(
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
    blocked = genomes.simulate_genomes_unordered(
        tree, transfer_to=mod.DrivenBy(str(driver), {"competent": 1.0, "normal": 0.0}), **kw)
    free = genomes.simulate_genomes_unordered(tree, transfer_to="uniform", **kw)
    assert not [e for e in blocked.events if e.kind == "transfer"]
    assert [e for e in free.events if e.kind == "transfer"]
    # a dropped event leaves the genomes untouched: the six crown families are simply inherited
    assert all(len(g) == 6 for g in blocked.genomes.values())


def test_both_couplings_compose(tmp_path):
    """The donor rate and the recipient weight are independent models and may be used together."""
    tree, tips, hot, driver = _flat_tree_and_driver(tmp_path, competent=4)
    res = genomes.simulate_genomes_unordered(
        tree, transfer=0.5 * mod.DrivenBy(str(driver), {"competent": 5.0, "normal": 0.0}),
        transfer_to=mod.DrivenBy(str(driver), {"competent": 0.0, "normal": 1.0}),
        initial_families=6, seed=5)
    donations = [e for e in res.events if e.kind == "transfer" and e.recipient is None]
    arrivals = [e for e in res.events if e.kind == "transfer" and e.recipient is not None]
    assert donations and arrivals
    assert all(e.lineage in hot for e in donations)          # only competent lineages donate
    assert all(e.lineage not in hot for e in arrivals)       # only non-competent lineages receive


# --- the guard: what is not wired ------------------------------------------------------------------

def test_transfer_to_rejects_a_rate():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="on its own, not a rate"):
        genomes.simulate_genomes_unordered(
            tree, transfer=0.1, transfer_to=1.0 * mod.DrivenBy("f.tsv", {"a": 2.0}),
            initial_families=1, seed=1)


def test_transfer_to_rejects_combining_distance_with_a_driven_weight():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="one recipient rule"):
        genomes.simulate_genomes_unordered(
            tree, transfer=0.1,
            transfer_to=(genomes.Distance(decay=1.0), mod.DrivenBy("f.tsv", {"a": 2.0})),
            initial_families=1, seed=1)


def test_transfer_to_rejects_an_unknown_rule():
    tree = simulate_species_tree(birth=1.0, total_time=1.0, seed=1).complete_tree
    with pytest.raises(ValueError, match="transfer_to must be"):
        genomes.simulate_genomes_unordered(tree, transfer=0.1, transfer_to="closest",
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
