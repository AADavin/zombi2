"""Replacement as a probability (issue #476): a transfer arriving where its family already has a copy
overwrites one resident with probability p and is added beside it otherwise, one Bernoulli draw per
transfer. At 0 and 1 the runs are byte-identical to the additive and replacing runs they replace."""

import collections
import hashlib

import pytest

from zombi2.genomes import family, group, simulate_genomes_family, simulate_genomes_ordered
from zombi2.species import simulate_species_tree


def _digest(r) -> str:
    return hashlib.sha256(repr([(e.time, e.kind, e.family, e.copy, e.parent, e.replaced)
                                for e in r.edges]).encode()).hexdigest()[:12]


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1).complete_tree


_RUN = dict(duplication=0.2, transfer=0.6, loss=0.3, origination=0.5, initial_families=30, seed=5)

#: the digests of these runs on main before #476, taken from an export of 442b17ff: additive and
#: replacing, on the three engines. A run at p = 0 or p = 1 must give exactly these.
_BEFORE = {("family", 0.0): "11ba08c736cd", ("family", 1.0): "36d71d2a8631",
           ("parallel", 0.0): "b44434f1e9aa", ("parallel", 1.0): "7144933d8212",
           ("ordered", 0.0): "efedbf271687", ("ordered", 1.0): "e1f7e2d0d2ca"}


def _engine(label):
    if label == "family":
        return simulate_genomes_family, {}
    if label == "parallel":
        return simulate_genomes_family, {"parallel": True}
    return simulate_genomes_ordered, {"transfer_extent": 2}


@pytest.mark.parametrize("label", ["family", "parallel", "ordered"])
def test_zero_and_one_are_the_additive_and_replacing_runs_of_before(tree, label):
    fn, extra = _engine(label)
    for p, as_bool in ((0.0, False), (1.0, True)):
        assert _digest(fn(tree, replacement=p, **_RUN, **extra)) == _BEFORE[(label, p)]
        assert _digest(fn(tree, replacement=as_bool, **_RUN, **extra)) == _BEFORE[(label, p)]


@pytest.mark.parametrize("label", ["family", "parallel", "ordered"])
def test_the_realised_share_of_replacing_transfers_matches_p(label):
    """Among the transfers that found a resident copy, the replacing ones are Binomial(N, p)."""
    fn, extra = _engine(label)
    t = simulate_species_tree(birth=1.0, death=0.0, n_extant=30, seed=2).complete_tree
    p = 0.3
    r = fn(t, duplication=0.1, transfer=0.5, loss=0.6, origination=1.0, initial_families=60,
           replacement=p, seed=3, **extra)
    c = r.replacements
    n, k = c["found_resident"], c["replaced"]
    assert n >= 200 and c["added_beside"] == n - k
    se = (p * (1 - p) / n) ** 0.5
    assert abs(k / n - p) < 4 * se, (k / n, p, n)
    # the replaced count is the count of replacing events in the log (one event is one draw; at the
    # ordered resolution one event is as many rows as families the segment carried)
    replacing_events = {(e.time, e.kind) for e in r.events if e.kind == "transfer_replacing"}
    assert len(replacing_events) == k
    assert r.summary()["transfers"] == c


def test_at_zero_nothing_is_counted_and_the_summary_says_nothing(tree):
    r = simulate_genomes_family(tree, **_RUN)
    assert r.replacements is None
    assert "transfers" not in r.summary()


def test_the_probability_is_checked(tree):
    with pytest.raises(ValueError, match="probability"):
        simulate_genomes_family(tree, replacement=1.5, **_RUN)
    with pytest.raises(ValueError, match="probability"):
        simulate_genomes_ordered(tree, replacement=-0.1, **_RUN)
    with pytest.raises(ValueError, match="replacement above 0"):
        simulate_genomes_ordered(tree, replacement=0.2, families=[family("A"), family("B")],
                                 groups=[group("op", ["A", "B"])], **_RUN)
    # p = 0 with groups is fine: nothing can be torn out of a unit
    simulate_genomes_ordered(tree, replacement=0.0, families=[family("A"), family("B")],
                             groups=[group("op", ["A", "B"])], **_RUN)


def test_the_flag_takes_a_probability_or_stands_alone_for_one(tmp_path, monkeypatch):
    from zombi2.cli import main as cli_main
    seen = {}
    def stub(args, parser):
        seen["p"] = args.replacement
        return 0
    monkeypatch.setattr(cli_main, "_RUN", {k: stub for k in cli_main._RUN})
    monkeypatch.chdir(tmp_path)
    for tokens, expected in ((["--replacement", "0.5"], 0.5), (["--replacement"], 1.0), ([], 0.0)):
        seen.clear()
        assert cli_main.main(["genomes", "out", *tokens, "--seed", "1"]) == 0
        assert seen["p"] == expected
    with pytest.raises(SystemExit):
        cli_main.main(["genomes", "out", "--replacement", "2", "--seed", "1"])


def test_the_counts_add_up_on_the_ordered_resolution_per_event(tree):
    r = simulate_genomes_ordered(tree, replacement=0.5, transfer_extent=3, **_RUN)
    c = r.replacements
    by_event = collections.Counter((e.time, e.kind) for e in r.events if e.kind == "transfer_replacing")
    assert len(by_event) == c["replaced"]
    assert max(by_event.values()) >= 1
