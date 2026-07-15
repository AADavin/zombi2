"""Per-species-branch event summary (``branch_events.tsv`` / ``--write branch_events``).

:func:`branch_events_table` aggregates the event log into one row per species-tree branch: the
count of each event that fired on it (D/T/L/O, plus inversion/transposition for ordered genomes),
with transfers split into ``transfer_out`` (donor) and ``transfer_in`` (recipient) and an
``is_extant`` flag derived from node times, so the extant-tree view is a filter on the table.

The aggregation-logic tests use the pure-Python engine (per-genome / ordered rates) so they run
with or without the Rust extension; the end-to-end CLI test exercises the default Rust path.
"""

import pytest

import zombi2 as z
from zombi2.genomes.events import EventLog, EventType
from zombi2.genomes.simulation import (BRANCH_EVENTS_HEADER, _extant_branches,
                                       branch_events_table)
from zombi2.cli import main
from zombi2.tree import read_newick


def _tree(n=25, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n, age=5.0, seed=seed)


def _parse(text):
    lines = text.strip().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]
    return header, rows


def _logcount(log, ev):
    return sum(1 for e in log if e.event is ev)


def test_header_and_one_row_per_branch():
    tree = _tree()
    g = z.simulate_genomes(tree, rates=z.PerGenomeRates(0.2, 0.2, 0.2, 0.3),
                           initial_families=20, seed=3)
    header, rows = _parse(branch_events_table(g.event_log, tree))
    assert "\t".join(header) == BRANCH_EVENTS_HEADER
    assert [r["branch"] for r in rows] == [n.name for n in tree.nodes_preorder()]  # every branch


def test_transfer_in_equals_transfer_out_equals_transfers():
    # every transfer is counted once as out (on the donor) and once as in (on the recipient)
    tree = _tree()
    g = z.simulate_genomes(tree, rates=z.PerGenomeRates(0.1, 0.5, 0.15, 0.3),
                           initial_families=20, seed=5)
    _, rows = _parse(branch_events_table(g.event_log, tree))
    tin = sum(int(r["transfer_in"]) for r in rows)
    tout = sum(int(r["transfer_out"]) for r in rows)
    assert tin == tout == _logcount(g.event_log, EventType.TRANSFER)
    assert tout > 0  # the seed actually produced transfers


def test_column_sums_match_event_log_and_total_excludes_transfer_in():
    tree = _tree()
    g = z.simulate_genomes(tree, rates=z.PerGenomeRates(0.25, 0.2, 0.25, 0.4),
                           initial_families=20, seed=9)
    _, rows = _parse(branch_events_table(g.event_log, tree))
    for col, ev in [("origination", EventType.ORIGINATION), ("duplication", EventType.DUPLICATION),
                    ("loss", EventType.LOSS), ("transfer_out", EventType.TRANSFER)]:
        assert sum(int(r[col]) for r in rows) == _logcount(g.event_log, ev)
    for r in rows:  # total = events that fired ON the branch (transfer_in fired on the donor)
        fired = ("origination", "duplication", "transfer_out", "loss", "inversion", "transposition")
        assert int(r["total"]) == sum(int(r[c]) for c in fired)


def test_ordered_model_counts_rearrangements():
    tree = _tree(n=20, seed=2)
    g = z.simulate_genomes(
        tree, rates=z.SharedRates(0.1, 0.0, 0.1, 0.0, inversion=0.4, transposition=0.4),
        initial_families=15, seed=4,
        genome_factory=lambda ids: z.OrderedGenome(ids, extension=0.7))
    _, rows = _parse(branch_events_table(g.event_log, tree))
    assert sum(int(r["inversion"]) for r in rows) == _logcount(g.event_log, EventType.INVERSION)
    assert (sum(int(r["transposition"]) for r in rows)
            == _logcount(g.event_log, EventType.TRANSPOSITION))
    assert sum(int(r["transposition"]) for r in rows) > 0


def test_extant_flag_derived_from_node_times():
    # b dies before the present (time 1.4 < total_age 2), so it is not on the extant tree; its
    # extant sibling a keeps the shared ancestor extant. Derived from times, not any stored flag.
    tree = read_newick("((a:1,b:0.4):1,c:2);")
    extant = _extant_branches(tree)
    assert {"a", "c"} <= extant
    assert "b" not in extant
    _, rows = _parse(branch_events_table(EventLog(), tree))  # flag is independent of the events
    by_name = {r["branch"]: r for r in rows}
    assert by_name["a"]["is_extant"] == "1" and by_name["c"]["is_extant"] == "1"
    assert by_name["b"]["is_extant"] == "0"


def test_ultrametric_tree_all_branches_extant():
    tree = _tree()  # backward/complete-sampling tree: every tip reaches the present
    _, rows = _parse(branch_events_table(EventLog(), tree))
    assert all(r["is_extant"] == "1" for r in rows)


@pytest.mark.skipif(not z.rust_available(), reason="zombi2_core (Rust extension) not built")
def test_cli_writes_branch_events(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--tips", "15", "--seed", "1", "-o", str(sp)])
    out = tmp_path / "g"
    rc = main(["genomes", "-t", str(sp / "species_tree.nwk"),
               "--dup", "0.2", "--trans", "0.3", "--loss", "0.25", "--orig", "0.3",
               "--seed", "7", "--write", "branch_events", "transfers", "-o", str(out)])
    assert rc == 0
    text = (out / "branch_events.tsv").read_text()
    assert text.splitlines()[0] == BRANCH_EVENTS_HEADER
    _, rows = _parse(text)
    # cross-check the transfer split against transfers.tsv
    n_transfers = len((out / "transfers.tsv").read_text().strip().splitlines()) - 1
    assert sum(int(r["transfer_out"]) for r in rows) == n_transfers
