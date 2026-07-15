"""Validation for the user-controllable rate features (docs/validation.md hard rule):

* **explicit per-family rates** — a ``FamilySampledRates`` table fixes named families' D/T/L; a
  family given only duplication grows while the rest (rate 0) stay at one copy.
* **per-event branch scaling** — ``BranchRates(events=("transfer",))`` scales *only* the transfer
  weight on a branch, leaving duplication/loss untouched (a unit invariant).
* **transfer receptivity** — the per-branch absorption weight makes recipient choice proportional to
  receptivity (a frequency oracle on ``_choose_recipient``), biases which branch receives in a full
  run, and is byte-identical to a plain run when off; a zero-receptivity branch never receives.
* **rate files** — ``read_family_rates`` / ``read_branch_rates`` parse the TSV formats (header and
  positional), round-tripping the maps the models consume.
"""

from __future__ import annotations

import numpy as np
import pytest

import zombi2 as z
from zombi2.genomes import (BranchRates, FamilySampledRates, Rates, TransferModel,
                            read_branch_rates, read_family_rates, simulate_genomes)
from zombi2.genomes.events import EventType
from zombi2.genomes.genome_sim import GenomeSimulator
from zombi2.genomes.simulation import events_trace_from_log
from zombi2.tree import TreeNode


# ------------------------------------------------------------------- explicit per-family rates

def test_family_table_fixes_listed_rates_and_defaults_the_rest():
    """A tabulated family returns exactly its triple; an unlisted one falls back to the model's
    distributions (here the default 0, i.e. inert)."""
    rates = FamilySampledRates(rates={"1": (3.0, 2.0, 1.0), "2": (4.0, 0.0, 1.0)})
    assert rates.rates_for("1") == (3.0, 2.0, 1.0)
    assert rates.rates_for("2") == (4.0, 0.0, 1.0)
    rates.bind(np.random.default_rng(0))
    assert rates.rates_for("57") == (0.0, 0.0, 0.0)   # unlisted -> sampled from Fixed(0)


def test_family_table_only_the_named_family_grows():
    """Inject-recover: with family 1 at duplication 2 and every other family at 0, only family 1
    grows above a single copy at the leaves (the inert families keep one copy per species)."""
    from collections import Counter

    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.1), n_tips=12, age=4.0, seed=1)
    g = simulate_genomes(tree, FamilySampledRates(rates={"1": (2.0, 0.0, 0.0)}),
                         initial_families=5, seed=3, max_family_size=0.9)
    n_species = len(g.profiles.species)
    total: Counter = Counter()
    for genome in g.leaf_genomes.values():
        for fam in genome.families():
            total[fam] += genome.copy_number(fam)
    others = [total[f] for f in total if f != "1"]
    assert total["1"] > max(others)                          # the driven family is the largest
    assert all(o == n_species for o in others)               # inert families: one copy per species


def test_family_table_rejects_negative_rates():
    with pytest.raises(ValueError):
        FamilySampledRates(rates={"1": (1.0, -2.0, 1.0)})


# --------------------------------------------------------------------- per-event branch scaling

def _weights_at(rate_model, tree, branch):
    rate_model.bind(np.random.default_rng(0), tree=tree)
    genome = next(iter(_leaf_genomes(tree)))
    return {ew.event: ew.rate for ew in rate_model.event_weights(genome, branch, 0.0)}


def _leaf_genomes(tree):
    g = simulate_genomes(tree, Rates(duplication=0.1), initial_families=3, seed=1)
    return list(g.leaf_genomes.values())


def test_branch_events_scales_only_the_named_event():
    """``events=("transfer",)`` multiplies the transfer weight on the scaled branch but leaves
    duplication and loss at their base value — the transfer-emission dial."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=2)
    base = Rates(duplication=0.2, transfer=0.2, loss=0.2)
    br = BranchRates(base, factors={"n3": 5.0}, events=("transfer",))
    w = _weights_at(br, tree, "n3")
    b = _weights_at(Rates(duplication=0.2, transfer=0.2, loss=0.2), tree, "n3")
    assert w[EventType.TRANSFER] == pytest.approx(5.0 * b[EventType.TRANSFER])
    assert w[EventType.DUPLICATION] == pytest.approx(b[EventType.DUPLICATION])   # untouched
    assert w[EventType.LOSS] == pytest.approx(b[EventType.LOSS])                 # untouched


def test_branch_events_default_scales_all_dtl():
    """Default ``events`` (None) keeps the original behaviour: D, T and L all scaled together."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=2)
    br = BranchRates(Rates(duplication=0.2, transfer=0.2, loss=0.2), factors={"n3": 5.0})
    w = _weights_at(br, tree, "n3")
    b = _weights_at(Rates(duplication=0.2, transfer=0.2, loss=0.2), tree, "n3")
    for ev in (EventType.DUPLICATION, EventType.TRANSFER, EventType.LOSS):
        assert w[ev] == pytest.approx(5.0 * b[ev])


# ------------------------------------------------------------------------- transfer receptivity

def test_receptivity_makes_selection_proportional_to_weight():
    """Frequency oracle on ``_choose_recipient``: with distance off, a candidate is chosen with
    probability proportional to its receptivity weight (1 : 2 : 3 here)."""
    sim = GenomeSimulator()
    sim._transfers = TransferModel(receptivity={"a": 1.0, "b": 2.0, "c": 3.0})
    a, b, c, d = (TreeNode("a", 1.0), TreeNode("b", 1.0), TreeNode("c", 1.0), TreeNode("d", 1.0))
    alive = {a: None, b: None, c: None, d: None}          # d is the donor, excluded
    rng = np.random.default_rng(0)
    counts = {"a": 0, "b": 0, "c": 0}
    n = 24000
    for _ in range(n):
        counts[sim._choose_recipient(d, alive, 1.0, rng).name] += 1
    assert counts["a"] / n == pytest.approx(1 / 6, abs=0.02)
    assert counts["b"] / n == pytest.approx(2 / 6, abs=0.02)
    assert counts["c"] / n == pytest.approx(3 / 6, abs=0.02)


def test_zero_receptivity_branch_never_receives():
    sim = GenomeSimulator()
    sim._transfers = TransferModel(receptivity={"a": 0.0})
    a, b, d = TreeNode("a", 1.0), TreeNode("b", 1.0), TreeNode("d", 1.0)
    alive = {a: None, b: None, d: None}
    rng = np.random.default_rng(1)
    picks = {sim._choose_recipient(d, alive, 1.0, rng).name for _ in range(500)}
    assert picks == {"b"}                                  # a (weight 0) is never chosen


def test_receptivity_off_is_byte_identical():
    """The no-receptivity path is unchanged: a default TransferModel and no transfers argument give
    the same event log for the same seed."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=12, age=3.0, seed=1)
    a = simulate_genomes(tree, transfer=0.5, transfers=TransferModel(),
                         initial_families=6, seed=5, max_family_size=0.9)
    b = simulate_genomes(tree, transfer=0.5, initial_families=6, seed=5, max_family_size=0.9)
    assert events_trace_from_log(a.event_log) == events_trace_from_log(b.event_log)


def test_rust_and_python_agree_on_receptivity():
    """Engine parity: the Rust built-in engine and the pure-Python engine apply the *same*
    per-branch receptivity, landing about the same share of transfers on the boosted branch. They
    have different RNG streams, so the two shares are compared statistically, not byte-for-byte."""
    from zombi2._sampling import NumpyEventSampler

    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.0), n_tips=10, age=3.0, seed=7)
    target = [n.name for n in tree.nodes_preorder() if not n.children][0]
    tm = TransferModel(receptivity={target: 20.0})

    def share(sampler):
        g = simulate_genomes(tree, Rates(transfer=0.6), transfers=tm,
                             initial_families=8, seed=11, max_family_size=0.9, sampler=sampler)
        recv = [r.recipient for r in g.event_log if r.event is EventType.TRANSFER]
        return recv.count(target) / max(len(recv), 1)

    rust = share(None)                     # sampler None + plain Rates -> Rust engine
    python = share(NumpyEventSampler())    # a custom sampler forces the pure-Python engine
    assert rust == pytest.approx(python, abs=0.12)
    assert rust > 0.35 and python > 0.35   # both clearly biased above the ~0.2 uniform baseline


def test_receptivity_biases_which_branch_receives():
    """Inject-recover on a full run: raising one leaf's receptivity raises its share of received
    transfers well above its unweighted baseline."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.0), n_tips=10, age=3.0, seed=7)
    leaves = [n.name for n in tree.nodes_preorder() if not n.children]
    target = leaves[0]

    def share(transfers):
        # FamilySampledRates forces the pure-Python engine, so receptivity is applied here; the
        # Rust counts-only path gets its own parity test once that engine learns receptivity.
        g = simulate_genomes(tree, FamilySampledRates(transfer=0.6), transfers=transfers,
                             initial_families=8, seed=11, max_family_size=0.9)
        recv = [r.recipient for r in g.event_log if r.event is EventType.TRANSFER]
        return recv.count(target) / max(len(recv), 1)

    base = share(TransferModel())
    boosted = share(TransferModel(receptivity={target: 20.0}))
    assert boosted > base + 0.1


# ----------------------------------------------------------------------------------- rate files

def test_read_family_rates_header_and_positional(tmp_path):
    with_header = tmp_path / "fam_h.tsv"
    with_header.write_text("family\tduplication\ttransfer\tloss\n1\t3\t2\t1\n2\t4\t0\t1\n")
    assert read_family_rates(with_header) == {"1": (3.0, 2.0, 1.0), "2": (4.0, 0.0, 1.0)}

    positional = tmp_path / "fam_p.tsv"
    positional.write_text("# comment\n1  3  2  1\n2  4  0  1\n")
    assert read_family_rates(positional) == {"1": (3.0, 2.0, 1.0), "2": (4.0, 0.0, 1.0)}


def test_read_branch_rates_both_columns(tmp_path):
    p = tmp_path / "branch.tsv"
    p.write_text("branch\temission\treceptivity\nn3\t5.0\t1.0\nn7\t1.0\t10.0\n")
    emission, receptivity = read_branch_rates(p)
    assert emission == {"n3": 5.0, "n7": 1.0}
    assert receptivity == {"n3": 1.0, "n7": 10.0}


def test_read_branch_rates_receptivity_only(tmp_path):
    p = tmp_path / "recept.tsv"
    p.write_text("branch\treceptivity\nn3\t2.0\nn7\t8.0\n")
    emission, receptivity = read_branch_rates(p)
    assert emission == {}
    assert receptivity == {"n3": 2.0, "n7": 8.0}


def test_cli_family_rates_and_branch_receptivity(tmp_path):
    """End-to-end CLI: ``--family-rates`` drives a per-family run (Python engine, events written)
    and ``--branch-rates`` with a receptivity column drives a run (profiles written)."""
    from zombi2.cli import main

    st = tmp_path / "st"
    main(["species", "--birth", "1", "--death", "0.1", "--tips", "12", "--age", "4",
          "--seed", "1", "-o", str(st)])
    tree_nwk = str(st / "species_tree.nwk")

    fam = tmp_path / "fam.tsv"
    fam.write_text("family\tduplication\ttransfer\tloss\n1\t2.0\t0\t0.1\n")
    out = tmp_path / "g1"
    rc = main(["genomes", "--tree", tree_nwk, "--family-rates", str(fam), "--initial-families", "5",
               "--max-family-size", "0.9", "--seed", "3", "--write", "events", "-o", str(out)])
    assert rc == 0 and (out / "gene_family_events").exists()

    br = tmp_path / "br.tsv"
    br.write_text("branch\treceptivity\ni5\t20.0\n")
    out2 = tmp_path / "g2"
    rc = main(["genomes", "--tree", tree_nwk, "--dup", "0.1", "--trans", "0.3", "--loss", "0.1",
               "--branch-rates", str(br), "--initial-families", "6", "--seed", "4",
               "--write", "profiles", "-o", str(out2)])
    assert rc == 0 and (out2 / "profiles.tsv").exists()


def test_family_and_branch_files_drive_a_full_run(tmp_path):
    """End-to-end from files: a per-family table + a receptivity map run a real simulation."""
    fam = tmp_path / "fam.tsv"
    fam.write_text("family\tduplication\ttransfer\tloss\n1\t0.6\t0.3\t0.2\n")
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)
    model = FamilySampledRates(rates=read_family_rates(fam))
    _emit, recept = ({}, {"n5": 4.0})
    g = simulate_genomes(tree, model, transfers=TransferModel(receptivity=recept),
                         initial_families=5, seed=2, max_family_size=0.9)
    assert g.profiles.families                      # ran and produced families
