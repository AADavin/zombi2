"""Intra-genome gene conversion — validation (the hard rule in docs/validation.md).

Gene conversion is a **core** model: ``SharedRates(conversion=...)`` emits the events and
``zombi2.ConversionModel`` sets the donor directionality. The *engine capability* it drives
(``EventType.CONVERSION``, ``UnorderedGenome.convert``, the simulator dispatch, the reconciliation)
lives in the core alongside it.

A conversion overwrites one gene copy of a family with a copy of *another* copy of the **same
family** in the same genome: non-reciprocal (the donor is unchanged), copy-number-neutral (a
replacement, not a duplication), and — repeated — homogenising: the converted copy descends from the
donor, so the two coalesce at the conversion time rather than at their duplication (concerted
evolution, and the gene-tree distortion the feature exists to produce).

The suite, strongest first (see docs/validation.md):

* **Reduction / determinism** — at rate 0 the model emits exactly the weights of a plain
  ``SharedRates`` (byte-identical log on the same engine); a fixed seed reproduces a run exactly.
* **Invariant** — ``convert`` conserves copy number and has the right genealogical shape; a hand-built
  log reconstructs to the converted copy coalescing with the donor at the conversion time.
* **Oracle** — the mean within-family coalescence depth of a stable two-copy family matches ``1/(2c)``.
* **Bias** — the directional knob tilts the donor toward the family's oldest lineage.
* **Core surface** — the public API carries it and it routes to the pure-Python engine.
* **Trace** — a conversion survives the ``events_trace.tsv`` round-trip and the compact-trace expansion.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import zombi2 as z
from zombi2.genomes.conversion import ConversionModel
from zombi2.genomes.events import EventRecord, EventType, GeneOp
from zombi2.genomes.genome import Gene, IdManager, UnorderedGenome
from zombi2.genomes.genome_sim import GenomeSimulator
from zombi2.genomes.rates import SharedRates
from zombi2.genomes.reconciliation import _node_tree, build_gene_trees, reconcile
from zombi2.genomes.simulation import events_trace_from_log, read_events_trace
from zombi2.tree import Tree, TreeNode


def _pair_tree(tau: float) -> Tree:
    """A two-tip tree of age ``tau`` (root at 0, two extant leaves at ``tau``). Each tip lineage is
    a single long branch on which one family can settle at exactly two copies — the clean setting
    for the coalescence-depth oracle."""
    root = TreeNode("root", 0.0)
    for name in ("n0", "n1"):
        root.add_child(TreeNode(name, float(tau)))
    return Tree(root, float(tau))


def _nodes_by_gid(root) -> dict:
    out: dict = {}

    def walk(n):
        out[n.gid] = n
        for c in n.children:
            walk(c)

    walk(root)
    return out


# --------------------------------------------------------------------------- reduction / determinism

def test_conversion_zero_is_byte_identical_to_no_conversion():
    """At ``conversion=0`` the model emits exactly the same weights as a plain ``SharedRates`` with
    no conversion argument, so — driven on the same (Python) engine with the same seed — the whole
    event log is byte-identical. The conversion machinery adds nothing until it is switched on."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=10, age=4.0, seed=1)
    with_conv = GenomeSimulator().simulate(
        tree, SharedRates(duplication=0.5, loss=0.1, origination=0.2, conversion=0.0),
        np.random.default_rng(7), initial_size=6)
    without = GenomeSimulator().simulate(
        tree, SharedRates(duplication=0.5, loss=0.1, origination=0.2),
        np.random.default_rng(7), initial_size=6)
    assert events_trace_from_log(with_conv.event_log) == events_trace_from_log(without.event_log)


def test_conversion_run_is_deterministic():
    """A fixed seed (and bias) reproduces a conversion run exactly."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=10, age=4.0, seed=2)

    def run():
        return z.simulate_genomes(
            tree, SharedRates(duplication=0.6, loss=0.1, conversion=1.0),
            conversions=ConversionModel(bias=0.5), initial_families=6, seed=13)

    assert events_trace_from_log(run().event_log) == events_trace_from_log(run().event_log)


# --------------------------------------------------------------------------------------- invariants

def test_convert_conserves_copy_number_and_shape():
    """Each conversion is net-zero on copy number and has the donor-bifurcation + recipient-loss
    shape; the overwritten copy is a distinct, pre-existing copy (not one of the donor's new ids)."""
    ids = IdManager()
    genome = UnorderedGenome(ids)
    for _ in range(4):
        genome._add(Gene(ids.new_gene(), "1", origin_order=ids.new_order()))
    rng = np.random.default_rng(0)
    for _ in range(50):
        before = genome.copy_number("1")
        donor_group, loss_group = genome.convert("1", rng, bias=0.0)
        assert genome.copy_number("1") == before  # replacement, not duplication
        assert [op.role for op in donor_group] == ["parent", "donor_copy", "converted_copy"]
        assert [op.role for op in loss_group] == ["converted_out"]
        assert all(op.family == "1" for op in donor_group + loss_group)
        new_gids = {donor_group[1].gid, donor_group[2].gid}
        assert loss_group[0].gid not in new_gids           # the recipient is not a donor product
        assert donor_group[0].gid != loss_group[0].gid     # donor and recipient are distinct copies


# ------------------------------------------------------------------------------------------- oracle

def test_reconstruction_places_converted_copy_on_donor_at_conversion_time():
    """White-box oracle: a hand-built log (originate, duplicate, then convert) reconstructs so the
    converted copy and the donor continuation coalesce at the conversion time (2.0), NOT at the
    duplication time (0.1), and the overwritten copy is a terminal loss at the conversion time."""
    total_age = 5.0
    records = [
        EventRecord(EventType.ORIGINATION, "root", 0.0, [GeneOp("g1", "1", "origin")]),
        EventRecord(EventType.DUPLICATION, "n0", 0.1,
                    [GeneOp("g1", "1", "parent"), GeneOp("g2", "1", "left"),
                     GeneOp("g3", "1", "right")]),
        # convert: donor g2 bifurcates into g4 (continuation) + g5 (converted); recipient g3 dies
        EventRecord(EventType.CONVERSION, "n0", 2.0,
                    [GeneOp("g2", "1", "parent"), GeneOp("g4", "1", "donor_copy"),
                     GeneOp("g5", "1", "converted_copy")], donor="n0", recipient="n0"),
        EventRecord(EventType.LOSS, "n0", 2.0, [GeneOp("g3", "1", "converted_out")]),
    ]
    gid2species = {"g4": "n0", "g5": "n0"}
    nodes = _nodes_by_gid(_node_tree(records, gid2species, total_age))

    assert nodes["g2"].kind is EventType.CONVERSION
    assert nodes["g2"].end == pytest.approx(2.0)                     # the donor bifurcation time
    assert {c.gid for c in nodes["g2"].children} == {"g4", "g5"}
    assert nodes["g4"].is_extant and nodes["g5"].is_extant
    # the two survivors coalesce at the conversion time, not the (earlier) duplication time
    assert nodes["g4"].birth == pytest.approx(2.0) == nodes["g5"].birth
    assert nodes["g2"].birth == pytest.approx(0.1)                   # donor lineage predates it
    # the overwritten copy's ancestry ends at the conversion time (a loss, pruned from the extant tree)
    assert nodes["g3"].is_loss and nodes["g3"].end == pytest.approx(2.0)

    # the public reconstruction keeps both survivors, labels the event 'C', and drops the loss
    _complete, extant = build_gene_trees(records, gid2species, total_age)
    assert "g4" in extant and "g5" in extant and "LOSS" not in extant
    recon = reconcile(records, gid2species, total_age)
    assert any(ev.event == "C" and ev.time == pytest.approx(2.0) for ev in recon.events)


@pytest.mark.parametrize("c", [0.5, 2.0])
def test_conversion_homogenizes_coalescence_depth_matches_theory(c):
    """Forward oracle. On a long lineage a family held at exactly two copies is converted as a
    Poisson process of rate ``c·n = 2c``; each conversion resets the pair's coalescence to that
    time, so the depth of the two extant copies (present minus their last conversion) is a
    (near-uncensored) exponential with mean ``1/(2c)``. Recovered across many two-tip replicates."""
    tau = 40.0
    n_reps = 300
    depths: list[float] = []
    for i in range(n_reps):
        g = z.simulate_genomes(_pair_tree(tau),
                               SharedRates(duplication=5.0, conversion=c),
                               initial_families=1, max_family_size=2, seed=1000 + i)
        records = g.gene_families.get("1", [])
        for tip in ("n0", "n1"):
            # the two copies on this tip coalesce at the last event that made a fresh pair there —
            # a duplication or (dominating, on a long branch) a conversion
            pair_times = [r.time for r in records if r.branch == tip
                          and r.event in (EventType.DUPLICATION, EventType.CONVERSION)]
            if pair_times:
                depths.append(tau - max(pair_times))
    assert len(depths) > 500
    assert np.mean(depths) == pytest.approx(1.0 / (2.0 * c), rel=0.15)


# --------------------------------------------------------------------------------------------- bias

def test_bias_directs_donor_toward_oldest():
    """The directional knob: bias 0 draws the donor uniformly (~1/3 of three candidates are the
    oldest), bias 1 always picks the oldest lineage, and it rises monotonically in between —
    P(oldest) = bias + (1-bias)/k for k candidates."""
    others = [Gene("a", "1", origin_order=1),   # oldest (smallest origin_order)
              Gene("b", "1", origin_order=2),
              Gene("c", "1", origin_order=3)]

    def oldest_fraction(bias, trials=5000):
        rng = np.random.default_rng(0)
        hits = sum(UnorderedGenome._choose_donor(others, rng, bias).origin_order == 1
                   for _ in range(trials))
        return hits / trials

    f0, f_half, f1 = oldest_fraction(0.0), oldest_fraction(0.5), oldest_fraction(1.0)
    assert f0 == pytest.approx(1 / 3, abs=0.05)
    assert f_half == pytest.approx(0.5 + 0.5 / 3, abs=0.05)
    assert f1 == 1.0
    assert f0 < f_half < f1


# --------------------------------------------------------------------------------------- core surface

def test_conversion_is_in_the_core_public_api():
    """Promoted to the core: ``ConversionModel`` is on the top-level namespace and ``SharedRates``
    takes a ``conversion`` rate, both constructed **without** an experimental (or any) warning."""
    assert z.ConversionModel is ConversionModel
    assert "ConversionModel" in z.__all__
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning (e.g. a leftover ExperimentalWarning) fails
        SharedRates(duplication=0.1, conversion=0.1)
        ConversionModel(bias=0.5)


def test_conversion_runs_on_the_python_engine():
    """``SharedRates(conversion>0)`` is deliberately **not** Rust-eligible — the counts-only Rust
    path tracks copy-number profiles, not gene-tree shape, so it would silently drop conversions —
    while a plain ``SharedRates`` still takes the Rust path (no regression)."""
    from zombi2 import _rust

    assert not _rust.eligible(SharedRates(duplication=0.2, conversion=0.5), UnorderedGenome, None)
    assert _rust.eligible(SharedRates(duplication=0.2), UnorderedGenome, None)
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.0), n_tips=6, age=4.0, seed=1)
    g = z.simulate_genomes(tree, SharedRates(duplication=0.6, conversion=1.5),
                           initial_families=4, seed=3)
    assert any(r.event is EventType.CONVERSION for r in g.event_log)


# -------------------------------------------------------------------------------------------- trace

def test_conversion_survives_events_trace_roundtrip():
    """A conversion (a 3-id ``C`` row + its paired ``L`` row) survives writing and re-reading
    ``events_trace.tsv``."""
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=8, age=4.0, seed=5)
    g = z.simulate_genomes(tree, SharedRates(duplication=0.5, conversion=1.0),
                           initial_families=5, seed=9)
    families = read_events_trace(events_trace_from_log(g.event_log))
    before = sum(1 for r in g.event_log if r.event is EventType.CONVERSION)
    after = sum(1 for recs in families.values() for r in recs if r.event is EventType.CONVERSION)
    assert before == after > 0


def test_expand_trace_handles_conversion():
    """Defensive: a compact trace (no speciation rows) carrying a conversion expands correctly — a
    conversion is a same-branch bifurcation like a duplication, and must not be misrouted as a
    transfer (whose recipient handling would mis-place, or crash on, the second child)."""
    from zombi2.genomes.reconciliation import expand_trace

    tree = _pair_tree(3.0)
    families = {"1": [
        EventRecord(EventType.ORIGINATION, "root", 0.0, [GeneOp("g1", "1", "origin")]),
        EventRecord(EventType.DUPLICATION, "root", 0.0,
                    [GeneOp("g1", "1", "p"), GeneOp("g2", "1", "l"), GeneOp("g3", "1", "r")]),
        EventRecord(EventType.CONVERSION, "n0", 1.0,
                    [GeneOp("g2", "1", "p"), GeneOp("g4", "1", "dc"), GeneOp("g5", "1", "cc")],
                    donor="n0", recipient="n0"),
        EventRecord(EventType.LOSS, "n0", 1.0, [GeneOp("g3", "1", "co")]),
    ]}
    expanded = expand_trace(families, tree)["1"]
    conv = [r for r in expanded if r.event is EventType.CONVERSION]
    assert conv, "conversion dropped or misrouted during compact-trace expansion"
    # both children of the expanded conversion sit on the conversion's own branch (n0), not elsewhere
    assert all(r.branch == "n0" and r.recipient == "n0" for r in conv)


# ----------------------------------------------------------------------------------------------- CLI

def test_cli_conversion_writes_conversion_events(tmp_path):
    """End-to-end: ``zombi2 genomes --conversion`` runs, takes the full (Python) path, and records
    conversion (``C``) events in the written event log."""
    from zombi2.cli import main

    st = tmp_path / "st"
    main(["species", "--birth", "1.0", "--death", "0.2", "--tips", "10", "--age", "4.0",
          "--seed", "4", "-o", str(st)])
    out = tmp_path / "gen"
    rc = main(["genomes", "--tree", str(st / "species_tree.nwk"),
               "--dup", "0.6", "--conversion", "1.5", "--conversion-bias", "1.0",
               "--initial-families", "6", "--seed", "3", "--write", "events", "-o", str(out)])
    assert rc == 0
    # events are written per family under gene_family_events/<fam>_events.tsv (columns:
    # time event branch donor recipient nodes); a conversion is the event code 'C'
    event_files = list((out / "gene_family_events").glob("*_events.tsv"))
    assert event_files
    codes = {line.split("\t")[1] for f in event_files for line in f.read_text().splitlines()[1:]}
    assert "C" in codes  # at least one conversion fired through the CLI
