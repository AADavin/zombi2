"""Joint runs and a family's own rates at the ordered resolution (issue #424).

A rate or an extent can read the gene content the run is building, and a declared family can have its
own duplication, transfer or loss, which applies to the segment an event covers (SPEC §6). The tests
use factors of 0 and 1 where they can, so the checks are exact, and each positive test first requires
at least one event.
"""

import collections

import pytest

import zombi2.genomes.ordered as ordered
from zombi2.genomes import family, simulate_genomes_ordered
from zombi2.genomes.links import Link
from zombi2.params import Extent, Fixed, LogNormal, PerCopy, PerLineage
from zombi2.params.conditioned import resolve_driver
from zombi2.species import simulate_species_tree

ONLY_WITHOUT = {"present": 0.0, "absent": 1.0}
#: a factor that leaves the rate positive either way, so the number a lineage reads moves with its
#: genome rather than sitting at zero
HALF_WITH = {"present": 0.5, "absent": 1.0}


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=1).complete_tree


@pytest.fixture(scope="module")
def clade(tree):
    """One of the root's two daughter lineages, where family A is planted."""
    return sorted(n for n in tree.nodes if tree.nodes[n].parent == tree.root)[0]


def _run(tree, **kw):
    base = dict(duplication=0.2, origination=0.0, initial_families=12, seed=5)
    return simulate_genomes_ordered(tree, **{**base, **kw})


def _before(tree, edge):
    """The instant before an event on its branch, so it reads the genome as it was."""
    return max(edge.time - 1e-9, tree.nodes[edge.lineage].birth_time)


def _presence(result, name, tree):
    return resolve_driver(result.presence(name), tree, step=None, level="genomes.ordered")


def _events(result, kind):
    """The edges of each event of ``kind``, grouped by (time, lineage): one segment is one event."""
    out = collections.defaultdict(list)
    for e in result.edges:
        if e.kind == kind:
            out[(e.time, e.lineage)].append(e)
    return out


def _guarded(f):
    return 0.0 if f == 1.0 else 1.0


def _always(f):
    return 1.0


# --- joint runs: a rate or an extent reads gene content ------------------------------------------

def test_a_run_rate_reads_a_family(tree, clade):
    g = _run(tree, joint=True, loss=PerCopy(0.4).scaled_by("genomes:A", ONLY_WITHOUT),
             families=[family("A", origin=clade)])
    a = _presence(g, "A", tree)
    lost = [e for e in g.edges if e.kind == "loss"]
    assert lost, "nothing was lost, so the test shows nothing"
    assert all(a.value(e.lineage, _before(tree, e)) == "absent" for e in lost)


@pytest.mark.parametrize("factor, anything_lost", [(_guarded, False), (_always, True)],
                         ids=["guarded", "unguarded"])
def test_a_run_rate_reads_a_module(tree, factor, anything_lost):
    g = _run(tree, joint=True, initial_families=0, loss=PerCopy(0.4).scaled_by("genomes:module:m", factor),
             families=[family(n, module="m") for n in ("x", "y", "z")])
    assert any(e.kind == "loss" for e in g.edges) == anything_lost


def test_an_extent_reads_a_family(tree, clade):
    """Where A is present a loss takes three genes, where it is absent one."""
    g = _run(tree, joint=True, loss=0.3,
             loss_extent=Extent(1.0).scaled_by("genomes:A", {"present": 3.0, "absent": 1.0}),
             families=[family("A", origin=clade)])
    a = _presence(g, "A", tree)
    sizes: dict[str, list[int]] = {"present": [], "absent": []}
    for edges in _events(g, "loss").values():
        sizes[a.value(edges[0].lineage, _before(tree, edges[0]))].append(len(edges))
    assert sizes["absent"] and sizes["present"]
    assert set(sizes["absent"]) == {1}
    assert 3 in sizes["present"]


def test_reading_gene_content_needs_joint(tree):
    with pytest.raises(ValueError, match="joint=True"):
        _run(tree, loss=PerCopy(0.4).scaled_by("genomes:A", ONLY_WITHOUT), families=[family("A")])


# --- a family's own rate, applied to the segment -------------------------------------------------

@pytest.mark.parametrize("kind", ["duplication", "transfer", "loss"])
def test_a_family_rate_is_its_own(tree, kind):
    """With the run's rate at 0 and only B writing one, every event of that kind is B's (extent 1)."""
    g = _run(tree, duplication=0.0, transfer=0.0, loss=0.0, families=[family("B", **{kind: 1.0})])
    b = g.family_names["B"]
    events = [e for e in g.edges if e.kind == kind]
    assert events, f"no {kind} happened, so the test shows nothing"
    assert all(e.family == b for e in events)


def test_every_segment_covers_a_gene_of_the_family(tree):
    """Extent 3, and only B has a loss rate: every loss event takes a B gene, and some also take a
    neighbour."""
    g = _run(tree, loss=0.0, loss_extent=Fixed(3), families=[family("B", loss=1.0)])
    b = g.family_names["B"]
    events = list(_events(g, "loss").values())
    assert events
    assert all(any(e.family == b for e in edges) for edges in events)
    assert any(any(e.family != b for e in edges) for edges in events)


def test_a_rate_of_zero_protects_a_family_only_from_events_that_start_on_it(tree):
    """With extent 1 a family whose own loss is 0 is never lost. With extent 4 a segment that starts
    on a neighbour can still take it."""
    def lost(**kw):
        g = _run(tree, loss=0.5, families=[family("steady", loss=0.0)], **kw)
        steady = g.family_names["steady"]
        return sum(1 for e in g.edges if e.kind == "loss" and e.family == steady)

    assert lost() == 0
    assert lost(loss_extent=Fixed(4)) > 0


def test_a_family_rate_reads_another_family(tree, clade):
    g = _run(tree, loss=0.0, joint=True, families=[
        family("A", origin=clade), family("B", loss=PerCopy(1.0).scaled_by("genomes:A", ONLY_WITHOUT))])
    a = _presence(g, "A", tree)
    b = g.family_names["B"]
    lost = [e for e in g.edges if e.kind == "loss" and e.family == b]
    assert lost
    assert all(a.value(e.lineage, _before(tree, e)) == "absent" for e in lost)


def test_a_family_rate_sits_beside_a_per_family_draw(tree, clade):
    """A draw on the run's duplication, and B's own loss read from A: the draw leaves B's loss alone."""
    g = _run(tree, duplication=PerCopy(0.2).varying_among("families", LogNormal(0.0, 0.5)), loss=0.0,
             joint=True, seed=2, families=[
                 family("A", origin=clade),
                 family("B", loss=PerCopy(0.3).scaled_by("genomes:A", ONLY_WITHOUT))])
    a = _presence(g, "A", tree)
    b = g.family_names["B"]
    lost = [e for e in g.edges if e.kind == "loss"]
    assert lost
    assert all(e.family == b for e in lost)
    assert all(a.value(e.lineage, _before(tree, e)) == "absent" for e in lost)


def test_a_family_rate_needs_a_per_copy_run_rate(tree):
    with pytest.raises(ValueError, match="counted per copy"):
        _run(tree, loss=PerLineage(0.5), families=[family("B", loss=1.0)])


# --- the shortcut the engine takes when it sums a family's own rates ------------------------------
# Every gene carries a rate, and a lineage's total is their sum. The engine adds up only the families
# that write their own rate, and takes every other gene at the run's rate, which is the lineage's
# whole weight less what those families hold. `_CHECK_OWN_SUMS` compares that against adding the
# families up one by one, which is the same number.

@pytest.fixture
def checked_sums(monkeypatch):
    monkeypatch.setattr(ordered, "_CHECK_OWN_SUMS", True)


@pytest.mark.parametrize("kw", [
    dict(loss=0.4, families=[family("B", loss=0.9)]),
    dict(loss=0.4, transfer=0.3, families=[family("B", duplication=0.7, loss=0.0),
                                           family("C", transfer=0.8)]),
    dict(loss=0.4, joint=True, families=[family("A"), family("B", loss=PerCopy(0.9).scaled_by(
        "genomes:A", ONLY_WITHOUT))]),
    dict(duplication=PerCopy(0.2).varying_among("families", LogNormal(0.0, 0.5)), loss=0.4,
         families=[family("B", loss=0.9), family("C", duplication=0.1)]),
], ids=["one own rate", "three own rates", "own rate reads a family", "own rate beside a draw"])
def test_the_summed_own_rates_match_the_family_by_family_sum(tree, checked_sums, kw):
    g = _run(tree, **{"origination": 0.2, **kw})
    assert sum(1 for e in g.edges if e.kind in ("duplication", "loss", "transfer")) > 10


def test_the_check_catches_a_wrong_sum(tree, checked_sums, monkeypatch):
    """A lineage weighed one gene too heavy is caught, so the comparison is not vacuous."""
    monkeypatch.setattr(ordered, "_genome_size", lambda genome: 1 + sum(len(c.genes) for c in genome))
    with pytest.raises(AssertionError, match="own-rate sum"):
        _run(tree, loss=0.5, families=[family("B", loss=0.2)])


# --- the rows the engine keeps between steps -----------------------------------------------------
# A step reads each living lineage on its own, and rebuilds only the lineages the step before
# changed. `_CHECK_ROWS` rebuilds every row at every step and compares the ones the engine kept
# against it: a lineage the engine did not mark has to read exactly what it read before.

@pytest.fixture
def checked_rows(monkeypatch):
    monkeypatch.setattr(ordered, "_CHECK_ROWS", True)


@pytest.mark.parametrize("kw", [
    dict(joint=True, loss=PerCopy(0.5).scaled_by("genomes:A", HALF_WITH),
         families=[family("A")]),
    dict(joint=True, loss=0.4, families=[family("A"), family("B", loss=PerCopy(0.9).scaled_by(
        "genomes:A", ONLY_WITHOUT))]),
    dict(joint=True, transfer=0.4, replacement=True, loss=0.2,
         families=[family("A", module="m"), family("B", module="m"),
                   family("C", module="m", loss=PerCopy(0.5).scaled_by("genomes:module:m",
                                                                       _guarded))]),
    dict(joint=True, loss=0.3, chromosomes=3, fission=0.1, fusion=0.1,
         chromosome_origination=0.05, chromosome_loss=0.05, translocation=0.05,
         families=[family("A"), family("B", loss=PerCopy(0.4).scaled_by("genomes:A",
                                                                        ONLY_WITHOUT))]),
    dict(loss=0.4, duplication=PerCopy(0.2).varying_among("families", LogNormal(0.0, 0.5)),
         families=[family("B", loss=0.9)]),
    dict(loss=PerCopy(0.4).changing_at({0.0: 1.0, 0.3: 4.0}), joint=True,
         families=[family("A"), family("B", loss=PerCopy(0.9).scaled_by("genomes:A",
                                                                        ONLY_WITHOUT))]),
], ids=["a run rate reads a family", "an own rate reads a family", "a module, with replacement",
        "chromosome events", "a per-family draw", "a rate on a schedule"])
def test_a_lineage_the_engine_did_not_mark_reads_what_it_read_before(tree, checked_rows, kw):
    g = _run(tree, **{"origination": 0.3, **kw})
    assert sum(1 for e in g.edges if e.kind in ("duplication", "loss", "transfer")) > 10


def test_the_check_catches_an_unmarked_lineage(tree, checked_rows, monkeypatch):
    """With the marking removed the rows go stale, and the check says so — so it is not vacuous."""
    monkeypatch.setattr(ordered._LineageRows, "touched",
                        lambda self, k, genome: self._count(k, genome))   # counts kept, row unmarked
    with pytest.raises(AssertionError, match="differs from a fresh one"):
        _run(tree, origination=0.3, loss=PerCopy(0.5).scaled_by("genomes:A", HALF_WITH),
             joint=True, families=[family("A")])


# --- links.tsv -----------------------------------------------------------------------------------

def test_the_links_of_an_ordered_run(tree):
    g = _run(tree, joint=True, loss=PerCopy(0.3).scaled_by("genomes:A", {"present": 0.5, "absent": 2.0}),
             loss_extent=Extent(1.0).scaled_by("genomes:A", {"present": 2.0, "absent": 1.0}),
             families=[family("A"), family("B", duplication=PerCopy(0.4).scaled_by(
                 "genomes:A", {"present": 2.0, "absent": 1.0}))])
    assert g.links == (
        Link("run", "loss", "genomes:A", "scaled_by", "present=0.5;absent=2.0"),
        Link("run", "loss_extent", "genomes:A", "scaled_by", "present=2.0;absent=1.0"),
        Link("B", "duplication", "genomes:A", "scaled_by", "present=2.0;absent=1.0"),
    )


def test_an_ordered_run_without_links_writes_the_header_alone(tree, tmp_path):
    _run(tree, loss=0.3).write(tmp_path)
    assert (tmp_path / "links.tsv").read_text(encoding="utf-8") == "family\ttarget\tdriver\tmodifier\tmapping\n"
