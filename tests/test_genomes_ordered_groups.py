"""Groups at the ordered resolution (issue #468): declared families that travel as one segment.

A copy of a group in a genome is a **unit**. The rules under test: a group is laid out whole; every
duplication, loss or transfer that touches one of its genes takes all of them; no rearrangement run
and no insertion ever falls between two genes of a unit, so on every node every unit is one
consecutive run in order; a transfer carries the unit whole to the recipient; and a run without
groups is the run it was.
"""

import collections
import hashlib
import pytest

from zombi2.genomes import family, group, simulate_genomes_ordered
from zombi2.genomes.ordered import Chromosome, Gene, _Units
from zombi2.species import simulate_species_tree


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1).complete_tree


#: the issue's own example, with every segmental and rearrangement event on
_RUN = dict(initial_families=200, duplication=0.1, transfer=0.4, loss=0.4, inversion=0.1,
            transposition=0.1, translocation=0.1, fission=0.05, fusion=0.05, chromosomes=2,
            origination=0.3, duplication_extent=3, loss_extent=2,
            families=[family("A"), family("B"), family("C")],
            groups=[group("operon1", ["A", "B", "C"], fates="group")])


def _units_on(result):
    """``{(node, chromosome index, unit): [positions]}`` of every unit on every node."""
    unit_of = {gene: unit for _g, unit, _f, gene in result.grouped_genes}
    out: dict = {}
    for node, genome in result.node_genomes.items():
        for ci, chrom in enumerate(genome):
            for i, g in enumerate(chrom.genes):
                u = unit_of.get(g.id)
                if u is not None:
                    out.setdefault((node, ci, u), []).append(i)
    return out


def _consecutive(chrom, idx) -> bool:
    n = len(chrom.genes)
    if idx == list(range(idx[0], idx[0] + len(idx))):
        return True
    return chrom.topology == "circular" and \
        sorted((i - idx[0]) % n for i in idx) == list(range(len(idx)))


# --- the three expectations the issue states -------------------------------------------------------

def test_the_group_is_laid_out_whole_and_in_order(tree):
    r = simulate_genomes_ordered(tree, seed=1, **_RUN)
    ids = [r.family_names[n] for n in "ABC"]
    seq = [g.family for chrom in r.initial_genome for g in chrom.genes]
    where = seq.index(ids[0])
    assert seq[where:where + 3] == ids
    assert r.groups == {"operon1": ("A", "B", "C")}
    assert group("x", ["A", "B"]).fates == "family"      # the default (#471)


def test_every_event_touching_a_member_takes_all_three(tree):
    r = simulate_genomes_ordered(tree, seed=1, **_RUN)
    ids = {r.family_names[n] for n in "ABC"}
    touched: dict = collections.defaultdict(set)
    for e in r.edges:
        if e.family in ids and e.kind in ("duplication", "loss", "transfer"):
            touched[(e.time, e.kind, e.lineage)].add(e.family)
    assert touched, "the run has to have hit the group for this to test anything"
    assert all(members == ids for members in touched.values())


def test_every_unit_is_one_consecutive_run_in_order_on_every_node(tree):
    r = simulate_genomes_ordered(tree, seed=1, **_RUN)
    fam_of = {v: k for k, v in r.family_names.items()}
    units = _units_on(r)
    assert len(units) > 10
    for (node, ci, _u), idx in units.items():
        chrom = r.node_genomes[node][ci]
        assert _consecutive(chrom, idx)
        names = "".join(fam_of[chrom.genes[i].family] for i in idx)
        strands = {chrom.genes[i].strand for i in idx}
        assert len(strands) == 1
        assert names == ("CBA" if strands == {-1} else "ABC")     # flipped whole, or not at all


# --- transfers carry the unit whole, duplications copy it whole --------------------------------------

def test_a_transferred_copy_is_a_new_whole_unit(tree):
    r = simulate_genomes_ordered(tree, seed=1, **_RUN)
    ids = {r.family_names[n] for n in "ABC"}
    unit_of = {gene: unit for _g, unit, _f, gene in r.grouped_genes}
    arrivals = [e for e in r.edges if e.kind == "transfer" and e.recipient is not None
                and e.family in ids]
    assert arrivals
    by_event = collections.defaultdict(list)
    for e in arrivals:
        by_event[(e.time, e.recipient)].append(e)
    for copies in by_event.values():
        assert {e.family for e in copies} == ids
        assert len({unit_of[e.copy] for e in copies}) == 1        # one new unit for the three
        donor_units = {unit_of[e.parent] for e in copies}
        assert len(donor_units) == 1 and donor_units != {unit_of[copies[0].copy]}


def test_groups_tsv_names_every_copy_by_group_unit_family_and_gene(tree, tmp_path):
    r = simulate_genomes_ordered(tree, seed=1, **_RUN)
    r.write(tmp_path)
    lines = (tmp_path / "groups.tsv").read_text().splitlines()
    assert lines[0] == "group\tunit\tfamily\tgene"
    assert [ln.split("\t")[:3] for ln in lines[1:4]] == [["operon1", "0", n] for n in "ABC"]
    assert len(lines) - 1 == len(r.grouped_genes)
    # a streamed run writes the same table
    simulate_genomes_ordered(tree, seed=1, stream_to=tmp_path / "stream", **_RUN)
    assert (tmp_path / "stream" / "groups.tsv").read_text() == (tmp_path / "groups.tsv").read_text()


# --- the rules on hand-built chromosomes -----------------------------------------------------------

def _hand():
    u = _Units({"op": ("A", "B", "C")})
    genes = [Gene(i, i, 1) for i in range(8)]
    unit = u.new_unit("op")
    for g in genes[3:6]:
        u.join(g, unit)
    return u, genes


def test_a_run_snaps_to_whole_units_and_counts_a_unit_as_one_gene():
    u, genes = _hand()
    ring = Chromosome(0, "circular", genes)
    assert u.snap(ring, 4, 1) == (3, 3)      # starts inside the unit: the unit, alone
    assert u.snap(ring, 5, 2) == (3, 4)      # the unit and one more gene
    assert u.snap(ring, 2, 2) == (2, 4)      # ends inside the unit: takes it whole
    assert u.snap(ring, 7, 3) == (7, 3)      # wraps position 0, no unit in the way
    line = Chromosome(1, "linear", list(genes))
    assert u.snap(line, 6, 5) == (6, 2)      # a linear chromosome still ends at its last gene


def test_an_insertion_and_a_cut_move_out_of_a_unit():
    u, genes = _hand()
    assert u.slot(genes, 4) == 6 and u.slot(genes, 3) == 3 and u.slot(genes, 6) == 6
    assert u.cut(genes, 4) == 6 and u.cut(genes, 1) == 1
    whole = [Gene(10 + i, i, 1) for i in range(3)]
    unit = u.new_unit("op")
    for g in whole:
        u.join(g, unit)
    assert u.cut(whole, 1) is None           # one unit is the whole chromosome: no cut exists


# --- declarations ---------------------------------------------------------------------------------

def test_a_group_is_checked_against_the_declared_families(tree):
    base = dict(initial_families=5, families=[family("A"), family("B")])
    with pytest.raises(ValueError, match="not a declared family"):
        simulate_genomes_ordered(tree, seed=1, groups=[group("op", ["A", "Z"])], **base)
    with pytest.raises(ValueError, match="two groups"):
        simulate_genomes_ordered(tree, seed=1, groups=[group("op", ["A", "B"]),
                                                       group("op2", ["B", "A"])], **base)
    with pytest.raises(ValueError, match="at least two"):
        group("op", ["A"])
    with pytest.raises(ValueError, match="carries an origin"):
        simulate_genomes_ordered(tree, seed=1, initial_families=5, groups=[group("op", ["A", "B"])],
                                 families=[family("A"), family("B", origin=(tree.root, None))])
    with pytest.raises(ValueError, match="replacement"):
        simulate_genomes_ordered(tree, seed=1, replacement=True, groups=[group("op", ["A", "B"])],
                                 **base)


def test_a_run_without_groups_is_the_run_it_was(tree):
    """Nothing in the plain path moved: the digest of a run with declared families, segments,
    rearrangements and chromosome events, taken before groups existed."""
    r = simulate_genomes_ordered(tree, seed=7, duplication=0.3, transfer=0.3, loss=0.3,
                                 origination=0.5, inversion=0.2, transposition=0.2, translocation=0.1,
                                 fission=0.1, fusion=0.1, chromosomes=2, duplication_extent=3,
                                 loss_extent=2, families=[family("A"), family("B", loss=0.1)],
                                 initial_families=40)
    h = hashlib.sha256()
    for node in sorted(r.node_genomes):
        h.update(repr(r.gene_order(node)).encode())
    assert r.groups == {} and r.grouped_genes == ()
    assert h.hexdigest() == _PLAIN_DIGEST


#: taken on main before #468, by the loop above
_PLAIN_DIGEST = "e6898c5426235a5bdf1d64e2ca8298b1769aea4db7444525e3b8787d6ca4b6a7"


# --- the three dials of #471 --------------------------------------------------------------------

#: the benchmark run of #471: a pack that travels six times in ten, not laid out together, a hub
#: and a decoy whose transfers prefer genomes that carry A but never travel with it
def _bench(_tree=None, **more):
    from zombi2.params import Recipients
    # its own tree: thirty tips with no extinction, so the pack is transferred often enough to count
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=30, seed=1).complete_tree
    return simulate_genomes_ordered(
        tree, initial_families=300, seed=11, joint=True, duplication=0.1, transfer=0.4, loss=0.4,
        inversion=0.1, transposition=0.1,
        families=[family("A"), family("B"), family("C"), family("H", transfer=40.0),
                  family("D", transfer_to=Recipients().weighted_by("genomes:A",
                                                                    {"present": 20.0, "absent": 1.0}))],
        groups=[group("pack", ["A", "B", "C"], together=0.6, adjacent=False)], **more)


def _shared(path):
    return [ln.split("\t") for ln in path.read_text().splitlines()[1:]]


def test_shared_transfers_lists_every_transfer_event_with_the_families_it_carried(tree, tmp_path):
    r = _bench(tree)
    r.write(tmp_path)
    rows = _shared(tmp_path / "shared_transfers.tsv")
    assert (tmp_path / "shared_transfers.tsv").read_text().splitlines()[0] == \
        "event\ttime\tdonor\trecipient\tfamilies"
    transfers = [e for e in r.events if e.kind.startswith("transfer")]
    assert sum(len(row[4].split(";")) for row in rows) == len(transfers)
    # the event id is the one genome_events.tsv carries on each of the event's rows
    ev = (tmp_path / "genome_events.tsv").read_text().splitlines()
    header = ev[0].split("\t")
    assert header[-1] == "event"
    by_id = {}
    for ln in ev[1:]:
        c = ln.split("\t")
        by_id.setdefault(c[-1], []).append(c)
    for event, time, donor, recipient, families in rows:
        own = by_id[event]
        assert {c[header.index("kind")] for c in own} <= {"transfer_additive", "transfer_replacing"}
        assert len(own) == len(families.split(";"))
        assert {c[header.index("time")] for c in own} == {time}


def test_together_is_the_share_of_a_member_s_transfers_that_carry_the_pack(tree, tmp_path):
    r = _bench(tree)
    r.write(tmp_path)
    rows = _shared(tmp_path / "shared_transfers.tsv")
    pack = {"A", "B", "C"}
    touching = [set(row[4].split(";")) for row in rows if pack & set(row[4].split(";"))]
    # the pack travelled: every member the donor still held went along, so two of three where the
    # third had been lost on that branch
    shared = [f for f in touching if len(f & pack) > 1]
    assert len(touching) >= 40
    assert 0.4 <= len(shared) / len(touching) <= 0.8       # 0.6, give or take the draw
    assert any(len(f & pack) == 1 for f in touching)        # and sometimes one member went alone
    # the decoy co-occurs with A and never shares an event with it
    assert not any("D" in f and "A" in f for f in touching)
    assert all(len(f) == 1 for f in (set(row[4].split(";")) for row in rows) if "D" in f)


def test_a_group_at_together_zero_and_not_adjacent_is_no_group_at_all(tree):
    """The null of the benchmark: the same seed with the group declared and with the families
    merely declared gives the same run, event for event."""
    base = dict(initial_families=120, seed=5, duplication=0.2, transfer=0.4, loss=0.3, inversion=0.1,
                transfer_extent=2, families=[family("A"), family("B"), family("C")])
    plain = simulate_genomes_ordered(tree, **base)
    null = simulate_genomes_ordered(tree, groups=[group("p", ["A", "B", "C"], together=0.0,
                                                        adjacent=False)], **base)
    assert [(e.time, e.kind, e.family) for e in plain.edges] == \
        [(e.time, e.kind, e.family) for e in null.edges]
    assert {n: plain.gene_order(n) for n in plain.node_genomes} == \
        {n: null.gene_order(n) for n in null.node_genomes}


def test_fates_family_lets_members_drift_apart_and_fates_group_does_not(tree):
    base = dict(initial_families=150, seed=3, duplication=0.3, loss=0.5,
                families=[family("A"), family("B"), family("C")])
    drift = simulate_genomes_ordered(tree, groups=[group("op", ["A", "B", "C"])], **base)
    held = simulate_genomes_ordered(tree, groups=[group("op", ["A", "B", "C"], fates="group")], **base)
    ids = {drift.family_names[n] for n in "ABC"}

    def partial_losses(r):
        by = {}
        for e in r.edges:
            if e.kind == "loss" and e.family in ids:
                by.setdefault((e.time, e.lineage), set()).add(e.family)
        return sum(1 for fams in by.values() if fams != ids)
    assert partial_losses(drift) > 0
    assert partial_losses(held) == 0
    # a duplicated member: the continuation stays in the unit, the tandem copy is not a member
    unit_of = {gene: unit for _g, unit, _f, gene in drift.grouped_genes}
    dup = [e.copy for e in drift.edges if e.kind == "duplication" and e.family in ids
           and e.parent in unit_of]              # per parent: the continuation, then the copy
    assert dup
    assert all(c in unit_of for c in dup[0::2])
    assert not any(c in unit_of for c in dup[1::2])


def test_the_streamed_benchmark_writes_the_same_tables(tree, tmp_path):
    _bench(tree).write(tmp_path / "mem")
    _bench(tree, stream_to=tmp_path / "stream")
    for f in ("shared_transfers.tsv", "genome_events.tsv", "groups.tsv"):
        assert (tmp_path / "mem" / f).read_text() == (tmp_path / "stream" / f).read_text()


def test_the_dials_are_checked():
    with pytest.raises(ValueError, match="probability"):
        group("p", ["A", "B"], together=1.5)
    with pytest.raises(ValueError, match="fates"):
        group("p", ["A", "B"], fates="unit")
