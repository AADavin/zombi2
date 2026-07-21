"""The written event positions must replay to the written genomes.

Ported from Krister Swenson's fork (``thekswenson/Zombi``, ``tests/test_geneorder_events.py``
``checkEventsAgainstGenomes``), and his argument is the one that makes this test worth having:
ZOMBI's output is a set of *files*, and the files are what someone inferring rearrangements
consumes. Checking the in-memory structures says nothing about whether the *written* coordinates
mean what they claim — a convention slip (0- vs 1-based, before- vs after-event, tandem copies
landing before rather than after) yields a plausible file that replays to the wrong genome, and
every in-memory test still passes.

So this test reads only what ``write`` put on disk — ``genome_event_positions.tsv``,
``rearrangements.tsv``, ``gene_order.tsv`` — replays the whole run forward from an empty root genome
with an independent implementation of each operation, and asserts the result matches
``gene_order.tsv`` at **every** node. It is a global replay, not a per-branch one, because a
transfer's arriving block comes from a contemporaneous donor mid-branch: the ``transfer_donor`` row
says what left, and the ``transfer_recipient`` row it pairs with says where it landed.

Genomes are compared as ``(family, strand)`` sequences: gene ids are re-minted at every event (the
per-segment id model), so they cannot be replayed, while families and orientations can.

Held to one chromosome per genome. The chromosome tier re-mints ids at speciation, fission and
fusion, so replaying a multi-chromosome run additionally needs ``chromosome_events.tsv`` to map a
parent's chromosome onto its daughters' — a separate concern from whether the coordinates are right.
"""

import pytest

from zombi2.genomes import simulate_genomes_ordered
from zombi2.species import simulate_species_tree


# --------------------------------------------------------------------------- #
# Reading the written files
# --------------------------------------------------------------------------- #
def _read_tsv(path):
    lines = path.read_text().splitlines()
    cols = lines[0].split("\t")
    return [dict(zip(cols, row.split("\t"))) for row in lines[1:] if row]


def _read_gene_order(path):
    """``gene_order.tsv`` -> ``{node: [(family, strand), ...]}`` in genome order."""
    genomes = {}
    for r in _read_tsv(path / "gene_order.tsv"):
        genomes.setdefault(int(r["species"]), []).append((int(r["family"]), int(r["strand"])))
    return genomes


def _read_steps(path):
    """Every written event that moves genes, merged into one time-ordered stream.

    Rows sharing a timestamp keep the order they were written (a replacing transfer's displacements
    precede its arrival), so the sort key carries each row's index within its own file.
    """
    steps = []
    for i, r in enumerate(_read_tsv(path / "genome_event_positions.tsv")):
        steps.append((float(r["time"]), 0, i, r))
    for i, r in enumerate(_read_tsv(path / "rearrangements.tsv")):
        steps.append((float(r["time"]), 1, i, r))
    steps.sort(key=lambda s: s[:3])
    return [r for *_, r in steps]


# --------------------------------------------------------------------------- #
# The replay — an independent implementation of each operation
# --------------------------------------------------------------------------- #
def _flip(segment):
    return [(fam, -strand) for fam, strand in reversed(segment)]


def _replay(steps, tree):
    """Run the written history forward from an empty root genome, returning ``{node: genome}``.

    Speciations are interleaved by the tree's own timing: when a branch ends, both daughters start
    from a copy of what it had.
    """
    live = {tree.root: []}
    ended = {}
    in_flight = {}                    # a transfer's donor row, waiting for its recipient row
    pending = sorted((n.end_time, n.id) for n in tree.nodes.values() if n.end_time is not None)
    p = 0

    def settle(upto):
        """Retire every branch that ends at or before ``upto``, seeding its daughters."""
        nonlocal p
        while p < len(pending) and pending[p][0] <= upto:
            _, nid = pending[p]
            p += 1
            genome = live.pop(nid)
            ended[nid] = genome
            for c in (tree.nodes[nid].children or ()):
                live[c] = list(genome)

    for r in steps:
        t, kind, lineage = float(r["time"]), r["kind"], int(r["lineage"])
        settle(t)
        # a branch that has already ended takes no more events; the engine never emits such a row
        assert lineage in live, f"event at {t} on branch {lineage}, which is not alive then"
        g = live[lineage]
        start, length = int(r["start"]), int(r["length"])

        if kind == "origination":
            g.insert(start, (int(r["family"]), +1))
        elif kind == "duplication":
            at = int(r["dest_position"])
            g[at:at] = g[start:start + length]
        elif kind == "loss":
            del g[start:start + length]
        elif kind == "transfer_donor":
            # the donor branch is unchanged; hold what left until its arriving row turns up
            in_flight[(t, int(r["donor"]), int(r["recipient"]))] = g[start:start + length]
        elif kind == "transfer_recipient":
            block = in_flight.pop((t, int(r["donor"]), int(r["recipient"])))
            assert len(block) == length, "the two rows of a transfer disagree on its length"
            g[start:start] = block                      # the block arrives at start
        elif kind == "inversion":
            g[start:start + length] = _flip(g[start:start + length])
        elif kind in ("transposition", "translocation"):
            segment = g[start:start + length]
            del g[start:start + length]                       # excised first, then placed
            at = int(r["dest_position"])
            g[at:at] = _flip(segment) if int(r["flipped"]) else segment
        else:
            raise AssertionError(f"unhandled event kind {kind!r}")

    settle(float("inf"))
    assert not in_flight, f"{len(in_flight)} transfer(s) left without an arriving row"
    return ended | live


def _run(tmp_path, *, seed, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.4, n_extant=12, seed=seed)
    params = dict(duplication=0.4, loss=0.3, origination=0.5, transfer=0.3, inversion=0.5,
                  transposition=0.4, chromosomes=1, initial_families=10, seed=seed)
    params.update(kw)
    r = simulate_genomes_ordered(sp, **params)
    r.write(tmp_path, outputs=("events", "gene_order", "rearrangements", "event_positions"))
    return r


# --------------------------------------------------------------------------- #
# The test
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_written_positions_replay_to_the_written_genomes(tmp_path, seed):
    r = _run(tmp_path, seed=seed)
    replayed = _replay(_read_steps(tmp_path), r.complete_tree)
    written = _read_gene_order(tmp_path)

    assert any(written.values()), "the fixture should produce genes to compare"
    for node, genome in written.items():
        assert replayed[node] == genome, f"replay diverges at node {node}"


def test_replacing_transfers_displace_before_they_arrive(tmp_path):
    # the ordering rule the file format promises: rows sharing a timestamp apply as written, so a
    # replacing transfer's losses land before its block does
    r = _run(tmp_path, seed=7, replacement=True, transfer=0.8)
    # the fixture must actually displace something, or the ordering rule goes untested
    arrivals = [i for i, p in enumerate(r.event_positions) if p.kind == "transfer_recipient"]
    assert arrivals, "fixture produced no transfers"
    displaced = [i for i, p in enumerate(r.event_positions)
                 if p.kind == "loss" and any(p.time == r.event_positions[j].time and i < j
                                             for j in arrivals)]
    assert displaced, "fixture produced no replacement displacements"

    replayed = _replay(_read_steps(tmp_path), r.complete_tree)
    for node, genome in _read_gene_order(tmp_path).items():
        assert replayed[node] == genome, f"replay diverges at node {node}"


def test_the_replay_is_sensitive_to_a_one_position_slip(tmp_path):
    # a negative control: if the test could not tell a correct coordinate from a wrong one, passing
    # would mean nothing. Nudge one written start by a single position and the replay must break.
    r = _run(tmp_path, seed=1)
    steps = _read_steps(tmp_path)
    victim = next(s for s in steps if s["kind"] == "loss" and int(s["start"]) > 0)
    victim["start"] = str(int(victim["start"]) - 1)

    replayed = _replay(steps, r.complete_tree)
    assert replayed != _read_gene_order(tmp_path)


def test_the_replay_is_sensitive_to_a_misplaced_arrival(tmp_path):
    # the same control for the half a transfer's second row is responsible for. Perturbing *one*
    # arrival need not show: genomes are compared as (family, strand), and a block landing between
    # two genes of the same family reads identically either side of the boundary. So perturb each
    # arrival in turn and require that the replay notices at least one.
    r = _run(tmp_path, seed=1)
    written = _read_gene_order(tmp_path)
    n_arrivals = sum(1 for s in _read_steps(tmp_path) if s["kind"] == "transfer_recipient")
    assert n_arrivals, "fixture produced no transfers"

    noticed = 0
    for i in range(n_arrivals):
        steps = _read_steps(tmp_path)                      # a fresh copy for each perturbation
        victim = [s for s in steps if s["kind"] == "transfer_recipient"][i]
        victim["start"] = str(int(victim["start"]) + 1)    # one position further along
        try:
            noticed += _replay(steps, r.complete_tree) != written
        except AssertionError:                             # a start past the end is also "noticed"
            noticed += 1
    assert noticed, f"none of the {n_arrivals} arrivals mattered — the replay is not reading them"


def test_every_gene_content_event_has_a_position(tmp_path):
    # the table is total over the events that change gene content: no silent gaps, and speciation
    # (which copies a genome wholesale) is correctly absent
    r = _run(tmp_path, seed=3)
    # a transfer writes one row per branch, so every genealogy row — donor side and recipient side —
    # finds its position under its own lineage once the two transfer kinds are folded together
    positioned = {(p.time, p.lineage, p.kind.split("_")[0]) for p in r.event_positions}

    kinds = set()
    for e in r.events:
        if e.kind == "speciation":            # a genome is copied wholesale: no position to record
            continue
        kinds.add(e.kind)
        assert (e.time, e.lineage, e.kind) in positioned, \
            f"{e.kind} at {e.time} on {e.lineage} has no position"
    assert kinds == {"origination", "duplication", "loss", "transfer"}, \
        f"the fixture should exercise every gene-content kind, got {sorted(kinds)}"
