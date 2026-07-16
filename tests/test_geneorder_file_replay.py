"""The written gene-order events must replay to the written genomes.

Ported from ``thekswenson/Zombi``'s ``tests/test_geneorder_events.py``
(``checkEventsAgainstGenomes``). Krister Swenson's point, and it is the right one: zombi2's
*output is a set of files*, and that output — not just the in-memory structures — is what a user
consumes. Someone inferring rearrangements reads ``Geneorder_events.tsv`` and replays it; this test
asserts that doing so actually reproduces the genomes zombi2 wrote.

The existing ``test_geneorder_events_output`` checks the file's *shape* (rows exist, breakpoints are
populated, the file is deterministic). It does **not** check that the logged breakpoints mean what
they claim — a coordinate-convention slip would produce a plausible file that replays wrong. This
test closes that gap end-to-end: read the events off disk, replay them onto the root genome with the
same primitives the simulator uses, and compare against the genomes on disk (``BED/<node>.bed``).

Scope: a content-conserving run (inversion + transposition), which is the gene-order /
rearrangement-inference case. Duplication / loss / transfer change gene content and are a natural
extension once their file-level semantics are pinned down.
"""

import numpy as np

from zombi2.cli import main
from zombi2.genomes.events import TargetParams
from zombi2.genomes.genome import IdManager
from zombi2.genomes.nucleotide_genome import NucleotideGenome, SegmentRegistry
from zombi2.tools.geneorder_export import read_node_orders
from zombi2.tree import read_newick

ROOT_LENGTH = 3000
GENES = [("g1", 0, 100), ("g2", 1000, 1100), ("g3", 2000, 2100)]
_GENES_TSV = "".join(f"{a}\t{b}\t{n}\n" for n, a, b in GENES)


# --------------------------------------------------------------------------- #
# Reading the written files
# --------------------------------------------------------------------------- #
def _read_events(path):
    """``Geneorder_events.tsv`` -> list of row dicts."""
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def _events_by_branch(rows):
    """``{branch: [row, ...]}``, each branch's events in the order they fired."""
    out = {}
    for r in rows:
        out.setdefault(r["branch"], []).append(r)
    for evs in out.values():
        evs.sort(key=lambda r: float(r["time"]))
    return out


# --------------------------------------------------------------------------- #
# The replay
# --------------------------------------------------------------------------- #
def _seed_root_genome():
    """A genome seeded exactly like the run's root: same length, same gene annotation."""
    reg = SegmentRegistry(pending_genes=[(a, b, n) for n, a, b in GENES])
    g = NucleotideGenome(IdManager(), root_length=ROOT_LENGTH, extension=0.9, registry=reg)
    g.originate(np.random.default_rng(0), TargetParams())  # deterministic gene/intergene tiling
    return g


def _gene_order(g):
    return [(s.gene_id, s.strand) for s in g._segments if s.gene_id is not None]


def _replay_event(g, row):
    """Apply one written event with the same primitive the simulator used."""
    kind = row["event"]
    if kind in ("S", "F", "O"):        # markers / the root seed origination carry no arc
        return
    start, length = int(row["start"]), int(row["length"])
    if kind == "I":
        g._apply_inversion(start, length)
    elif kind == "P":
        g._apply_transposition(start, length, int(row["dest"]))
    else:
        raise AssertionError(f"this replay only covers I/P (content-conserving); saw {kind!r}")


def _path_to_root(node):
    chain = []
    while node is not None:
        chain.append(node)
        node = node.parent
    return list(reversed(chain))


def _replay_to(node, by_branch):
    """Replay every event on the root->node path onto a fresh seed genome."""
    g = _seed_root_genome()
    for anc in _path_to_root(node):
        for row in by_branch.get(anc.name, []):
            _replay_event(g, row)
    return g


# --------------------------------------------------------------------------- #
# The test
# --------------------------------------------------------------------------- #
def _run(tmp_path, seed):
    sp = tmp_path / "S"
    assert main(["species", "--birth", "1", "--death", "0.3", "--tips", "6", "--age", "3",
                 "--seed", str(seed), "-o", str(sp)]) == 0
    genes = tmp_path / "genes.tsv"
    genes.write_text(_GENES_TSV)
    g = tmp_path / "G"
    assert main(["genomes", "-t", str(sp / "species_tree.nwk"), "--genome-model", "nucleotide",
                 "--genes", str(genes), "--root-length", str(ROOT_LENGTH),
                 "--inversion", "0.02", "--transposition", "0.01",
                 "--write", "geneorder", "bed", "--seed", str(seed), "-o", str(g)]) == 0
    return g


def test_written_events_replay_to_the_written_genomes(tmp_path):
    out = _run(tmp_path, seed=11)
    by_branch = _events_by_branch(_read_events(out / "Geneorder_events.tsv"))
    written = read_node_orders(str(out))            # the genomes zombi2 wrote (BED/<node>.bed)
    tree = read_newick((out / "species_tree.nwk").read_text())

    checked = 0
    for node in tree.nodes_preorder():
        if node.name not in written:
            continue
        replayed = _gene_order(_replay_to(node, by_branch))
        assert replayed == written[node.name], (
            f"replaying Geneorder_events.tsv does not reproduce BED/{node.name}.bed\n"
            f"  replayed: {replayed}\n  written:  {written[node.name]}")
        checked += 1
    assert checked >= 2, "expected to check the root and at least one descendant"


def test_replay_is_exercised_by_real_rearrangements(tmp_path):
    # guard: the run must actually contain the events whose replay we claim to verify
    out = _run(tmp_path, seed=11)
    kinds = {r["event"] for r in _read_events(out / "Geneorder_events.tsv")}
    assert "I" in kinds and "P" in kinds
