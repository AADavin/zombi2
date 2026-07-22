"""The gene-genealogy event log — the shared source of truth every resolution writes.

An :class:`Event` records one moment in a gene family's history; the per-family gene trees are
*derived* from a run's events (see :mod:`.gene_trees`), identically whether the genome was an
unordered multiset or an ordered set of chromosomes. Position and orientation are **not** here —
they live in the genome snapshots and the rearrangement log — because an event is about gene
*identity and descent*, which is resolution-blind. So this module is imported by both the unordered
core and the ordered engine, and neither owns it.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Event:
    """A recorded genome event — the true history every per-family gene tree is derived from. Gene
    ids are **per segment** (the ZOMBI1 model): every event ends a gene and starts fresh ids for its
    descendants, so an id belongs to exactly one species branch and every node's genome has its own
    ids. ``lineage`` is the species-tree node the event fired on; ``time`` is when (crown-forward,
    the species tree's clock). By kind:

    - ``"origination"`` — ``copy`` is a founding gene of a fresh family (``parent`` ``None``): a root.
    - ``"duplication"`` — the gene ``parent`` ends; ``copy`` is one of its **two** descendants (the
      continuation and the new copy — two rows, same ``parent``), both on ``lineage``.
    - ``"transfer"`` — the donor gene ``parent`` ends; ``copy`` is one of its two descendants: the
      continuation on the donor ``lineage``, or the transferred copy on the ``recipient`` lineage (a
      horizontal edge). Two rows, same ``parent``, and **both name ``donor``** — the branch the
      material left. Without it the arriving row said only where the copy landed (twice over:
      ``lineage`` and ``recipient`` are the same branch there), so reading who donated to whom meant
      pairing the two rows on ``(time, parent)``. A transfer is an edge; each of its rows names both
      ends.
    - ``"speciation"`` — the gene ``parent`` ends at a split; ``copy`` is its descendant in daughter
      species ``lineage`` (one row per daughter — two, same ``parent``).
    - ``"loss"`` — the gene ``copy`` ends with no descendant (``parent`` ``None``).
    """

    time: float
    kind: str  # "origination" | "duplication" | "loss" | "transfer" | "speciation"
    lineage: int  # the species-tree node id where it fired (for a transfer: the donor lineage)
    family: int
    copy: int  # the copy born (origination / duplication / transfer) or removed (loss)
    parent: int | None = None  # duplication & transfer: the source copy (which survives)
    recipient: int | None = None  # transfer only: the species lineage the new copy is born on
    #: transfer only: the species lineage the material left. Set on **both** rows, so either names
    #: the whole edge; on the donor's own row it repeats ``lineage``, which is the price of that.
    donor: int | None = None


_COLS = ("time", "kind", "lineage", "family", "copy", "parent", "recipient", "donor")

#: Columns holding a species-tree node, written as ``n<id>``. ``parent`` is a *gene copy*, not a
#: lineage, so it stays a bare id.
_NODE_COLS = frozenset({"lineage", "recipient", "donor"})


def node_label(node_id: int | None) -> str:
    """A species-tree node as every ZOMBI2 table writes it: ``n<id>``, the same token the Newick
    uses and the species and trait tables already used. Empty for ``None``."""
    return "" if node_id is None else f"n{node_id}"


def node_from_label(cell: str) -> int:
    """The inverse of :func:`node_label`. A bare integer is accepted too, so a log written before the
    node columns carried their ``n`` still replays."""
    return int(cell[1:] if cell[:1] == "n" else cell)


def events_tsv(events: list[Event]) -> str:
    """The event log as TSV — one row per event; empty cells for the fields a kind does not use."""
    def cell(e: Event, col: str) -> str:
        v = getattr(e, col)
        if v is None:
            return ""
        return node_label(v) if col in _NODE_COLS else str(v)

    rows = ["\t".join(cell(e, c) for c in _COLS) for e in events]
    return "\n".join(["\t".join(_COLS), *rows]) + "\n"


def events_from_tsv(text: str) -> list[Event]:
    """Parse the TSV :func:`events_tsv` writes back into a ``list[Event]`` — the deserializer twin, so
    a written ``genome_events.tsv`` can be replayed (a downstream level's gene trees are derived from
    the log by :func:`~zombi2.genomes.gene_trees.gene_trees_from_events`). ``time`` is a float, the
    id columns are ints, and the optional ``parent`` / ``recipient`` are ints or ``None`` (empty)."""
    lines = text.splitlines()
    if not lines:
        raise ValueError("empty genome event log — is the file empty?")
    header = lines[0].split("\t")
    if tuple(header[:len(_COLS)]) == _COLS:
        # the ordered resolution writes the same genealogy with each event's position beside it, and
        # its rearrangements in the same table. The extra columns are ignored here and the
        # rearrangement rows skipped: this reader is about identity and descent, which they do not
        # touch. (A genome level that needs the positions reads them itself.)
        return _parse(lines, header)
    if tuple(header) != _COLS:
        # the nucleotide resolution writes its own, wider log to the same filename — a likely
        # mistake worth naming, since the two look alike until you read the columns
        hint = ("; this looks like a --resolution nucleotide log, whose events are keyed by "
                "ancestral interval rather than gene family. Read one with "
                "zombi2.genomes.nucleotide.read_nucleotide_genomes, which the sequence level uses "
                "when its handoff is a nucleotide run"
                if "source" in header and "family" not in header else "")
        raise ValueError(f"unexpected genome-event columns {header}; expected {list(_COLS)}{hint}")
    return _parse(lines, header)


#: the kinds this log records about gene identity and descent. A wider table may also carry the
#: ancestry-**neutral** rearrangements; they end no gene lineage, so they are not events here.
_GENEALOGY = frozenset({"origination", "duplication", "loss", "transfer", "speciation"})


def _parse(lines: list[str], header: list[str]) -> list[Event]:
    """Read the rows by column **name**, so a table carrying more than the seven canonical columns
    parses unchanged and only the genealogy rows come back."""
    at = {c: i for i, c in enumerate(header)}
    events: list[Event] = []
    for lineno, raw in enumerate(lines[1:], 2):
        if not raw:                                     # tolerate a trailing blank line
            continue
        cells = raw.split("\t")
        if len(cells) != len(header):
            raise ValueError(f"genome event log line {lineno}: expected {len(header)} columns, "
                             f"got {len(cells)}")
        if cells[at["kind"]] not in _GENEALOGY:
            continue
        get = lambda c: cells[at[c]]                    # noqa: E731
        events.append(Event(
            time=float(get("time")), kind=get("kind"), lineage=node_from_label(get("lineage")),
            family=int(get("family")), copy=int(get("copy")),
            parent=int(get("parent")) if get("parent") else None,
            recipient=node_from_label(get("recipient")) if get("recipient") else None,
            donor=node_from_label(get("donor")) if get("donor") else None))
    return events
