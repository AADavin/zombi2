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
      horizontal edge). Two rows, same ``parent``.
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


_COLS = ("time", "kind", "lineage", "family", "copy", "parent", "recipient")


def events_tsv(events: list[Event]) -> str:
    """The event log as TSV — one row per event; empty cells for the fields a kind does not use."""
    rows = ["\t".join("" if (v := getattr(e, c)) is None else str(v) for c in _COLS) for e in events]
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
    if tuple(header) != _COLS:
        # the nucleotide resolution writes its own, wider log to the same filename — a likely
        # mistake worth naming, since the two look alike until you read the columns
        hint = ("; this looks like a --resolution nucleotide log, whose events are keyed by "
                "ancestral interval rather than gene family — sequences replays the unordered or "
                "ordered log" if "source" in header and "family" not in header else "")
        raise ValueError(f"unexpected genome-event columns {header}; expected {list(_COLS)}{hint}")
    events: list[Event] = []
    for lineno, raw in enumerate(lines[1:], 2):
        if not raw:                                     # tolerate a trailing blank line
            continue
        cells = raw.split("\t")
        if len(cells) != len(_COLS):
            raise ValueError(f"genome event log line {lineno}: expected {len(_COLS)} columns, "
                             f"got {len(cells)}")
        time, kind, lineage, family, copy, parent, recipient = cells
        events.append(Event(
            time=float(time), kind=kind, lineage=int(lineage), family=int(family), copy=int(copy),
            parent=int(parent) if parent else None, recipient=int(recipient) if recipient else None))
    return events
