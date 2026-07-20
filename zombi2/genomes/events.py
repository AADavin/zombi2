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


def events_tsv(events: list[Event]) -> str:
    """The event log as TSV — one row per event; empty cells for the fields a kind does not use."""
    cols = ("time", "kind", "lineage", "family", "copy", "parent", "recipient")
    rows = ["\t".join("" if (v := getattr(e, c)) is None else str(v) for c in cols) for e in events]
    return "\n".join(["\t".join(cols), *rows]) + "\n"
