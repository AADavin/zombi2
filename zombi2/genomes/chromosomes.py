"""The chromosome genealogy — shared by every genome resolution that has chromosomes.

A chromosome carries identity: its id is re-minted at every event that reshapes it, and each such edge
is recorded as a :class:`ChromosomeEvent`. The resulting edge list *is* the **chromosome network** —
the middle tier between the species tree and the gene trees. It is resolution-agnostic (the ids are
just chromosome ids, whether a chromosome holds ordered gene tokens or nucleotide blocks), so it lives
here, one home, imported by both the ordered and the nucleotide engines.
"""

from __future__ import annotations

from dataclasses import dataclass

from .events import node_label


@dataclass(frozen=True)
class ChromosomeEvent:
    """One edge of the **chromosome genealogy** — a chromosome lineage's birth, split, merge, or
    death, fired on species branch ``lineage`` at ``time``. ``parents`` → ``children`` are chromosome
    ids and the arity names the event: ``"origination"`` (``()`` → one child: an initial or de-novo
    replicon, a **root**), ``"speciation"`` and ``"fission"`` (one parent → two children, a
    **bifurcation**), ``"fusion"`` (two parents → one child, the **reticulation** — in-degree 2, what
    makes this a network and not a tree), ``"loss"`` (one parent → ``()``, a **leaf**). The edge list
    is the network's ground truth; a graph serialisation is derived from it (never eNewick — a
    multi-rooted, reticulating graph is not a tree)."""

    time: float
    kind: str  # "origination" | "speciation" | "fission" | "fusion" | "loss"
    lineage: int
    parents: tuple[int, ...]
    children: tuple[int, ...]


_COLS = ("time", "kind", "lineage", "parents", "children")


def chromosome_events_tsv(chromosome_events: list[ChromosomeEvent]) -> str:
    """The chromosome network as TSV — one row per edge, the ids of a multi-ended side joined by
    ``;`` (a fusion has two parents, a fission two children). Both resolutions write this file, so
    the writer lives with the record it writes."""
    rows = [f"{e.time}\t{e.kind}\t{node_label(e.lineage)}\t{';'.join(map(str, e.parents))}\t"
            f"{';'.join(map(str, e.children))}" for e in chromosome_events]
    return "\n".join(["\t".join(_COLS), *rows]) + "\n"


__all__ = ["ChromosomeEvent", "chromosome_events_tsv"]
