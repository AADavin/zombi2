"""The chromosome genealogy — shared by every genome resolution that has chromosomes.

A chromosome carries identity: its id is re-minted at every event that reshapes it, and each such edge
is recorded as a :class:`ChromosomeEvent`. The resulting edge list *is* the **chromosome network** —
the middle tier between the species tree and the gene trees. It is resolution-agnostic (the ids are
just chromosome ids, whether a chromosome holds ordered gene tokens or nucleotide blocks), so it lives
here, one home, imported by both the ordered and the nucleotide engines.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ChromosomeEvent:
    """One edge of the **chromosome genealogy** — a chromosome lineage's birth, split, merge, or
    death, fired on species branch ``lineage`` at ``time``. ``parents`` → ``children`` are chromosome
    ids and the arity names the event: ``"origination"`` (``()`` → one child: a seed or de-novo
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


__all__ = ["ChromosomeEvent"]
