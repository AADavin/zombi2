"""Sequence composition as a conditioning driver — how much of this lineage's sequence is *these*
letters?

One statistic, two front doors. ``result.gc()`` is the named one, G+C over a nucleotide run;
``result.composition(letters)`` is the general one, any letters of the run's own alphabet, which is
how an amino-acid frequency is asked for. Both return a `Composition`, and the pooling,
the direction rule and the interpolation below are shared.

This is what makes Sequences a driver and not only a target (SPEC §3, Traits–Sequences). It points
**forwards only**: a sequence is grown along the gene trees the genome level produced, so a genome
reading it back would condition a run on its own output — Genomes–Sequences can be joined and never
conditioned. `Composition.refuses` says so by name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ..params.driver import DriverTrajectory, interpolated_segments
from ..tree import node_from_label, node_label

if TYPE_CHECKING:                      # `gc()` imports this module, so a run-time import back is a cycle
    from . import SequencesResult

#: the engines that may not read this driver, named as each level names itself
#: (`zombi2.params.modifiers.Modifier.implemented_for`). All three are the genome — a sequence is
#: downstream of every one of them. ``species`` and ``joint`` cannot read it either but never reach
#: here: species takes no ``scaled_by`` at all, and joint refuses every object source up front.
_UPSTREAM = frozenset({"genomes.family", "genomes.ordered", "genomes.nucleotide"})


@dataclass(frozen=True)
class Composition:
    """The share of a lineage's sequence that is one of ``letters``, over time.

    Built by ``result.gc()`` or ``result.composition(letters)`` and handed to
    `~zombi2.params.rate.Rate.scaled_by` like any grown driver::

        seqs = simulate_sequences(g, model=hky85(2.0), length=300, seed=1)
        simulate_discrete(tree, states=["mesophile", "thermophile"], start="mesophile",
                          switch=PerLineage(0.2).scaled_by(seqs.gc(),
                                                           Curve(lambda x: 20.0 ** (x - 0.5))),
                          seed=2)

        proteins = simulate_sequences(g, model=lg(), length=300, seed=1)
        proteins.composition("KR")          # the basic residues — an amino-acid frequency
        proteins.composition("AVLIMFWP")    # the hydrophobic set

    ``letters`` must be letters of the run's own ``alphabet``, so a set that could never occur is
    refused rather than driving every lineage at 0.0.

    **Pooled over the lineage, not one family**, because that is what a composition means. One
    family's is refused rather than offered: it is undefined wherever that family is absent, and a
    driver has to answer for every branch the target walks.

    One value per species node — the lineage's whole complement at the end of its branch — with the
    path between two nodes taken as the straight line, cut at ``step`` (`interpolated_segments`,
    where that approximation's caveats are written).
    """

    result: "SequencesResult"
    letters: str

    def __post_init__(self) -> None:
        alphabet = self.result.alphabet
        if not self.letters:
            raise ValueError("composition() needs at least one letter to count, e.g. 'KR'")
        if not alphabet:
            raise ValueError(
                "this result records no alphabet, so the letters cannot be checked against the one "
                "its sequences are written in. It was not built by simulate_sequences().")
        stray = sorted(set(self.letters) - set(alphabet))
        if stray:
            raise ValueError(
                f"composition({self.letters!r}) names {stray}, which are not in this run's alphabet "
                f"({alphabet}). They occur nowhere, so the driver would read 0.0 on every lineage.")

    def refuses(self, level: str) -> str | None:
        """Why ``level`` may not read this driver, or ``None`` when it may — the hook
        `~zombi2.params.driver.refuse_wrong_direction` calls. It lives on the driver because direction
        is a fact about the pair, and this is the half of the pair that knows it."""
        if level not in _UPSTREAM:
            return None
        return (
            "the genome level cannot be driven by a sequence: a sequence is grown along the gene trees "
            "this level produces, so reading it back conditions a run on its own output. SPEC §3 lets "
            "Genomes and Sequences be joined, never conditioned, and no joint engine for the pair "
            "exists. A composition drives what comes after a sequence — a trait, or a further "
            "sequence run.")

    def _node_values(self, tree) -> dict[int, float]:
        """``{species node: the share of its sequence that is`` ``letters``\\ ``}`` for every node of
        ``tree``.

        Sequences are keyed ``n7_g12``, so the lineage is the half before the ``_g`` — the join every
        ZOMBI2 file uses. A lineage carrying no sequence at all takes its parent's value: 0.0 would
        drive that branch as though its whole complement had turned over."""
        counts: dict[int, list[int]] = {}
        for by_gene in (self.result.alignments, self.result.ancestral):
            for sequences in by_gene.values():
                for label, seq in sequences.items():
                    node = _species_of(label, tree)
                    c = counts.setdefault(node, [0, 0])
                    c[0] += sum(seq.count(x) for x in self.letters)
                    c[1] += len(seq)
        values: dict[int, float] = {}
        stack = [tree.root]
        while stack:                             # pre-order, so a lineage can read its parent's value
            i = stack.pop()
            node = tree.nodes[i]
            hits, total = counts.get(i, (0, 0))
            if total:
                values[i] = hits / total
            elif node.parent is not None:
                values[i] = values[node.parent]
            else:
                raise ValueError(
                    "this sequence run put no sequence on the root lineage, so there is no composition "
                    "to start from — check that the genome run it replayed had families in it.")
            if node.children is not None:
                stack.extend(node.children)
        return values

    def as_driver_trajectory(self, tree, *, step: float | None = None) -> DriverTrajectory:
        """The per-lineage trajectory, for `zombi2.params.driver.resolve_driver`. ``step`` is the
        resolution the path between two nodes is read at; ``None`` takes
        `~zombi2.params.driver.CONTINUOUS_DRIVER_FRACTION` of the tree's height."""
        return DriverTrajectory(interpolated_segments(tree, self._node_values(tree), step))

    def __repr__(self) -> str:
        return f"composition({self.letters!r})"


def _species_of(label: str, tree) -> int:
    """The species lineage a sequence's ``n7_g12`` key names, checked against ``tree`` — the same
    same-tree check the genome drivers make, against the labels rather than a stored tree."""
    head = label.rsplit("_g", 1)[0]
    try:
        node = node_from_label(head)
    except ValueError:
        raise ValueError(f"sequence {label!r} is not named <lineage>_g<copy>, so the lineage it "
                         f"belongs to cannot be read off it") from None
    if node not in tree.nodes:
        raise ValueError(
            f"the sequence run and the run reading it are on different species trees: the sequences "
            f"name lineage {node_label(node)}, which is not in this tree. Grow both on the same "
            f"complete tree.")
    return node
