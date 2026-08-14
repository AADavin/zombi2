"""Indels at the sequence level — which sites a lineage carries, before any of them has a letter.

The **genome** owns indels at the nucleotide resolution, where a genome is a coordinate space of base
pairs and an indel has a position: it can fall inside a gene, move a GFF coordinate, change how long
an assembled chromosome is. The family and ordered resolutions have no such space — a genome there is
gene families and their copies — so nothing at that level can say how long a gene's sequence is. The
sequence level can, and this is where an indel lives for them. One word, one meaning, at whichever
level owns the sites (``docs/design/indels.md``).

**Two passes, not one.** The obvious way to do this is to walk a branch as a process, interleaving
indels and substitutions, and it is the wrong way: it splits every branch at every event, shifts the
site indices at every node, breaks the once-per-family rate classes, and turns each sub-interval into
a branch length that never recurs, so the transition-matrix cache stops paying. None of that is
necessary. Draw the **geometry** first — which columns exist and who carries them — and then let the
existing engine evolve those columns down the tree exactly as it evolves any fixed-width alignment.

That separation is exact rather than convenient. An inserted run is founded from the model's
stationary frequencies, and a reversible model started at ``π`` is still at ``π`` after any branch
length; everything below the insertion depends only on the state where it happened. So evolving the
inserted columns from the root — through lineages that do not carry them, where the letters are
masked and never read — gives the same distribution as founding them at the moment they arrive. The
indel history and the substitution history are independent given the tree.

**Column order.** Every column is an id, and the alignment's order is one list of them. A lineage
carries a subsequence. A deletion drops ids from that lineage's list; an insertion mints fresh ids
and splices them into both its own list and the shared order, immediately after the id they follow.
Two insertions that landed in the same place in different lineages therefore stay two runs, never
one — merging them would assert homology between material that was never related — and an insertion
inside an earlier insertion nests without any special case, because it anchors to an id like any
other.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .clock import Clock  # noqa: F401  — the `clock:` annotation names it
from ..genomes.gene_trees import GeneNode


@dataclass(frozen=True)
class IndelHistory:
    """What the geometry pass produces: how wide the alignment is, and who carries which column.

    ``present[id(node)]`` is a boolean over the ``width`` columns — the mask that turns a node's
    evolved states into its row, gaps where it carries nothing. Keyed by object identity for the
    same reason the states are: a gene-tree node has no unique id of its own."""

    width: int
    present: dict[int, np.ndarray]
    #: How many insertions and deletions actually fired, for the run's counters.
    insertions: int
    deletions: int


def draw_indel_history(root, length: int, *, insertion: float, deletion: float,
                       insertion_extent, deletion_extent, rate_base: float,
                       clock: "Clock | None", rng: np.random.Generator) -> "IndelHistory | None":
    """Walk the gene tree rooted at ``root`` drawing indels, and return who ends up carrying what.

    ``None`` when both rates are zero, and then **not one draw is taken** — a run that did not ask
    for indels is bit-identical to one from before they existed.

    Rates are **relative to substitution**: the events on a branch are Poisson with mean
    ``rate × bl × sites``, where ``bl`` is that branch's length in substitutions per site and
    ``sites`` is how long the sequence is when the branch begins. So ``deletion=0.05`` reads "five
    deletions for every hundred substitutions a site expects", the clock applies to indels for free
    (a lineage that evolves fast gains and loses sites fast), and the number means the same thing on
    a tree of any height — which a per-unit-time rate does not.

    An extent is drawn from its own distribution rather than a mean, so ``Fixed(1)`` really is one
    site (`~zombi2.params.parameter.Extent`). A deletion is clipped to what is there; a deletion
    that would empty the sequence is refused, as one that would empty a chromosome is at the
    nucleotide resolution."""
    if insertion <= 0.0 and deletion <= 0.0:
        return None

    # Where each run sits is recorded as "immediately after this id" rather than spliced into a
    # shared list as it happens: splicing needs the anchor's index, and looking that up per event
    # made the whole pass quadratic in the number of insertions. The order is built once at the end.
    after: dict[int, list[int]] = {}
    head: list[int] = []                           # runs that landed before every original column
    carried: dict[int, list[int]] = {}
    next_id = length
    n_ins = n_del = 0

    # Iterative pre-order, the same walk and the same branch length `evolve_gene_tree` uses, so a
    # lineage's clock reaches its indels as well as its substitutions.
    # The root's parent time is `None` rather than its origination: the founding sequence IS the
    # root's column set, so its stem carries substitutions but no indels. Giving the stem indels
    # would mean the family started at a length nobody asked for.
    stack: "list[tuple[GeneNode, float | None, list[int]]]" = [(root, None, list(range(length)))]
    while stack:
        node, parent_time, parent_seq = stack.pop()
        seq = list(parent_seq)
        if parent_time is not None:
            bl = (rate_base * (node.time - parent_time) if clock is None
                  else clock.branch_length(rate_base, node.species, parent_time, node.time))
            if bl > 0.0 and seq:
                total = (insertion + deletion) * bl * len(seq)
                for _ in range(int(rng.poisson(total))):
                    if not seq:
                        break
                    if rng.random() * (insertion + deletion) < insertion:
                        size = max(1, int(insertion_extent.sample(rng)))
                        at = int(rng.integers(len(seq) + 1))
                        fresh = list(range(next_id, next_id + size))
                        next_id += size
                        # anchored to the id it follows, so a lineage's own order is preserved, two
                        # runs in the same place stay two, and a run inside an earlier run nests
                        (after.setdefault(seq[at - 1], []) if at else head).extend(fresh)
                        seq[at:at] = fresh
                        n_ins += 1
                    else:
                        size = max(1, int(deletion_extent.sample(rng)))
                        if size >= len(seq):           # never empty the sequence
                            continue
                        at = int(rng.integers(len(seq) - size + 1))
                        del seq[at:at + size]
                        n_del += 1
        carried[id(node)] = seq
        for child in reversed(node.children or ()):
            stack.append((child, node.time, seq))

    order: list[int] = []                          # every column, in the alignment's order
    work = list(reversed(head + list(range(length))))
    while work:                                    # depth-first: a run follows the id it anchors to
        cid = work.pop()
        order.append(cid)
        runs = after.get(cid)
        if runs:
            work.extend(reversed(runs))

    where = {cid: i for i, cid in enumerate(order)}
    present: dict[int, np.ndarray] = {}
    for key, seq in carried.items():
        mask = np.zeros(len(order), dtype=bool)
        if seq:
            mask[[where[c] for c in seq]] = True
        present[key] = mask
    return IndelHistory(len(order), present, n_ins, n_del)


__all__ = ["IndelHistory", "draw_indel_history"]
