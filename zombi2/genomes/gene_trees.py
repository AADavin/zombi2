"""Gene trees — the true genealogy of each gene family, evolving inside the complete species tree.

A gene tree **is** the history the run generated. (It is not a *reconciliation* — that word belongs
to inference, recovering an embedding of an *observed* gene tree against the *extant* species tree,
and lives in the tools; here we simulated the embedding, so we simply record it.) One gene tree gives
both:

- the **complete** gene tree — every copy-lineage, the lost ones and those in extinct species included;
- the **extant** gene tree — pruned to the copies that survive at the extant tips, degree-two nodes
  suppressed — mirroring the species result's ``complete_tree`` / ``extant_tree``.

Every node is annotated from the start with the **species branch** it sits on, the **event** that made
it (``origination`` · ``duplication`` · ``transfer`` · ``speciation`` · ``loss`` · ``extant`` ·
``extinct`` · ``unsampled``), its **time**, and its gene **copy**.

Derivation, per family, is a recursive descent of the complete species tree: a copy is threaded along
its species branch — a duplication or transfer bifurcates it (the source continues, a new copy starts,
on the same branch for a duplication or on the recipient branch for a transfer), a loss caps it — and
at the end of a branch the copy either bifurcates at the **speciation** into both daughter species (a
copy id persists across a split) or ends at a tip (extant / extinct / unsampled).
"""

from __future__ import annotations

import collections
import contextlib
import sys
from dataclasses import dataclass, field
from functools import cached_property

# Gene trees are as deep as duplication ladders and the species tree together allow — unbounded in
# principle, well past Python's default limit. Lift it for the recursive build / walk, then restore.
_DEPTH = 100_000


@contextlib.contextmanager
def _deep_recursion():
    old = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old, _DEPTH))
    try:
        yield
    finally:
        sys.setrecursionlimit(old)


@dataclass
class GeneNode:
    """One node in a gene tree — an event in a family's true history.

    ``kind`` — ``origination`` (the founding copy; a unary root) · ``duplication`` · ``transfer``
    (the two horizontal-edge children sit on different species branches) · ``speciation`` (a species
    split the gene was present at) · ``loss`` (a copy lost mid-branch) · ``extant`` / ``extinct`` /
    ``unsampled`` (a copy that reached a species tip of that fate).
    ``species`` — the species-tree node id it sits on. ``time`` — when (crown-forward).
    ``copy`` — the gene copy id along this segment.
    """

    kind: str
    species: int
    time: float
    copy: int
    children: list["GeneNode"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


@dataclass
class GeneTree:
    """One gene family's true genealogy. ``complete`` is the whole tree (lost and extinct-species
    lineages included); ``extant`` is it pruned to the copies surviving at the extant tips (degree-two
    nodes suppressed), or ``None`` if the family left no extant copy. ``to_newick`` serialises either."""

    family: int
    complete: GeneNode

    @cached_property
    def extant(self) -> GeneNode | None:
        with _deep_recursion():
            return _prune_to_extant(self.complete)

    def to_newick(self, which: str = "extant", *, annotate: bool = True) -> str | None:
        """Newick of the ``"extant"`` (default) or ``"complete"`` tree; ``None`` if it is empty.
        Leaves are ``g<copy>_n<species>``; with ``annotate`` internal nodes carry ``<kind>_n<species>``;
        branch lengths are time differences."""
        root = self.extant if which == "extant" else self.complete
        if root is None:
            return None
        return _to_newick(root, annotate) + ";"


def gene_trees_from_events(events: list, tree) -> dict[int, GeneTree]:
    """Derive ``{family id: GeneTree}`` from the event log inside the complete ``tree``."""
    originations: dict[int, tuple[int, int, float]] = {}          # family -> (copy, lineage, time)
    child_ev: dict = collections.defaultdict(list)               # (parent, lineage) -> [(time, kind, new, recip)]
    loss_ev: dict[tuple[int, int], float] = {}                    # (copy, lineage) -> time
    for e in events:
        if e.kind == "origination":
            originations[e.family] = (e.copy, e.lineage, e.time)
        elif e.kind == "duplication":
            child_ev[(e.parent, e.lineage)].append((e.time, "duplication", e.copy, None))
        elif e.kind == "transfer":
            child_ev[(e.parent, e.lineage)].append((e.time, "transfer", e.copy, e.recipient))
        elif e.kind == "loss":
            loss_ev[(e.copy, e.lineage)] = e.time
    for key in child_ev:
        child_ev[key].sort()
    builder = _Builder(tree.nodes, child_ev, loss_ev)
    with _deep_recursion():
        return {fam: GeneTree(fam, builder.origination(copy, lineage, t))
                for fam, (copy, lineage, t) in originations.items()}


class _Builder:
    """Threads one family's copies through the species tree into a gene tree (see module docstring)."""

    def __init__(self, nodes: dict, child_ev: dict, loss_ev: dict):
        self.nodes = nodes
        self.child_ev = child_ev
        self.loss_ev = loss_ev

    def origination(self, copy: int, lineage: int, time: float) -> GeneNode:
        """The family's founding copy: a unary ``origination`` root over its threaded lineage."""
        root = GeneNode("origination", lineage, time, copy)
        root.children = [self._segment(copy, lineage, time)]
        return root

    def _segment(self, copy: int, lineage: int, t0: float) -> GeneNode:
        """The gene-subtree for ``copy`` entering species branch ``lineage`` at ``t0``."""
        timeline = [(t, kind, new, recip)
                    for (t, kind, new, recip) in self.child_ev.get((copy, lineage), ())
                    if t >= t0]
        loss = self.loss_ev.get((copy, lineage))
        if loss is not None and loss >= t0:
            timeline.append((loss, "loss", None, None))
        timeline.sort(key=lambda x: x[0])
        return self._thread(copy, lineage, timeline, 0)

    def _thread(self, copy: int, lineage: int, timeline: list, i: int) -> GeneNode:
        """Walk ``copy``'s events on ``lineage`` in time order; then close at the branch end."""
        if i < len(timeline):
            t, kind, new, recip = timeline[i]
            if kind == "loss":
                return GeneNode("loss", lineage, t, copy)
            node = GeneNode(kind, lineage, t, copy)
            cont = self._thread(copy, lineage, timeline, i + 1)              # the source copy continues
            branch = lineage if kind == "duplication" else recip            # transfer lands on the recipient
            node.children = [cont, self._segment(new, branch, t)]
            return node
        return self._branch_end(copy, lineage)

    def _branch_end(self, copy: int, lineage: int) -> GeneNode:
        """No more events: bifurcate at the speciation into both daughters, or end at the tip."""
        s = self.nodes[lineage]
        if s.children is None:                                              # a species tip
            return GeneNode(s.fate, lineage, s.end_time, copy)             # extant / extinct / unsampled
        node = GeneNode("speciation", lineage, s.end_time, copy)
        node.children = [self._segment(copy, child, s.end_time) for child in s.children]
        return node


def _prune_to_extant(node: GeneNode) -> GeneNode | None:
    """A fresh tree keeping only lineages that reach an ``extant`` tip; degree-two nodes suppressed."""
    if node.is_leaf:
        return GeneNode(node.kind, node.species, node.time, node.copy) if node.kind == "extant" else None
    kept = [k for k in (_prune_to_extant(c) for c in node.children) if k is not None]
    if not kept:
        return None
    if len(kept) == 1:
        return kept[0]                                                     # suppress the degree-two node
    out = GeneNode(node.kind, node.species, node.time, node.copy)
    out.children = kept
    return out


def _to_newick(root: GeneNode, annotate: bool) -> str:
    """Serialise iteratively — gene trees run deeper than CPython's C-stack recursion guard, which
    ``setrecursionlimit`` cannot lift, so recursion would crash on deep (high-turnover) trees."""
    # explicit post-order stack; each frame is [node, parent_time, next_child, collected_child_strings]
    stack: list[list] = [[root, None, 0, []]]
    result = ""
    while stack:
        frame = stack[-1]
        node, parent_time, ci, parts = frame
        if ci < len(node.children):
            frame[2] = ci + 1
            stack.append([node.children[ci], node.time, 0, []])
            continue
        bl = "" if parent_time is None else f":{node.time - parent_time:.6g}"
        if node.is_leaf:
            s = f"g{node.copy}_n{node.species}{bl}"
        else:
            label = f"{node.kind}_n{node.species}" if annotate else ""
            s = f"({','.join(parts)}){label}{bl}"
        stack.pop()
        if stack:
            stack[-1][3].append(s)
        else:
            result = s
    return result
