"""Gene trees — the true genealogy of each gene family, evolving inside the complete species tree.

A gene tree **is** the history the run generated. (It is not a *reconciliation* — that word belongs
to inference, recovering an embedding of an *observed* gene tree against the *extant* species tree,
and lives in the tools; here we simulated the embedding, so we simply record it.) One gene tree gives
both:

- the **complete** gene tree — every gene lineage, the lost ones and those in extinct species included;
- the **extant** gene tree — pruned to the genes that survive at the extant tips, degree-two nodes
  suppressed — mirroring the species result's ``complete_tree`` / ``extant_tree``.

Gene ids are **per segment** (the ZOMBI1 model): every event ends a gene and starts fresh ids for its
descendants, so an id belongs to one species branch and every gene-tree tip is unique. That makes each
gene a single node whose children are recorded directly in the event log, so building the tree is a
plain parent→children graph read: a gene's node ``kind`` is the event that *ended* it —
``duplication`` · ``transfer`` · ``speciation`` (two children each) internally, or ``loss`` (dead) /
``extant`` · ``extinct`` · ``unsampled`` (a tip of that species fate) at a leaf. The founding gene of a
family (its origination) is simply the root. Each node also knows its species branch, its gene id, and
its end time (branch lengths are time differences).
"""

from __future__ import annotations

import collections
from dataclasses import dataclass, field


@dataclass
class GeneNode:
    """One gene in a gene tree. ``kind`` is the event that ended it — ``duplication`` · ``transfer``
    (the two children sit on different species branches) · ``speciation`` internally, or ``loss`` /
    ``extant`` · ``extinct`` · ``unsampled`` at a leaf. ``species`` is the species-tree node the gene
    lived on; ``time`` is when it ended (crown-forward); ``copy`` is the gene id."""

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
    lineages included); ``extant`` is it pruned to the genes surviving at the extant tips (degree-two
    nodes suppressed), or ``None`` if the family left no extant gene. ``to_newick`` serialises either.

    ``origination`` is when the family was founded — the exact time of its origination event, or the
    root lineage's start for a family seeded by ``initial_families``. A :class:`GeneNode` records when
    it *ended*, so this is the one time the tree cannot derive: it is where the root's branch begins."""

    family: int
    complete: GeneNode
    origination: float

    @property
    def extant(self) -> GeneNode | None:
        if not hasattr(self, "_extant"):
            self._extant = _prune_to_extant(self.complete)
        return self._extant

    def to_newick(self, which: str = "extant", *, annotate: bool = True) -> str | None:
        """Newick of the ``"extant"`` (default) or ``"complete"`` tree; ``None`` if it is empty.
        Leaves are ``g<id>``; with ``annotate`` internal nodes carry ``<kind>_n<species>``; branch
        lengths are time differences.

        The root carries one too, running from ``origination`` to where the root gene ended — the
        stem of the family, real time in which that founding gene existed. On the extant tree the
        root may be a node whose ancestors were suppressed; its branch still starts at ``origination``
        and so absorbs them, exactly as the species tree's extant root absorbs its own."""
        root = self.extant if which == "extant" else self.complete
        if root is None:
            return None
        return _to_newick(root, annotate, self.origination) + ";"


def write_gene_trees(gene_trees: dict[int, "GeneTree"], directory) -> None:
    """Write ``gene_tree_fam<family>_complete.nwk`` and ``…_extant.nwk``, one pair per family, into
    ``directory``. Every resolution writes them the same way, so the writer lives here rather than
    three times over. A family with no surviving copy has no extant tree and writes no ``_extant``
    file — its ``_complete`` one still records the lineages that died."""
    import pathlib

    d = pathlib.Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    for fam, gt in sorted(gene_trees.items()):
        (d / f"gene_tree_fam{fam}_complete.nwk").write_text(gt.to_newick("complete") + "\n")
        extant = gt.to_newick("extant")
        if extant is not None:
            (d / f"gene_tree_fam{fam}_extant.nwk").write_text(extant + "\n")


def gene_trees_from_events(events: list, tree) -> dict[int, GeneTree]:
    """Derive ``{family id: GeneTree}`` from the event log inside the complete ``tree``. Each event
    records a gene ending and its descendants beginning, so this is a direct parent→children read."""
    birth: dict[int, tuple[int, int]] = {}                # gene -> (family, species it lived on)
    children: dict[int, list[int]] = collections.defaultdict(list)   # gene -> descendant genes
    ended_by: dict[int, str] = {}                         # gene -> the event kind that ended it
    end_time: dict[int, float] = {}                       # gene -> when it ended
    origin_time: dict[int, float] = {}                    # founding gene -> when its family began
    lost: set[int] = set()
    roots: list[int] = []
    for e in events:
        if e.kind == "origination":
            birth[e.copy] = (e.family, e.lineage)
            origin_time[e.copy] = e.time                  # where the root's branch begins
            roots.append(e.copy)
        elif e.kind in ("duplication", "transfer", "speciation"):
            birth[e.copy] = (e.family, e.lineage)
            children[e.parent].append(e.copy)
            ended_by[e.parent] = e.kind
            end_time[e.parent] = e.time                   # the parent ends here (its children's birth)
        elif e.kind == "loss":
            lost.add(e.copy)
            end_time[e.copy] = e.time
    return {birth[root][0]: GeneTree(birth[root][0],
                                     _build(root, birth, children, ended_by, end_time, lost, tree.nodes),
                                     origin_time[root])
            for root in roots}


def _build(root, birth, children, ended_by, end_time, lost, nodes) -> GeneNode:
    """Assemble one family's tree bottom-up (iterative — gene trees run past the recursion limit)."""
    order: list[int] = []                                 # a pre-order: every gene before its descendants
    stack = [root]
    while stack:
        g = stack.pop()
        order.append(g)
        stack.extend(children.get(g, ()))
    built: dict[int, GeneNode] = {}
    for g in reversed(order):                             # descendants first, so a node's children exist
        _, species = birth[g]
        kids = children.get(g)
        if kids:                                          # internal: ended by dup / transfer / speciation
            node = GeneNode(ended_by[g], species, end_time[g], g)
            node.children = [built[c] for c in kids]
        elif g in lost:                                   # a dead leaf
            node = GeneNode("loss", species, end_time[g], g)
        else:                                             # alive at its species' tip
            s = nodes[species]
            node = GeneNode(s.fate, species, s.end_time, g)
        built[g] = node
    return built[root]


def _prune_to_extant(root: GeneNode) -> GeneNode | None:
    """A fresh tree keeping only lineages reaching an ``extant`` tip; degree-two nodes suppressed."""
    order: list[GeneNode] = []
    stack = [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(n.children)
    pruned: dict[int, GeneNode | None] = {}
    for n in reversed(order):
        if n.is_leaf:
            pruned[id(n)] = GeneNode(n.kind, n.species, n.time, n.copy) if n.kind == "extant" else None
            continue
        kept = [k for k in (pruned[id(c)] for c in n.children) if k is not None]
        if not kept:
            pruned[id(n)] = None
        elif len(kept) == 1:
            pruned[id(n)] = kept[0]                        # suppress the degree-two node
        else:
            node = GeneNode(n.kind, n.species, n.time, n.copy)
            node.children = kept
            pruned[id(n)] = node
    return pruned[id(root)]


def _to_newick(root: GeneNode, annotate: bool, origination: float) -> str:
    """Serialise iteratively (gene trees run deeper than CPython's C-stack recursion guard). The root
    is seeded with ``origination`` as its parent time, so it gets a branch length like every other
    node instead of the bare label that would drop the family's stem."""
    stack: list[list] = [[root, origination, 0, []]]       # [node, parent_time, next_child, child_strings]
    result = ""
    while stack:
        frame = stack[-1]
        node, parent_time, ci, parts = frame
        if ci < len(node.children):
            frame[2] = ci + 1
            stack.append([node.children[ci], node.time, 0, []])
            continue
        bl = f":{node.time - parent_time:.6g}"             # the root's parent time is `origination`
        if node.is_leaf:
            s = f"g{node.copy}{bl}"
        else:
            label = f"{node.kind}_n{node.species}" if annotate else ""
            s = f"({','.join(parts)}){label}{bl}"
        stack.pop()
        if stack:
            stack[-1][3].append(s)
        else:
            result = s
    return result
