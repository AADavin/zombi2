"""Gene-tree reconstruction from the event log (post-processing, like ZOMBI 1).

Because every event re-mints lineage ids (see :mod:`zombi2.genome`), the per-family
event log is a complete genealogy: each event maps one incoming lineage id to its
outgoing ids (``genes[0]`` -> ``genes[1:]``). We rebuild a node-per-segment tree from
those edges and emit:

* the **complete** tree — every lineage, including losses;
* the **extant (pruned)** tree — only lineages ancestral to a surviving copy, with
  degree-two nodes suppressed.

Speciation events give internal bifurcations; duplications/transfers give the extra
copies; losses are terminal "LOSS" leaves in the complete tree; surviving copies are
extant leaves, labelled ``<species>_<gid>`` (species taken from the final genome of each
extant leaf).
"""

from __future__ import annotations

from .events import EventType


class _Node:
    __slots__ = ("gid", "birth", "end", "kind", "children", "species", "is_loss", "is_extant")

    def __init__(self, gid, birth):
        self.gid = gid
        self.birth = birth
        self.end = birth
        self.kind = None
        self.children: list["_Node"] = []
        self.species: str | None = None
        self.is_loss = False
        self.is_extant = False


_INTERNAL = (EventType.DUPLICATION, EventType.TRANSFER, EventType.SPECIATION)


def build_gene_trees(records, gid2species, total_age):
    """Return ``(complete_newick, extant_newick)`` for one family; extant may be ``None``."""
    records = sorted(records, key=lambda r: r.time)

    children: dict[str, list[str]] = {}
    end_time: dict[str, float] = {}
    kind: dict[str, EventType] = {}
    birth: dict[str, float] = {}
    root = None

    for r in records:
        ev = r.event
        if ev is EventType.ORIGINATION:
            root = r.genes[0].gid
            birth.setdefault(root, r.time)
        elif ev in _INTERNAL:
            frm = r.genes[0].gid
            tos = [op.gid for op in r.genes[1:]]
            children[frm] = tos
            end_time[frm] = r.time
            kind[frm] = ev
            for c in tos:
                birth[c] = r.time
        elif ev is EventType.LOSS:
            frm = r.genes[0].gid
            children[frm] = []
            end_time[frm] = r.time
            kind[frm] = ev

    if root is None:
        return None, None

    def build(gid: str) -> _Node:
        node = _Node(gid, birth.get(gid, 0.0))
        if gid in children:  # terminated by an event
            node.end = end_time[gid]
            node.kind = kind[gid]
            if kind[gid] is EventType.LOSS:
                node.is_loss = True
            else:
                node.children = [build(c) for c in children[gid]]
        else:  # alive at the present -> extant leaf
            node.end = total_age
            node.species = gid2species.get(gid)
            node.is_extant = node.species is not None
        return node

    root_node = build(root)
    complete = _to_newick(root_node) + ";"
    pruned = _prune(root_node)
    extant = (_to_newick(pruned) + ";") if pruned is not None else None
    return complete, extant


def _bl(node: _Node) -> float:
    return max(0.0, node.end - node.birth)


def _to_newick(node: _Node) -> str:
    if not node.children:
        if node.is_loss:
            name = f"LOSS_{node.gid}"
        elif node.is_extant:
            name = f"{node.species}_{node.gid}"
        else:
            name = node.gid
        return f"{name}:{_bl(node):.6g}"
    inner = ",".join(_to_newick(c) for c in node.children)
    return f"({inner}){node.gid}:{_bl(node):.6g}"


def _prune(node: _Node) -> _Node | None:
    """Keep only lineages leading to an extant leaf; suppress degree-two nodes."""
    if not node.children:
        if not node.is_extant:
            return None
        leaf = _Node(node.gid, node.birth)
        leaf.end, leaf.species, leaf.is_extant = node.end, node.species, True
        return leaf

    kept = [k for k in (_prune(c) for c in node.children) if k is not None]
    if not kept:
        return None
    if len(kept) == 1:  # suppress: absorb this node's branch into the survivor
        survivor = kept[0]
        survivor.birth = node.birth
        return survivor

    inner = _Node(node.gid, node.birth)
    inner.end, inner.kind, inner.children = node.end, node.kind, kept
    return inner
