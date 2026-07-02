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

from collections import namedtuple

from .events import EventType


class _Node:
    __slots__ = ("gid", "birth", "end", "kind", "children", "species", "is_loss", "is_extant",
                 "branch", "recipient", "donor_lost")

    def __init__(self, gid, birth):
        self.gid = gid
        self.birth = birth
        self.end = birth
        self.kind = None
        self.children: list["_Node"] = []
        self.species: str | None = None
        self.is_loss = False
        self.is_extant = False
        self.branch: str | None = None       # species branch/node where the event happened
        self.recipient: str | None = None     # transfer recipient species branch
        self.donor_lost = False               # transfer kept because only the copy survives


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


# ============ reconciliation: embed the extant gene tree in the species tree ============

#: One row of the reconciliation events table. ``event`` in {"S","D","T","L"}; ``species`` is
#: the species node/branch it maps to; ``recipient`` is the transfer recipient (else None);
#: ``gene`` is the gene-lineage id for D/T (None for inferred S/L).
ReconEvent = namedtuple("ReconEvent", ["event", "species", "recipient", "time", "gene"])


class _RNode:
    """A node of the reconciled tree: an event mapped to a species location."""
    __slots__ = ("event", "species", "recipient", "gene", "children", "time")

    def __init__(self, event, species, time, gene=None, recipient=None):
        self.event = event          # "S" | "D" | "T" | "L" | "Leaf"
        self.species = species
        self.recipient = recipient
        self.gene = gene
        self.time = time
        self.children: list["_RNode"] = []


def _copy_internal(node: _Node) -> _Node:
    new = _Node(node.gid, node.birth)
    new.end, new.kind, new.branch, new.recipient = node.end, node.kind, node.branch, node.recipient
    return new


def _prune_for_recon(node: _Node) -> _Node | None:
    """Prune to extant lineages, keeping **visible** transfers.

    Speciations / duplications with one surviving child are suppressed (degree-2) — their
    missing side is re-inferred as a loss during the species-tree walk. A transfer is kept
    whenever its *transferred copy* survives (observable as a horizontal move), even if the
    donor continuation dies; a transfer where only the donor continuation survives is
    invisible and suppressed.
    """
    if not node.children:
        if not node.is_extant:
            return None
        leaf = _Node(node.gid, node.birth)
        leaf.end, leaf.species, leaf.is_extant = node.end, node.species, True
        return leaf

    kids = [_prune_for_recon(c) for c in node.children]
    if node.kind is EventType.TRANSFER:  # children = [donor_continuation, transferred_copy]
        donor_c, transf_c = kids[0], kids[1]
        if transf_c is not None and donor_c is not None:
            keep = _copy_internal(node)
            keep.children = [donor_c, transf_c]
            return keep
        if transf_c is not None:                      # donor side died -> still a visible transfer
            keep = _copy_internal(node)
            keep.children = [transf_c]
            keep.donor_lost = True
            return keep
        if donor_c is not None:                       # transferred copy died -> transfer invisible
            donor_c.birth = node.birth
            return donor_c
        return None

    kept = [k for k in kids if k is not None]
    if not kept:
        return None
    if len(kept) == 1:                                # suppressed speciation/duplication
        kept[0].birth = node.birth
        return kept[0]
    keep = _copy_internal(node)
    keep.children = kept
    return keep


def reconcile(records, gid2species, species_tree):
    """Reconcile one family's extant gene tree against ``species_tree``.

    Returns ``(reconciled_newick, events)``: the extant gene tree embedded in the species
    tree, with a LOSS inferred at every species split the lineage skips, and ``events`` a
    list of :class:`ReconEvent` (S/D/T/L). ``(None, [])`` if nothing survives.

    Internal node labels are ``species|EVENT`` (transfers ``donor|T>recipient``), extant
    tips ``species|gid``, inferred losses ``LOSS|species``.
    """
    records = sorted(records, key=lambda r: r.time)
    children: dict[str, list[str]] = {}
    end_time: dict[str, float] = {}
    kind: dict[str, EventType] = {}
    birth: dict[str, float] = {}
    branch: dict[str, str] = {}
    recipient: dict[str, str] = {}
    root = origin_branch = None
    origin_time = 0.0

    for r in records:
        ev = r.event
        if ev is EventType.ORIGINATION:
            root, origin_branch, origin_time = r.genes[0].gid, r.branch, r.time
            birth.setdefault(root, r.time)
        elif ev in _INTERNAL:
            frm = r.genes[0].gid
            children[frm] = [op.gid for op in r.genes[1:]]
            end_time[frm], kind[frm], branch[frm] = r.time, ev, r.branch
            if ev is EventType.TRANSFER:
                recipient[frm] = r.recipient
            for c in children[frm]:
                birth[c] = r.time
        elif ev is EventType.LOSS:
            frm = r.genes[0].gid
            children[frm], end_time[frm], kind[frm] = [], r.time, ev

    if root is None:
        return None, []

    def build(gid: str) -> _Node:
        node = _Node(gid, birth.get(gid, 0.0))
        if gid in children:
            node.end, node.kind = end_time[gid], kind[gid]
            node.branch, node.recipient = branch.get(gid), recipient.get(gid)
            if kind[gid] is EventType.LOSS:
                node.is_loss = True
            else:
                node.children = [build(c) for c in children[gid]]
        else:
            node.species = gid2species.get(gid)
            node.is_extant = node.species is not None
        return node

    pruned = _prune_for_recon(build(root))
    if pruned is None:
        return None, []

    sp = {n.name: n for n in species_tree.nodes_preorder()}

    def other_child(parent, child):
        a, b = parent.children
        return b if a is child else a

    def descend_losses(top_name, bottom_name):
        """(speciation_node, lost_sibling) pairs passed descending species top -> bottom."""
        top, node, out = sp[top_name], sp[bottom_name], []
        while node is not top:
            p = node.parent
            out.append((p, other_child(p, node)))
            node = p
        return out  # bottom-up

    def child_species(g):
        return g.species if g.is_extant else g.branch

    def species_child_toward(spec, target_name):
        node = sp[target_name]
        while node.parent is not spec:
            node = node.parent
        return node

    def wrap(rnode, entry_name, target_name):
        """Insert the speciation+loss nodes skipped descending species entry -> target."""
        for spec, lost in descend_losses(entry_name, target_name):  # bottom-up
            s = _RNode("S", spec.name, spec.time)
            s.children = [rnode, _RNode("L", lost.name, spec.time)]
            rnode = s
        return rnode

    def build_recon(g, entry_name):
        if g.is_extant:
            return wrap(_RNode("Leaf", g.species, species_tree.total_age, gene=g.gid),
                        entry_name, g.species)
        loc = g.branch
        if g.kind is EventType.SPECIATION:
            r = _RNode("S", loc, sp[loc].time, gene=g.gid)
            r.children = [build_recon(c, species_child_toward(sp[loc], child_species(c)).name)
                          for c in g.children]
        elif g.kind is EventType.DUPLICATION:
            r = _RNode("D", loc, g.end, gene=g.gid)
            r.children = [build_recon(c, loc) for c in g.children]
        else:  # TRANSFER
            r = _RNode("T", loc, g.end, gene=g.gid, recipient=g.recipient)
            if g.donor_lost:  # only the transferred copy survives; donor side is a loss
                r.children = [build_recon(g.children[0], g.recipient), _RNode("L", loc, g.end)]
            else:
                r.children = [build_recon(g.children[0], loc),
                              build_recon(g.children[1], g.recipient)]
        return wrap(r, entry_name, loc)

    rtree = build_recon(pruned, origin_branch)

    events: list[ReconEvent] = []

    def collect(n):
        if n.event in ("S", "D", "T", "L"):
            events.append(ReconEvent(n.event, n.species, n.recipient, n.time, n.gene))
        for c in n.children:
            collect(c)
    collect(rtree)

    return _recon_newick(rtree, origin_time) + ";", events


def _recon_newick(node: _RNode, parent_time: float) -> str:
    bl = max(0.0, node.time - parent_time)
    if not node.children:
        name = f"LOSS|{node.species}" if node.event == "L" else f"{node.species}|{node.gene}"
        return f"{name}:{bl:.6g}"
    inner = ",".join(_recon_newick(c, node.time) for c in node.children)
    label = f"{node.species}|T>{node.recipient}" if node.event == "T" else f"{node.species}|{node.event}"
    return f"({inner}){label}:{bl:.6g}"
