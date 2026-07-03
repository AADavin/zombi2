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


def build_gene_trees(records, gid2species, total_age, annotate_species=False):
    """Return ``(complete_newick, extant_newick)`` for one family; extant may be ``None``.

    With ``annotate_species=True`` each internal gene node is labelled ``<gid>|<species-branch>``
    (the species branch the event happened on) instead of just ``<gid>``.
    """
    records = sorted(records, key=lambda r: r.time)

    children: dict[str, list[str]] = {}
    end_time: dict[str, float] = {}
    kind: dict[str, EventType] = {}
    birth: dict[str, float] = {}
    branch: dict[str, str] = {}
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
            branch[frm] = r.branch
            for c in tos:
                birth[c] = r.time
        elif ev is EventType.LOSS:
            frm = r.genes[0].gid
            children[frm] = []
            end_time[frm] = r.time
            kind[frm] = ev
            branch[frm] = r.branch

    if root is None:
        return None, None

    def build(gid: str) -> _Node:
        node = _Node(gid, birth.get(gid, 0.0))
        node.branch = branch.get(gid)
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
    complete = _to_newick(root_node, annotate_species) + ";"
    pruned = _prune(root_node)
    extant = (_to_newick(pruned, annotate_species) + ";") if pruned is not None else None
    return complete, extant


def _bl(node: _Node) -> float:
    return max(0.0, node.end - node.birth)


def _to_newick(node: _Node, annotate: bool = False) -> str:
    if not node.children:
        if node.is_loss:
            name = f"LOSS_{node.gid}"
        elif node.is_extant:
            name = f"{node.species}_{node.gid}"
        else:
            name = node.gid
        return f"{name}:{_bl(node):.6g}"
    inner = ",".join(_to_newick(c, annotate) for c in node.children)
    label = f"{node.gid}|{node.branch}" if annotate and node.branch else node.gid
    return f"({inner}){label}:{_bl(node):.6g}"


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
    inner.end, inner.kind, inner.children, inner.branch = node.end, node.kind, kept, node.branch
    return inner


# ============ reconciliation: annotate the gene tree with its species mapping ============

#: One row of the reconciliation events table. ``event`` in {"S","D","T","L"}; ``species`` is
#: the species branch/node the event maps to; ``recipient`` is the transfer recipient (else
#: None); ``time`` is the event time; ``gene`` the gene-lineage id.
ReconEvent = namedtuple("ReconEvent", ["event", "species", "recipient", "time", "gene"])

#: A family's reconciliation. ``complete`` / ``extant`` are annotated Newick strings (or None);
#: ``events`` is the list of :class:`ReconEvent` read off the complete tree.
Reconciliation = namedtuple("Reconciliation", ["complete", "extant", "events"])

_EV_CHAR = {EventType.DUPLICATION: "D", EventType.TRANSFER: "T", EventType.SPECIATION: "S"}


def _prune_recon(node: "_Node"):
    """Extant lineages only: drop losses, suppress degree-2 nodes, keep species branch/recipient."""
    if not node.children:
        if not node.is_extant:
            return None
        leaf = _Node(node.gid, node.birth)
        leaf.end, leaf.species, leaf.is_extant = node.end, node.species, True
        return leaf
    kept = [k for k in (_prune_recon(c) for c in node.children) if k is not None]
    if not kept:
        return None
    if len(kept) == 1:
        kept[0].birth = node.birth
        return kept[0]
    inner = _Node(node.gid, node.birth)
    inner.end, inner.kind = node.end, node.kind
    inner.branch, inner.recipient, inner.children = node.branch, node.recipient, kept
    return inner


def _recon_newick(node: "_Node") -> str:
    if not node.children:
        if node.is_loss:
            name = f"LOSS|{node.branch}"
        elif node.is_extant:
            name = f"{node.species}|{node.gid}"
        else:  # non-extant dead-end (e.g. a ghost / extinct tip)
            name = f"{node.branch}|{node.gid}"
        return f"{name}:{_bl(node):.6g}"
    inner = ",".join(_recon_newick(c) for c in node.children)
    if node.kind is EventType.TRANSFER:
        label = f"{node.branch}|T>{node.recipient}"
    else:
        label = f"{node.branch}|{_EV_CHAR.get(node.kind, '?')}"
    return f"({inner}){label}:{_bl(node):.6g}"


def reconcile(records, gid2species, total_age) -> "Reconciliation":
    """Reconcile one family's gene tree against the species tree.

    The simulator records the true species branch of every event, so reconciliation is exact
    **annotation** — no LCA/parsimony inference. Returns a :class:`Reconciliation` with two
    trees:

    * ``complete`` — the **complete** gene tree reconciled: every event, including the real
      ``LOSS|branch`` leaves (the ground-truth history);
    * ``extant`` — the **extant** (pruned) gene tree reconciled: only the observable lineages
      (the cherries), with no losses (degree-2 nodes suppressed);

    plus ``events``, the list of :class:`ReconEvent` (S/D/T/L) read off the complete tree. Node
    labels: ``branch|EVENT`` internal (``donor|T>recipient`` for transfers), ``species|gid``
    extant tips, ``LOSS|branch`` losses. ``Reconciliation(None, None, [])`` if nothing originated.
    """
    records = sorted(records, key=lambda r: r.time)
    children: dict[str, list[str]] = {}
    end_time: dict[str, float] = {}
    kind: dict[str, EventType] = {}
    birth: dict[str, float] = {}
    branch: dict[str, str] = {}
    recipient: dict[str, str] = {}
    root = None

    for r in records:
        ev = r.event
        if ev is EventType.ORIGINATION:
            root = r.genes[0].gid
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
            children[frm], end_time[frm], kind[frm], branch[frm] = [], r.time, ev, r.branch

    if root is None:
        return Reconciliation(None, None, [])

    def build(gid: str) -> "_Node":
        node = _Node(gid, birth.get(gid, 0.0))
        if gid in children:
            node.end, node.kind = end_time[gid], kind[gid]
            node.branch, node.recipient = branch.get(gid), recipient.get(gid)
            if kind[gid] is EventType.LOSS:
                node.is_loss = True
            else:
                node.children = [build(c) for c in children[gid]]
        else:
            node.end = total_age
            node.species = gid2species.get(gid)
            node.is_extant = node.species is not None
        return node

    full = build(root)
    events: list[ReconEvent] = []

    def collect(n: "_Node") -> None:
        if n.is_loss:
            events.append(ReconEvent("L", n.branch, None, n.end, n.gid))
        elif n.children:
            events.append(ReconEvent(_EV_CHAR.get(n.kind, "?"), n.branch, n.recipient, n.end, n.gid))
        for c in n.children:
            collect(c)
    collect(full)

    pruned = _prune_recon(full)
    return Reconciliation(
        complete=_recon_newick(full) + ";",
        extant=(_recon_newick(pruned) + ";") if pruned is not None else None,
        events=events,
    )
