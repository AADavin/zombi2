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

import contextlib
import sys
from collections import namedtuple

from zombi2.genomes.events import EventRecord, EventType, GeneOp


@contextlib.contextmanager
def _deep_recursion():
    """Temporarily lift the interpreter's recursion limit for the gene-tree walks.

    Every walk here (building the node tree, pruning it, serialising it to Newick) recurses on
    *gene-tree depth*, which is unbounded: a long duplication ladder is a routine outcome of a
    high-duplication run, and the default ~1000-frame limit turns it into a ``RecursionError``
    from an ordinary ``gene_trees()`` / ``reconciliations()`` call. Same treatment, and the same
    reason, as :func:`expand_trace`. Restores the previous limit on the way out, and nests safely.
    """
    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, 1_000_000))
    try:
        yield
    finally:
        sys.setrecursionlimit(old_limit)


class _Node:
    __slots__ = ("gid", "birth", "end", "kind", "children", "species", "is_loss", "is_extant",
                 "branch", "recipient", "donor_lost", "is_pseudo")

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
        self.is_pseudo = False                # pseudogenization: a unary gene->intergene state flip


# Events that terminate an incoming lineage and open child lineages (a from->to bifurcation). A
# CONVERSION is a donor bifurcation just like a DUPLICATION (both children stay on the same species
# branch — it is intra-genome, unlike TRANSFER), so it reconstructs through the same machinery.
_INTERNAL = (EventType.DUPLICATION, EventType.TRANSFER, EventType.SPECIATION, EventType.CONVERSION)


def _node_tree(records, gid2species, total_age) -> _Node | None:
    """Build the complete gene-lineage node tree from one family's event records.

    Returns the root :class:`_Node` (each node carrying ``birth``/``end`` times, the species
    ``branch`` its terminating event fired on, and — for survivors — the leaf ``species``), or
    ``None`` when nothing originated. Shared by :func:`build_gene_trees` and the
    sequence-evolution scaler, which both rebuild this genealogy from the log.
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
        elif ev in _INTERNAL or ev is EventType.PSEUDOGENIZATION:
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
        return None

    def build(gid: str) -> _Node:
        node = _Node(gid, birth.get(gid, 0.0))
        node.branch = branch.get(gid)
        if gid in children:  # terminated by an event
            node.end = end_time[gid]
            node.kind = kind[gid]
            if kind[gid] is EventType.LOSS:
                node.is_loss = True
            else:
                node.is_pseudo = kind[gid] is EventType.PSEUDOGENIZATION
                node.children = [build(c) for c in children[gid]]
        else:  # alive at the present -> extant leaf
            node.end = total_age
            node.species = gid2species.get(gid)
            node.is_extant = node.species is not None
        return node

    with _deep_recursion():          # gene-tree depth is unbounded (duplication ladders)
        return build(root)


def build_gene_trees(records, gid2species, total_age, annotate_species=False):
    """Return ``(complete_newick, extant_newick)`` for one family; extant may be ``None``.

    With ``annotate_species=True`` each internal gene node is labelled ``<gid>|<species-branch>``
    (the species branch the event happened on) instead of just ``<gid>``.
    """
    root_node = _node_tree(records, gid2species, total_age)
    if root_node is None:
        return None, None
    with _deep_recursion():          # _to_newick / _prune recurse on the same unbounded depth
        complete = _to_newick(root_node, annotate_species) + ";"
        pruned = _prune(root_node)
        extant = (_to_newick(pruned, annotate_species) + ";") if pruned is not None else None
    return complete, extant


def extant_species_from_records(families, species_tree) -> dict:
    """Reconstruct ``{gid: extant_species}`` for surviving gene lineages, from event records only.

    The live simulator reads this off its ``leaf_genomes``; when replaying from a written
    ``events_trace.tsv`` there are none, so recover it from the genealogy. A lineage that is
    never the *parent* of a later event survives to the present, and the species it lives in is
    the branch it was **born into**, fixed by the event that created it:

    * speciation on species ``X`` → its children go to ``X``'s daughter species, in tree order;
    * transfer → the first child (donor copy) stays on the donor branch, the second (the
      transferred copy) lands on the ``recipient``;
    * duplication / pseudogenization → the children stay on the event's own branch.

    Only lineages born into an **extant leaf** species are extant (others sit on extinct or
    internal branches and are pruned away), so exactly those are returned — matching
    :meth:`~zombi2.Genomes._gid_to_species`.
    """
    node_by_name = {n.name: n for n in species_tree.nodes_preorder()}
    out: dict = {}
    for records in families.values():
        parents: set = set()
        born: dict = {}          # gid -> species branch it was born into
        origin_gid = origin_branch = None
        for r in sorted(records, key=lambda rec: rec.time):
            gids = [op.gid for op in r.genes]
            if r.event is EventType.ORIGINATION:
                origin_gid, origin_branch = gids[0], r.branch
                continue
            parents.add(gids[0])
            children = gids[1:]
            if not children:
                continue
            if r.event is EventType.SPECIATION:
                daughters = node_by_name[r.branch].children
                for i, c in enumerate(children):
                    born[c] = daughters[i].name if i < len(daughters) else r.branch
            elif r.event is EventType.TRANSFER:
                born[children[0]] = r.branch                 # donor copy stays
                if len(children) > 1:
                    born[children[1]] = r.recipient          # transferred copy -> recipient
            else:                                            # duplication / pseudogenization
                for c in children:
                    born[c] = r.branch
        if origin_gid is not None and origin_gid not in parents:
            born.setdefault(origin_gid, origin_branch)       # originated and never split -> survivor
        for gid in set(born) - parents:
            node = node_by_name.get(born.get(gid))
            if node is not None and not node.children and node.is_extant:
                out[gid] = node.name
    return out


# ============ expand a compact (speciation-free) trace back into a full record list ============

def expand_trace(families, species_tree):
    """Replay a compact event trace against ``species_tree`` into a full record list.

    The ``output="trace"`` engine keeps gene ids across speciations and emits **no** speciation
    records, so a lineage id ``g`` is shared by every branch it spreads into; ``(g, branch)`` is
    then a unique lineage-instance key. This walks the species tree from each origination,
    re-inserting the speciation bifurcations the tree implies and assigning a **fresh** id to
    every lineage instance — reproducing exactly the O/S/D/T/L structure the full-log engine
    would have emitted (with reminted ids). The result feeds the ordinary reconstruction
    (:func:`build_gene_trees`, :func:`reconcile`, :func:`extant_species_from_records`, and the
    sequence-evolution scaler) unchanged.

    ``families`` is ``{family: [EventRecord]}`` of O/D/T/L records; returns the same shape with
    speciations inserted.
    """
    node_by_name = {n.name: n for n in species_tree.nodes_preorder()}
    return {fam: _expand_one(fam, recs, node_by_name) for fam, recs in families.items()}


def _expand_one(family, records, node_by_name):
    origins: list[tuple[str, str, float]] = []       # (gid, branch, time)
    event_at: dict[tuple[str, str], object] = {}     # (parent_gid, branch) -> the D/T/L record
    for r in records:
        if r.event is EventType.ORIGINATION:
            origins.append((r.genes[0].gid, r.branch, r.time))
        else:
            event_at[(r.genes[0].gid, r.branch)] = r

    emitted: list[EventRecord] = []
    counter = 0

    def fresh() -> str:
        nonlocal counter
        counter += 1
        return f"{family}.{counter}"

    def walk(gid: str, node, birth_time: float, *, origin: bool = False) -> str:
        """Emit records for the lineage instance (``gid`` on species ``node``) and return the
        fresh id assigned to it."""
        my = fresh()
        if origin:
            emitted.append(EventRecord(EventType.ORIGINATION, node.name, birth_time,
                                       [GeneOp(my, family, "origin")]))
        ev = event_at.get((gid, node.name))
        if ev is not None:                                   # consumed on this branch
            if ev.event is EventType.LOSS:
                emitted.append(EventRecord(EventType.LOSS, node.name, ev.time,
                                           [GeneOp(my, family, "lost")]))
            elif ev.event is EventType.DUPLICATION:
                a = walk(ev.genes[1].gid, node, ev.time)
                b = walk(ev.genes[2].gid, node, ev.time)
                emitted.append(EventRecord(EventType.DUPLICATION, node.name, ev.time,
                    [GeneOp(my, family, "parent"), GeneOp(a, family, "left"),
                     GeneOp(b, family, "right")]))
            elif ev.event is EventType.CONVERSION:           # donor bifurcation, both on this branch
                a = walk(ev.genes[1].gid, node, ev.time)
                b = walk(ev.genes[2].gid, node, ev.time)
                emitted.append(EventRecord(EventType.CONVERSION, node.name, ev.time,
                    [GeneOp(my, family, "parent"), GeneOp(a, family, "donor_copy"),
                     GeneOp(b, family, "converted_copy")],
                    donor=node.name, recipient=node.name))
            else:                                            # TRANSFER
                cont = walk(ev.genes[1].gid, node, ev.time)
                tc = walk(ev.genes[2].gid, node_by_name[ev.recipient], ev.time)
                emitted.append(EventRecord(EventType.TRANSFER, node.name, ev.time,
                    [GeneOp(my, family, "parent"), GeneOp(cont, family, "donor_copy"),
                     GeneOp(tc, family, "transfer_copy")],
                    donor=node.name, recipient=ev.recipient))
        else:                                                # flows to the end of the branch
            kids = node.children
            if len(kids) >= 2:                               # speciation at node.time
                children_ids = [walk(gid, c, node.time) for c in kids]
                emitted.append(EventRecord(EventType.SPECIATION, node.name, node.time,
                    [GeneOp(my, family, "parent"),
                     *(GeneOp(c, family, "child") for c in children_ids)]))
            elif len(kids) == 1:                             # degree-2 pass-through (non-binary)
                child_id = walk(gid, kids[0], node.time)
                emitted.append(EventRecord(EventType.SPECIATION, node.name, node.time,
                    [GeneOp(my, family, "parent"), GeneOp(child_id, family, "child")]))
            # leaf: survivor / extinct tip — no record; extant-ness recovered by the caller
            # (extant_species_from_records) from the emitted speciation edges.
        return my

    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, 1_000_000))         # species-tree depth can be large
    try:
        for gid, branch, t in origins:
            node = node_by_name.get(branch)
            if node is not None:
                walk(gid, node, t, origin=True)
    finally:
        sys.setrecursionlimit(old_limit)
    return emitted


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
    if node.is_pseudo:  # a gene->intergene state flip is always marked
        label = f"{node.gid}|G@{node.branch}" if annotate and node.branch else f"{node.gid}|G"
    else:
        label = f"{node.gid}|{node.branch}" if annotate and node.branch else node.gid
    return f"({inner}){label}:{_bl(node):.6g}"


def _prune(node: _Node) -> _Node | None:
    """Keep only lineages leading to an extant leaf; suppress degree-two nodes.

    A pseudogenization node is *not* suppressed even though it is degree-one — the gene->intergene
    state flip must stay visible in the extant tree.
    """
    if not node.children:
        if not node.is_extant:
            return None
        leaf = _Node(node.gid, node.birth)
        leaf.end, leaf.species, leaf.is_extant = node.end, node.species, True
        return leaf

    kept = [k for k in (_prune(c) for c in node.children) if k is not None]
    if not kept:
        return None
    if len(kept) == 1 and not node.is_pseudo:  # suppress: absorb this node's branch into the survivor
        survivor = kept[0]
        survivor.birth = node.birth
        return survivor

    inner = _Node(node.gid, node.birth)
    inner.end, inner.kind, inner.children, inner.branch = node.end, node.kind, kept, node.branch
    inner.is_pseudo = node.is_pseudo
    return inner


# ============ reconciliation: annotate the gene tree with its species mapping ============

#: One row of the reconciliation events table. ``event`` in {"S","D","T","L"}; ``species`` is
#: the species branch/node the event maps to; ``recipient`` is the transfer recipient (else
#: None); ``time`` is the event time; ``gene`` the gene-lineage id.
ReconEvent = namedtuple("ReconEvent", ["event", "species", "recipient", "time", "gene"])

#: A family's reconciliation. ``complete`` / ``extant`` are annotated Newick strings (or None);
#: ``events`` is the list of :class:`ReconEvent` read off the complete tree.
Reconciliation = namedtuple("Reconciliation", ["complete", "extant", "events"])

_EV_CHAR = {EventType.DUPLICATION: "D", EventType.TRANSFER: "T", EventType.SPECIATION: "S",
            EventType.PSEUDOGENIZATION: "G", EventType.CONVERSION: "C"}


def _prune_recon(node: "_Node"):
    """Extant lineages only: drop losses, suppress degree-2 nodes, keep species branch/recipient.

    A pseudogenization node (degree-one) is kept so the state flip stays in the extant tree.
    """
    if not node.children:
        if not node.is_extant:
            return None
        leaf = _Node(node.gid, node.birth)
        leaf.end, leaf.species, leaf.is_extant = node.end, node.species, True
        return leaf
    kept = [k for k in (_prune_recon(c) for c in node.children) if k is not None]
    if not kept:
        return None
    if len(kept) == 1 and not node.is_pseudo:
        kept[0].birth = node.birth
        return kept[0]
    inner = _Node(node.gid, node.birth)
    inner.end, inner.kind, inner.is_pseudo = node.end, node.kind, node.is_pseudo
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
    if node.is_pseudo:
        label = f"{node.branch}|G"
    elif node.kind is EventType.TRANSFER:
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
        elif ev in _INTERNAL or ev is EventType.PSEUDOGENIZATION:
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
                node.is_pseudo = kind[gid] is EventType.PSEUDOGENIZATION
                node.children = [build(c) for c in children[gid]]
        else:
            node.end = total_age
            node.species = gid2species.get(gid)
            node.is_extant = node.species is not None
        return node

    events: list[ReconEvent] = []

    def collect(n: "_Node") -> None:
        if n.is_loss:
            events.append(ReconEvent("L", n.branch, None, n.end, n.gid))
        elif n.children:
            events.append(ReconEvent(_EV_CHAR.get(n.kind, "?"), n.branch, n.recipient, n.end, n.gid))
        for c in n.children:
            collect(c)

    # build / collect / _prune_recon / _recon_newick all recurse on gene-tree depth
    with _deep_recursion():
        full = build(root)
        collect(full)
        pruned = _prune_recon(full)
        return Reconciliation(
            complete=_recon_newick(full) + ";",
            extant=(_recon_newick(pruned) + ";") if pruned is not None else None,
            events=events,
        )
