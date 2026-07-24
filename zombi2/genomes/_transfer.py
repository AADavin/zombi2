"""Transfer mechanics shared across genome resolutions — the ``transfer_to`` weighting.

A transfer's *rate* is an ordinary rate; what is special is **who receives** once it fires. That
mechanic is the same whether the genome is an unordered multiset or an ordered set of chromosomes,
so it lives here, imported by every resolution. ``transfer_to`` is the **choice slot** of SPEC §5 —
the numbers in it are per-candidate weights, normalised across the contemporaneous lineages, so they
change neither how fast nor how many transfers happen, only **who** receives. Three rules:

- ``"uniform"`` — every contemporaneous lineage gets equal weight;
- :class:`Distance` — weight by relatedness (closer relatives likelier), which needs the tree's mean
  root-to-tip time to stay scale-free;
- :class:`~zombi2.rates.modifiers.DrivenBy` — weight by **another level**: candidate ``k``'s weight is
  the mapping of the driver's value on lineage ``k`` at this instant (a trait that makes a lineage
  competent to take DNA up). Wired for the unordered resolution only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..rates.mapping import Between
from ..rates.modifiers import DrivenBy
from ..species import _weighted_index


@dataclass(frozen=True)
class Distance:
    """A ``transfer_to`` weighting by relatedness: a recipient at patristic distance ``d`` from the
    donor gets weight ``exp(-decay × d / depth)``, where ``depth`` is the tree's mean root-to-tip
    time — so ``decay`` is **scale-free** (in units of tree depth), meaning the same across trees of
    different absolute timescales. ``transfer_to="distance"`` is ``Distance(decay=1.0)``."""

    decay: float = 1.0

    def __post_init__(self) -> None:
        if isinstance(self.decay, bool) or not isinstance(self.decay, (int, float)) \
                or not math.isfinite(self.decay) or self.decay < 0:
            raise ValueError(f"Distance decay must be a finite non-negative number, got {self.decay!r}")


@dataclass
class Clades:
    """A ``transfer_to`` weighting by **named clades** — the topological, *donor-conditioned* sibling of
    :class:`Distance`. Each group is a clade of the species tree, and a
    :class:`~zombi2.rates.mapping.Between` kernel weights a candidate recipient by the **pair** (donor's
    clade, recipient's clade), so a transfer can be steered to run *between* two clades rather than
    within them — which the per-recipient weight of a :class:`~zombi2.rates.modifiers.DrivenBy` cannot
    express::

        transfer_to = Clades({"A": ["n12", "n27"], "B": 40},
                             Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0))

    A clade is named either by a **set of tips** (a list — the clade is the subtree below their MRCA) or
    by a single **node id** (an int, or an ``"n<id>"`` label — the clade is that node's whole subtree).
    Groups must be disjoint; a lineage in none of them is in the implicit group ``"rest"``, usable as a
    kernel key. Membership is read from the **tree** (a clade is a fact about the tree, not another
    level), so this is a topological rule like ``"distance"``, resolved once per run — **not** a
    ``DrivenBy`` coupling and needing no driver file."""

    groups: dict
    between: object

    def __post_init__(self) -> None:
        if not isinstance(self.groups, dict) or not self.groups:
            raise ValueError(
                "Clades needs a non-empty {label: clade} dict, where a clade is a list of tips (its "
                "MRCA's subtree) or a single node id — e.g. Clades({'A': ['n1', 'n2'], 'B': 40}, ...)")
        for label in self.groups:
            if not isinstance(label, str) or not label.strip():
                raise ValueError(f"clade labels must be non-empty strings, got {label!r}")
            if label == "rest":
                raise ValueError(
                    "'rest' is reserved for lineages in no named clade — name your clade something else")
        if isinstance(self.between, dict):
            self.between = Between(self.between)
        if not isinstance(self.between, Between):
            raise ValueError(
                "Clades takes a Between kernel (or a plain {(from, to): weight} dict) as its second "
                f"argument — the per-pair recipient weights — got {self.between!r}")
        unknown = self.between.groups() - (set(self.groups) | {"rest"})
        if unknown:
            raise ValueError(
                f"the Between kernel names groups {sorted(unknown)} that are not defined clades; "
                f"defined clades are {sorted(self.groups)} (plus the implicit 'rest')")


def _resolve_node(tree, spec) -> int:
    """A lineage reference — an ``int`` node id or an ``"n<id>"`` / ``"<id>"`` string — to a node id in
    ``tree`` (the form ``species_events.tsv`` and ``to_newick`` write). Raises if it is neither, or not
    a node of this tree."""
    if isinstance(spec, bool):
        raise ValueError(f"a clade reference must be a node id or 'n<id>' label, got {spec!r}")
    if isinstance(spec, int):
        nid = spec
    elif isinstance(spec, str):
        s = spec.strip()
        s = s[1:] if s[:1] == "n" else s
        try:
            nid = int(s)
        except ValueError:
            raise ValueError(
                f"clade reference {spec!r} is not a node id — name a lineage by its integer id or its "
                f"'n<id>' label (as species_events.tsv / to_newick write it)") from None
    else:
        raise ValueError(f"a clade reference must be a node id (int) or an 'n<id>' label, got {spec!r}")
    if nid not in tree.nodes:
        raise ValueError(f"clade reference n{nid} is not a lineage of this tree")
    return nid


def _mrca(tree, node_ids) -> int:
    """The most recent common ancestor of ``node_ids`` in ``tree`` — the deepest node on every one of
    their root-ward paths. A single id is its own MRCA."""
    def ancestors(i):
        chain = []
        while i is not None:
            chain.append(i)
            i = tree.nodes[i].parent
        return chain

    common = ancestors(node_ids[0])  # deepest-first: self, parent, …, root
    for other in node_ids[1:]:
        others = set(ancestors(other))
        common = [a for a in common if a in others]  # keep depth order, drop non-shared
    return common[0]


def _subtree(tree, root_id) -> set:
    """Every node in the clade below ``root_id`` (inclusive), extinct and internal lineages included."""
    out, stack = set(), [root_id]
    while stack:
        i = stack.pop()
        out.add(i)
        kids = tree.nodes[i].children
        if kids is not None:
            stack.extend(kids)
    return out


def resolve_groups(tree, groups) -> dict:
    """Paint every node of the complete ``tree`` with its clade label — ``{node_id: label}``, ``"rest"``
    for a lineage in no named clade. A clade named by a list of tips is the subtree below their MRCA; a
    clade named by a single node id is that node's subtree. Clades must be **disjoint** (an overlap —
    one clade nested in another — is refused). Computed once per run; membership is constant along a
    branch, so unlike a driver it adds no Gillespie breakpoints."""
    group_of = {i: "rest" for i in tree.nodes}
    claimed: dict[int, str] = {}
    for label, spec in groups.items():
        root = _mrca(tree, [_resolve_node(tree, t) for t in spec]) if isinstance(spec, (list, tuple)) \
            else _resolve_node(tree, spec)
        for i in _subtree(tree, root):
            if i in claimed and claimed[i] != label:
                raise ValueError(
                    f"clades {claimed[i]!r} and {label!r} overlap at n{i}; groups must be disjoint — "
                    f"is one clade nested inside the other?")
            claimed[i] = label
            group_of[i] = label
    return group_of


def recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth, to_traj=None, groups=None):
    """Pick a recipient lineage index (into ``alive``) from the candidate indices ``cand`` by the
    ``transfer_to`` rule: ``"uniform"`` gives every contemporaneous lineage equal weight; a
    :class:`Distance` weights by relatedness (closer relatives likelier); a :class:`Clades` weights by
    the kernel on (donor's clade, candidate's clade), read from the precomputed ``groups`` map; a
    :class:`~zombi2.rates.modifiers.DrivenBy` weights by the driver's value on each candidate, read
    from ``to_traj`` (the trajectory the engine resolved for that source) — and, with a
    :class:`~zombi2.rates.mapping.Between` mapping, by the donor's value too.

    Returns ``None`` — "nobody can receive" — when a driven weighting gives **every** candidate a
    weight of 0. The caller must then make the event a **no-op**: leaving it unrecorded is exactly the
    model in which the transfer rate itself drops to zero while no eligible recipient exists, because
    rejecting an event whose acceptance depends only on the current state is Poisson thinning, and a
    rejected event changes nothing (see :func:`~zombi2.genomes._do_transfer`)."""
    if transfer_to == "uniform":
        return cand[int(rng.integers(len(cand)))]
    if isinstance(transfer_to, Clades):
        # topological, donor-conditioned: candidate k's weight is the kernel on (donor's clade, k's
        # clade), read from the precomputed membership map. A weight of 0 means "cannot receive".
        g_d = groups[donor]
        weights = [transfer_to.between.weight(g_d, groups[alive[k]]) for k in cand]
        total = sum(weights)
        if total <= 0.0:
            return None
        return cand[_weighted_index(rng, weights, total)]
    if isinstance(transfer_to, DrivenBy):
        # the choice slot: candidate k's weight is the mapping of the driver on lineage k right now,
        # normalised over the candidates. A weight of 0 means "cannot receive". A Between mapping is
        # donor-conditioned — the weight reads the driver on the DONOR too — so a trait can steer
        # transfer between guilds exactly as Clades does between clades.
        if isinstance(transfer_to.mapping, Between):
            g_d = to_traj.value(donor, t)
            weights = [transfer_to.mapping.weight(g_d, to_traj.value(alive[k], t)) for k in cand]
        else:
            weights = [transfer_to.mapping.multiplier(to_traj.value(alive[k], t)) for k in cand]
        total = sum(weights)
        if total <= 0.0:
            return None
        return cand[_weighted_index(rng, weights, total)]
    # Distance: patristic distance d(donor, x) = 2·(t − t_mrca); scale-free in the tree depth. Mark
    # the donor's ancestor end-times once, then climb each candidate to its first marked ancestor.
    anc = {}
    p = tree.nodes[donor].parent
    while p is not None:
        anc[p] = tree.nodes[p].end_time
        p = tree.nodes[p].parent
    dists = []
    for k in cand:
        x = alive[k]
        if x == donor:
            dists.append(0.0)  # self (only reachable under self_transfer): closest
            continue
        q = x
        while q not in anc:
            q = tree.nodes[q].parent
        dists.append(2.0 * (t - anc[q]))
    dmin = min(dists)
    weights = [math.exp(-transfer_to.decay * (d - dmin) / depth) for d in dists]  # dmin: softmax-stable
    return cand[_weighted_index(rng, weights, sum(weights))]


def mean_root_to_tip(tree) -> float:
    """The tree's mean root-to-tip time — the timescale that makes :class:`Distance` decay scale-free.
    Over the extant tips (all leaves if none survive); 1.0 for a degenerate zero-height tree."""
    root_t = tree.nodes[tree.root].birth_time
    tips = tree.extant() or tree.leaves()
    depth = sum(n.end_time - root_t for n in tips) / len(tips)
    return depth if depth > 0 else 1.0
