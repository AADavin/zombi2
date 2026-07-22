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


def recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth, to_traj=None):
    """Pick a recipient lineage index (into ``alive``) from the candidate indices ``cand`` by the
    ``transfer_to`` rule: ``"uniform"`` gives every contemporaneous lineage equal weight; a
    :class:`Distance` weights by relatedness (closer relatives likelier); a
    :class:`~zombi2.rates.modifiers.DrivenBy` weights by the driver's value on each candidate, read
    from ``to_traj`` (the trajectory the engine resolved for that source).

    Returns ``None`` — "nobody can receive" — when a driven weighting gives **every** candidate a
    weight of 0. The caller must then make the event a **no-op**: leaving it unrecorded is exactly the
    model in which the transfer rate itself drops to zero while no eligible recipient exists, because
    rejecting an event whose acceptance depends only on the current state is Poisson thinning, and a
    rejected event changes nothing (see :func:`~zombi2.genomes._do_transfer`)."""
    if transfer_to == "uniform":
        return cand[int(rng.integers(len(cand)))]
    if isinstance(transfer_to, DrivenBy):
        # the choice slot: candidate k's weight is the mapping of the driver on lineage k right now,
        # normalised over the candidates. A weight of 0 means "cannot receive".
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
