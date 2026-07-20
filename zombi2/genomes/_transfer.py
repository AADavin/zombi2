"""Transfer mechanics shared across genome resolutions — the ``transfer_to`` weighting.

A transfer's *rate* is an ordinary rate; what is special is **who receives** once it fires. That
mechanic is the same whether the genome is an unordered multiset or an ordered set of chromosomes,
so it lives here, imported by every resolution. ``transfer_to="uniform"`` gives every
contemporaneous lineage equal weight; :class:`Distance` weights by relatedness (closer relatives
likelier), which needs the tree's mean root-to-tip time to stay scale-free.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

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


def recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth) -> int:
    """Pick a recipient lineage index (into ``alive``) from the candidate indices ``cand`` by the
    ``transfer_to`` rule: ``"uniform"`` gives every contemporaneous lineage equal weight; a
    :class:`Distance` weights by relatedness (closer relatives likelier)."""
    if transfer_to == "uniform":
        return cand[int(rng.integers(len(cand)))]
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
