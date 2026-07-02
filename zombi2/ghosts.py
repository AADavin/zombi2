"""Un-pruning: graft dead (ghost) lineages onto a reconstructed species tree.

The backward :func:`~zombi2.simulate_species_tree` yields only lineages with sampled extant
descendants. This module adds the dead ones back using the exact conditional law of the
complete tree given the reconstructed one: along each reconstructed edge, ghost lineages
attach as an inhomogeneous Poisson process with intensity ``λ(t)·E(t)``, where ``E`` is the
probability of leaving no sampled descendant (Nee, May & Harvey 1994; Stadler 2009;
Lambert & Stadler 2013). Each ghost roots a birth–death subtree conditioned on no sampled
descendant. See ``docs/ghost_lineages.md``.

Supports the constant-rate :class:`~zombi2.BirthDeath`/:class:`~zombi2.Yule` and the
:class:`~zombi2.EpisodicBirthDeath` (time-varying rates, incomplete sampling ``ρ<1``). With
``ρ=1`` ghosts are lineages that go extinct before the present; with ``ρ<1`` they also include
lineages alive today but unsampled. Once grafted, the forward gene-family simulator uses them
automatically as transfer donors/recipients ("transfer from the dead").
"""

from __future__ import annotations

import bisect

import numpy as np

from .species_model import BirthDeath, EpisodicBirthDeath
from .tree import Tree, TreeNode


class _ProcessView:
    """Model-agnostic access to the rates the grafting needs, as functions of age ``τ`` (time
    before the present): ``E(τ)`` (no-sampled-descendant prob), ``rates(τ) -> (λ, μ)``, the
    sampling fraction ``ρ``, and thinning bounds ``birth_bound`` (max λ) / ``rate_bound``
    (max λ+μ)."""

    __slots__ = ("E", "rates", "rho", "birth_bound", "rate_bound")

    def __init__(self, model, total_age):
        if isinstance(model, EpisodicBirthDeath):
            model.validate()
            if getattr(model, "_cache_A", None) != total_age:
                model._prepare(total_age)
            ages, egrid = model._ages, model._E
            shifts, births, deaths = model.shifts, model.birth, model.death

            def E(tau):
                return float(np.interp(tau, ages, egrid))

            def rates(tau):
                i = bisect.bisect_right(shifts, tau)
                return births[i], deaths[i]

            self.E = E
            self.rates = rates
            self.rho = model.rho
            self.birth_bound = max(births)
            self.rate_bound = max(b + d for b, d in zip(births, deaths))
        elif isinstance(model, BirthDeath):  # Yule is a subclass -> supported (μ=0 -> no ghosts)
            model.validate()
            b, d = model.birth, model.death
            self.E = model.extinction_prob
            self.rates = lambda tau: (b, d)
            self.rho = 1.0  # constant-rate model assumes complete sampling
            self.birth_bound = b
            self.rate_bound = b + d
        else:
            raise NotImplementedError(
                f"add_ghost_lineages supports BirthDeath/Yule and EpisodicBirthDeath, "
                f"not {type(model).__name__}"
            )


def _next_event(bt, total_age, view, rng):
    """First event of a lineage born at absolute time ``bt``, via thinning (exact for the
    piecewise-constant episodic rates). Returns ``("present"|"birth"|"death", time)``."""
    t = bt
    bound = view.rate_bound
    while True:
        t += rng.exponential(1.0 / bound)
        if t >= total_age:
            return "present", total_age
        lam, mu = view.rates(total_age - t)
        r = lam + mu
        if r > 0.0 and rng.random() < r / bound:
            return ("birth", t) if rng.random() < lam / r else ("death", t)


def _grow_ghost(t_start, total_age, view, rng, max_size):
    """Grow one birth–death subtree forward from ``t_start``, conditioned (by rejection) on
    leaving no sampled descendant. Returns the subtree root, or ``None`` to reject."""
    root = TreeNode(name="", time=t_start, is_extant=False)
    root._birth = t_start
    stack = [root]
    size = 1
    while stack:
        node = stack.pop()
        bt = node._birth
        del node._birth
        kind, et = _next_event(bt, total_age, view, rng)
        node.time = et
        if kind == "present":
            if rng.random() < view.rho:  # sampled tip -> violates the conditioning
                return None
            # else: an unsampled extant ghost (ρ<1); is_extant already False
        elif kind == "birth":
            for _ in range(2):
                child = TreeNode(name="", time=et, is_extant=False)
                child._birth = et
                node.add_child(child)
                stack.append(child)
            size += 2
            if size > max_size:  # runaway -> almost surely leaves a sampled tip
                return None
        # "death": a dead-before-present leaf; nothing more to grow
    return root


def add_ghost_lineages(
    tree: Tree,
    model,
    *,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_subtree_size: int = 4096,
    max_attempts: int = 100_000,
) -> Tree:
    """Graft dead (ghost) lineages onto a reconstructed ``tree`` **in place**, and return it.

    ``model`` should be the same species-tree model used to build the tree
    (:class:`~zombi2.BirthDeath`/:class:`~zombi2.Yule`, or
    :class:`~zombi2.EpisodicBirthDeath`). Ghost lineages attach along each edge as a Poisson
    process with intensity ``λ(t)·E(t)`` (τ = time before present) and each roots a
    birth–death subtree conditioned on leaving no sampled descendant; every ghost node gets
    ``is_extant=False`` and a ``ghost_*`` name. The sampled (extant) leaves are left untouched,
    so pruning back to them recovers the original tree.

    Parameters
    ----------
    max_subtree_size:
        Reject (and retry) a ghost subtree that exceeds this many lineages — a cheap guard
        against runaway supercritical growth (bias is negligible).
    max_attempts:
        Safety cap on rejection retries per ghost.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    view = _ProcessView(model, tree.total_age)
    total_age = tree.total_age
    birth_bound = view.birth_bound
    if birth_bound <= 0.0:  # no speciation -> no ghosts possible
        return tree

    # snapshot the reconstructed edges before mutating the tree
    edges = [n for n in tree.nodes_preorder() if n.parent is not None]

    for node in edges:
        parent = node.parent  # current upper endpoint of the (shrinking) segment
        t = parent.time
        t1 = node.time
        # inhomogeneous Poisson via thinning: propose at rate birth_bound, accept w.p.
        # λ(τ)·E(τ) / birth_bound (≤ 1 since λ ≤ birth_bound and E ≤ 1)
        while True:
            t += rng.exponential(1.0 / birth_bound)
            if t >= t1:
                break
            tau = total_age - t
            lam, _ = view.rates(tau)
            if rng.random() >= lam * view.E(tau) / birth_bound:
                continue  # thinned out
            ghost = None
            for _ in range(max_attempts):
                ghost = _grow_ghost(t, total_age, view, rng, max_subtree_size)
                if ghost is not None:
                    break
            if ghost is None:
                continue  # gave up on this attachment (extremely rare)
            # splice a binary junction between `parent` and `node` at time t
            junction = TreeNode(name="", time=t, is_extant=False)
            parent.children[parent.children.index(node)] = junction
            junction.parent = parent
            junction.add_child(node)   # the reconstructed lineage continues
            junction.add_child(ghost)  # the dead sibling
            parent = junction          # next ghost on this edge sits above `node`, below M

    # name every new (ghost / junction) node uniquely; leave reconstructed names intact
    k = 0
    for n in tree.nodes_preorder():
        if not n.name:
            n.name = f"ghost_{k}"
            k += 1
    return tree
