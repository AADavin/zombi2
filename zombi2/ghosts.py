"""Un-pruning: graft dead (ghost) lineages onto a reconstructed species tree.

The backward :func:`~zombi2.simulate_species_tree` yields only lineages with sampled extant
descendants. This module adds the dead ones back using the exact conditional law of the
complete tree given the reconstructed one: along each reconstructed edge, ghost lineages
attach as an inhomogeneous Poisson process with intensity ``λ(t)·E(t)``, where ``E`` is the
probability of leaving no sampled descendant (Nee, May & Harvey 1994; Stadler 2009;
Lambert & Stadler 2013). Each ghost roots a birth–death subtree conditioned on no sampled
descendant. See ``docs/ghost_lineages.md``.

v1 supports the constant-rate :class:`~zombi2.BirthDeath` (complete sampling, ρ=1), so ghosts
are lineages that go extinct before the present. Once grafted, the forward gene-family
simulator uses them automatically as transfer donors/recipients ("transfer from the dead").
"""

from __future__ import annotations

import math

import numpy as np

from .species_model import BirthDeath
from .tree import Tree, TreeNode


def _grow_ghost(t_start, total_age, lam, mu, rho, rng, max_size):
    """Grow one birth–death subtree forward from ``t_start`` to the present, conditioned on
    leaving no sampled descendant (by rejection). Returns the subtree root, or ``None`` if the
    attempt should be rejected (a tip got sampled, or it exceeded ``max_size``)."""
    rate = lam + mu
    p_birth = lam / rate
    root = TreeNode(name="", time=t_start, is_extant=False)
    root._birth = t_start
    stack = [root]
    size = 1
    while stack:
        node = stack.pop()
        bt = node._birth
        del node._birth
        et = bt + rng.exponential(1.0 / rate)
        if et >= total_age:  # lineage reaches the present
            node.time = total_age
            if rng.random() < rho:  # sampled -> violates the conditioning
                return None
            # else: an unsampled extant ghost tip (is_extant already False)
        elif rng.random() < p_birth:  # speciation
            node.time = et
            for _ in range(2):
                child = TreeNode(name="", time=et, is_extant=False)
                child._birth = et
                node.add_child(child)
                stack.append(child)
            size += 2
            if size > max_size:  # runaway supercritical -> almost surely leaves a sampled tip
                return None
        else:  # extinction (a dead-before-present leaf)
            node.time = et
    return root


def add_ghost_lineages(
    tree: Tree,
    model: BirthDeath,
    *,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_subtree_size: int = 4096,
    max_attempts: int = 100_000,
) -> Tree:
    """Graft dead (ghost) lineages onto a reconstructed ``tree`` **in place**, and return it.

    ``model`` should be the same constant-rate :class:`~zombi2.BirthDeath` used to build the
    tree. Ghost lineages attach along each edge as a Poisson process with intensity
    ``birth · E(τ)`` (τ = time before present) and each roots a birth–death subtree conditioned
    on leaving no sampled descendant; every ghost node gets ``is_extant=False`` and a
    ``ghost_*`` name. The sampled (extant) leaves are left untouched, so pruning back to them
    recovers the original tree.

    Parameters
    ----------
    max_subtree_size:
        Reject (and retry) a ghost subtree that exceeds this many lineages — a cheap guard
        against runaway supercritical growth (bias is negligible).
    max_attempts:
        Safety cap on rejection retries per ghost.
    """
    if not isinstance(model, BirthDeath):  # Yule is a BirthDeath subclass -> supported
        raise NotImplementedError(
            "add_ghost_lineages currently supports the constant-rate BirthDeath/Yule model; "
            "episodic rates and incomplete sampling (ρ<1) are planned (see docs/ghost_lineages.md)"
        )
    model.validate()
    if rng is None:
        rng = np.random.default_rng(seed)

    lam, mu, rho = model.birth, model.death, 1.0
    total_age = tree.total_age

    # snapshot the reconstructed edges before mutating the tree
    edges = [n for n in tree.nodes_preorder() if n.parent is not None]

    for node in edges:
        parent = node.parent  # current upper endpoint of the (shrinking) segment
        t = parent.time
        t1 = node.time
        # inhomogeneous Poisson via thinning: propose at rate `lam`, accept w.p. E(τ) (≤1)
        while True:
            t += rng.exponential(1.0 / lam)
            if t >= t1:
                break
            if rng.random() >= model.extinction_prob(total_age - t):
                continue  # thinned out
            # place a ghost at time t: grow its conditioned subtree
            ghost = None
            for _ in range(max_attempts):
                ghost = _grow_ghost(t, total_age, lam, mu, rho, rng, max_subtree_size)
                if ghost is not None:
                    break
            if ghost is None:
                continue  # gave up on this attachment (extremely rare)
            # splice a binary junction M between `parent` and `node` at time t
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
