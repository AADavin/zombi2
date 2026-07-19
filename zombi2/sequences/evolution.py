"""Evolve a sequence down one gene tree.

A gene tree is a **timetree** — its branch lengths are time (``child.time - parent.time``). A
sequence living inside that gene converts time to *substitutions/site* through the substitution
**rate**: a branch spanning ``Δt`` accrues ``rate · Δt`` substitutions/site. Under the **strict
clock** (this slice) the rate is one number for the whole tree, so the branch length in subs/site is
just ``rate_base · Δt``. (The relaxed-clock family — a per-branch ``clock`` modifier riding the
species tree — and across-site ``+Γ`` are later slices; they scale this same branch length.)

The engine draws the root sequence from the model's stationary frequencies, then walks the tree from
root to tips: a child's sequence is sampled site-by-site from ``P(bl)[parent_state]``, where
``P(bl) = exp(Q·bl)`` (the reversible eigendecomposition in :mod:`.substitution_models`). Only the branch
*endpoints* are sampled — this gives the sequence at every node (the observable tip alignment and the
ancestral reconstructions) but not the individual substitution events, which are a later opt-in
``record=`` slice. Everything is vectorised over sites; a zero-length branch copies its parent.

The walk is **iterative** (an explicit stack): gene trees run deeper than CPython's C-stack recursion
guard on high-turnover families — the same reason :func:`~zombi2.genomes.gene_trees._to_newick` is
iterative — so recursion would crash on deep trees.
"""

from __future__ import annotations

import numpy as np

from .substitution_models import SubstitutionModel


def evolve_gene_tree(root, model: SubstitutionModel, length: int, rate_base: float,
                     rng: np.random.Generator) -> dict[int, np.ndarray]:
    """Evolve a sequence of ``length`` sites down the gene tree rooted at ``root`` (a
    :class:`~zombi2.genomes.gene_trees.GeneNode`).

    Returns ``{id(node): states}`` for **every** node — integer state arrays over the model's
    alphabet, keyed by object identity (gene-tree nodes carry no unique id, and identity is unique and
    stable for the run). The caller decodes and labels only the nodes it keeps. Deterministic given
    ``rng``; the branch length in substitutions/site is ``rate_base · (child.time - parent.time)``.
    """
    pi = model.stationary
    k = model.k
    root_states = rng.choice(k, size=length, p=pi).astype(np.int8)
    out: dict[int, np.ndarray] = {}
    pcache: dict[float, np.ndarray] = {}

    # Iterative pre-order. Each stack frame carries the parent's end time and states so a node's own
    # states are sampled when it is popped (strict pre-order rng consumption); children are pushed
    # reversed so they pop in forward order. ``parent_time is None`` marks the root.
    stack: list[tuple[object, float | None, np.ndarray | None]] = [(root, None, None)]
    while stack:
        node, parent_time, parent_states = stack.pop()
        if parent_states is None:
            states = root_states
        else:
            bl = rate_base * (node.time - parent_time)
            states = parent_states if bl <= 0.0 else _sample(parent_states, _p_for(pcache, model, bl), rng)
        out[id(node)] = states
        for child in reversed(node.children):
            stack.append((child, node.time, states))
    return out


def _p_for(cache: dict, model: SubstitutionModel, bl: float) -> np.ndarray:
    """``P(bl)``, cached by branch length rounded to 12 decimals (identical lengths reuse one matrix)."""
    key = round(float(bl), 12)
    P = cache.get(key)
    if P is None:
        P = model.p_matrix(key)
        cache[key] = P
    return P


def _sample(parent_states: np.ndarray, P: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Draw each site's child state from ``P[parent_state]`` (vectorised over sites)."""
    cum = P.cumsum(1)
    # P rows are clipped, not renormalised, so a row's final cumulative can land a hair below 1.0
    # (~1e-15). Pin it to 1.0 so a maximal draw can't slip past every threshold and make ``argmax``
    # silently return state 0. ``cum`` is a fresh array, so mutating it is safe.
    cum[:, -1] = 1.0
    r = rng.random(parent_states.shape[0])
    return (r[:, None] < cum[parent_states]).argmax(1).astype(np.int8)


__all__ = ["evolve_gene_tree"]
