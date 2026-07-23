"""Evolve a sequence down one gene tree.

A gene tree is a **timetree** — its branch lengths are time (``child.time - parent.time``). A
sequence living inside that gene converts time to *substitutions/site* through the substitution
**rate**: a branch spanning ``Δt`` accrues ``rate · Δt`` substitutions/site. Under the **strict
clock** (this slice) the rate is one number for the whole tree, so the branch length in subs/site is
just ``rate_base · Δt``. (The relaxed-clock family — a per-branch ``clock`` modifier riding the
species tree — and across-site ``+Γ`` are later slices; they scale this same branch length.)

The engine draws the founding sequence from the model's stationary frequencies **at the family's
origination**, then walks the tree from root to tips: a child's sequence is sampled site-by-site from
``P(bl)[parent_state]``, where
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
                     clock: "dict[int, float] | None", rng: np.random.Generator,
                     origination: float,
                     founding: "np.ndarray | None" = None) -> tuple[dict[int, np.ndarray], np.ndarray]:
    """Evolve a sequence of ``length`` sites down the gene tree rooted at ``root`` (a
    :class:`~zombi2.genomes.gene_trees.GeneNode`), starting at ``origination``.

    Returns ``({id(node): states}, founding_states)``. The first is **every** node — integer state
    arrays over the model's alphabet, keyed by object identity (gene-tree nodes carry no unique id,
    and identity is unique and stable for the run); the caller decodes and labels only the nodes it
    keeps. The second is the sequence the family began with. Deterministic given ``rng``.

    The branch ending at a node lies on that node's species branch, so its length in substitutions/site
    is ``rate_base · clock[node.species] · (node.time - parent.time)`` — the lineage clock (``clock``,
    one value per species branch, shared across families) rescales it. ``clock=None`` (the strict
    clock) uses factor 1 everywhere.

    The **root is an ordinary node** here, its parent time being ``origination``: a family exists from
    the moment it originates, and its founding gene evolves across the stem before whatever event ends
    it. Drawing the root's own sequence from the stationary frequencies instead would leave that
    stretch of the gene's life un-evolved, and give the phylogram a root branch nothing happened on.

    ``founding`` supplies the sequence the family began with (integer states, ``length`` long) instead
    of drawing it from the stationary frequencies — how a run seeded from a real ``fasta=`` starts each
    block from the supplied DNA. It still evolves across the stem; at rate 0 it survives unchanged,
    which is what makes the assembled root genome equal the input.
    """
    pi = model.stationary
    k = model.k
    if founding is None:
        founding_states = rng.choice(k, size=length, p=pi).astype(np.int8)
    else:
        founding_states = np.asarray(founding, dtype=np.int8)
        if founding_states.shape != (length,):
            raise ValueError(f"founding sequence is {founding_states.shape}, expected ({length},)")
    out: dict[int, np.ndarray] = {}
    pcache: dict[float, np.ndarray] = {}

    # Iterative pre-order. Each stack frame carries the parent's end time and states so a node's own
    # states are sampled when it is popped (strict pre-order rng consumption); children are pushed
    # reversed so they pop in forward order. The root's "parent" is the origination.
    stack: list[tuple[object, float, np.ndarray]] = [(root, origination, founding_states)]
    while stack:
        node, parent_time, parent_states = stack.pop()
        factor = 1.0 if clock is None else clock.get(node.species, 1.0)
        bl = rate_base * factor * (node.time - parent_time)
        states = parent_states if bl <= 0.0 else _sample(parent_states, _p_for(pcache, model, bl), rng)
        out[id(node)] = states
        for child in reversed(node.children):
            stack.append((child, node.time, states))
    return out, founding_states


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
