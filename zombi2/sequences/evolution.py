"""Evolve a sequence down one gene tree.

A gene tree is a **timetree** — its branch lengths are time (``child.time - parent.time``). A
sequence living inside that gene converts time to *substitutions/site* through the substitution
**rate**: a branch spanning ``Δt`` accrues ``rate · Δt`` substitutions/site. Under the **strict
clock** the rate is one number for the whole tree, so the branch length in subs/site is just
``rate_base · Δt``. Any per-lineage variation — a relaxed clock drawn per species branch, a trait
driving the rate, or both — comes from the run's `Clock`, which converts a
stretch of species branch into subs/site and is asked for every branch length here.

**Across-site rate variation** scales the same branch length again, per **site**. A model decorated
with `SubstitutionModel.across_sites()` carries rate
classes, and this draws each site's class **once per family** — that is what across-site variation
means, as against a rate that changes along a branch: a slow site is slow down the whole tree. The
sites of one class then share a branch length and so share a cached transition matrix, which is why
the classes are discretised rather than one rate per site. A class at rate 0 (``+I``) is skipped
entirely, so those sites keep their founding state forever. A model without classes takes not one
extra draw, so a run that did not ask for this is bit-identical to what it was.

The engine draws the founding sequence from the model's stationary frequencies **at the family's
origination**, then walks the tree from root to tips: a child's sequence is sampled site-by-site from
``P(bl)[parent_state]``, where
``P(bl) = exp(Q·bl)`` (the reversible eigendecomposition in `substitution_models`). Only the branch
*endpoints* are sampled — this gives the sequence at every node (the observable tip alignment and the
ancestral reconstructions) but not the individual substitution events, which are not recorded at all —
there is no ``record=`` parameter and no substitution event log. Everything is vectorised over sites; a zero-length branch copies its parent.

The walk is **iterative** (an explicit stack): gene trees run deeper than CPython's C-stack recursion
guard on high-turnover families — the same reason `_to_newick()` is
iterative — so recursion would crash on deep trees.
"""

from __future__ import annotations

import numpy as np

from .clock import Clock  # noqa: F401  — the `clock:` annotation names it; kept out of __all__
from ..genomes.gene_trees import GeneNode
from .substitution_models import SubstitutionModel


def evolve_gene_tree(root, model: SubstitutionModel, length: int, rate_base: float,
                     clock: "Clock | None", rng: np.random.Generator,
                     origination: float,
                     founding: "np.ndarray | None" = None,
                     cdf_cache: "dict[tuple[int, float], np.ndarray] | None" = None,
                     models: "dict[int, SubstitutionModel] | None" = None,
                     record=None,
                     ) -> tuple[dict[int, np.ndarray], np.ndarray]:
    """Evolve a sequence of ``length`` sites down the gene tree rooted at ``root`` (a
    `GeneNode`), starting at ``origination``.

    Returns ``({id(node): states}, founding_states)``. The first is **every** node — integer state
    arrays over the model's alphabet, keyed by object identity (gene-tree nodes carry no unique id,
    and identity is unique and stable for the run); the caller decodes and labels only the nodes it
    keeps. The second is the sequence the family began with. Deterministic given ``rng``.

    The branch ending at a node lies on that node's species branch, so its length in substitutions/site
    is what the run's `Clock` makes of that stretch — the drawn per-lineage
    factor (shared across families), times a driver integrated over the stretch where the rate reads a
    trait. ``clock=None`` is the strict, undriven rate: ``rate_base · (node.time - parent.time)``.
    Asking the clock for the stretch rather than for a factor is what makes a driver that switches
    *inside* this branch come out exact — a factor per branch could not express it.

    The **root is an ordinary node** here, its parent time being ``origination``: a family exists from
    the moment it originates, and its founding gene evolves across the stem before whatever event ends
    it. Drawing the root's own sequence from the stationary frequencies instead would leave that
    stretch of the gene's life un-evolved, and give the phylogram a root branch nothing happened on.

    ``founding`` supplies the sequence the family began with (integer states, ``length`` long) instead
    of drawing it from the stationary frequencies — how a run founded from a real ``fasta=`` starts each
    block from the supplied DNA. It still evolves across the stem; at rate 0 it survives unchanged,
    which is what makes the assembled root genome equal the input.

    ``cdf_cache`` memoises the per-branch-length transition CDF (see `_cdf_for()`). It is keyed by
    the model's identity **and** the branch length, so one cache holds every model the run uses and
    the caller passes the same one everywhere. Branch lengths recur massively across
    blocks — a block passing straight through a species branch reuses that branch's length — so a
    run-wide cache computes a few hundred matrices where a per-block cache recomputes tens of thousands.
    ``None`` makes a fresh local cache (a standalone call is then self-contained). Across-site rate
    variation needed **no** change here: a class-scaled branch length is just another length of the
    same model, so the classes multiply the number of cache keys and nothing else.

    ``models`` gives each **species branch** its own model, so a clade can evolve under a different
    matrix rather than only at a different speed. Branch lengths do not change meaning: every model
    is normalised to one expected substitution per site per unit branch length, so ``rate · Δt`` says
    the same thing on every branch and the phylogram stays exact — with one caveat worth stating,
    because it is the transient such a study is usually about. That normalisation holds **at
    stationarity**. A lineage that has just entered a clade whose model has different equilibrium
    frequencies is not yet at them, so while its composition is still relaxing it accrues somewhat
    fewer substitutions than its nominal branch length claims. That is the ordinary price of a
    non-homogeneous model rather than a defect. ``None`` is the plain path: not one extra lookup and
    not one extra draw, so a run written before this is bit-identical.
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
    cache = {} if cdf_cache is None else cdf_cache

    # Which rate class each site belongs to, drawn once for the family and fixed for the whole tree.
    # `None` for a model with no across-site variation — and then not a single draw is taken, which
    # is what keeps every run written before this bit-identical.
    groups = _site_classes(model, length, rng)

    # Iterative pre-order. Each stack frame carries the parent's end time and states so a node's own
    # states are sampled when it is popped (strict pre-order rng consumption); children are pushed
    # reversed so they pop in forward order. The root's "parent" is the origination.
    stack: list[tuple[GeneNode, float, np.ndarray]] = [(root, origination, founding_states)]
    while stack:
        node, parent_time, parent_states = stack.pop()
        bl = (rate_base * (node.time - parent_time) if clock is None
              else clock.branch_length(rate_base, node.species, parent_time, node.time))
        # The model is a property of the SPECIES branch this stretch of gene tree sits on, so it is
        # read here beside the branch length and from the same `node.species`. One dict lookup per
        # branch; the matrices themselves are still built once per (model, length) and cached, and a
        # branch length belongs to one species branch and so to one model, so a per-clade run builds
        # the same number of matrices a one-model run does.
        m = model if models is None else models[node.species]
        if bl <= 0.0:
            states = parent_states
        elif record is not None:
            # A recorded branch walks the path rather than jumping to its end: same process, same
            # distribution at the end, and the only one of the two that can say what happened on
            # the way (`_record`).
            from ._record import walk as _walk
            states, rows = _walk(parent_states, m.Q, bl, rng, present=record.mask(node))
            record.add(node, parent_time, node.time - parent_time, rows, m.alphabet)
        elif groups is None:
            states = _sample(parent_states, _cdf_for(cache, m, bl), rng)
        else:
            # One class at a time, ascending, each at its own scaled branch length — so a class is
            # one `_sample` over its own sites rather than one per site, and its scaled length is
            # just another key in the same per-model cache. A rate-0 class is never sampled, so
            # invariant sites keep the parent's state (and hence the founding state) untouched.
            states = parent_states.copy()
            for rate, idx in zip(m.site_rates, groups):
                if rate > 0.0 and idx.size:
                    states[idx] = _sample(parent_states[idx], _cdf_for(cache, m, bl * rate), rng)
        out[id(node)] = states
        for child in reversed(node.children or ()):
            stack.append((child, node.time, states))
    return out, founding_states


def _site_classes(model: SubstitutionModel, length: int,
                  rng: np.random.Generator) -> "list[np.ndarray] | None":
    """Assign each of ``length`` sites to one of the model's rate classes — the site indices per
    class, in class order — or ``None`` when the model has no across-site variation.

    Drawn from the cumulated class shares by the same threshold idiom `_cdf_for()` uses, with the
    last cumulative pinned to 1.0 so a maximal draw cannot fall off the end; `np.searchsorted` rather
    than ``rng.choice(p=)`` because a whole-genome block would otherwise build a ``(length, classes)``
    temporary to pick from.

    ``None`` short-circuits the caller onto the path it took before this existed, and — the point —
    consumes no randomness at all, so a model with one class draws exactly what it always drew."""
    if len(model.site_rates) == 1:
        return None
    cum = np.cumsum(model.site_shares)
    cum[-1] = 1.0
    classes = np.searchsorted(cum, rng.random(length), side="right")
    return [np.flatnonzero(classes == c) for c in range(len(cum))]


def _cdf_for(cache: dict, model: SubstitutionModel, bl: float) -> np.ndarray:
    """The transition **CDF** over branch length ``bl`` — ``P(bl)`` cumulated along each row, cached by
    the model's identity **and** the branch length rounded to 12 decimals, so one cache holds every
    model a run uses and the same model over an identical length reuses one matrix. Caching the *cumulated*
    matrix, not ``P`` itself, means `_sample()` never recomputes the cumsum: it is the same for
    every node on a branch of this length, so it is built once here and reused for all of them."""
    key = (id(model), round(float(bl), 12))
    cum = cache.get(key)
    if cum is None:
        cum = model.p_matrix(key[1]).cumsum(1)
        # P rows are clipped, not renormalised, so a row's final cumulative can land a hair below 1.0
        # (~1e-15). Pin it to 1.0 so a maximal draw can't slip past every threshold and make the
        # ``argmax`` in _sample silently return state 0.
        cum[:, -1] = 1.0
        cache[key] = cum
    return cum


def _sample(parent_states: np.ndarray, cum: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Draw each site's child state from the transition CDF ``cum[parent_state]`` (vectorised over
    sites). ``cum`` is the row-cumulated, 1.0-pinned transition matrix from `_cdf_for()`."""
    r = rng.random(parent_states.shape)
    return (r[:, None] < cum[parent_states]).argmax(1).astype(np.int8)


__all__ = ["evolve_gene_tree"]
