"""Traits — a value riding the species tree (level 4).

A trait is a **value that rides the tree** — a body size, a habitat, a presence/absence — evolved
along the branches of a fixed tree. The result records the value at
**every** node (``node_values``, so the ancestral states are exact, not inferred) and, like the other
levels, an **event log** (``events``): a **discrete** trait mirrors the genome level exactly — its
transitions are timestamped events, the source of truth, and the per-branch stochastic character map
(``history``) is derived from them; a **continuous** trait diffuses with no along-branch events, so
its log holds the ``initial`` row and the jumps at speciation nodes
(just the ``initial`` row without ``at_speciation=``),
and ``node_values`` carries the diffusion. What keeps traits inside the one framework is that the
*ways* a value evolves reuse the same ``scope(base).verb(...)`` rate grammar (SPEC §5).

This is the **continuous** trait level — ``simulate_continuous`` — and its variants are the
same diffusion wearing different knobs, not a class each (SPEC §4):

- **Brownian motion**, the native process: over a branch the value moves by ``Normal(0, σ²·dt)``, so
  node-by-node in preorder it reproduces the exact tip law (Felsenstein 1985): the extant tips are
  multivariate-normal with variance ``σ² ×`` (root-to-tip depth) and covariance ``σ² ×`` (shared
  path length). ``rate`` is the variance-rate σ².
- **Ornstein–Uhlenbeck**: add ``reverts_to`` (the optimum θ) and ``pull`` (the strength α) and the
  diffusion is pulled toward θ — stabilizing selection. The exact per-branch transition is normal
  with mean ``θ + (x−θ)·e^{−α·dt}`` and variance ``σ²/(2α)·(1−e^{−2α·dt})``.
- **Early burst / ACDC**: give ``rate`` a skyline (``rate = PerLineage(σ²).changing_at({0: 1, 5: 0.2})``)
  and σ² changes through time — the *same* verb that gives the species tree its skyline.
  The per-branch variance is then the exact integral ``∫ σ²(t) dt`` over the branch.
- **Variable-rates BM** ("ClaDS for traits"): give ``rate`` a ``Drift`` law
  (``rate = PerLineage(σ²).varying_among('lineages', Drift(LogNormal(0.0, 0.3)))``) and σ² drifts branch-to-branch — each lineage inherits
  its parent's σ² times a lognormal kick at the split — the *same* law that drifts
  the species rate (ClaDS) and the autocorrelated clock, one level over. (``reverts_to`` / ``pull`` are
  OU function arguments that revert the trait *value*, **not** a modifier — a modifier acts on a
  *rate*.)
- **Diversity-dependent** (ecological limits): give ``rate`` the ``TotalDiversity`` driver
  (``rate = PerLineage(σ²).scaled_by(TotalDiversity(cap=100))``) and σ² slows as the clade fills — scaled by
  ``(1 − standing_diversity/cap)`` as the tree's lineages-through-time grows — the *same* driver
  that slows species diversification, read here off the fixed tree (one-way, tree → trait).
- **Driven by another level**: scale ``rate`` by a driver
  (``rate = PerLineage(σ²).scaled_by(habitat, {"aquatic": 3.0, "terrestrial": 1.0})``) and σ² reads a value
  grown first on this same tree: a second trait, or a genome's ``presence`` / ``completion`` — the
  *same* verb that drives a genome rate.
  One trait driving another is conditioning like any other (SPEC §3): the driver can be finished
  before the driven level starts, so it is two ordinary runs in order, handed over as the grown result or
  as its written ``trait_events.tsv``. A discrete driver switches *mid-branch*, so the per-branch
  variance is the integral across those pieces, not one sample per branch. The **discrete** engine
  takes it too: write ``switch`` as a rate expression and a trait's switch rate is driven the same way.

``rate`` thus takes the whole vocabulary — ``changing_at``, an inherited value, ``TotalDiversity``,
``scaled_by`` — like any other rate, and the verbs chain
(``PerLineage(σ²).changing_at({…}).varying_among('lineages', Drift(LogNormal(0.0, …)))``).

``rate`` is *per lineage*: each lineage carries its own independent diffusion, never pooled across the
tree — the engine evaluates the rate one lineage at a time (``lineages=1``), where the event levels
sum a per-unit rate over everything alive at once. The σ² modifiers **compose with**
``reverts_to`` / ``pull``: an early-bursting, drifting or driven σ² under an OU pull gives a
per-branch variance of ``∫ e^{−2α(t₁−s)}·σ²(s) ds``, the diffusion integral weighted by how much of
each moment's noise the pull has not yet erased by the end of the branch.

The **discrete** twin is ``simulate_discrete`` — a state switching along the tree (the Mk model),
simulated *exactly* by the Gillespie algorithm along every branch. Its ``events`` log (each transition
timestamped, on a lineage, ``from_state → to_state``) is the source of truth, exactly as at the genome
level; ``history`` (each node's ``(state, duration)`` segments) is the derived stochastic character
map. ``switch`` gives the rates (symmetric shortcut, ``{"a->b": rate}`` dict, or a ``k×k`` matrix), and
a switch rate may be written with ``scaled_by`` — the trait switching faster where another level says
so. The
**threshold** model (``liability=`` / ``threshold=``) reads a discrete state off a continuous Brownian
liability; the crossings are un-timed, so it carries no event log or map.

Also built: correlated traits (the ``correlation=`` overlay), jumps at speciation
(``at_speciation=``, in either engine, and under ``correlation=`` drawn through the same overlay),
and multi-optimum OU (``regimes=``, which takes the jumps too). ``reverts_to`` / ``pull`` also
apply to a correlated set, as **multivariate OU restricted to a diagonal drift** — each trait
reverts to its own optimum at its own strength, and the correlation rides in the diffusion rather
than in the reversion; a full drift matrix, where one trait's deviation pulls another, is refused
by name. SSE is **not** a trait model — it is trait↔species *joint*, Chapter 8 (Dependent runs),
`zombi2.joint`, which grows the trait with the tree: `discrete` for BiSSE and MuSSE, `continuous`
for QuaSSE. A diffusing driver moves at every instant, so that one run slices.
"""

from __future__ import annotations

from .continuous import ContinuousTrait, continuous, simulate_continuous
from .continuous import IMPLEMENTED_MODIFIERS  # noqa: F401  (re-exported for the CLI, not in __all__)
from .discrete import DiscreteTrait, discrete, simulate_discrete, simulate_traits
from .result import Change, TraitsResult

__all__ = ["simulate_continuous", "simulate_discrete", "simulate_traits", "TraitsResult",
           "Change", "DiscreteTrait", "discrete", "ContinuousTrait", "continuous"]
