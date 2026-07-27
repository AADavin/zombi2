"""Traits — a value riding the species tree (level 4).

A trait is a **value that rides the tree** — a body size, a habitat, a presence/absence — evolved
along the branches of a fixed tree. The result records the value at
**every** node (``node_values``, so the ancestral states are exact, not inferred) and, like the other
levels, an **event log** (``events``): a **discrete** trait mirrors the genome level exactly — its
transitions are timestamped events, the source of truth, and the per-branch stochastic character map
(``history``) is derived from them; a **continuous** trait diffuses with no along-branch events, so
its log holds only the jumps at speciation nodes (empty without ``at_speciation=``),
and ``node_values`` carries the diffusion. What keeps traits inside the one framework is that the
*ways* a value evolves reuse the same ``scope(base) × modifiers`` rate grammar (SPEC §5).

This is the **continuous** trait level — ``simulate_continuous`` — and its three variants are the
same diffusion wearing different knobs, not three classes (SPEC §4):

- **Brownian motion**, the native process: over a branch the value moves by ``Normal(0, σ²·dt)``, so
  node-by-node in preorder it reproduces the exact tip law (Felsenstein 1985): the extant tips are
  multivariate-normal with variance ``σ² ×`` (root-to-tip depth) and covariance ``σ² ×`` (shared
  path length). ``rate`` is the variance-rate σ².
- **Ornstein–Uhlenbeck**: add ``reverts_to`` (the optimum θ) and ``pull`` (the strength α) and the
  diffusion is pulled toward θ — stabilizing selection. The exact per-branch transition is normal
  with mean ``θ + (x−θ)·e^{−α·dt}`` and variance ``σ²/(2α)·(1−e^{−2α·dt})``. (These are the same two
  knobs the CIR clock grows one level over — a shared vocabulary, not shared code.)
- **Early burst / ACDC**: give ``rate`` a ``OnTime`` skyline (``rate = σ² * mod.OnTime({0: 1, 5: 0.2})``)
  and σ² changes through time — the *same* ``OnTime`` modifier that gives the species tree its skyline.
  The per-branch variance is then the exact integral ``∫ σ²(t) dt`` over the branch.
- **Variable-rates BM** ("ClaDS for traits"): give ``rate`` an ``FromParent`` modifier
  (``rate = σ² * mod.FromParent(spread=0.3)``) and σ² drifts branch-to-branch — each lineage inherits
  its parent's σ² times a lognormal kick at the split — the *same* ``FromParent`` modifier that drifts
  the species rate (ClaDS) and the autocorrelated clock, one level over. (``reverts_to`` / ``pull`` are
  OU function arguments that revert the trait *value*, **not** a modifier — a rate modifier reverts a
  *rate*, which is the sequences level's CIR clock, a different mechanism.)
- **Diversity-dependent** (ecological limits): give ``rate`` a ``OnTotalDiversity`` modifier
  (``rate = σ² * mod.OnTotalDiversity(cap=100)``) and σ² slows as the clade fills — scaled by
  ``(1 − standing_diversity/cap)`` as the tree's lineages-through-time grows — the *same* ``OnTotalDiversity``
  modifier that slows species diversification, read here off the fixed tree (one-way, tree → trait).

``rate`` thus takes the whole modifier vocabulary — ``OnTime``, ``FromParent``, ``OnTotalDiversity`` — like any
other rate, and they compose (``σ² * OnTime({…}) * FromParent(spread=…)``).

``rate`` is *per lineage*: each lineage carries its own independent diffusion, never pooled across the
tree — the engine evaluates the rate one lineage at a time (``lineages=1``), where the event levels
sum a per-unit rate over everything alive at once. (OU with a time-varying σ² — the two knob-sets at
once — is not available; use one or the other.)

The **discrete** twin is ``simulate_discrete`` — a state switching along the tree (the Mk model),
simulated *exactly* by the Gillespie algorithm along every branch. Its ``events`` log (each transition
timestamped, on a lineage, ``from_state → to_state``) is the source of truth, exactly as at the genome
level; ``history`` (each node's ``(state, duration)`` segments) is the derived stochastic character
map. ``switch`` gives the rates (symmetric shortcut, ``{"a->b": rate}`` dict, or a ``k×k`` matrix). The
**threshold** model (``liability=`` / ``threshold=``) reads a discrete state off a continuous Brownian
liability; the crossings are un-timed, so it carries no event log or map.

Also built: correlated traits (the ``correlation=`` overlay), jumps at speciation (``at_speciation=``),
and multi-optimum OU (``regimes=``). SSE (BiSSE/MuSSE/QuaSSE) is
**not** a trait model — it is trait↔species *joint*, Part III.
"""

from __future__ import annotations

from .continuous import simulate_continuous
from .continuous import WIRED_MODIFIERS  # noqa: F401  (re-exported for the CLI, not in __all__)
from .discrete import DiscreteTrait, discrete, simulate_discrete
from .result import Change, TraitsResult

__all__ = ["simulate_continuous", "simulate_discrete", "TraitsResult", "Change",
           "DiscreteTrait", "discrete"]
