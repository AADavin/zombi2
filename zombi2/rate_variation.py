"""Relaxed molecular clocks — turn a chronogram (timetree) into a phylogram.

Our species trees and gene trees are **timetrees**: every branch length is an amount of
*time*. What sequence evolution actually accumulates along a branch is a number of
*substitutions per site*, which is time multiplied by an evolutionary **rate** that varies
across the tree. A *relaxed molecular clock* is the model of that variation. Applying one
rescales every branch from time into substitutions and yields a **phylogram**.

This module is a coherent family of clocks with one interface (:class:`Clock`):

* :class:`StrictClock` — a single rate everywhere: the baseline with *no* rate variation.
* Uncorrelated clocks — each branch draws an **i.i.d.** rate multiplier, so a branch's rate
  tells you nothing about its neighbours' (the PhyloBayes-style relaxed clocks):
  :class:`UncorrelatedLogNormalClock`, :class:`UncorrelatedGammaClock`,
  :class:`WhiteNoiseClock`.
* Autocorrelated clocks — a branch's rate is anchored to its parent's, so nearby lineages
  evolve at similar rates: :class:`AutocorrelatedLogNormalClock` (a geometric Brownian /
  Thorne–Kishino walk), :class:`CIRClock` (a Cox–Ingersoll–Ross diffusion), and
  :class:`RateVariation` (the discrete-bin, nearest-neighbour GTDB model).

Every clock turns any :class:`~zombi2.tree.Tree` — a species tree, or a gene tree loaded via
:func:`~zombi2.tree.read_newick` — into a :class:`RateScaledTree` (the phylogram) via
``scale(tree, seed=...)``, and is deterministic given a seed. The same clocks drive the
shared **lineage** clock inside :class:`~zombi2.SequenceEvolution`, which multiplies them by
a per-family speed and integrates them along reconciled gene trees.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .tree import Tree, TreeNode


@dataclass
class RateScaledTree:
    """The result of applying a :class:`Clock` to a timetree.

    ``branch_lengths`` maps each node to its substitution branch length; ``branch_rate`` is
    the (time-averaged) rate multiplier applied to each node's branch; ``segments`` is the
    per-branch list of within-branch pieces (a single piece for a constant-per-branch clock,
    several for a within-branch clock such as :class:`CIRClock`). Each piece is
    ``(rate, duration)`` for the general clocks, but ``(bin_index, duration)`` for
    :class:`RateVariation`, whose ``end_bin`` (the bin index at the end of each branch) is the
    only clock to populate that field; it is empty for every other clock. :meth:`to_newick`
    emits the phylogram.
    """

    tree: Tree
    branch_lengths: dict
    end_bin: dict = field(default_factory=dict)
    segments: dict = field(default_factory=dict)
    branch_rate: dict = field(default_factory=dict)

    def to_newick(self, include_internal_names: bool = True) -> str:
        def rec(node: TreeNode) -> str:
            if node.children:
                inner = ",".join(rec(c) for c in node.children)
                label = node.name if include_internal_names else ""
                s = f"({inner}){label}"
            else:
                s = node.name
            if node.parent is not None:
                s += f":{self.branch_lengths[node]:.10g}"
            return s

        return rec(self.tree.root) + ";"


class Clock:
    """Common interface for a relaxed molecular clock: chronogram → phylogram.

    A clock draws a substitution-rate multiplier for every branch of a timetree and rescales
    branch lengths from *time* into *expected substitutions per site*. Concrete clocks differ
    only in how the per-branch rate is drawn (see the module docstring). All share two entry
    points:

    ``scale(tree, seed=...)`` → :class:`RateScaledTree`
        the phylogram (a rescaled copy of ``tree``); call ``.to_newick()`` to emit it.
    ``lineage_segments(tree, rng)`` → ``(segments, avg)``
        the per-branch rate as piecewise-constant ``(rate, t0, t1)`` intervals in absolute
        time (keyed by branch/node **name**), plus the time-averaged rate per branch. This is
        what :class:`~zombi2.SequenceEvolution` integrates when the clock is used as the
        shared lineage clock. Constant-per-branch clocks return one interval per branch;
        within-branch clocks (:class:`CIRClock`, :class:`RateVariation`) return several.

    A clock is deterministic given ``seed`` (or an explicit ``numpy`` ``Generator``).
    """

    #: rate reported for the (zero-length) root branch and, for autocorrelated clocks, the
    #: value the process starts from at the root.
    root_rate: float = 1.0

    def _branch_rate(self, node: TreeNode, parent_rate: float, rng) -> float:
        """Draw the constant rate multiplier for the branch above ``node`` (subclass hook)."""
        raise NotImplementedError

    def lineage_segments(self, tree: Tree, rng):
        """Per-branch rate as absolute-time ``(rate, t0, t1)`` intervals, plus per-branch avg.

        The default is one interval per branch at a constant rate drawn by :meth:`_branch_rate`;
        within-branch clocks override this.
        """
        segments: dict = {}
        avg: dict = {}
        for node in tree.nodes_preorder():
            if node.parent is None:
                segments[node.name] = []
                avg[node.name] = self.root_rate
                continue
            r = self._branch_rate(node, avg[node.parent.name], rng)
            segments[node.name] = [(r, node.parent.time, node.time)]
            avg[node.name] = r
        return segments, avg

    def scale(self, tree: Tree, rng: np.random.Generator | None = None,
              seed: int | None = None) -> RateScaledTree:
        """Overlay the clock on ``tree`` and return the resulting phylogram."""
        if rng is None:
            rng = np.random.default_rng(seed)
        segments, avg = self.lineage_segments(tree, rng)

        branch_lengths: dict = {}
        branch_rate: dict = {}
        segs: dict = {}
        for node in tree.nodes_preorder():
            pieces = segments[node.name]
            branch_lengths[node] = sum(rate * (hi - lo) for rate, lo, hi in pieces)
            branch_rate[node] = avg[node.name]     # the time-averaged rate (exact; no round-trip)
            segs[node] = [(rate, hi - lo) for rate, lo, hi in pieces]
        return RateScaledTree(tree, branch_lengths, {}, segs, branch_rate)


class StrictClock(Clock):
    """The strict clock: a single rate ``rate`` on every branch (no rate variation).

    The phylogram is the chronogram uniformly stretched by ``rate`` — relative branch
    proportions are unchanged. This is the baseline every relaxed clock relaxes.
    """

    def __init__(self, rate: float = 1.0):
        if rate <= 0:
            raise ValueError(f"rate must be > 0, got {rate}")
        self.rate = float(rate)
        self.root_rate = self.rate

    def _branch_rate(self, node: TreeNode, parent_rate: float, rng) -> float:
        return self.rate


class UncorrelatedLogNormalClock(Clock):
    """Uncorrelated lognormal relaxed clock (Drummond et al. 2006).

    Each branch draws an **independent** rate multiplier from a lognormal distribution with
    mean ``mean``: ``rate = mean · exp(𝒩(−σ²/2, σ))``, so ``E[rate] = mean`` for any ``σ``.
    Larger ``sigma`` means more rate heterogeneity; ``sigma = 0`` is the strict clock. Because
    the draws are i.i.d., a branch's rate is uninformative about its neighbours'.
    """

    def __init__(self, sigma: float, mean: float = 1.0):
        if sigma < 0:
            raise ValueError(f"sigma must be >= 0, got {sigma}")
        if mean <= 0:
            raise ValueError(f"mean must be > 0, got {mean}")
        self.sigma = float(sigma)
        self.mean = float(mean)
        self.root_rate = self.mean

    def _branch_rate(self, node: TreeNode, parent_rate: float, rng) -> float:
        if self.sigma == 0.0:
            return self.mean
        return float(self.mean * math.exp(rng.normal(-0.5 * self.sigma ** 2, self.sigma)))


class UncorrelatedGammaClock(Clock):
    """Uncorrelated gamma relaxed clock (Drummond et al. 2006; PhyloBayes ``-ugam``).

    Each branch draws an **independent** rate from a gamma distribution with mean ``mean`` and
    variance ``mean²/shape``: ``rate ~ Gamma(shape, mean/shape)``. The single ``shape``
    parameter controls dispersion — large ``shape`` concentrates rates near ``mean`` (→ strict
    clock), small ``shape`` spreads them widely.
    """

    def __init__(self, shape: float, mean: float = 1.0):
        if shape <= 0:
            raise ValueError(f"shape must be > 0, got {shape}")
        if mean <= 0:
            raise ValueError(f"mean must be > 0, got {mean}")
        self.shape = float(shape)
        self.mean = float(mean)
        self.root_rate = self.mean

    def _branch_rate(self, node: TreeNode, parent_rate: float, rng) -> float:
        return float(rng.gamma(self.shape, self.mean / self.shape))


class WhiteNoiseClock(Clock):
    """White-noise relaxed clock (PhyloBayes ``-wn``).

    An uncorrelated clock in which the branch multiplier is the **integral of a white-noise
    rate** over the branch: it is gamma-distributed with mean ``mean`` and variance
    ``mean²·σ²/Δt`` inversely proportional to the branch duration ``Δt``. Long branches
    average the noise away (rate → ``mean``); short branches are highly variable. This
    branch-length dependence is what distinguishes it from :class:`UncorrelatedGammaClock`.
    ``sigma = 0`` is the strict clock.
    """

    def __init__(self, sigma: float, mean: float = 1.0):
        if sigma < 0:
            raise ValueError(f"sigma must be >= 0, got {sigma}")
        if mean <= 0:
            raise ValueError(f"mean must be > 0, got {mean}")
        self.sigma = float(sigma)
        self.mean = float(mean)
        self.root_rate = self.mean

    def _branch_rate(self, node: TreeNode, parent_rate: float, rng) -> float:
        dt = node.branch_length()
        if self.sigma == 0.0 or dt <= 0.0:
            return self.mean
        shape = dt / (self.sigma ** 2)                 # k
        scale = self.mean * (self.sigma ** 2) / dt     # θ  → mean=kθ, var=kθ²=mean²σ²/Δt
        return float(rng.gamma(shape, scale))


class AutocorrelatedLogNormalClock(Clock):
    """Autocorrelated lognormal relaxed clock (Thorne–Kishino–Painter 1998).

    The rate evolves down the tree as a geometric random walk anchored to the parent,
    ``R_child = R_parent · exp(𝒩(0, σ·√ℓ))``, where ``ℓ`` is the branch length in time. A
    child's rate is centred on its parent's, so nearby lineages have similar rates — the clock
    is *autocorrelated*. ``sigma = 0`` freezes the walk into a strict clock at ``root_rate``.
    This is the shared-lineage clock that :class:`~zombi2.SequenceEvolution`'s ``branch_sigma``
    selects.
    """

    def __init__(self, sigma: float, root_rate: float = 1.0):
        if sigma < 0:
            raise ValueError(f"sigma must be >= 0, got {sigma}")
        if root_rate <= 0:
            raise ValueError(f"root_rate must be > 0, got {root_rate}")
        self.sigma = float(sigma)
        self.root_rate = float(root_rate)

    def _branch_rate(self, node: TreeNode, parent_rate: float, rng) -> float:
        scale = self.sigma * math.sqrt(max(node.branch_length(), 0.0))
        drift = math.exp(rng.normal(0.0, scale)) if scale > 0 else 1.0
        return parent_rate * drift


class CIRClock(Clock):
    """Autocorrelated Cox–Ingersoll–Ross relaxed clock (Lepage et al. 2007).

    The instantaneous rate follows a mean-reverting CIR diffusion,

    .. math:: dr = \\theta\\,(\\mu - r)\\,dt + \\sigma\\sqrt{r}\\;dW,

    which stays strictly positive and pulls back toward the long-run mean ``mean`` (``μ``) at
    speed ``theta`` (``θ``); ``sigma`` (``σ``) sets the volatility. The path is simulated
    *within* each branch by full-truncation Euler–Maruyama on sub-steps of length at most
    ``max_step``, and a branch's substitution length is the integral of that piecewise path.
    Like every autocorrelated clock a child starts from where its parent ended, so adjacent
    branches have correlated rates; unlike the lognormal walk the rate also varies *within* a
    branch. With ``sigma = 0`` the diffusion is deterministic and, started at ``mean``, gives a
    strict clock.
    """

    def __init__(self, theta: float, sigma: float, mean: float = 1.0,
                 root_rate: float | None = None, max_step: float = 0.05):
        if theta < 0:
            raise ValueError(f"theta must be >= 0, got {theta}")
        if sigma < 0:
            raise ValueError(f"sigma must be >= 0, got {sigma}")
        if mean <= 0:
            raise ValueError(f"mean must be > 0, got {mean}")
        if max_step <= 0:
            raise ValueError(f"max_step must be > 0, got {max_step}")
        self.theta = float(theta)
        self.sigma = float(sigma)
        self.mean = float(mean)
        self.max_step = float(max_step)
        self.root_rate = self.mean if root_rate is None else float(root_rate)
        if self.root_rate <= 0:
            raise ValueError(f"root_rate must be > 0, got {self.root_rate}")

    def _simulate_branch(self, r0: float, t0: float, t1: float, rng):
        """Return ``(pieces, end_rate)`` for a branch simulated from rate ``r0``."""
        dur = t1 - t0
        if dur <= 0:
            return [], r0
        n = max(1, int(math.ceil(dur / self.max_step)))
        dt = dur / n
        sqrt_dt = math.sqrt(dt)
        pieces = []
        r = r0
        t = t0
        for _ in range(n):
            rp = max(r, 0.0)
            pieces.append((rp, t, t + dt))             # left value integrates the sub-step
            r = rp + self.theta * (self.mean - rp) * dt + self.sigma * math.sqrt(rp) * sqrt_dt * rng.normal()
            t += dt
        return pieces, max(r, 0.0)

    def lineage_segments(self, tree: Tree, rng):
        segments: dict = {}
        avg: dict = {}
        end: dict = {}
        for node in tree.nodes_preorder():
            if node.parent is None:
                segments[node.name] = []
                avg[node.name] = end[node.name] = self.root_rate
                continue
            pieces, end_rate = self._simulate_branch(end[node.parent.name],
                                                     node.parent.time, node.time, rng)
            segments[node.name] = pieces
            end[node.name] = end_rate
            length = node.branch_length()
            subst = sum(rate * (hi - lo) for rate, lo, hi in pieces)
            avg[node.name] = subst / length if length > 0 else end_rate
        return segments, avg


class RateVariation(Clock):
    """Autocorrelated rate variation over ordered rate bins (nearest-neighbour walk).

    The discrete-bin, within-branch clock used in the GTDB archaea study. An **ordered** set of
    rate **bins** is laid down and a continuous-time Markov process runs *along the phylogeny*,
    stepping only to an **adjacent bin** (index ± 1); because the rate changes gradually a
    single branch may be split into several **segments** in neighbouring bins, and its
    substitution length is ``Σ (segment_duration × bin_rate)``.

    Parameters
    ----------
    bins:
        **Ordered** rate multipliers, e.g. ``[0.25, 0.5, 1.0, 2.0, 4.0]`` (slow → fast),
        all > 0. Ordering matters: the process only moves between neighbouring bins.
    switch_rate:
        Rate at which an interior bin switches to a neighbour (per unit time). ``0`` freezes
        the process in its starting bin (a strict clock).
    up_bias:
        Probability that a switch goes to the *next-higher* bin (default ``0.5``, a
        symmetric walk). At the boundaries only the available direction is possible, so a
        symmetric walk has a uniform stationary distribution over bins.
    start:
        Index of the initial bin at the root (default: the middle bin).
    """

    def __init__(self, bins, switch_rate: float, up_bias: float = 0.5, start: int | None = None):
        self.bins = [float(b) for b in bins]
        if not self.bins or any(b <= 0 for b in self.bins):
            raise ValueError("bins must be a non-empty ordered list of positive rate multipliers")
        if switch_rate < 0:
            raise ValueError("switch_rate must be >= 0")
        if not (0.0 <= up_bias <= 1.0):
            raise ValueError("up_bias must be in [0, 1]")
        self.switch_rate = float(switch_rate)
        self.up_bias = float(up_bias)
        self.start = len(self.bins) // 2 if start is None else int(start)
        if not (0 <= self.start < len(self.bins)):
            raise ValueError(f"start must be a valid bin index in [0, {len(self.bins) - 1}]")

    def _rates(self, i: int) -> tuple[float, float]:
        """(up_rate, down_rate) out of bin ``i`` — a neighbour is unavailable at a boundary."""
        up = self.switch_rate * self.up_bias if i < len(self.bins) - 1 else 0.0
        down = self.switch_rate * (1.0 - self.up_bias) if i > 0 else 0.0
        return up, down

    def _simulate_branch(self, start_bin: int, duration: float, rng):
        """Return (segments, end_bin) for a branch of the given duration."""
        segments = []
        elapsed = 0.0
        current = start_bin
        while True:
            up, down = self._rates(current)
            total = up + down
            if total <= 0.0:  # frozen (switch_rate=0 or a single bin)
                segments.append((current, duration - elapsed))
                return segments, current
            dt = rng.exponential(1.0 / total)
            if elapsed + dt >= duration:
                segments.append((current, duration - elapsed))
                return segments, current
            segments.append((current, dt))
            elapsed += dt
            current = current + 1 if rng.random() < up / total else current - 1

    def scale(self, tree: Tree, rng: np.random.Generator | None = None,
              seed: int | None = None) -> RateScaledTree:
        """Overlay rate variation on ``tree`` and return the resulting phylogram."""
        if rng is None:
            rng = np.random.default_rng(seed)

        branch_lengths: dict = {}
        end_bin: dict = {}
        segments: dict = {}
        branch_rate: dict = {}

        root = tree.root
        end_bin[root] = self.start   # bin state at the root node
        branch_lengths[root] = 0.0
        segments[root] = []
        branch_rate[root] = self.bins[self.start]

        for node in tree.nodes_preorder():
            if node.parent is None:
                continue
            segs, eb = self._simulate_branch(end_bin[node.parent], node.branch_length(), rng)
            branch_lengths[node] = sum(self.bins[b] * d for b, d in segs)
            end_bin[node] = eb
            segments[node] = segs
            dur = node.branch_length()
            branch_rate[node] = branch_lengths[node] / dur if dur > 0 else self.bins[eb]

        return RateScaledTree(tree, branch_lengths, end_bin, segments, branch_rate)

    def lineage_segments(self, tree: Tree, rng):
        """Convert the within-branch bin path into absolute-time ``(rate, t0, t1)`` intervals."""
        scaled = self.scale(tree, rng=rng)
        segments: dict = {}
        avg: dict = {}
        for node in tree.nodes_preorder():
            if node.parent is None:
                segments[node.name] = []
                avg[node.name] = self.root_rate
                continue
            start = node.parent.time
            pieces = []
            for b, dur in scaled.segments[node]:
                pieces.append((self.bins[b], start, start + dur))
                start += dur
            segments[node.name] = pieces
            length = node.branch_length()
            avg[node.name] = (scaled.branch_lengths[node] / length if length > 0
                              else self.bins[scaled.end_bin[node]])
        return segments, avg
