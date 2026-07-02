"""Substitution-rate variation on a timetree (an autocorrelated relaxed clock).

Our species trees and gene trees are **timetrees** (branch lengths are time). To obtain
the branch lengths one would see from *sequence evolution*, we overlay a rate that varies
across the tree. This class implements the discrete-bin model used in the GTDB archaea
study:

* there is an **ordered** set of rate **bins** — multipliers, some above 1 (fast) and some
  below (slow);
* a continuous-time Markov process runs **along the phylogeny**, and it can only step to an
  **adjacent bin** (index ± 1). Because the rate can only change gradually, nearby lineages
  have similar rates — the process is *autocorrelated*;
* the current bin is inherited by both descendants at every node, so a single branch may be
  split into several **segments** in different (adjacent) bins; its substitution length is
  ``Σ (segment_duration × bin_rate)``.

The result is a *phylogram* (substitution lengths) built from the *chronogram* (times). It
applies to any :class:`~zombi2.tree.Tree` — a species tree, or a gene tree loaded via
:func:`~zombi2.tree.read_newick`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .tree import Tree, TreeNode


@dataclass
class RateScaledTree:
    """The result of applying :class:`RateVariation` to a timetree.

    ``branch_lengths`` maps each node to its substitution branch length; ``end_bin`` is the
    bin index at the end of each node's branch; ``segments`` is the per-branch list of
    ``(bin_index, duration)`` pieces. :meth:`to_newick` emits the phylogram.
    """

    tree: Tree
    branch_lengths: dict
    end_bin: dict
    segments: dict

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


class RateVariation:
    """Autocorrelated rate variation over ordered rate bins (nearest-neighbour walk).

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

        root = tree.root
        end_bin[root] = self.start   # bin state at the root node
        branch_lengths[root] = 0.0
        segments[root] = []

        for node in tree.nodes_preorder():
            if node.parent is None:
                continue
            segs, eb = self._simulate_branch(end_bin[node.parent], node.branch_length(), rng)
            branch_lengths[node] = sum(self.bins[b] * d for b, d in segs)
            end_bin[node] = eb
            segments[node] = segs

        return RateScaledTree(tree, branch_lengths, end_bin, segments)
