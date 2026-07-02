"""Substitution-rate variation on a timetree (a Markov-modulated relaxed clock).

Our species trees and gene trees are **timetrees** (branch lengths are time). To obtain
the branch lengths one would see from *sequence evolution*, we overlay a rate that varies
across the tree. This class implements the discrete-bin, Markov-switching model used in the
GTDB archaea study:

* there is a set of rate **bins** — multipliers, some above 1 (fast) and some below (slow);
* a continuous-time Markov process runs **along the phylogeny**, switching bins at a
  constant rate; the current bin is inherited by both descendants at every node;
* a single branch may therefore be split into several **segments** in different bins; its
  substitution length is ``Σ (segment_duration × bin_rate)``.

The result is a *phylogram* (substitution lengths) built from the *chronogram* (times).
It applies to any :class:`~zombi2.tree.Tree` — a species tree, or a gene tree loaded via
:func:`~zombi2.tree.read_newick`.
"""

from __future__ import annotations

import math
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
    """Markov-modulated rate variation with discrete rate bins.

    Parameters
    ----------
    bins:
        Rate multipliers (e.g. ``[0.5, 1.0, 2.0]``), all > 0.
    switch_rate:
        Rate of the continuous-time process switching bins (per unit time). ``0`` means a
        single bin is used along the whole tree (a strict clock with one multiplier).
    weights:
        Probabilities of the bins — the stationary distribution the process switches to
        (and the root's initial bin). Defaults to uniform.
    """

    def __init__(self, bins, switch_rate: float, weights=None):
        self.bins = [float(b) for b in bins]
        if not self.bins or any(b <= 0 for b in self.bins):
            raise ValueError("bins must be a non-empty list of positive rate multipliers")
        if switch_rate < 0:
            raise ValueError("switch_rate must be >= 0")
        self.switch_rate = float(switch_rate)
        if weights is None:
            weights = [1.0 / len(self.bins)] * len(self.bins)
        if len(weights) != len(self.bins) or any(w < 0 for w in weights) or sum(weights) <= 0:
            raise ValueError("weights must be non-negative, one per bin, and sum to > 0")
        w = np.asarray(weights, dtype=float)
        self.weights = (w / w.sum()).tolist()
        self._cumw = np.cumsum(self.weights)

    def _draw_bin(self, rng) -> int:
        return int(np.searchsorted(self._cumw, rng.random(), side="right"))

    def _simulate_branch(self, start_bin: int, duration: float, rng):
        """Return (segments, end_bin) for a branch of the given duration."""
        segments = []
        elapsed = 0.0
        current = start_bin
        while True:
            if self.switch_rate <= 0.0:
                segments.append((current, duration - elapsed))
                return segments, current
            dt = rng.exponential(1.0 / self.switch_rate)
            if elapsed + dt >= duration:
                segments.append((current, duration - elapsed))
                return segments, current
            segments.append((current, dt))
            elapsed += dt
            current = self._draw_bin(rng)

    def scale(self, tree: Tree, rng: np.random.Generator | None = None,
              seed: int | None = None) -> RateScaledTree:
        """Overlay rate variation on ``tree`` and return the resulting phylogram."""
        if rng is None:
            rng = np.random.default_rng(seed)

        branch_lengths: dict = {}
        end_bin: dict = {}
        segments: dict = {}

        root = tree.root
        end_bin[root] = self._draw_bin(rng)   # the bin state at the root node
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
