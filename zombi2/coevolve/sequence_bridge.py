"""Bridge: grammar couplings onto the sequence engine (the diamond's bottom tier).

Sequences evolve by molecular clocks (per-lineage substitution *rate*) and codon models (selection,
dN/dS ω) — a different engine from the genome rate machinery, so the sequence tier needs its own
bridge, not the :class:`~zombi2.coevolve.rate_bridge.CouplingModifier`.

This module handles the ``substitution_speed`` target-variable — a coupling on how fast a lineage's
sequences evolve. :class:`DriverClock` is a :class:`~zombi2.sequences.clocks.Clock` whose per-branch
rate is set by a grammar :class:`~zombi2.coevolve.grammar.Response` applied to a
:class:`~zombi2.coevolve.grammar.DriverSignal`, so it drops into
:class:`~zombi2.SequenceEvolution`'s ``lineage_segments`` contract unchanged. This realizes the
**T→Σ** (a trait sets substitution speed) and, with a gene-content driver, part of **G→Σ** edges.

The ``selection`` (ω) target-variable — a coupling on dN/dS via the codon models — is a separate,
heavier piece (a codon matrix per rate class); it is not in this module yet. Neither the clock nor
the codon machinery lives in ``zombi2.genomes.rates``, so this tier is independent of the rate rename.
See ``docs/design/coevolve-grammar.md`` §5.
"""

from __future__ import annotations

from zombi2.coevolve.grammar import DriverSignal, Response
from zombi2.sequences.clocks import Clock
from zombi2.tree import Tree


class DriverClock(Clock):
    """A molecular clock whose per-branch substitution rate is set by a grammar coupling on
    ``sequences.substitution_speed``.

    The rate on lineage ``b`` at time ``t`` is ``base_rate · response.rate_multiplier(driver_value)``,
    where the driver value is ``driver.value(b, t)``. Each branch is sub-segmented at the driver's
    interior change points (:meth:`DriverSignal.refresh_times`), so the rate tracks a within-branch
    driver change exactly. A null response (``Scalar(0)``) reduces this to a strict clock at
    ``base_rate``.

    Deterministic given the (already-simulated) driver — ``lineage_segments`` ignores its ``rng``.
    Satisfies the :class:`~zombi2.sequences.clocks.Clock` contract, so it is used exactly like any
    other clock (``.scale(tree)``, or as the shared clock in :class:`~zombi2.SequenceEvolution`).
    """

    def __init__(self, driver: DriverSignal, response: Response, *, base_rate: float = 1.0):
        if base_rate <= 0:
            raise ValueError(f"base_rate must be > 0, got {base_rate}")
        self.driver = driver
        self.response = response
        self.base_rate = float(base_rate)
        self.root_rate = self.base_rate

    def lineage_segments(self, tree: Tree, rng):
        segments: dict = {}
        avg: dict = {}
        for node in tree.nodes_preorder():
            if node.parent is None:
                segments[node.name] = []
                avg[node.name] = self.root_rate
                continue
            b0, b1 = node.parent.time, node.time
            cuts = sorted(t for (t, br) in self.driver.refresh_times(b0, b1)
                          if br == node.name and b0 < t < b1)
            bounds = [b0, *cuts, b1]
            segs = []
            rate_time = 0.0
            for s0, s1 in zip(bounds[:-1], bounds[1:]):
                if s1 <= s0:
                    continue
                r = self.base_rate * self.response.rate_multiplier(self.driver.value(node.name, s0))
                segs.append((r, s0, s1))
                rate_time += r * (s1 - s0)
            segments[node.name] = segs
            span = b1 - b0
            # A zero-length branch contributes no substitution length; report its instantaneous
            # driver-scaled rate (not a bare base_rate) so branch_rate matches sibling branches.
            avg[node.name] = (rate_time / span) if span > 0 else (
                self.base_rate * self.response.rate_multiplier(self.driver.value(node.name, b0)))
        return segments, avg
