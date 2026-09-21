"""A continuous trait's path **within** a branch, drawn after the run.

ZOMBI2 records a continuous trait only at the nodes. A reader that needs the trait *along* a branch
— a rate scaled by it, an optimum set by it — used to take the straight line between the two node
values. That line is the mean of the bridge between them, not a path: a real path wanders either
side of it, and under a non-linear mapping the excursions do not average out, so a smaller reading
step does not remove the bias (issue #454).

This module draws the path instead, conditioned on both node values. It runs **after** the trait is
simulated, off its own random stream, so the engine's draw order does not change at all: the same
seed gives the same ``node_values`` it always gave, path or no path.

The law. On a branch spanning ``[t0, t1]`` with the trait at ``a`` when the branch starts and ``b``
when it ends, write ``V(p, q)`` for the variance the trait accrues from ``p`` to ``q`` — exactly
what `~zombi2.traits.continuous._accrued_variance` computes, which is ``∫ σ²`` under Brownian motion
and the pull-weighted integral ``∫ e^{−2α(q−r)} σ²(r) dr`` under Ornstein–Uhlenbeck. With
``s = t − t0``, ``u = t1 − t``, ``m_t = θ + (a − θ)e^{−αs}``, ``m_T = θ + (a − θ)e^{−α(t1−t0)}`` and
``C = e^{−αu}·V(t0, t)``, the trait at ``t`` given both ends is normal with

    mean = m_t + C/V(t0, t1) · (b − m_T)
    var  = V(t0, t) − C²/V(t0, t1)

A whole branch is sampled by walking the stretch midpoints left to right: each stretch applies the
same two-point formula with ``a`` the value just drawn and ``t0`` the previous midpoint, one normal
per stretch. That is the exact joint law of the path at those times. Building it on the accrued
variance rather than on ``σ²·(t1−t0)`` is what keeps it exact under a ``changing_at`` schedule, a
`~zombi2.params.law.Drift` among lineages, and a diversity-dependent σ²: the same rate, LTT, driver
trajectories and inherited factor the branch itself used are carried here and integrated again.

Two shapes get no drawn path, and keep the straight line. A trait grown with ``regimes=`` or with a
driven optimum (``reverts_to=set_by(...)``) has a θ that is piecewise constant along the branch, so
the two-point formula above does not hold over the whole branch; conditioning stretch by stretch
across the optimum's own breakpoints is a separate piece of work. And a multivariate trait
(``correlation=``) is refused as a driver anyway. `BridgeLaw` is simply not built in those cases,
and the reader falls back to the line.

The path is also **written**, to ``trait_path.tsv``, so a run driven from files reads the same path
a run driven in memory does. `WrittenPath` reads it back and answers at any resolution.
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:                      # only for the annotation; importing it here would be a cycle
    from .continuous import _LTT


@dataclass(frozen=True)
class BridgeLaw:
    """What a continuous trait run needs to hand a reader so the reader can redraw the path.

    Everything here is what the branch walk itself used: ``rate`` (the σ² object, with its
    modifiers), ``pull`` (α, ``0.0`` for Brownian motion), ``theta`` (the optimum, ignored when
    ``pull`` is 0), ``inherited`` (each lineage's variable-rates factor, constant along its branch),
    ``ltt`` (the standing-diversity curve when σ² reads it) and ``trajs`` (the driver trajectories
    when σ² is itself driven).

    The branch endpoints are not here: they come from the run's node values and its event log, which
    a reader already has (`zombi2.params.conditioned.branch_starts`)."""

    rate: object
    pull: float
    theta: float
    inherited: dict
    ltt: object
    trajs: dict

    def variance(self, node_id: int, p: float, q: float) -> float:
        """``V(p, q)`` on lineage ``node_id`` — the variance the trait accrues from ``p`` to ``q``."""
        from .continuous import _accrued_variance

        if q <= p:
            return 0.0
        return _accrued_variance(self.rate, p, q, inherited=float(self.inherited.get(node_id, 1.0)),
                                 ltt=cast("_LTT | None", self.ltt), trajs=self.trajs,
                                 node_id=node_id, pull=self.pull)

    def sample(self, node_id: int, t0: float, t1: float, a: float, b: float,
               times, rng) -> list[tuple[float, float]]:
        """The trait at each time in ``times`` on lineage ``node_id``, drawn as the bridge from
        ``a`` at ``t0`` to ``b`` at ``t1``.

        ``times`` are the stretch midpoints, in order, strictly inside the branch. One normal is
        drawn per time, each conditioned on the value drawn before it and on ``b``. ``rng`` is this
        branch's own generator (`zombi2.rng.branch_stream`).

        Returns one ``(value, variance)`` pair per time: the trait there, and ``V`` from the
        previous point to this one, which is what a reader needs to refine between two of these
        points later (`WrittenPath`).

        A stretch whose accrued variance is zero — σ² is zero there, or the branch has no length
        left — takes the conditional mean with nothing added, which is the right answer for a
        process that cannot move."""
        alpha = float(self.pull)
        theta = float(self.theta)
        out: list[tuple[float, float]] = []
        p = t0
        x = float(a)
        for t in times:
            v_t = self.variance(node_id, p, t)
            v_full = self.variance(node_id, p, t1)
            if not (v_full > 0.0):
                out.append((x, v_t))
                p = t
                continue
            if alpha > 0.0:
                m_t = theta + (x - theta) * math.exp(-alpha * (t - p))
                m_full = theta + (x - theta) * math.exp(-alpha * (t1 - p))
                c = math.exp(-alpha * (t1 - t)) * v_t
            else:
                m_t = m_full = x
                c = v_t
            mean = m_t + c / v_full * (float(b) - m_full)
            var = v_t - c * c / v_full
            if var > 0.0:
                x = mean + float(rng.normal(0.0, math.sqrt(var)))
            else:
                x = mean                       # the two ends already pin this point
            out.append((x, v_t))
            p = t
        return out


@dataclass(frozen=True)
class WrittenPath:
    """The path read back from ``trait_path.tsv`` — the same reading interface `BridgeLaw` offers,
    built from written points rather than from the law that drew them.

    **The written path is the truth.** A run writes its path at one step; a reader may want another.
    Rather than redraw — which would give a different path, and two readers of one run that disagree
    — a reader derives what it needs from the two written points that surround the time it asks for,
    applying the same two-point bridge conditioned on those two values. So a reader at the written
    step gets the written values exactly (the bridge between two points, read at one of them, is
    that point), a coarser reader reads values consistent with them, and a finer one refines between
    them.

    ``points`` is ``{node: [(time, value, variance), …]}`` in time order, where ``variance`` is what
    the trait accrued since the previous point on that branch. The first point of a branch is its
    left endpoint and the last is the node value, so a branch is self-contained: nothing outside
    this file is needed to read anywhere on it.

    Within one written gap the process is treated as Brownian — σ² constant, the pull negligible —
    so the accrued variance grows linearly across the gap. That is exact for Brownian motion, and
    for everything else the error is bounded by the written step rather than by the branch, which is
    the whole difference #454 is about."""

    points: dict

    def sample(self, node_id: int, t0: float, t1: float, a: float, b: float,
               times, rng) -> list[tuple[float, float]]:
        """The trait at each time in ``times``, read off the written path. ``a``, ``b``, ``t0`` and
        ``t1`` are ignored: the written points carry the branch's endpoints themselves."""
        pts = self.points.get(node_id)
        out: list[tuple[float, float]] = []
        if not pts:
            return [(float(b), 0.0) for _ in times]
        times_w = [t for t, _, _ in pts]
        held_t, held_x = None, 0.0        # the last value drawn inside the gap being refined
        for t in times:
            k = bisect.bisect_right(times_w, t)
            if k == 0:                                   # before the first written point
                out.append((pts[0][1], 0.0))
                continue
            if k >= len(pts):                            # at or after the last written point
                out.append((pts[-1][1], 0.0))
                continue
            p, x_p, _ = pts[k - 1]
            q, x_q, v_gap = pts[k]
            span = q - p
            # Several asked-for times can fall in ONE written gap, when the reader's step is finer
            # than the step the path was written at. Each must be conditioned on the value drawn
            # just before it as well as on the gap's right end, exactly as `BridgeLaw.sample` walks
            # a branch. Conditioning each of them on the gap's two ends alone would draw them
            # independently, and independent draws at neighbouring times are not a path: they
            # scatter instead of wandering.
            if held_t is not None and p < held_t < t:
                p, x_p = held_t, held_x
                span = q - p
            f = 0.0 if span <= 0.0 else (t - p) / span
            v_span = v_gap if span <= 0.0 or q - pts[k - 1][0] <= 0.0 else \
                v_gap * span / (q - pts[k - 1][0])       # the gap's variance, its share of the span
            mean = x_p + f * (x_q - x_p)
            var = v_span * f * (1.0 - f)
            if var > 0.0:
                x = mean + float(rng.normal(0.0, math.sqrt(var)))
            else:
                x = mean
            out.append((x, v_span * f))
            held_t, held_x = t, x
        return out

    def variance(self, node_id: int, p: float, q: float) -> float:
        """``V(p, q)`` summed over the written gaps between ``p`` and ``q``, each gap's share taken
        in proportion to the time of it that falls inside."""
        pts = self.points.get(node_id)
        if not pts or q <= p:
            return 0.0
        total = 0.0
        for k in range(1, len(pts)):
            a_t, _, _ = pts[k - 1]
            b_t, _, v = pts[k]
            lo, hi = max(a_t, p), min(b_t, q)
            if hi > lo and b_t > a_t:
                total += v * (hi - lo) / (b_t - a_t)
        return total


__all__ = ["BridgeLaw", "WrittenPath"]
