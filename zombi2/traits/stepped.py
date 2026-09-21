"""Traits — the stepped engine: two traits that feed back into each other with no exact solution.

**experimental** (SPEC §9). Python-first, no CLI flag yet.

Two traits drive each other, so neither can be grown first: the pair is joint (SPEC §2–3). Where
both are discrete, the pair is one Markov chain over the product of their states and
`zombi2.traits.several.simulate_traits` solves it exactly. Where a **continuous** trait is one of
the two, there is no such solution, and this engine holds the traits still over short steps instead.

Two models run here:

- **A curved optimum.** Two continuous traits, each one's optimum a curved function of the other
  (``reverts_to=set_by("traits:<name>", f, step=...)`` on both).
- **A continuous trait with a discrete one.** The discrete trait paints the continuous trait's
  optimum (``regimes="traits:<name>"`` with ``reverts_to={state: θ}``), and the discrete trait's
  switch rate reads the continuous one (``scaled_by("traits:<name>", …, step=…)``).

A third model — two traits setting each other's σ² — is refused. It is not built.

**The walk.** Every branch is cut into slices of at most ``step``, the last one short, and the grid
is reset at every node so a slice never spans a split. In one slice:

1. every live driver is read at the **start** of the slice and held there. This is the only
   approximation in the run.
2. each discrete trait is advanced across the slice by the plain Gillespie walk under its frozen
   generator. Within the slice its own switches are exact, mid-slice ones included.
3. the continuous traits are advanced across the segments step 2 produced, one exact
   Ornstein–Uhlenbeck draw per segment under the optimum that segment holds. With no discrete trait
   a slice is one segment.

Step 3 is why a discrete driver of an optimum is **not** frozen: the optimum switches exactly where
the discrete trait switches, and the OU transition is exact between two switches. Only the discrete
trait's switch rate is frozen. `zombi2.joint._traits_sequences` walks the same shape.

A waiting time that runs past a slice boundary is not carried over. A Gillespie wait is memoryless,
so redrawing it in the next slice gives the same distribution.

``step`` rides on the **connection**, never on the run (`zombi2._runtime.slicing`), and it is what
bounds how stale a frozen driver gets. A run that freezes nothing needs none: where the only reading
is ``regimes=``, the optimum switches exactly where the discrete trait switches, so each branch is
walked in one slice and the whole run is exact.

**The error is first order in** ``step``: halve ``step`` and the error halves. Halving ``step``
also **changes the numbers** for a given seed, because it changes how many draws are taken. The two
runs are two samples from two nearby laws, not one answer refined. So convergence is checked
**across seeds**, never within one: run a set of seeds at ``step``, ``step/2`` and ``step/4``, pick
a summary the model is about (the mean tip value, the tip variance, the correlation between the
traits at the tips, the fraction of tips in a discrete state), and compare the shift between
consecutive steps against the seed-to-seed standard error. Chapter 8 writes the recipe out.

**The draw order is fixed**, so that adding a trait that draws nothing perturbs nothing: the root
start values in the order the traits were written; the nodes in `zombi2.traits._shared._preorder`;
per node the jumps at the split in written order, then the branch; per slice each discrete trait's
Gillespie walk in written order, then one ``rng.standard_normal(k)`` per segment for all ``k``
continuous traits at once, pushed through the square root of their covariance. The full ``k``-vector
is drawn every segment even where a trait's variance is zero, so the number of variates per segment
does not depend on the parameters.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, cast

from .._runtime.slicing import check_step, step_of
from ..params.connection import Driven
from ..params.driver import OnTime
from ..params.evaluate import describe
from ..params.parameter import Rate, as_rate
from ..params.scope import PerLineage
from ..rng import stream

from ._shared import _preorder
from .continuous import (ContinuousTrait, _accrued_variance, _at_speciation_jump_sd,
                         _driven_optimum)
from .discrete import DiscreteTrait, _driven_entries, _driven_q, _gillespie, _merge_runs
from .result import Change, TraitsResult

#: How far either side of its starting value an optimum curve is probed, in standard deviations of
#: the trait being read. Two either side covers the values that trait plausibly reaches.
_PROBE_SPREAD = 2.0


@dataclass
class _Cont:
    """One continuous trait, settled: what the walk reads and nothing else."""

    index: int                              # where it was written, which fixes its place in the draws
    name: str
    rate: Rate                              # σ², per lineage
    start: float
    alpha: float                            # the OU strength; 0.0 is Brownian motion
    jump_sd: float                          # the width of the jump at a split; 0.0 for none
    theta: float | None = None              # a fixed optimum
    theta_of: Callable[[object], float] | None = None   # value -> optimum, read off a trait
    reads: object | None = None             # the lookup key that curve's driver is threaded under
    regime_of: int | None = None            # the discrete trait whose state paints the optimum
    regime_theta: list = field(default_factory=list)   # that trait's optima, in its own state order


@dataclass
class _Disc:
    """One discrete trait, settled."""

    index: int
    name: str
    states: list
    entries: list                           # the off-diagonal switch rates, left as specs
    spec_start: object                      # the root state as written; None draws one
    shift: float                            # the chance of a hop at a split
    start: int = 0                          # filled once the run's stream exists


def _settle_continuous(spec: ContinuousTrait, index: int, by_name: dict) -> _Cont:
    """A `~zombi2.traits.continuous.ContinuousTrait` spec turned into what the walk reads.

    ``by_name`` maps ``"traits:<name>"`` to that trait's position in the run. The refusals here
    belong to this engine: a σ² reading another live trait is the model Adrián set aside, and the
    other σ² modifiers read a finished tree — its standing diversity, or a driver already grown on
    it — which is `~zombi2.traits.continuous.simulate_continuous`'s ground, not this one's."""
    r = as_rate(spec.rate, default_scope=PerLineage)
    assert r.scope is not None          # `as_rate` filled the level's default just above
    if r.scope is not PerLineage:
        raise ValueError(
            f"trait {spec.name!r}: rate is the variance-rate σ² per lineage — write PerLineage(...), "
            f"or a bare number, which is per lineage here; got a {r.scope.__name__} scope.")
    for m in r.modifiers:
        if isinstance(m, Driven) and isinstance(m.driver, str) and m.driver in by_name:
            raise ValueError(
                f"trait {spec.name!r}: σ² reads {m.driver!r}, a trait growing beside it. Two traits "
                f"that set each other's σ² are not built. The driven value this engine takes is the "
                f"optimum: write reverts_to=set_by({m.driver!r}, f, step=...) for a curved optimum, "
                f"or regimes={m.driver!r} with reverts_to={{state: θ}} for a discrete driver. To "
                f"scale σ² by a trait grown EARLIER, pass that run's result rather than its name, "
                f"which is conditioning and runs in traits.simulate_continuous.")
        if not isinstance(m, OnTime):
            raise ValueError(
                f"trait {spec.name!r}: σ² carries {describe(m)}, and this engine takes a bare "
                f"variance-rate or a changing_at(...) skyline. The others read a finished tree — "
                f"its standing diversity, or a driver already grown on it — which is "
                f"traits.simulate_continuous(tree, rate=...).")
    r.check_one_base(f"trait {spec.name!r}: rate")

    out = _Cont(index=index, name=cast(str, spec.name), rate=r, start=float(spec.start),
                alpha=0.0 if spec.pull is None else float(spec.pull),
                jump_sd=_at_speciation_jump_sd(spec.at_speciation))
    if spec.regimes is not None:
        if spec.regimes not in by_name:
            raise ValueError(
                f"trait {spec.name!r} takes its regimes from {spec.regimes!r}, which is not a trait "
                f"in this run. The traits here are {sorted(by_name)}.")
        out.regime_of = by_name[spec.regimes]
    elif isinstance(spec.reverts_to, Driven):
        out.theta_of = _driven_optimum(spec.reverts_to)
        out.reads = spec.reverts_to.key
    elif spec.reverts_to is not None:
        out.theta = float(cast(float, spec.reverts_to))
    return out


def _settle_discrete(spec: DiscreteTrait, index: int) -> _Disc:
    """A `~zombi2.traits.discrete.DiscreteTrait` spec turned into what the walk reads. The root
    state is left as written: drawing it needs the run's stream, and the order the draws come in is
    part of what this engine promises."""
    states = list(spec.states)
    shift = 0.0 if spec.at_speciation is None else float(cast(float, spec.at_speciation))
    if spec.at_speciation is not None and (isinstance(spec.at_speciation, bool)
            or not isinstance(spec.at_speciation, (int, float)) or not 0.0 <= shift <= 1.0):
        raise ValueError(f"at_speciation must be a probability in [0, 1] (the shift chance), got "
                         f"{spec.at_speciation!r}")
    return _Disc(index=index, name=cast(str, spec.name), states=states,
                 entries=_driven_entries(states, spec.switch), spec_start=spec.start, shift=shift)


def _wire_regimes(plans: list, specs: list) -> None:
    """Point each regime-driven optimum at the discrete trait that paints it, and check the table
    covers that trait's whole alphabet. An optimum missing a state would be read as nothing at the
    moment the trait switched into it."""
    for plan in plans:
        if not isinstance(plan, _Cont) or plan.regime_of is None:
            continue
        painter = plans[plan.regime_of]
        if not isinstance(painter, _Disc):
            raise ValueError(
                f"trait {plan.name!r} takes its regimes from {painter.name!r}, which is a continuous "
                f"trait. Regimes are painted by a discrete trait; an optimum read off a continuous "
                f"one is reverts_to=set_by(...).")
        table = cast(dict, specs[plan.index].reverts_to)
        missing = [s for s in painter.states if s not in table]
        if missing:
            raise ValueError(
                f"trait {plan.name!r}: reverts_to is missing an optimum for regime state(s) "
                f"{sorted(missing)} of {painter.name!r}.")
        plan.regime_theta = [float(table[s]) for s in painter.states]


# --- the straight-line case, which has an exact answer and is refused here -------------------------

def _sigma2_at(plan: _Cont) -> float:
    """The trait's σ² as one number, for a message. A skyline has no single one, and the value it
    opens at is better in a message than nothing."""
    return plan.rate.effective(lineages=1, time=0.0)


def _probe_span(plan: _Cont, tallest: float) -> tuple[float, float]:
    """The stretch of values ``plan`` plausibly reaches: its starting value, two standard deviations
    either side. Under a pull that spread is the stationary ``σ²/(2α)``; with no pull it is what a
    Brownian trait has by the tallest tip. A trait that cannot move is given a unit span, so the
    probes still straddle its starting value."""
    total = _accrued_variance(plan.rate, 0.0, tallest, pull=0.0)
    sd = math.sqrt(total / (2.0 * plan.alpha * tallest)) if plan.alpha > 0.0 else math.sqrt(total)
    if not (sd > 0.0) or not math.isfinite(sd):
        sd = 1.0
    return plan.start - _PROBE_SPREAD * sd, plan.start + _PROBE_SPREAD * sd


def _affine_fit(fn, lo: float, hi: float) -> tuple[float, float] | None:
    """The line ``(intercept, slope)`` a curve follows over ``[lo, hi]``, or ``None`` where it
    follows none. Five probes: the line through the outer two, then every probe checked against it
    to 1e-9 relative. Five is the fewest that sees a bend at the quarter points as well as at the
    middle."""
    if hi <= lo:
        return None
    xs = [lo + (hi - lo) * j / 4.0 for j in range(5)]
    ys = [float(fn(x)) for x in xs]
    slope = (ys[4] - ys[0]) / (xs[4] - xs[0])
    intercept = ys[0] - slope * xs[0]
    for x, y in zip(xs, ys):
        if abs(y - (intercept + slope * x)) > 1e-9 * max(1.0, abs(y)):
            return None
    return intercept, slope


def _g(v: float) -> str:
    """A number as a reader would paste it back."""
    return f"{v:.10g}"


def _refuse_affine(plans: list, tallest: float) -> None:
    """Refuse a pair of straight optima, and name the exact call to write instead.

    Two continuous traits whose optima are straight lines in each other are one multivariate
    Ornstein–Uhlenbeck with a full drift matrix, which `simulate_continuous` solves exactly, one
    draw per branch. Running it here would approximate a model that has an answer.

    The run is **refused** rather than sent there. Five probes of a Python callable are a guess
    about the whole line, and changing engines on that guess would change the draw order and what
    ``step`` means: the same call would give different numbers for a property nobody asked about.

    Writing ``θ_x(y) = a_x + b_x·y`` into ``dx = −α_x(x − θ_x(y))dt`` gives the drift matrix
    ``P_xx = α_x`` and ``P_xy = −α_x·b_x``, and the pair's own optima are where the two lines
    cross: ``θ*_x = (a_x + b_x·a_y)/(1 − b_x·b_y)``. Two lines that never cross leave the pair with
    no stationary mean, and then only the matrix can be written out."""
    conts = [p for p in plans if isinstance(p, _Cont)]
    if len(conts) != len(plans) or len(conts) != 2:
        return                              # a discrete trait is here: no exact engine covers it
    lines = []
    for p in conts:
        if p.theta_of is not None:
            other = next(q for q in conts if q is not p)
            fit = _affine_fit(p.theta_of, *_probe_span(other, tallest))
            if fit is None:
                return
            lines.append(fit)
        elif p.theta is not None:
            lines.append((p.theta, 0.0))    # a fixed optimum is a flat line in the other trait
        else:
            return                          # Brownian motion: no optimum, so no exact engine either
    x, y = conts
    (ax, bx), (ay, by) = lines
    det = 1.0 - bx * by
    if abs(det) < 1e-12:
        optima = ("# the two lines never cross, so the pair has no stationary optimum:\n"
                  "                        # only the matrix is determined, and reverts_to is yours "
                  "to pick.")
    else:
        optima = (f"reverts_to={{{x.name!r}: {_g((ax + bx * ay) / det)}, "
                  f"{y.name!r}: {_g((ay + by * ax) / det)}}},")
    raise ValueError(
        f"each optimum here is a straight line in the other trait, and that model is solved "
        f"exactly: the two traits are one multivariate Ornstein–Uhlenbeck with a full drift "
        f"matrix, one draw per branch, no step and no error. Write it as one call:\n"
        f"    simulate_continuous(tree, start={{{x.name!r}: {_g(x.start)}, "
        f"{y.name!r}: {_g(y.start)}}},\n"
        f"                        rate={{{x.name!r}: {_g(_sigma2_at(x))}, "
        f"{y.name!r}: {_g(_sigma2_at(y))}}},\n"
        f"                        {optima}\n"
        f"                        pull={{({x.name!r}, {x.name!r}): {_g(x.alpha)}, "
        f"({x.name!r}, {y.name!r}): {_g(-x.alpha * bx)},\n"
        f"                              ({y.name!r}, {y.name!r}): {_g(y.alpha)}, "
        f"({y.name!r}, {x.name!r}): {_g(-y.alpha * by)}}})\n"
        f"This run is refused rather than sent there: five probes of a callable are a guess about "
        f"the whole line, and changing engines on that guess would change the draw order and what "
        f"step means.")


# --- the walk -------------------------------------------------------------------------------------

def _boundaries(t0: float, t1: float, step: float | None) -> list[float]:
    """The slice boundaries across one branch: whole steps of ``step``, the last slice short.

    ``step`` is ``None`` when nothing in the run is frozen, and the branch is then one slice: with
    no driver held still there is nothing for a boundary to release, and the walk is exact.

    The count is computed rather than accumulated, so a branch that is an exact multiple of ``step``
    cannot gain a sliver of a slice from round-off, and how many draws a branch takes is a fact
    about the numbers written rather than about the order they were added in."""
    span = t1 - t0
    if span <= 0.0:
        return []
    if step is None:
        return [t0, t1]
    n = max(1, int(math.ceil(span / step - 1e-9)))
    return [t0 + j * step for j in range(n)] + [t1]


def _live_reads(spec, by_name: dict) -> list:
    """Every `~zombi2.params.connection.Driven` on ``spec`` that reads a trait growing beside it.

    A read of a name that is not in the run is left alone here. `simulate_traits` checks the names,
    and one gate saying it is better than two saying it differently."""
    if isinstance(spec, ContinuousTrait):
        mods = list(as_rate(spec.rate, default_scope=PerLineage).modifiers)
        if isinstance(spec.reverts_to, Driven):
            mods.append(spec.reverts_to)
    else:
        mods = [m for e in _driven_entries(list(spec.states), spec.switch) for m in e[2].modifiers]
    return [m for m in mods
            if isinstance(m, Driven) and isinstance(m.driver, str) and m.driver in by_name]


def _theta(plan: _Cont, frozen: dict, pieces: dict, at: float):
    """The optimum a continuous trait heads for over one segment, or ``None`` for Brownian motion.

    A **curved** optimum reads its driver frozen at the top of the slice, which is the run's one
    approximation. A **regime** optimum reads the discrete trait's real state on this segment, which
    is exact: the segments are cut where that trait switched, so the optimum is constant across
    one."""
    if plan.regime_of is not None:
        spans = pieces[plan.regime_of]
        for state, lo, hi in spans:
            if lo <= at < hi:
                return plan.regime_theta[state]
        return plan.regime_theta[spans[-1][0]]
    if plan.theta_of is not None:
        return plan.theta_of(frozen[plan.reads])
    return plan.theta


def simulate(tree, specs, *, seed=None, progress=False, check_affine: bool = True):
    """Walk two traits that feed back into each other, slice by slice, and return
    ``{name: TraitsResult}``.

    `zombi2.traits.several.simulate_traits` is the front door; this is the engine behind it, for the
    case where one of the two traits is continuous. The module docstring holds the method, the draw
    order, and what ``step`` costs.

    ``check_affine`` is the gate that refuses a pair of straight optima and names the exact call
    instead (`_refuse_affine`). It is on for every run. The one test that compares this engine
    against the exact one turns it off, because a straight optimum is the only case where the two
    engines describe the same law, and so the only case where they can be compared at all.

    Deterministic given ``seed``. Halving ``step`` changes the numbers for a given seed, by design:
    the module docstring says how to check convergence.
    """
    tallest = max(n.end_time for n in tree.nodes.values())
    by_name = {f"traits:{s.name}": j for j, s in enumerate(specs)}
    plans: list = [_settle_continuous(s, j, by_name) if isinstance(s, ContinuousTrait)
                   else _settle_discrete(s, j) for j, s in enumerate(specs)]
    _wire_regimes(plans, specs)

    # `step` bounds how stale a **frozen** driver gets, and the frozen drivers are exactly the
    # connections written with a modifier. A run whose only reading is regimes= freezes nothing —
    # the optimum switches where the discrete trait switches, which the segments give exactly — so
    # it takes no step and walks each branch in one slice.
    live = [m for s in specs for m in _live_reads(s, by_name)]
    step = None if not live else check_step(
        step_of(live, what="a trait here",
                how='set_by("traits:<name>", f, step=0.01) on the optimum, or '
                    'scaled_by("traits:<name>", f, step=0.01) on a switch rate'),
        tallest)
    if check_affine:
        _refuse_affine(plans, tallest)

    conts = [p for p in plans if isinstance(p, _Cont)]
    discs = [p for p in plans if isinstance(p, _Disc)]
    k = len(conts)
    # where a trait's value is threaded for a rate to read it. A connection that declares a step
    # keys on (name, step, path), so the key is taken from the connection rather than rebuilt here.
    reads_from: dict[object, int] = {m.key: by_name[cast(str, m.driver)] for m in live}

    rng, seed = stream("traits", seed)
    for plan in plans:                      # the root's starting values, in written order
        if isinstance(plan, _Disc):
            idx = {s: i for i, s in enumerate(plan.states)}
            if plan.spec_start is None:
                plan.start = int(rng.integers(len(plan.states)))
            elif plan.spec_start in idx:
                plan.start = idx[plan.spec_start]
            else:
                raise ValueError(f"start must be one of states={plan.states} (or None for a "
                                 f"uniform draw), got {plan.spec_start!r}")

    root = tree.nodes[tree.root]
    start_values = [p.start for p in plans]
    values: dict[int, list] = {}            # each node's end value, one per trait, in written order
    log: list[list[Change]] = [
        [Change(root.birth_time, "initial", tree.root, None,
                p.start if isinstance(p, _Cont) else p.states[p.start])] for p in plans]

    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        cur = list(start_values) if node.parent is None else list(values[node.parent])
        if node.parent is not None:
            # the jumps at the split, in written order: each trait hops on its own, exactly as it
            # would in a run of its own, and the pair lands wherever the hops leave it
            for plan in plans:
                if isinstance(plan, _Cont):
                    if plan.jump_sd > 0.0:
                        was = cur[plan.index]
                        cur[plan.index] = was + float(rng.normal(0.0, plan.jump_sd))
                        log[plan.index].append(
                            Change(node.birth_time, "on_speciation", i, was, cur[plan.index]))
                elif plan.shift > 0.0 and float(rng.random()) < plan.shift:
                    was = cur[plan.index]
                    j = int(rng.integers(len(plan.states) - 1))   # to a uniform *other* state
                    cur[plan.index] = j if j < was else j + 1
                    log[plan.index].append(Change(node.birth_time, "on_speciation", i,
                                                  plan.states[was], plan.states[cur[plan.index]]))

        branch: dict[int, list] = {p.index: [] for p in discs}
        bounds = _boundaries(node.birth_time, node.end_time, step)
        for lo, hi in zip(bounds, bounds[1:]):
            # 1. every live driver, read at the start of the slice and held across it
            frozen = {key: (cur[j] if isinstance(plans[j], _Cont) else plans[j].states[cur[j]])
                      for key, j in reads_from.items()}
            # 2. each discrete trait across the whole slice, under its frozen generator
            cuts = {lo, hi}
            pieces: dict[int, list] = {}
            for plan in discs:
                Q = _driven_q(plan.entries, len(plan.states), frozen, lo)
                end, segs = _gillespie(cur[plan.index], hi - lo, Q, rng)
                at, spans = lo, []
                for state, dur in segs:
                    spans.append((state, at, min(at + dur, hi)))
                    at = min(at + dur, hi)
                    cuts.add(at)
                pieces[plan.index] = spans
                branch[plan.index].extend(segs)
                cur[plan.index] = end
            # 3. the continuous traits across the segments step 2 produced, one exact OU draw each
            edges = sorted(cuts)
            for a, b in zip(edges, edges[1:]):
                if b <= a:
                    continue
                z = rng.standard_normal(k)  # one vector for every continuous trait, always drawn
                moved = []
                for c, plan in enumerate(conts):
                    theta = _theta(plan, frozen, pieces, a)
                    x = cur[plan.index]
                    var = _accrued_variance(plan.rate, a, b, pull=plan.alpha)
                    if theta is None:       # Brownian motion: nothing to revert to
                        mean = x
                    else:
                        mean = theta + (x - theta) * math.exp(-plan.alpha * (b - a))
                    # the square root of the traits' covariance. Nothing here writes a correlation
                    # between them, so that matrix is diagonal and its root is the standard
                    # deviations — one line to widen the day a correlation= is written.
                    moved.append((mean + math.sqrt(var) * float(z[c])) if var > 0.0 else mean)
                for plan, v in zip(conts, moved):
                    cur[plan.index] = v

        for plan in discs:                  # the switches between the segments are the events
            t = node.birth_time
            merged = _merge_runs(branch[plan.index])
            for (s1, d1), (s2, _d) in zip(merged, merged[1:]):
                t += d1
                log[plan.index].append(Change(t, "on_branch", i, plan.states[s1], plan.states[s2]))
        values[i] = cur

    out = {}
    for plan in plans:
        log[plan.index].sort(key=lambda c: c.time)
        if isinstance(plan, _Cont):
            node_values: dict[int, object] = {i: float(v[plan.index]) for i, v in values.items()}
            out[plan.name] = TraitsResult(tree, node_values, log[plan.index], seed)
        else:
            node_values = {i: plan.states[v[plan.index]] for i, v in values.items()}
            out[plan.name] = TraitsResult(tree, node_values, log[plan.index], seed, kind="discrete")
    return out
