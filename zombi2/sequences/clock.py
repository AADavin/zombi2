"""The sequences level's per-species-branch rate variation — the **clock** (SPEC §5, §7).

SPEC §7 reserves the word *clock* for one thing: the sequences level's by-lineage substitution-rate
modifier. That is exactly and only what this module is about — how the substitution rate varies from
lineage to lineage, and nothing else. There are two ways it can vary, and they compose:

- **drawn** — a per-lineage draw gives every species branch one i.i.d. multiplier (the uncorrelated,
  "relaxed" clock), an inherited value lets the multiplier drift parent→child down the species tree (the
  autocorrelated clock). One draw per *species* branch, shared by every gene family passing through
  it: a species that runs hot runs hot for all of its genes.
- **read off a driver** — ``ScaledBy`` reads a **trait** grown first and maps its state to a factor,
  so a lineage's habitat or lifestyle sets how fast its sequences evolve. SPEC §3 allows the pair
  Traits→Sequences to be *conditioned* (the trait can be grown first and held fixed), which is what
  makes this an ordinary modifier on an ordinary run rather than a joint engine.

Both end up as a factor on the same axis — the species branch — so `Clock` carries them together and
answers the one question the engine asks: **how long is this stretch of branch, in substitutions per
site?** Every place the level converts time into substitutions goes through
`Clock.branch_length()`, so the alignment and the phylograms cannot disagree about the tree the
sequences were drawn along.

**A driver that switches mid-branch is integrated, not sampled.** A discrete trait switches partway
along a branch, so the driver factor is piecewise-constant in time rather than constant per branch,
and a stretch of branch from ``t0`` to ``t1`` is worth ``base × clock × ∫ g(s) ds`` over that
stretch, with ``g`` the driver factor. Reading ``g`` once per branch would credit the whole branch to
whichever state happened to be read — the same error `zombi2.traits.continuous` refuses for a driven
σ² and the genome engines refuse for their Gillespie horizons. The integral is exact here because
``g`` is piecewise-constant: this module precomputes each branch's breakpoints and the running
integral at each of them once, and `Clock.branch_length()` is then two bisected lookups.

Precomputing to plain floats is also what lets a driven run go **parallel**: the resolved clock is
shipped to every worker through the process-pool initializer, and a `~zombi2.rates.mapping.Curve`
mapping is usually a lambda, which does not pickle. The trajectories and the mappings are read at
build time and dropped; only numbers cross the process boundary.
"""

from __future__ import annotations

import math

import bisect

from ..rates.modifiers import INHERITED


def _all_species(gene_trees) -> list[int]:
    """The sorted set of species-branch ids the gene trees touch — the lineages the clock is drawn
    over. Collected from the gene trees rather than the species tree, so the draws depend only on the
    branches genes actually pass through; every branch that needs a clock value has its species
    branch present as some node's ``species`` (a branch no gene crossed keeps the factor 1.0)."""
    ids: set[int] = set()
    for gt in gene_trees.values():
        stack = [gt.complete]
        while stack:
            n = stack.pop()
            ids.add(n.species)
            stack.extend(n.children)
    return sorted(ids)


def _preorder(tree) -> list[int]:
    """Species-tree node ids, parent before child — the order the autocorrelated clock descends."""
    order: list[int] = []
    stack = [tree.root]
    while stack:
        i = stack.pop()
        order.append(i)
        kids = tree.nodes[i].children
        if kids is not None:
            stack.extend(kids)
    return order


def _draw_clock(clock_mods, species_tree, gene_trees, rng) -> "dict[int, float]":
    """The drawn part of the clock: one factor per species branch, drawn once and shared by every
    family (a hot species runs hot for all its genes). A `DRAWN`
    modifier draws each branch i.i.d.; an `INHERITED` one drifts the factor
    parent→child down the species tree (autocorrelated). ``{}`` for a strict clock — no draw is
    taken, and no randomness consumed.

    ``clock_mods`` is what the rate carries per lineage, from `Rate.carried_modifiers`, and
    **every** one of them is drawn and multiplied in. The two kinds cannot be mixed on one rate
    (`check_one_memory`), so the walk below is one shape or
    the other, never both.

    Factored out so the serial loop and the parallel engine draw it the same way — each just hands in
    its own ``rng`` (the serial run's shared generator, or a stream spawned for the clock so the
    parallel engine is worker-count invariant)."""
    if not clock_mods:
        return {}
    if clock_mods[0].reads[0] == INHERITED:
        walk: dict[int, tuple[float, ...]] = {}
        for i in _preorder(species_tree):                       # parent before child
            p = species_tree.nodes[i].parent
            walk[i] = (tuple(m.initial() for m in clock_mods) if p is None
                       else tuple(m.descend(v, rng) for m, v in zip(clock_mods, walk[p])))
        return {i: math.prod(v) for i, v in walk.items()}
    return {sid: math.prod(tuple(m.draw(rng) for m in clock_mods))
            for sid in _all_species(gene_trees)}


class Clock:
    """The resolved per-species-branch rate variation of one run — what `resolve_clock()` returns.

    It holds the drawn clock factor per species branch and, when the substitution rate carries a
    ``ScaledBy``, that branch's driver factor as a piecewise-constant function of time already
    integrated into a running total. Both are plain numbers: the object is read-only, is shipped
    to the parallel engine's workers, and holds no trajectory, mapping or tree.

    ``branch_length()`` is the whole interface. Everything the level computes in substitutions per
    site — the branch a sequence is sampled across, the gene phylogram, the species phylogram —
    asks it, so those three cannot drift apart."""

    __slots__ = ("_factor", "_breaks", "_level", "_cum")

    def __init__(self, factor: "dict[int, float]",
                 breaks: "dict[int, list[float]] | None" = None,
                 level: "dict[int, list[float]] | None" = None,
                 cum: "dict[int, list[float]] | None" = None) -> None:
        self._factor = factor       # {species branch: the drawn clock factor}; empty ⇒ strict
        # The driver track, or None throughout when the rate carries no Driven. Per species branch:
        # the times the driver factor changes (the first is the branch's birth), the factor on each
        # stretch, and the integral of the factor from the branch's birth to each breakpoint.
        self._breaks = breaks
        self._level = level
        self._cum = cum

    def branch_length(self, rate_base: float, species: int, start: float, end: float) -> float:
        """The stretch of species branch ``species`` from ``start`` to ``end``, in **substitutions per
        site** — ``base × clock × ∫ driver``.

        ``rate_base`` is an argument rather than a stored field because a nucleotide run scales it per
        block (spacer evolves at ``intergene_speed`` times the gene rate), so there is no single base
        for the run.

        With no driver this is ``rate_base × clock × Δt``, in that association order, which is the
        expression the level computed before drivers existed — so an undriven run at a given seed
        produces the same bytes it always did. A branch no gene passed through, and every branch under
        a strict clock, has factor 1.0."""
        f = self._factor.get(species, 1.0)
        if self._cum is None:
            return rate_base * f * (end - start)
        return rate_base * f * (self._integral(species, end) - self._integral(species, start))

    def _integral(self, species: int, t: float) -> float:
        """``∫`` of the driver factor along species branch ``species``, from its birth to ``t``.

        Piecewise-linear in ``t`` because the factor is piecewise-constant, so the value at ``t`` is
        the running total at the last breakpoint at or before it, plus that stretch's factor times how
        far into the stretch ``t`` sits. A time outside the branch extends the first or last stretch
        rather than clamping, which is what the driver itself does with its state (`DriverTrajectory`
        is right-continuous and holds the last state) — so a query a hair past a branch end from
        floating-point drift costs nothing rather than silently shortening the branch."""
        breaks = self._breaks[species]                          # type: ignore[index]
        i = bisect.bisect_right(breaks, t) - 1
        if i < 0:
            i = 0
        return self._cum[species][i] + self._level[species][i] * (t - breaks[i])  # type: ignore[index]


def resolve_clock(clock_mods, driven, species_tree, gene_trees, rng) -> "Clock | None":
    """Build the run's `Clock` from the substitution rate's modifiers, or ``None`` when the rate
    carries neither a lineage clock nor a driver.

    ``clock_mods`` is what the rate carries per lineage (possibly none), and ``driven`` the list
    of ``(ScaledBy modifier, DriverTrajectory)`` pairs the caller already resolved. Returning ``None``
    for a strict, undriven rate keeps ``clock is None`` meaning what it has always meant downstream:
    the fast path with no lookups at all.

    The draw comes **first** and is untouched, so a relaxed-clock run at a given seed draws exactly
    what it drew before drivers existed; building the driver track consumes no randomness.

    The track is built by walking each species branch from its birth, asking each driver where it next
    switches (`DriverTrajectory.next_change`) and taking the earliest, exactly as the genome and trait
    engines advance their own integrals. Several ``ScaledBy`` on one rate are allowed and their factors
    multiply, which is what SPEC §5 says modifiers do; a branch the drivers never switch on gets a
    single stretch, and the arithmetic degenerates to a constant factor times Δt."""
    factor = _draw_clock(clock_mods, species_tree, gene_trees, rng)
    if not driven:
        return Clock(factor) if factor else None
    breaks: dict[int, list[float]] = {}
    level: dict[int, list[float]] = {}
    cum: dict[int, list[float]] = {}
    for i, node in species_tree.nodes.items():
        end = node.end_time
        ts: list[float] = []
        gs: list[float] = []
        t = node.birth_time
        while True:
            g = 1.0
            for m, traj in driven:
                g *= m.factor(drivers={m.key: traj.value(i, t)})
            ts.append(t)
            gs.append(g)
            nxt = min(traj.next_change(i, t) for _, traj in driven)
            if not (nxt < end):     # also stops at inf: the driver is constant to the branch's end
                break
            t = nxt
        running = [0.0]
        for k in range(1, len(ts)):
            running.append(running[k - 1] + gs[k - 1] * (ts[k] - ts[k - 1]))
        breaks[i], level[i], cum[i] = ts, gs, running
    return Clock(factor, breaks, level, cum)


__all__ = ["Clock", "resolve_clock"]
