"""The parallel engine for `simulate_genomes_family()` — one gene family at a time.

The default engine is a single Gillespie over the whole species tree, because a transfer at time ``t``
couples two lineages alive at ``t``. But families never mix: no event ever spans two families (a
transfer moves a copy *within* one family; caps, replacement and the recipient rule are all
within-family or read only the tree), so the global loop is a **superposition** of independent
per-family processes. This engine runs them as such — reached only when ``parallel`` is truthy.

Two passes. **Pass 1** (serial, cheap) runs the per-lineage origination Poisson over the tree to
enumerate every family and the exact point it originates — origination reads only the tree (and, when
skyline, time), never genome content, so it splits off cleanly. **Pass 2** (parallel) evolves each
family's duplication / transfer / loss from its origination down the tree, one worker process per
family, each under its own spawned RNG stream — so the result is identical for any worker count. A
family "roams" across lineages: it is inherited by both daughters at a speciation and carried to a new
lineage by a transfer, but it is always self-contained, which is what makes the decomposition exact.

Pass 2's per-family evolution is `simulate_one_family()` — a standalone primitive that, given a
prepared `FamilyContext` (the run's tree, rates and contemporaneous-lineage schedule, built once
by `prepare_family_context()`), evolves one family from any origination point down the whole tree.
The in-memory and streamed engines are both thin maps over it; the workers only ship that context once
and read it from module state.

The realisation differs from the serial reference engine for a given seed (a different, equally valid
draw — the "A" decision). Everything the family resolution accepts runs here: duplication / transfer /
loss / origination, every recipient rule, a skyline ``changing_at``, the family cap, a per-family draw
heterogeneity, ``self_transfer``, ``replacement``, named families — and a
**conditioned** rate, which does not couple families either: a ``Driven`` driver was grown before
this run and is an input to it, so a lineage's factor is the same number whichever family is asking.
The workers thread the driver trajectories with the rest of the context.

**One model does couple families, and this engine refuses it by name.** A ``joint=True`` run whose
rate reads live gene content — ``scaled_by("genomes:<family>", …)``, `resolve_live_drivers` — makes
every family's rate depend on what the rest of the genome is doing right now, which is exactly the
assumption above. The refusal is in `zombi2.genomes.family`, it names ``parallel`` / ``stream_to``,
and it fires before any work starts; such a run is serial. That is the cost of the decomposition,
not a gap in it.
"""

from __future__ import annotations

import bisect
import heapq
import math
import os
import pathlib
import shutil
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from .._runtime.parallel import guard_pool_workers, pool_errors, resolve_workers
from .._runtime.progress import progress_bar
from ..params.evaluate import values_at_birth
from ..rng import seed_sequence
from ._live import enter, retire, weighted_index
from ._transfer import mean_root_to_tip, recipient_index
from .events import EVENTS_HEADER, GeneEdge, event_rows, gene_label, node_label
from .._runtime.outputs import fresh_dirs
from .gene_trees import gene_trees_from_edges, write_gene_trees


def _unsupported_reason(dup, tra, los, org, transfer_to) -> str | None:
    """A one-line reason the parallel engine cannot run this configuration (so it falls back to the
    serial loop, loudly), or ``None`` when it can.

    Nothing is left in this list. **Conditioning does not couple families**, which is what makes the
    decomposition survive it: a ``Driven`` rate reads a driver that was grown *before* this run and
    is an input to it, so a lineage's factor at time ``t`` is the same number whichever family is
    asking, and no family can influence another through it. The workers thread the driver
    trajectories, and a ``Clades`` recipient rule threads the same way — clade membership is a fact
    about the tree, painted once. The function stays because the next model that *does* couple
    families (an epistatic rate, a genome-wide budget) belongs here rather than in a silent wrong
    answer."""
    return None


# --- Pass 1: enumerate every family and where it originates (serial, cheap) ------------------------

def _enumerate_families(tree, org, initial_families, families_named, placed, rng, trajs=None,
                        driven=False, lineage_factor=None):
    """``[(family_id, birth_lineage, birth_time), …]`` for every family, and ``{name: family_id}``.

    Initial and named families originate at the origin, then the ones ``origins=`` ``placed`` at a
    chosen point; the rest are drawn by the per-lineage origination Poisson walked over the tree
    schedule — a mini-Gillespie with only origination live, which is exact because origination reads
    only the number of living lineages and the time (and, when ``driven``, the driver on each of them
    — still no genome content, so the split is unaffected; the rate is then summed per lineage and
    the birth lineage drawn with the same weights). ``lineage_factor`` is origination's drawn factor
    per species branch, when the rate carries one, and it makes the walk per-lineage the same way a
    driver does: a branch that originates fast is likelier to be the one a family is born on.

    A placed family is one entry in this list like any other, which is the whole of what this engine
    needs to know about `resolve_origins`: Pass 2 already evolves a family from *any* origination
    point. The ids come before the drawn ones and in the same order the serial loop mints them, so
    the two engines agree on which family is the one you placed."""
    root = tree.nodes[tree.root]
    t0 = root.birth_time
    families: list[tuple[int, int, float]] = []
    fid = 0
    for _ in range(initial_families):
        families.append((fid, tree.root, t0)); fid += 1
    named: dict[str, int] = {}
    for name in families_named:
        named[name] = fid
        families.append((fid, tree.root, t0)); fid += 1
    for start_time, lineage in placed:
        families.append((fid, lineage, start_time)); fid += 1

    # the origination-only Gillespie: the live-lineage set follows the tree schedule exactly as the
    # global loop's does; when origination fires, a fresh family is born on a uniform living lineage.
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)
    alive: list[int] = []
    dummy: list = []
    pos: dict[int, int] = {}
    enter(alive, dummy, pos, tree.root, None)
    t = t0
    si = 0
    while si < len(schedule):
        k_alive = len(alive)
        weights = None
        if driven and lineage_factor is not None:
            weights = [org.effective(copies=0, lineages=1, time=t,
                                     drivers={key: trajs[key].value(alive[k], t) for key in trajs},
                                     carried_factor=lineage_factor[alive[k]])
                       for k in range(k_alive)]
        elif driven:
            weights = [org.effective(copies=0, lineages=1, time=t,
                                     drivers={key: trajs[key].value(alive[k], t) for key in trajs})
                       for k in range(k_alive)]
        elif lineage_factor is not None:
            weights = [org.effective(copies=0, lineages=1, time=t,
                                     carried_factor=lineage_factor[alive[k]])
                       for k in range(k_alive)]
        if weights is not None:
            rate = sum(weights)
        else:
            rate = org.effective(copies=0, lineages=k_alive, time=t)
        next_species = schedule[si][0]
        horizon = min(next_species, org.next_change(t))
        if driven:      # the driver switching on a branch changes the rate: stop and re-evaluate
            horizon = min(horizon, min((trajs[key].next_change(alive[k], t) for key in trajs
                                        for k in range(k_alive)), default=math.inf))
        if rate > 0.0:
            t_ev = t + float(rng.exponential(1.0 / rate))
            if t_ev < horizon:
                t = t_ev
                # origination is per lineage: uniform, or with the driver's weights when conditioned
                k = int(rng.integers(k_alive)) if weights is None else weighted_index(rng, weights, rate)
                families.append((fid, alive[k], t)); fid += 1
                continue
        if horizon == next_species:
            t = next_species
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                retire(alive, dummy, pos, pos[i])
                node = tree.nodes[i]
                if node.children:
                    for c in node.children:
                        enter(alive, dummy, pos, c, None)
                si += 1
        else:
            t = horizon
    return families, named


# --- Pass 2: evolve one family down the tree (a worker) -------------------------------------------

# Each family owns a copy-id range ``[fid << SHIFT, (fid+1) << SHIFT)``, so a worker mints
# **globally-unique** ids as ``base + local`` with no coordination and no post-hoc offset — the merge
# and the streamed shards simply concatenate. 2^30 leaves room for a billion copies in one family
# (never approached); ids are Python ints, so a large ``fid`` never overflows.
_COPY_ID_SHIFT = 30


def _copy_base(fid: int) -> int:
    return fid << _COPY_ID_SHIFT


@dataclass(frozen=True)
class FamilyContext:
    """Everything constant across a run's families — the read-only context `simulate_one_family()`
    evolves each family against. Built once by `prepare_family_context()` and, in the parallel
    engine, shipped to each worker a single time via the pool initializer (never re-pickled per family).
    Holds the tree, the D/T/L rates and their per-family a per-family draw slots, the transfer settings, the
    family cap, and the precomputed contemporaneous-lineage schedule (sorted birth / death times whose
    two pointers give the set alive at any instant, plus the times the "≥ 2 lineages" transfer gate
    flips). Origination is *not* here: a family's origination point is an input to
    `simulate_one_family()`, drawn once for the whole run in Pass 1."""

    tree: object
    dup: object
    tra: object
    los: object
    transfer_to: object
    replacement: bool
    self_transfer: bool
    cap: "int | None"
    depth: float
    fam_by: dict
    birth_times: list
    birth_nodes: list
    death_times: list
    death_nodes: list
    cross2: list
    #: driver key -> DriverTrajectory, for the rates a Driven conditions. Resolved once in
    #: `simulate_genomes_family()`, before the engines split, and shipped to each worker with the rest
    #: of the context — a trajectory is per-lineage segment lists, so it pickles like any other data.
    trajs: dict
    #: the transfer_to driver, read only at the instant a transfer fires (it changes no rate, so it is
    #: deliberately not in ``trajs`` and adds no horizon breakpoint)
    to_traj: object
    #: lineage -> its named clade, for a Clades recipient rule; ``None`` otherwise
    group_of: object
    #: which of duplication / transfer / loss carry a Driven — the per-lineage path is taken only
    #: for those, so an unconditioned run keeps the pooled arithmetic exactly as it was
    driven: dict
    #: ``{target: {node: factor}}`` — each rate's drawn multiplier per species branch, drawn once in
    #: the parent (which is what keeps the run worker-count invariant) and shipped to every worker
    #: with the rest of the context. Plain floats, so it pickles like the trajectories. ``{}`` unless
    #: some rate varies among lineages, and then a rate carrying none reads 1.0 everywhere.
    lin_mult: dict
    #: which of duplication / transfer / loss vary among lineages — the same role ``driven`` plays,
    #: for the other reason a rate is read lineage by lineage
    varying: dict


def prepare_family_context(tree, *, dup, tra, los, transfer_to, replacement, self_transfer,
                           cap, trajs=None, to_traj=None, group_of=None,
                           driven=None, lin_mult=None, varying=None) -> FamilyContext:
    """Precompute the per-run `FamilyContext` — the schedule and rate metadata every family
    reuses — so `simulate_one_family()` can evolve any family from any origination point without
    recomputing it. ``dup`` / ``tra`` / ``los`` are resolved `Rate`s (per
    copy)."""
    depth = mean_root_to_tip(tree)
    # the per-family modifiers each rate carries (origination is excluded upstream) — drawn once per
    # family and multiplied onto that rate, exactly the serial engine's fam_mult placement.
    fam_by = {"duplication": tuple(m for m, _ in dup.carried_modifiers(unit="families")),
              "transfer": tuple(m for m, _ in tra.carried_modifiers(unit="families")),
              "loss": tuple(m for m, _ in los.carried_modifiers(unit="families"))}
    # the contemporaneous-lineage machinery: sorted birth / death times (two pointers give the set alive
    # at any t) and the times can_xfer (≥ 2 lineages alive) flips.
    births = sorted((tree.nodes[i].birth_time, i) for i in tree.nodes)
    deaths = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)
    return FamilyContext(
        tree=tree, dup=dup, tra=tra, los=los, transfer_to=transfer_to, replacement=replacement,
        self_transfer=self_transfer, cap=cap, depth=depth, fam_by=fam_by,
        birth_times=[t for t, _ in births], birth_nodes=[i for _, i in births],
        death_times=[t for t, _ in deaths], death_nodes=[i for _, i in deaths],
        cross2=_cross2_times(tree),
        trajs=trajs or {}, to_traj=to_traj, group_of=group_of,
        driven=driven or {"duplication": False, "transfer": False, "loss": False},
        lin_mult=lin_mult or {},
        varying=varying or {"duplication": False, "transfer": False, "loss": False})


# The shared context and stream config, shipped once per worker by the initializer (never re-pickled per
# family). In-process (inline) runs set them directly; ``_STREAM`` is ``None`` off the streaming path.
_CTX: "FamilyContext | None" = None
_STREAM: "dict | None" = None


def _init_worker(ctx, stream=None) -> None:
    global _CTX, _STREAM
    _CTX = ctx
    _STREAM = stream


def _family_mults(rng, fam_by):
    """The family's rate multipliers, drawn once from its own stream (so they are worker-invariant).
    The draw order is fixed — duplication, then transfer, then loss, and within a rate the order its
    modifiers were written — so it is reproducible. ``1.0`` where a rate carries none. One per-family draw
    object read by two rates is one draw shared between them, which is how a family-wide tempo is
    said."""
    shared: dict[int, float] = {}    # across this family's rates: one object, one number
    out = {}
    for key in ("duplication", "transfer", "loss"):
        out[key] = math.prod(values_at_birth(fam_by.get(key, ()), rng, shared))
    return out


def simulate_one_family(ctx, *, family, lineage, time, rng, copy_id_base=0):
    """Evolve **one** gene family from a given origination point down the whole tree, and return its
    ``(events, node_genomes, multipliers)`` — the compact event log, this family's copies at every node
    it reaches, and the rate multipliers it drew (``{target: multiplier}``, 1.0 where a rate draws none).

    This is the per-family primitive the parallel engine is built on. The family resolution's D/T/L process is a
    superposition of independent per-family processes (no event ever spans two families), so a family
    born as ``family`` in lineage ``lineage`` at time ``time`` can be evolved on its own, against the
    contemporaneous-lineage schedule prepared once in ``ctx`` (a `FamilyContext`; see
    `prepare_family_context()`). A transfer may still hand a copy to *any* lineage alive at the
    instant — the family "roams" — which is why the whole tree's schedule is needed, not just
    ``lineage``'s subtree.

    ``rng`` is this family's own generator (spawn one stream per family for a worker-count-invariant
    run); ``copy_id_base`` offsets the minted copy ids so several families' logs concatenate without a
    rewrite (the parallel engine passes ``family << 30``). Mirrors the serial loop's inner event
    handling, scoped to this family's footprint (the lineages it occupies)."""
    from .family import GeneCopy, _at_cap, _duplicate, _lose_at   # package helpers; no cycle

    tree, dup, tra, los = ctx.tree, ctx.dup, ctx.tra, ctx.los
    transfer_to, replacement, self_transfer = ctx.transfer_to, ctx.replacement, ctx.self_transfer
    cap, depth, fam_by = ctx.cap, ctx.depth, ctx.fam_by
    birth_times, birth_nodes = ctx.birth_times, ctx.birth_nodes
    death_times, death_nodes = ctx.death_times, ctx.death_nodes
    cross2 = ctx.cross2
    trajs, to_traj, group_of, driven = ctx.trajs, ctx.to_traj, ctx.group_of, ctx.driven
    any_driven = bool(trajs)
    # The other reason a rate is read lineage by lineage: a factor drawn for each species branch,
    # shared by every family passing through it. Drawn in the parent and shipped in `ctx`, so every
    # family reads the same branch factors whatever the worker count.
    lin_mult, varying = ctx.lin_mult, ctx.varying
    any_lineage = bool(lin_mult)

    mult = _family_mults(rng, fam_by)
    m_dup, m_tra, m_los = mult["duplication"], mult["transfer"], mult["loss"]

    def weigh(key, rate, family_mult, values, t):
        """One rate on every lineage the family occupies, times the family's own multiplier.

        Three shapes rather than one context built per lineage, because this runs once per rate per
        event: a driven rate reads the branch's driver, a rate varying among lineages reads the
        branch's drawn factor, and a rate doing both reads both. The driven-only shape is the
        expression it always was, so a conditioned run is unchanged to the last bit. ``alive`` and
        ``gen`` are the loop's own arrays, mutated in place, so they are read rather than passed."""
        table = lin_mult[key] if varying[key] else None
        if table is None:                          # driven only — why `weigh` was called at all
            return [rate.effective(copies=len(gen[k]), lineages=1, time=t, drivers=values[k])
                    * family_mult for k in range(len(alive))]
        if not driven[key]:
            return [rate.effective(copies=len(gen[k]), lineages=1, time=t,
                                   carried_factor=table[alive[k]])
                    * family_mult for k in range(len(alive))]
        return [rate.effective(copies=len(gen[k]), lineages=1, time=t, drivers=values[k],
                               carried_factor=table[alive[k]])
                * family_mult for k in range(len(alive))]

    events: list[GeneEdge] = []
    node_genomes: dict[int, list] = {}
    base = copy_id_base
    counter = 0

    def new_copy(fam: int) -> "GeneCopy":
        nonlocal counter
        if counter >= (1 << _COPY_ID_SHIFT):                # a billion copies in one family — unreachable
            raise OverflowError(f"family {family} exceeded {1 << _COPY_ID_SHIFT} copies; raise _COPY_ID_SHIFT")
        gc = GeneCopy(base + counter, fam)
        counter += 1
        return gc

    # the occupied set — the family's footprint — kept as the same (alive, gen, pos) parallel arrays the
    # global loop uses, but holding only lineages where this family is present. A min-heap of their end
    # times gives the next structural event; stale entries (a lineage the family left) are skipped lazily.
    alive: list[int] = []
    gen: list[list] = []
    pos: dict[int, int] = {}
    heap: list[tuple[float, int]] = []

    founding = new_copy(family)                             # the founding gene (local id 0)
    events.append(GeneEdge(time, "origination", lineage, family, founding.id))
    enter(alive, gen, pos, lineage, [founding])
    heapq.heappush(heap, (tree.nodes[lineage].end_time, lineage))

    # the contemporaneous lineage set (all lineages alive now), maintained by two pointers as t rises —
    # needed to choose a transfer recipient and to know whether a recipient can exist at all.
    bi = bisect.bisect_right(birth_times, time)
    di = bisect.bisect_right(death_times, time)
    contemp = set(birth_nodes[:bi]) - set(death_nodes[:di])

    def advance_contemp(t):
        nonlocal bi, di
        while bi < len(birth_times) and birth_times[bi] <= t:
            contemp.add(birth_nodes[bi]); bi += 1
        while di < len(death_times) and death_times[di] <= t:
            contemp.discard(death_nodes[di]); di += 1

    def next_valid_end():
        while heap and heap[0][1] not in pos:              # drop lineages the family already left
            heapq.heappop(heap)
        return heap[0][0] if heap else math.inf

    def next_cross2(t):
        j = bisect.bisect_right(cross2, t)
        return cross2[j] if j < len(cross2) else math.inf

    t = time
    total = 1
    while alive:
        advance_contemp(t)
        k_alive = len(contemp)
        can_xfer = total > 0 and (k_alive >= 2 or self_transfer)
        # A driven rate is per lineage: it reads the driver on the branch, so it cannot be pooled over
        # the family's copies. Summed over the lineages this family occupies — its footprint, not the
        # whole tree — and the same weights then draw which of them the event lands on, or the rate
        # would say one thing and the picking another. The undriven rates stay pooled, so a run with
        # no conditioning does exactly the arithmetic it did before.
        w_dup = w_los = w_tra = None
        if (any_driven or any_lineage) and alive:
            values = ([{key: trajs[key].value(alive[k], t) for key in trajs}
                       for k in range(len(alive))] if any_driven else None)
            if driven["duplication"] or varying["duplication"]:
                w_dup = weigh("duplication", dup, m_dup, values, t)
            if driven["loss"] or varying["loss"]:
                w_los = weigh("loss", los, m_los, values, t)
            if (driven["transfer"] or varying["transfer"]) and can_xfer:
                w_tra = weigh("transfer", tra, m_tra, values, t)
        r_dup = (sum(w_dup) if w_dup is not None else
                 dup.effective(copies=total, lineages=1, time=t) * m_dup if total else 0.0)
        r_los = (sum(w_los) if w_los is not None else
                 los.effective(copies=total, lineages=1, time=t) * m_los if total else 0.0)
        r_tra = (sum(w_tra) if w_tra is not None else
                 tra.effective(copies=total, lineages=1, time=t) * m_tra if can_xfer else 0.0)
        rate_total = r_dup + r_los + r_tra

        next_struct = next_valid_end()
        horizon = min(next_struct, dup.next_change(t), los.next_change(t), tra.next_change(t),
                      next_cross2(t))
        if any_driven:   # a driven rate also changes when the driver switches on an occupied branch
            horizon = min(horizon, min((trajs[key].next_change(alive[k], t) for key in trajs
                                        for k in range(len(alive))), default=math.inf))

        if rate_total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / rate_total))
            if t_ev < horizon:
                t = t_ev
                r = float(rng.random()) * rate_total
                if r < r_dup:
                    k, j = _pick_in(rng, gen, total, w_dup, r_dup)
                    if not _at_cap(gen[k], gen[k][j].family, cap):
                        _duplicate(gen[k], j, tree.nodes[alive[k]], t, events, new_copy)
                        total += 1
                elif r < r_dup + r_los:
                    k, j = _pick_in(rng, gen, total, w_los, r_los)
                    _lose_at(gen[k], j, tree.nodes[alive[k]], t, events)
                    total -= 1
                    if not gen[k]:                          # the family left this lineage
                        retire(alive, gen, pos, k)
                else:
                    # the event fires at t_ev, which may be past a non-occupied lineage's birth/death
                    # (those are not in this family's horizon), so refresh the contemporaries to the
                    # event time before choosing a recipient — otherwise a dead or unborn lineage could
                    # be picked, and entering one with an end time behind t would rewind the clock.
                    advance_contemp(t)
                    kd = None if w_tra is None else weighted_index(rng, w_tra, r_tra)
                    total += _family_transfer(rng, tree, contemp, alive, gen, pos, heap, total, t,
                                              events, new_copy, transfer_to, replacement,
                                              self_transfer, depth, cap, to_traj, group_of, kd)
                continue

        # advance to the horizon: a structural event (an occupied lineage ends) or a rate breakpoint
        if next_struct < math.inf and horizon == next_struct:
            t = next_struct
            while heap and heap[0][0] == t and heap[0][1] in pos:
                _, i = heapq.heappop(heap)
                if i not in pos:                            # stale duplicate
                    continue
                g = gen[pos[i]]
                node_genomes[i] = list(g)                   # finalise this family's copies at node i
                retire(alive, gen, pos, pos[i])
                node = tree.nodes[i]
                if node.children and g:         # speciation: re-id each copy into daughters
                    total -= len(g)
                    per_daughter = []
                    for ch in node.children:
                        child = [new_copy(old.family) for old in g]
                        per_daughter.append([GeneEdge(t, "speciation", ch, old.family, nc.id,
                                                   parent=old.id) for old, nc in zip(g, child)])
                        enter(alive, gen, pos, ch, child)
                        heapq.heappush(heap, (tree.nodes[ch].end_time, ch))
                        total += len(child)
                    for rows in zip(*per_daughter):     # one gene, one event: its two edges together
                        events.extend(rows)
                else:                                       # a tip / extinction / empty: copies end here
                    total -= len(g)
        else:
            t = horizon                                     # a skyline or transfer-window breakpoint

    return events, node_genomes, mult


def _evolve_one(family, lineage, birth_time, seedseq):
    """Worker adapter: evolve one family from the shared per-run context under its own spawned RNG
    stream (so the run is identical for any worker count). The pool ships the context once via the
    initializer; this reads it from module state and hands it to `simulate_one_family()`."""
    return simulate_one_family(_CTX, family=family, lineage=lineage, time=birth_time,
                               rng=np.random.default_rng(seedseq), copy_id_base=_copy_base(family))


def _evolve_family(task):
    """Collect-mode worker: evolve one family and hand its log back for the in-memory merge."""
    fid, birth_lineage, birth_time, seedseq = task
    events, node_genomes, multipliers = _evolve_one(fid, birth_lineage, birth_time, seedseq)
    return fid, events, node_genomes, multipliers


def _pick_in(rng, gen, total, weights, weight_total):
    """``(lineage index, copy index)`` for the lineage an event lands on and the copy within it.

    ``weights is None`` — the undriven case — is the uniform pick over the family's whole copy pool,
    unchanged. Driven, the lineage is drawn with the same per-lineage weights the rate was summed
    with, and then a copy uniformly inside it, which is what makes the pick agree with the rate."""
    from .family import _pick_copy

    if weights is None:
        return _pick_copy(rng, gen, total)
    k = weighted_index(rng, weights, weight_total)
    if not gen[k]:                       # only via weighted_index's r == total float guard
        k = max(range(len(gen)), key=lambda i: weights[i])
    return k, int(rng.integers(len(gen[k])))


def _family_transfer(rng, tree, contemp, alive, gen, pos, heap, total, t, events, new_copy,
                     transfer_to, replacement, self_transfer, depth, cap,
                     to_traj=None, group_of=None, donor_k=None) -> int:
    """One transfer for the current family. Mirrors `_do_transfer()` exactly (donor
    continuation re-ids, optional homologous replacement, the cap thinning) but over this family's
    footprint: the recipient is a contemporaneous lineage picked by ``transfer_to``, and a recipient
    the family had not reached is entered into the footprint. Returns the change in copy count
    (+1 additive, 0 replacement/no-op).

    The donor copy is a uniform pick across the family's copies, except under a **driven** transfer
    rate, where the caller has already drawn the donor lineage with the rate's own weights and passes
    it as ``donor_k`` — a driven transfer weights who donates. ``to_traj`` / ``group_of`` are the
    recipient-side weightings, read here at the instant the transfer fires."""
    from .family import _at_cap, _pick_copy

    if donor_k is None:
        kd, jd = _pick_copy(rng, gen, total)               # a uniform donor copy across the family
    else:
        kd, jd = donor_k, int(rng.integers(len(gen[donor_k])))
    donor = alive[kd]
    src = gen[kd][jd]
    fam = src.family

    contemp_list = sorted(contemp)                         # deterministic order for the recipient pick
    cand = [k for k in range(len(contemp_list)) if self_transfer or contemp_list[k] != donor]
    if not cand:
        return 0
    kr = recipient_index(rng, tree, contemp_list, cand, donor, t, transfer_to, depth,
                         to_traj, group_of)
    if kr is None:                                         # nobody can receive (weighting thinning)
        return 0
    recipient = contemp_list[kr]

    rg = gen[pos[recipient]] if recipient in pos else None
    if rg is not None and _at_cap(rg, fam, cap):           # recipient full of this family: no-op
        return 0

    cont, xfer = new_copy(fam), new_copy(fam)
    gen[kd][jd] = cont                                     # the donor gene ends; continuation re-ids
    delta = 1
    if rg is None:                                         # the family reaches a new lineage
        enter(alive, gen, pos, recipient, [])
        heapq.heappush(heap, (tree.nodes[recipient].end_time, recipient))
        rg = gen[pos[recipient]]
    replaced = None
    if replacement:
        residents = [p for p, cpy in enumerate(rg) if cpy.family == fam and cpy.id != cont.id]
        if residents:                                      # homologous overwrite (else additive)
            p = residents[int(rng.integers(len(residents)))]
            victim = rg[p]
            rg[p] = rg[-1]; rg.pop()
            replaced = victim.id      # named on the transfer, so the log writes them as one event
            events.append(GeneEdge(t, "loss", recipient, fam, replaced))
            delta = 0
    rg.append(xfer)
    events.append(GeneEdge(t, "transfer", donor, fam, cont.id, parent=src.id, donor=donor,
                        replaced=replaced))
    events.append(GeneEdge(t, "transfer", recipient, fam, xfer.id, parent=src.id, recipient=recipient,
                        donor=donor, replaced=replaced))
    return delta


# --- output vocabulary for a streamed run (the files match FamilyGenomesResult.write) -------------------

#: the outputs a streamed run can produce and their top-level filenames — the same names the in-memory
#: ``FamilyGenomesResult.write`` uses. Gene trees are the exception: one Newick pair per family under a
#: ``gene_trees/`` subdirectory, so a million families do not land as two million files in the run root.
_STREAM_OUTPUTS = ("events", "profiles", "genomes", "initial_genome", "gene_trees", "species_tree",
                   "links", "family_multipliers", "lineage_multipliers")
_STREAM_FILENAMES = {"events": "genome_events.tsv", "profiles": "profiles.tsv",
                     "genomes": "genomes.tsv", "initial_genome": "initial_genome.tsv",
                     "species_tree": "species_complete.nwk", "links": "links.tsv",
                     "family_multipliers": "family_multipliers.tsv",
                     "lineage_multipliers": "lineage_multipliers.tsv",
                     # the files only an ordered run writes (`zombi2.genomes.ordered`)
                     "gene_order": "gene_order.tsv", "chromosome_events": "chromosome_events.tsv",
                     "summary": "genome_summary.json"}
_DEFAULT_STREAM_OUTPUTS = _STREAM_OUTPUTS

#: families per streamed chunk — **fixed**, independent of the worker count, so a chunk is a contiguous
#: family-id range and the shards concatenate in a deterministic order: a streamed run's files are
#: byte-identical for any number of workers. A million families is ~4000 chunks, plenty to fill the cores.
_STREAM_CHUNK = 256


@dataclass(frozen=True)
class StreamedRun:
    """A genome run written **straight to disk** — what ``stream_to=`` returns, for a scale where a
    whole result would not fit in memory. At the family resolution the files are written family by
    family; at the ordered resolution they are written as the run goes. Thin by design: the outputs
    *are* the files and the disk is the handoff (the sequences level reads them back), so this carries
    where they are and how big the run was, not the run itself. ``n_events`` counts gene-tree edges,
    and ``n_families`` the families the run began."""

    directory: str
    seed: "int | None"
    n_families: int
    n_events: int
    outputs: tuple

    def path(self, output: str) -> str:
        """The path of a written top-level file — e.g. ``path("events")`` → ``…/genome_events.tsv``.
        Gene trees are not a single file; they live one pair per family under ``gene_trees/``."""
        if output not in _STREAM_FILENAMES:
            raise KeyError(f"{output!r} is not a top-level streamed file (gene trees are under "
                           f"gene_trees/); files are {sorted(_STREAM_FILENAMES)}")
        return os.path.join(self.directory, _STREAM_FILENAMES[output])


def _stream_chunk(task):
    """Streaming worker: evolve a contiguous chunk of families and write as it goes — a per-chunk shard
    for each row output (events / genomes / profiles) and one Newick pair per family for the gene trees.
    Nothing run-sized is held; the parent concatenates the shards afterwards. Returns
    ``(chunk_index, n_families, n_events)``."""
    from .multipliers import FAMILY_TARGETS, multipliers_row

    chunk_index, family_list = task
    tree, s = _CTX.tree, _STREAM
    out_dir, outputs, extant_ids, shard_dir = s["out_dir"], s["outputs"], s["extant_ids"], s["shard_dir"]
    want = {name: name in outputs for name in ("events", "genomes", "profiles", "gene_trees")}
    # a run whose rates do not vary among families has no multiplier to write: its table is the header
    want["family_multipliers"] = "family_multipliers" in outputs and any(_CTX.fam_by.values())
    trees_dir = os.path.join(out_dir, "gene_trees")

    files = {name: open(os.path.join(shard_dir, f"{name}_{chunk_index}.tsv"), "w", encoding="utf-8")
             for name in ("events", "genomes", "profiles", "family_multipliers") if want[name]}
    names = tree.labels()   # e<id> for a lineage that died; once per chunk, not once per family
    n_events = 0
    try:
        for (fid, lineage, birth_time, seedseq) in family_list:
            events, node_genomes, multipliers = _evolve_one(fid, lineage, birth_time, seedseq)
            n_events += len(events)
            if want["family_multipliers"]:
                files["family_multipliers"].write(multipliers_row(fid, multipliers, FAMILY_TARGETS) + "\n")
            if want["events"]:
                f = files["events"]
                for row in event_rows(events, names):
                    f.write(row); f.write("\n")
            if want["genomes"]:
                f = files["genomes"]
                for node_id, copies in node_genomes.items():
                    label = names[node_id]
                    for cp in copies:
                        f.write(f"{label}\t{cp.family}\t{gene_label(cp.id)}\n")
            if want["profiles"]:
                counts = [len(node_genomes.get(sp, ())) for sp in extant_ids]
                if any(counts):                             # a family absent from every extant tip: no row
                    files["profiles"].write(f"{fid}\t" + "\t".join(map(str, counts)) + "\n")
            if want["gene_trees"]:
                write_gene_trees(gene_trees_from_edges(events, tree), trees_dir, names)
    finally:
        for f in files.values():
            f.close()
    return chunk_index, len(family_list), n_events


# --- the public entry: two passes, then either an in-memory merge or a stream to disk -------------

def run_parallel_family(tree, *, dup, tra, los, org, transfer_to, replacement, self_transfer,
                        initial_families, family_names, modules, cap, seed, parallel,
                        progress, placed=(), stream_to=None, outputs=None,
                        trajs=None, to_traj=None, group_of=None, driven=None, lin_by=None):
    """Run the per-family engine. Returns a `FamilyGenomesResult` (the in-memory
    merge), or a `StreamedRun` when ``stream_to`` is a directory — each family written straight
    to disk, for a scale a whole result would not hold. It can still return ``None`` (a loud fallback to the serial
    loop) for the in-memory path, but nothing currently asks for one — conditioning does not couple
    families, so a driven rate runs here as well; a streamed run **raises** rather than silently
    pulling the whole thing back into memory.

    Copy ids are global from the start (``fid << SHIFT`` + local), so the merge and the streamed shards
    both just concatenate — no id rewrite, no run-sized bottleneck beyond the (serial) in-memory merge
    the streaming path exists to avoid."""
    from .family import GeneCopy, FamilyGenomesResult
    from .multipliers import FAMILY_TARGETS, draw_lineage_multipliers, lineage_multipliers_of

    reason = _unsupported_reason(dup, tra, los, org, transfer_to)
    if reason is not None:
        if stream_to is not None:
            raise ValueError(
                f"a streamed run cannot handle this: {reason}. Streaming needs the per-family engine "
                "(families written independently); run without the driver, or drop stream_to to fall "
                "back to the serial in-memory engine.")
        print(f"note: --parallel not applied — {reason}; running serially instead")
        return None
    if stream_to is not None:
        outputs = tuple(outputs) if outputs is not None else _DEFAULT_STREAM_OUTPUTS
        unknown = [o for o in outputs if o not in _STREAM_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown stream outputs {unknown}; choose from {list(_STREAM_OUTPUTS)}")

    workers = guard_pool_workers(resolve_workers(parallel))
    driven = driven or {}
    lin_by = lin_by or {}
    varying = {key: bool(lin_by.get(key)) for key in ("duplication", "transfer", "loss")}

    # Pass 0: the per-lineage multipliers, drawn in the parent from a stream of their own, before
    # the families are enumerated — origination may carry one, and the enumeration needs it. Drawing
    # here is what keeps the run worker-count invariant, exactly as the sequences level draws its
    # clock in the parent. A run whose rates do not vary among lineages spawns nothing, so its
    # family streams are the ones it always had.
    root_ss = seed_sequence("genomes", seed)[0]
    lin_mult: dict = {}
    if any(lin_by.values()):
        lin_mult = draw_lineage_multipliers(lin_by, tree, np.random.default_rng(root_ss.spawn(1)[0]))

    ctx = prepare_family_context(
        tree, dup=dup, tra=tra, los=los, transfer_to=transfer_to, replacement=replacement,
        self_transfer=self_transfer, cap=cap,
        trajs=trajs, to_traj=to_traj, group_of=group_of, driven=driven,
        lin_mult=lin_mult, varying=varying)

    # Pass 1: who originates, and where. One reserved stream for it; one per family after.
    families_meta, named = _enumerate_families(
        tree, org, initial_families, family_names, placed,
        np.random.default_rng(root_ss.spawn(1)[0]),
        trajs=trajs, driven=driven.get("origination", False),
        lineage_factor=lin_mult["origination"] if lin_by.get("origination") else None)
    n_families = len(families_meta)
    family_seeds = root_ss.spawn(n_families) if n_families else []
    per_family = [(fid, lin, bt, family_seeds[k]) for k, (fid, lin, bt) in enumerate(families_meta)]

    if stream_to is not None:
        return _run_streaming(tree, ctx, per_family, n_families, workers, seed, initial_families,
                              family_names, str(stream_to), outputs, progress)

    # In-memory: evolve each family, then merge. Inline for a small run (the pool's spawn + IPC would
    # cost more than it saves); one process per family otherwise. Same streams either way.
    results = []
    bar = progress_bar(max(1, n_families), "genomes", unit="family", enabled=progress)
    if workers > 1 and n_families >= 2:
        w = min(workers, n_families)
        with pool_errors(), ProcessPoolExecutor(max_workers=w, initializer=_init_worker,
                                                initargs=(ctx,)) as ex:
            for out in ex.map(_evolve_family, per_family, chunksize=max(1, n_families // (w * 8))):
                results.append(out); bar.update()
    else:
        _init_worker(ctx)
        for task in per_family:
            results.append(_evolve_family(task)); bar.update()
    bar.close()

    # The merge is now a concatenation: copy ids are already globally unique, so the per-family logs and
    # node snapshots stitch together with no rewrite. Every node appears, empty where no family reached it.
    events: list[GeneEdge] = []
    genomes: dict[int, list] = {i: [] for i in tree.nodes}
    for _fid, fam_events, node_genomes, _multipliers in results:
        events.extend(fam_events)
        for node_id, copies in node_genomes.items():
            genomes[node_id].extend(copies)
    events.sort(key=lambda e: e.time)                      # a chronological log, like the serial one
    genomes_final = {i: tuple(g) for i, g in genomes.items()}
    # each family's drawn multipliers, kept when some rate varies among families, as the serial engine does
    multipliers = ({fid: {target: float(m[target]) for target in FAMILY_TARGETS}
                    for fid, _e, _g, m in sorted(results, key=lambda r: r[0])}
                   if any(ctx.fam_by.values()) else {})

    # the genome the run started with: every initial and named family's founding gene (its base id),
    # before the stem — the snapshot the serial engine takes as `initial_genome`.
    n_seeded = initial_families + len(family_names)
    initial_genome = tuple(GeneCopy(_copy_base(fid), fid) for fid in range(n_seeded))
    return FamilyGenomesResult(tree, genomes_final, events, seed, named, dict(modules or {}),
                               initial_genome,
                               ctx.cap if hasattr(ctx, 'cap') else None,
                               family_multipliers=multipliers,
                               lineage_multipliers=(lineage_multipliers_of(lin_mult, tree.labels())
                                                    if lin_mult else {}))


def _run_streaming(tree, ctx, per_family, n_families, workers, seed, initial_families, family_names,
                   out_dir, outputs, progress):
    """The streaming half of `run_parallel_family()`: fixed contiguous chunks written to
    per-chunk shards, concatenated in chunk order (so the files are byte-identical for any worker
    count), then the shards removed. Returns a `StreamedRun`."""
    os.makedirs(out_dir, exist_ok=True)
    # `gene_trees/` is emptied first, as every other writer of a per-family directory does: a
    # streamed run writes each family as it goes rather than through `.write()`, so it skipped this
    # and a re-run with fewer families left the previous run's trees sitting beside the new ones.
    fresh_dirs(pathlib.Path(out_dir), ("gene_trees",), flat=False)
    shard_dir = os.path.join(out_dir, "_shards")
    os.makedirs(shard_dir, exist_ok=True)
    extant_ids = sorted(tree.extant_leaves())
    stream_cfg = {"out_dir": out_dir, "outputs": set(outputs), "extant_ids": extant_ids,
                  "shard_dir": shard_dir}
    # The tree the run evolved along, beside its outputs — every one of them is indexed by this
    # tree's node labels, so without it the directory is not a dataset anyone (or `read_run`) can
    # reopen. The in-memory `.write` learned this first; a streamed run needs it more, being the one
    # whose only handoff *is* the directory.
    if "species_tree" in outputs:
        with open(os.path.join(out_dir, "species_complete.nwk"), "w", encoding="utf-8") as f:
            f.write(tree.to_newick() + "\n")
    # This engine refuses every run that reads its own gene content, so a streamed run has no link.
    # Its table is the header alone, written so the directory holds the same files `.write` gives.
    if "links" in outputs:
        from .links import links_tsv
        with open(os.path.join(out_dir, "links.tsv"), "w", encoding="utf-8") as f:
            f.write(links_tsv(()))
    # The per-lineage multipliers belong to the run rather than to any family — they were drawn in
    # the parent, before the families were enumerated — so they are written here, whole, rather than
    # stitched from per-chunk shards the way the per-family table is.
    if "lineage_multipliers" in outputs:
        from .multipliers import lineage_multipliers_of, lineage_multipliers_tsv
        labels = tree.labels()
        table = lineage_multipliers_of(ctx.lin_mult, labels) if ctx.lin_mult else {}
        with open(os.path.join(out_dir, "lineage_multipliers.tsv"), "w", encoding="utf-8") as f:
            f.write(lineage_multipliers_tsv(table, labels))

    chunks = [per_family[i:i + _STREAM_CHUNK] for i in range(0, n_families, _STREAM_CHUNK)]
    tasks = list(enumerate(chunks))
    n_chunks = len(chunks)
    total_events = 0
    bar = progress_bar(max(1, n_chunks), "genomes", unit="chunk", enabled=progress)
    if workers > 1 and n_chunks >= 2:
        w = min(workers, n_chunks)
        with pool_errors("stream_to="), ProcessPoolExecutor(max_workers=w, initializer=_init_worker,
                                                            initargs=(ctx, stream_cfg)) as ex:
            for (_ci, _nfam, nev) in ex.map(_stream_chunk, tasks):
                total_events += nev; bar.update()
    else:
        _init_worker(ctx, stream_cfg)
        for task in tasks:
            _ci, _nfam, nev = _stream_chunk(task)
            total_events += nev; bar.update()
    bar.close()

    _finalize_stream(out_dir, shard_dir, outputs, extant_ids, n_chunks, initial_families, family_names)
    return StreamedRun(out_dir, seed, n_families, total_events, tuple(outputs))


def _finalize_stream(out_dir, shard_dir, outputs, extant_ids, n_chunks, initial_families, family_names):
    """Stitch the per-chunk shards into the run's files — the header once, then every shard in chunk
    order (pure I/O, never a run-sized allocation) — write ``initial_genome.tsv`` from the seeded
    families' base ids, and drop the shard directory."""
    from .multipliers import FAMILY_TARGETS, multipliers_header

    headers = {"events": EVENTS_HEADER,
               "genomes": "lineage\tfamily\tcopy",
               "family_multipliers": multipliers_header(FAMILY_TARGETS),
               # extant tips only, so every column is n<id>: a profile never names a dead lineage
               "profiles": "family\t" + "\t".join(node_label(s) for s in extant_ids)}
    for name, header in headers.items():
        if name in outputs:
            with open(os.path.join(out_dir, _STREAM_FILENAMES[name]), "w", encoding="utf-8") as out:
                out.write(header + "\n")
                for ci in range(n_chunks):
                    shard = os.path.join(shard_dir, f"{name}_{ci}.tsv")
                    if os.path.exists(shard):
                        with open(shard, encoding="utf-8") as sf:
                            shutil.copyfileobj(sf, out)
    if "initial_genome" in outputs:
        n_seeded = initial_families + len(family_names)
        with open(os.path.join(out_dir, _STREAM_FILENAMES["initial_genome"]), "w", encoding="utf-8") as out:
            out.write("family\tcopy\n")
            for fid in range(n_seeded):
                out.write(f"{fid}\t{gene_label(_copy_base(fid))}\n")
    shutil.rmtree(shard_dir, ignore_errors=True)


def _cross2_times(tree) -> list[float]:
    """The times at which the number of living lineages crosses 2 — where a transfer's "a recipient can
    exist" gate (``k_alive ≥ 2``) flips. Usually just the first speciation (1→2). Included in each
    family's horizon so a long waiting time cannot straddle the moment transfer becomes possible."""
    deltas = sorted([(tree.nodes[i].birth_time, 1) for i in tree.nodes]
                    + [(tree.nodes[i].end_time, -1) for i in tree.nodes])
    out: list[float] = []
    count = 0
    prev = False
    idx = 0
    while idx < len(deltas):
        t = deltas[idx][0]
        while idx < len(deltas) and deltas[idx][0] == t:   # apply every change at this instant together
            count += deltas[idx][1]; idx += 1
        can = count >= 2
        if can != prev:
            out.append(t); prev = can
    return out
