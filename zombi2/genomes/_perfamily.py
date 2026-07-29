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
loss / origination, every recipient rule, skyline ``OnTime``, the family cap, ``ByFamily`` /
``family_speed`` heterogeneity, ``self_transfer``, ``replacement``, named families — and a
**conditioned** rate, which does not couple families either: a ``DrivenBy`` driver was grown before
this run and is an input to it, so a lineage's factor is the same number whichever family is asking.
The workers thread the driver trajectories with the rest of the context. The engine still *has* a
loud fallback, for the next model that genuinely couples families.
"""

from __future__ import annotations

import bisect
import heapq
import math
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from .._runtime.parallel import guard_pool_workers, resolve_workers
from .._runtime.progress import progress_bar
from ..rates.modifiers import ByFamily
from ._live import enter, retire, weighted_index
from ._transfer import mean_root_to_tip, recipient_index
from .events import EVENTS_HEADER, Event, event_rows, gene_label, node_label
from .gene_trees import gene_trees_from_events, write_gene_trees


def _unsupported_reason(dup, tra, los, org, transfer_to) -> str | None:
    """A one-line reason the parallel engine cannot run this configuration (so it falls back to the
    serial loop, loudly), or ``None`` when it can.

    Nothing is left in this list. **Conditioning does not couple families**, which is what makes the
    decomposition survive it: a ``DrivenBy`` rate reads a driver that was grown *before* this run and
    is an input to it, so a lineage's factor at time ``t`` is the same number whichever family is
    asking, and no family can influence another through it. The workers thread the driver
    trajectories, and a ``Clades`` recipient rule threads the same way — clade membership is a fact
    about the tree, painted once. The function stays because the next model that *does* couple
    families (an epistatic rate, a genome-wide budget) belongs here rather than in a silent wrong
    answer."""
    return None


# --- Pass 1: enumerate every family and where it originates (serial, cheap) ------------------------

def _enumerate_families(tree, org, initial_families, families_named, rng, trajs=None, driven=False):
    """``[(family_id, birth_lineage, birth_time), …]`` for every family, and ``{name: family_id}``.

    Initial and named families originate at the origin; the rest are drawn by the per-lineage
    origination Poisson walked over the tree schedule — a mini-Gillespie with only origination live,
    which is exact because origination reads only the number of living lineages and the time (and,
    when ``driven``, the driver on each of them — still no genome content, so the split is unaffected;
    the rate is then summed per lineage and the birth lineage drawn with the same weights)."""
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
        if driven:
            weights = [org.effective(copies=0, lineages=1, time=t,
                                     drivers={key: trajs[key].value(alive[k], t) for key in trajs})
                       for k in range(k_alive)]
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
                if node.children is not None:
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
    Holds the tree, the D/T/L rates and their per-family ``ByFamily`` slots, the transfer settings, the
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
    family_speed: object
    fam_by: dict
    birth_times: list
    birth_nodes: list
    death_times: list
    death_nodes: list
    cross2: list
    #: driver key -> DriverTrajectory, for the rates a DrivenBy conditions. Resolved once in
    #: `simulate_genomes_family()`, before the engines split, and shipped to each worker with the rest
    #: of the context — a trajectory is per-lineage segment lists, so it pickles like any other data.
    trajs: dict
    #: the transfer_to driver, read only at the instant a transfer fires (it changes no rate, so it is
    #: deliberately not in ``trajs`` and adds no horizon breakpoint)
    to_traj: object
    #: lineage -> its named clade, for a Clades recipient rule; ``None`` otherwise
    group_of: object
    #: which of duplication / transfer / loss carry a DrivenBy — the per-lineage path is taken only
    #: for those, so an unconditioned run keeps the pooled arithmetic exactly as it was
    driven: dict


def prepare_family_context(tree, *, dup, tra, los, transfer_to, replacement, self_transfer,
                           cap, family_speed, trajs=None, to_traj=None, group_of=None,
                           driven=None) -> FamilyContext:
    """Precompute the per-run `FamilyContext` — the schedule and rate metadata every family
    reuses — so `simulate_one_family()` can evolve any family from any origination point without
    recomputing it. ``dup`` / ``tra`` / ``los`` are resolved `Rate`s (per
    copy)."""
    depth = mean_root_to_tip(tree)
    # which single rate each per-family ByFamily slot sits on (origination is excluded upstream) — drawn
    # once per family and multiplied onto that rate, exactly the serial engine's fam_mult placement.
    fam_by = {"duplication": next((m for m in dup.modifiers if isinstance(m, ByFamily)), None),
              "transfer": next((m for m in tra.modifiers if isinstance(m, ByFamily)), None),
              "loss": next((m for m in los.modifiers if isinstance(m, ByFamily)), None)}
    # the contemporaneous-lineage machinery: sorted birth / death times (two pointers give the set alive
    # at any t) and the times can_xfer (≥ 2 lineages alive) flips.
    births = sorted((tree.nodes[i].birth_time, i) for i in tree.nodes)
    deaths = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)
    return FamilyContext(
        tree=tree, dup=dup, tra=tra, los=los, transfer_to=transfer_to, replacement=replacement,
        self_transfer=self_transfer, cap=cap, depth=depth, family_speed=family_speed, fam_by=fam_by,
        birth_times=[t for t, _ in births], birth_nodes=[i for _, i in births],
        death_times=[t for t, _ in deaths], death_nodes=[i for _, i in deaths],
        cross2=_cross2_times(tree),
        trajs=trajs or {}, to_traj=to_traj, group_of=group_of,
        driven=driven or {"duplication": False, "transfer": False, "loss": False})


# The shared context and stream config, shipped once per worker by the initializer (never re-pickled per
# family). In-process (inline) runs set them directly; ``_STREAM`` is ``None`` off the streaming path.
_CTX: "FamilyContext | None" = None
_STREAM: "dict | None" = None


def _init_worker(ctx, stream=None) -> None:
    global _CTX, _STREAM
    _CTX = ctx
    _STREAM = stream


def _family_mults(rng, family_speed, fam_by):
    """The family's rate multipliers, drawn once from its own stream (so they are worker-invariant):
    ``family_speed`` scales every rate together (one draw), a per-rate ``ByFamily`` varies that rate on
    its own. The draw order is fixed — speed, then duplication / transfer / loss — so it is
    reproducible. ``1.0`` where a slot carries neither."""
    speed = family_speed.draw(rng) if family_speed is not None else 1.0
    out = {}
    for key in ("duplication", "transfer", "loss"):
        m = fam_by.get(key)
        out[key] = speed * (m.draw(rng) if m is not None else 1.0)
    return out


def simulate_one_family(ctx, *, family, lineage, time, rng, copy_id_base=0):
    """Evolve **one** gene family from a given origination point down the whole tree, and return its
    ``(events, node_genomes)`` — the compact event log and this family's copies at every node it reaches.

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
    cap, depth, family_speed, fam_by = ctx.cap, ctx.depth, ctx.family_speed, ctx.fam_by
    birth_times, birth_nodes = ctx.birth_times, ctx.birth_nodes
    death_times, death_nodes = ctx.death_times, ctx.death_nodes
    cross2 = ctx.cross2
    trajs, to_traj, group_of, driven = ctx.trajs, ctx.to_traj, ctx.group_of, ctx.driven
    any_driven = bool(trajs)

    mult = _family_mults(rng, family_speed, fam_by)
    m_dup, m_tra, m_los = mult["duplication"], mult["transfer"], mult["loss"]

    events: list[Event] = []
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
    events.append(Event(time, "origination", lineage, family, founding.id))
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
        if any_driven and alive:
            values = [{key: trajs[key].value(alive[k], t) for key in trajs}
                      for k in range(len(alive))]
            if driven["duplication"]:
                w_dup = [dup.effective(copies=len(gen[k]), lineages=1, time=t, drivers=values[k])
                         * m_dup for k in range(len(alive))]
            if driven["loss"]:
                w_los = [los.effective(copies=len(gen[k]), lineages=1, time=t, drivers=values[k])
                         * m_los for k in range(len(alive))]
            if driven["transfer"] and can_xfer:
                w_tra = [tra.effective(copies=len(gen[k]), lineages=1, time=t, drivers=values[k])
                         * m_tra for k in range(len(alive))]
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
                if node.children is not None and g:         # speciation: re-id each copy into daughters
                    total -= len(g)
                    for ch in node.children:
                        child = [new_copy(old.family) for old in g]
                        for old, nc in zip(g, child):
                            events.append(Event(t, "speciation", ch, old.family, nc.id, parent=old.id))
                        enter(alive, gen, pos, ch, child)
                        heapq.heappush(heap, (tree.nodes[ch].end_time, ch))
                        total += len(child)
                else:                                       # a tip / extinction / empty: copies end here
                    total -= len(g)
        else:
            t = horizon                                     # a skyline or transfer-window breakpoint

    return events, node_genomes


def _evolve_one(family, lineage, birth_time, seedseq):
    """Worker adapter: evolve one family from the shared per-run context under its own spawned RNG
    stream (so the run is identical for any worker count). The pool ships the context once via the
    initializer; this reads it from module state and hands it to `simulate_one_family()`."""
    return simulate_one_family(_CTX, family=family, lineage=lineage, time=birth_time,
                               rng=np.random.default_rng(seedseq), copy_id_base=_copy_base(family))


def _evolve_family(task):
    """Collect-mode worker: evolve one family and hand its log back for the in-memory merge."""
    fid, birth_lineage, birth_time, seedseq = task
    events, node_genomes = _evolve_one(fid, birth_lineage, birth_time, seedseq)
    return fid, events, node_genomes


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
    if replacement:
        residents = [p for p, cpy in enumerate(rg) if cpy.family == fam and cpy.id != cont.id]
        if residents:                                      # homologous overwrite (else additive)
            p = residents[int(rng.integers(len(residents)))]
            victim = rg[p]
            rg[p] = rg[-1]; rg.pop()
            events.append(Event(t, "loss", recipient, fam, victim.id))
            delta = 0
    rg.append(xfer)
    events.append(Event(t, "transfer", donor, fam, cont.id, parent=src.id, donor=donor))
    events.append(Event(t, "transfer", recipient, fam, xfer.id, parent=src.id, recipient=recipient,
                        donor=donor))
    return delta


# --- output vocabulary for a streamed run (the files match FamilyGenomesResult.write) -------------------

#: the outputs a streamed run can produce and their top-level filenames — the same names the in-memory
#: ``FamilyGenomesResult.write`` uses. Gene trees are the exception: one Newick pair per family under a
#: ``gene_trees/`` subdirectory, so a million families do not land as two million files in the run root.
_STREAM_OUTPUTS = ("events", "profiles", "genomes", "initial_genome", "gene_trees")
_STREAM_FILENAMES = {"events": "genome_events.tsv", "profiles": "profiles.tsv",
                     "genomes": "genomes.tsv", "initial_genome": "initial_genome.tsv"}
_DEFAULT_STREAM_OUTPUTS = _STREAM_OUTPUTS

#: families per streamed chunk — **fixed**, independent of the worker count, so a chunk is a contiguous
#: family-id range and the shards concatenate in a deterministic order: a streamed run's files are
#: byte-identical for any number of workers. A million families is ~4000 chunks, plenty to fill the cores.
_STREAM_CHUNK = 256


@dataclass(frozen=True)
class StreamedRun:
    """A genome run written **straight to disk**, family by family — what ``stream_to=`` returns, for a
    scale where a whole `FamilyGenomesResult` would not fit in memory. Thin by design:
    the outputs *are* the files and the disk is the handoff (the sequences level reads them back), so
    this carries where they are and how big the run was, not the run itself."""

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
    chunk_index, family_list = task
    tree, s = _CTX.tree, _STREAM
    out_dir, outputs, extant_ids, shard_dir = s["out_dir"], s["outputs"], s["extant_ids"], s["shard_dir"]
    want = {name: name in outputs for name in ("events", "genomes", "profiles", "gene_trees")}
    trees_dir = os.path.join(out_dir, "gene_trees")

    files = {name: open(os.path.join(shard_dir, f"{name}_{chunk_index}.tsv"), "w", encoding="utf-8")
             for name in ("events", "genomes", "profiles") if want[name]}
    names = tree.labels()   # e<id> for a lineage that died; once per chunk, not once per family
    n_events = 0
    try:
        for (fid, lineage, birth_time, seedseq) in family_list:
            events, node_genomes = _evolve_one(fid, lineage, birth_time, seedseq)
            n_events += len(events)
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
                write_gene_trees(gene_trees_from_events(events, tree), trees_dir, names)
    finally:
        for f in files.values():
            f.close()
    return chunk_index, len(family_list), n_events


# --- the public entry: two passes, then either an in-memory merge or a stream to disk -------------

def run_parallel_family(tree, *, dup, tra, los, org, transfer_to, replacement, self_transfer,
                        initial_families, family_names, family_speed, cap, seed, parallel,
                        progress, stream_to=None, outputs=None,
                        trajs=None, to_traj=None, group_of=None, driven=None):
    """Run the per-family engine. Returns a `FamilyGenomesResult` (the in-memory
    merge), or a `StreamedRun` when ``stream_to`` is a directory — each family written straight
    to disk, for a scale a whole result would not hold. Returns ``None`` (a loud fallback to the serial
    loop) only for the in-memory path: a driven rate has no per-family engine, and a streamed run
    **raises** rather than silently pulling the whole thing back into memory.

    Copy ids are global from the start (``fid << SHIFT`` + local), so the merge and the streamed shards
    both just concatenate — no id rewrite, no run-sized bottleneck beyond the (serial) in-memory merge
    the streaming path exists to avoid."""
    from .family import GeneCopy, FamilyGenomesResult

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
    ctx = prepare_family_context(
        tree, dup=dup, tra=tra, los=los, transfer_to=transfer_to, replacement=replacement,
        self_transfer=self_transfer, cap=cap, family_speed=family_speed,
        trajs=trajs, to_traj=to_traj, group_of=group_of, driven=driven)

    # Pass 1: who originates, and where. One reserved stream for it; one per family after.
    root_ss = np.random.SeedSequence(seed)
    families_meta, named = _enumerate_families(
        tree, org, initial_families, family_names, np.random.default_rng(root_ss.spawn(1)[0]),
        trajs=trajs, driven=driven.get("origination", False))
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
        with ProcessPoolExecutor(max_workers=w, initializer=_init_worker, initargs=(ctx,)) as ex:
            for out in ex.map(_evolve_family, per_family, chunksize=max(1, n_families // (w * 8))):
                results.append(out); bar.update()
    else:
        _init_worker(ctx)
        for task in per_family:
            results.append(_evolve_family(task)); bar.update()
    bar.close()

    # The merge is now a concatenation: copy ids are already globally unique, so the per-family logs and
    # node snapshots stitch together with no rewrite. Every node appears, empty where no family reached it.
    events: list[Event] = []
    genomes: dict[int, list] = {i: [] for i in tree.nodes}
    for _fid, fam_events, node_genomes in results:
        events.extend(fam_events)
        for node_id, copies in node_genomes.items():
            genomes[node_id].extend(copies)
    events.sort(key=lambda e: e.time)                      # a chronological log, like the serial one
    genomes_final = {i: tuple(g) for i, g in genomes.items()}

    # the genome the run started with: every initial and named family's founding gene (its base id),
    # before the stem — the snapshot the serial engine takes as `initial_genome`.
    n_seeded = initial_families + len(family_names)
    initial_genome = tuple(GeneCopy(_copy_base(fid), fid) for fid in range(n_seeded))
    return FamilyGenomesResult(tree, genomes_final, events, seed, named, initial_genome)


def _run_streaming(tree, ctx, per_family, n_families, workers, seed, initial_families, family_names,
                   out_dir, outputs, progress):
    """The streaming half of `run_parallel_family()`: fixed contiguous chunks written to
    per-chunk shards, concatenated in chunk order (so the files are byte-identical for any worker
    count), then the shards removed. Returns a `StreamedRun`."""
    os.makedirs(out_dir, exist_ok=True)
    shard_dir = os.path.join(out_dir, "_shards")
    os.makedirs(shard_dir, exist_ok=True)
    extant_ids = sorted(n.id for n in tree.extant())
    stream_cfg = {"out_dir": out_dir, "outputs": set(outputs), "extant_ids": extant_ids,
                  "shard_dir": shard_dir}

    chunks = [per_family[i:i + _STREAM_CHUNK] for i in range(0, n_families, _STREAM_CHUNK)]
    tasks = list(enumerate(chunks))
    n_chunks = len(chunks)
    total_events = 0
    bar = progress_bar(max(1, n_chunks), "genomes", unit="chunk", enabled=progress)
    if workers > 1 and n_chunks >= 2:
        w = min(workers, n_chunks)
        with ProcessPoolExecutor(max_workers=w, initializer=_init_worker,
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
    headers = {"events": EVENTS_HEADER,
               "genomes": "lineage\tfamily\tcopy",
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
