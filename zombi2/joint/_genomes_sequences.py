"""A genome and a gene's sequence, each driving the other, on a tree the run is handed.

The last cell of the map (the manual's Dependent runs chapter). The genome decides which sequences
exist — when a gene is copied, when a copy moves lineage, when a copy dies — and the sequences decide
how fast the genome changes, because a genome rate reads their composition.

Neither half can go first. To evolve the sequences you need to know which copies exist and when. To
know which copies exist you need the loss rate, which reads the sequences.

**The walk.** This is the genome race of `_genomes_traits`, on the same given tree, with one thing
added: every live copy of the **declared** gene carries its own sequence. Species time is sliced.
Inside a slice the composition each lineage's rates read is held where it was; at the boundary every
live copy's sequence is advanced and the composition is read again.

A genome event does one of three things to a sequence: clone one, end one, or draw a fresh one. All
three fall out of the edges the event appends. An edge with a parent means that parent ended and this
copy descends from it, so the parent's sequence is carried to the instant, kept, and cloned. A loss
edge means the copy ended. An origination edge means a new family, which is never the declared one.

That the parent is carried **to the instant** before cloning is the point. Cloning the sequence as it
stood at the top of the slice would have the two copies each redraw the stretch they actually shared.
"""

from __future__ import annotations

import numpy as np

from .._runtime.grown import Grown
from ..genomes import FamilyGenome, GeneCopy, GeneEdge, GeneFamily
from ..genomes.family import _duplicate, _lose_at, _originate
from ..params.connection import Driven
from ..params.driver import OnTime
from ..params.evaluate import describe
from ..params.scope import PerLineage, PerSite
from ._runaway import GENOME_COUNT as _GENOME_COUNT


def grow(rng, tree, genome: FamilyGenome, gene, gene_keys, step: float, record=None) -> Grown:
    """Race the genome on ``tree`` while the declared gene's sequences evolve beside it.

    ``gene`` is the `~zombi2.sequences.GeneSpec`; ``gene_keys`` are the lookup keys the genome rates
    read its composition under. Returns a `~zombi2._runtime.grown.Grown` carrying both levels."""
    from ..genomes._live import enter, retire, weighted_index
    from ..genomes._transfer import mean_root_to_tip
    from ..genomes.family import _FamilyCounts, _do_transfer, resolve_families
    from ..params.parameter import as_rate
    from ..params.scope import PerCopy
    from ..sequences._record import walk as _walk
    from ..sequences.evolution import _cdf_for, _sample

    model, sites = gene.model, gene.length
    letters = tuple(model.alphabet.index(c) for c in gene.offers.letters)
    absent = gene.offers.absent
    founder = model if gene.start is None else gene.start
    sub = as_rate(1.0 if gene.substitution is None else gene.substitution, default_scope=PerSite)
    if sub.modifiers:
        raise ValueError(
            f"gene {gene.name!r}'s substitution rate carries "
            f"{', '.join(dict.fromkeys(describe(m) for m in sub.modifiers))}. Here the gene is what "
            f"the genome reads, so its own rate is a plain number: PerSite(1.0).")
    assert sub.base is not None                  # `as_rate` fills the level's default
    rate_base = float(sub.base)

    dup = as_rate(genome.duplication, default_scope=PerCopy, label="duplication")
    los = as_rate(genome.loss, default_scope=PerCopy, label="loss")
    tra = as_rate(genome.transfer, default_scope=PerCopy, label="transfer")
    org = as_rate(genome.origination, default_scope=PerLineage, label="origination")
    for label, rate in (("duplication", dup), ("loss", los), ("transfer", tra),
                        ("origination", org)):
        for m in rate.modifiers:
            if not isinstance(m, (Driven, OnTime)):
                raise ValueError(
                    f"{label} carries {describe(m)}, which this engine does not thread. Here a "
                    f"genome rate takes changing_at and scaled_by — the verb that reads the gene "
                    f"simulated beside it.")
    declared, _modules, _planted = resolve_families(
        [GeneFamily(n) for n in genome.family_names], tree)

    counter = {"copy": 0, "family": 0}

    def new_copy(fam):
        c = GeneCopy(counter["copy"], fam)
        counter["copy"] += 1
        return c

    def new_family():
        f = counter["family"]
        counter["family"] += 1
        return f

    root = tree.nodes[tree.root]
    t = root.birth_time
    alive: list[int] = []
    gen: list[list] = []
    pos: dict[int, int] = {}
    genomes_out: dict[int, tuple] = {}
    events: list[GeneEdge] = []
    enter(alive, gen, pos, root.id, [])
    for _ in range(genome.initial_families):
        _originate(gen[0], root, t, events, new_copy, new_family)
    named: dict[str, int] = {}
    for spec in declared:
        fid = new_family()
        named[spec.name] = fid
        c = new_copy(fid)
        gen[0].append(c)
        events.append(GeneEdge(t, "origination", root.id, fid, c.id))
    if gene.name not in named:
        raise ValueError(
            f"gene {gene.name!r} names no family this genome spec declared. Declare it there — "
            f"genomes.genome(..., families=[family({gene.name!r})]).")
    watched = named[gene.name]

    # the sequence side: only the declared gene's copies carry one
    live_seq: dict[int, np.ndarray] = {}
    carried_to: dict[int, float] = {}
    accrued: dict[int, float] = {}
    kept: dict[int, np.ndarray] = {}
    length_of: dict[int, float] = {}
    founding = rng.choice(model.k, size=sites, p=founder.stationary).astype(np.int8)
    for c in gen[0]:
        if c.family == watched:
            live_seq[c.id] = founding
            carried_to[c.id] = t
            accrued[c.id] = 0.0
    cache: dict[tuple[int, float], np.ndarray] = {}

    def carry(copy_id: int, lineage: int, until: float) -> None:
        """Advance one copy's sequence to ``until``. The one place substitutions happen."""
        was = carried_to[copy_id]
        if until <= was:
            return
        bl = rate_base * (until - was)
        if bl > 0.0 and record is not None:
            live_seq[copy_id], rows = _walk(live_seq[copy_id], model.Q, bl, rng)
            record.add(_Where(lineage, copy_id), was, until - was, rows, model.alphabet)
        elif bl > 0.0:
            live_seq[copy_id] = _sample(live_seq[copy_id], _cdf_for(cache, model, bl), rng)
        accrued[copy_id] += bl
        carried_to[copy_id] = until

    def settle(copy_id: int) -> None:
        """This copy's branch has ended: keep what it holds and stop carrying it."""
        kept[copy_id] = live_seq.pop(copy_id)
        length_of[copy_id] = accrued.pop(copy_id)
        carried_to.pop(copy_id)

    def follow(mark: int, when: float) -> None:
        """Apply the edges an event just appended to the sequences that event touched."""
        for e in events[mark:]:
            if e.parent is not None and e.parent in live_seq:
                carry(e.parent, e.lineage, when)
                live_seq[e.copy] = live_seq[e.parent].copy()
                carried_to[e.copy] = when
                accrued[e.copy] = 0.0
            elif e.kind == "loss" and e.copy in live_seq:
                carry(e.copy, e.lineage, when)
                settle(e.copy)
        for e in events[mark:]:                    # every parent has handed its sequence on by now
            if e.parent is not None and e.parent in live_seq:
                settle(e.parent)

    counts = _FamilyCounts(gen)
    depth = mean_root_to_tip(tree)
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)
    si = 0
    total_copies = len(gen[0])
    initial_genome = tuple(gen[0])
    slice_end = t + step
    shares = _shares(gen, alive, live_seq, watched, letters)

    while si < len(schedule):
        n_alive = len(alive)
        drivers: list[dict] = []
        for k in range(n_alive):
            d: dict = {_GENOME_COUNT: len(gen[k])}
            for name, fid in named.items():
                d[f"genomes:{name}"] = "present" if counts.holds(k, fid) else "absent"
            for key in gene_keys:
                d[key] = shares.get(alive[k], absent)
            drivers.append(d)
        can_xfer = total_copies > 0 and n_alive >= 2
        w_dup = [dup.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                 for k in range(n_alive)]
        w_los = [los.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                 for k in range(n_alive)]
        w_org = [org.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                 for k in range(n_alive)]
        w_tra = ([tra.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                  for k in range(n_alive)] if can_xfer else [0.0] * n_alive)
        r_dup, r_los, r_org, r_tra = (sum(w) for w in (w_dup, w_los, w_org, w_tra))
        total = r_dup + r_los + r_org + r_tra

        horizon = min(schedule[si][0], slice_end, dup.next_change(t), los.next_change(t),
                      org.next_change(t), tra.next_change(t))
        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:
                t = t_ev
                mark = len(events)
                r = float(rng.random()) * total
                if r < r_dup:
                    k = weighted_index(rng, w_dup, r_dup)
                    j = int(rng.integers(len(gen[k])))
                    fam = gen[k][j].family
                    if not counts.at_cap(k, fam, genome.max_family_size):
                        _duplicate(gen[k], j, tree.nodes[alive[k]], t, events, new_copy)
                        counts.added(k, fam); total_copies += 1
                elif r < r_dup + r_los:
                    k = weighted_index(rng, w_los, r_los)
                    j = int(rng.integers(len(gen[k])))
                    counts.removed(k, gen[k][j].family)
                    _lose_at(gen[k], j, tree.nodes[alive[k]], t, events)
                    total_copies -= 1
                elif r < r_dup + r_los + r_org:
                    k = weighted_index(rng, w_org, r_org)
                    _originate(gen[k], tree.nodes[alive[k]], t, events, new_copy, new_family)
                    counts.added(k, gen[k][-1].family); total_copies += 1
                else:
                    kd = weighted_index(rng, w_tra, r_tra)
                    jd = int(rng.integers(len(gen[kd])))
                    delta, _kr = _do_transfer(rng, tree, alive, gen, counts, kd, jd, t, events,
                                              new_copy, "uniform", False, False, depth, None,
                                              genome.max_family_size, None)
                    total_copies += delta
                follow(mark, t)
                continue

        if horizon == slice_end and slice_end < schedule[si][0]:
            for k in range(len(alive)):                # the boundary: release the frozen driver
                for c in gen[k]:
                    if c.id in live_seq:
                        carry(c.id, alive[k], slice_end)
            t = slice_end
            slice_end += step
            shares = _shares(gen, alive, live_seq, watched, letters)
            continue
        if horizon == schedule[si][0]:
            t = schedule[si][0]
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                k_out = pos[i]
                g = gen[k_out]
                for c in g:                            # a node's sequence is what it holds at its end
                    if c.id in live_seq:
                        carry(c.id, i, t)
                genomes_out[i] = tuple(g)
                total_copies -= len(g)
                retire(alive, gen, pos, k_out)
                inherited = counts.retired(k_out)
                node = tree.nodes[i]
                mark = len(events)
                if node.children:
                    per_daughter = []
                    for c in node.children:
                        child, rows = [], []
                        for old in g:
                            nc = new_copy(old.family)
                            child.append(nc)
                            rows.append(GeneEdge(t, "speciation", c, old.family, nc.id,
                                                 parent=old.id))
                        per_daughter.append(rows)
                        enter(alive, gen, pos, c, child)
                        counts.entered_like(inherited)
                        total_copies += len(child)
                    for pair in zip(*per_daughter):
                        events.extend(pair)
                    follow(mark, t)
                else:                                  # a tip: every copy it holds ends here
                    for c in g:
                        if c.id in live_seq:
                            settle(c.id)
                si += 1
            continue
        t = horizon

    for copy_id in list(live_seq):                     # anything still in flight ends at the present
        settle(copy_id)
    return Grown(tree, genomes=genomes_out, genome_events=events, genome_names=named,
                 genome_initial=initial_genome,
                 sequences={gene.name: (kept, founding, length_of)})


def _shares(gen, alive, live_seq, watched, letters) -> dict[int, float]:
    """Each lineage's share of the counted letters, pooled over the copies of the watched family it
    is carrying right now. A lineage with none is simply not a key; the caller falls back to the
    ``absent`` the gene declared."""
    out: dict[int, float] = {}
    for k, lineage in enumerate(alive):
        hits = total = 0
        for c in gen[k]:
            if c.family == watched and c.id in live_seq:
                states = live_seq[c.id]
                hits += int(np.isin(states, letters).sum())
                total += states.size
        if total:
            out[lineage] = hits / total
    return out


class _Where:
    """What `~zombi2.sequences._record.Recorder` reads off a gene-tree node, for a race that has no
    node yet: the species lineage it is on, and the copy it is."""

    __slots__ = ("species", "copy")

    def __init__(self, species, copy):
        self.species = species
        self.copy = copy
