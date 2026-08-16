"""Species with gene content, grown together.

The tree comes out of the run. Birth and death read a summary of each lineage's live genome
while that genome evolves by duplication, loss and origination on the growing tree."""

from __future__ import annotations

import math


from .._runtime.draw import weighted_index as _weighted_index
from .._runtime.grown import Grown
from ..genomes import FamilyGenome, GeneCopy, GeneEdge
from ..genomes.family import _duplicate, _lose_at, _originate, _pick_copy
from ..params.connection import Driven
from ..species import Event as SpeciesEvent
from ..tree import Node, Tree
from ._runaway import GENOME_COUNT as _GENOME_COUNT
from ._runaway import runaway as _runaway

def grow(rng, birth_rate, death_rate, spec: FamilyGenome, driver_names, n_extant,
                       total_time, max_lineages=None):
    """Grow a forward birth-death tree whose birth/death read the genome's **live gene content**, while
    the genome (duplication/loss/origination) evolves on that same growing tree. The species race and
    the genome's own D/L/O race run in one Gillespie over a shared living set. Returns
    ``(tree, species_events, genomes_out, genome_events, family_names)``."""
    dup, los, org = spec._resolve()

    # The same guard the joint TRAIT path has applied one function up: a mapping that can never fire
    # leaves every lineage at the default factor, so the run is the UNDRIVEN model while reporting that
    # gene content drove it. Both gene-content drivers have an alphabet known before the race starts —
    # a named family is present or absent, and a count is a number — so the check is exhaustive here
    # too. Without it a typo'd `{"presnt": 3.0}` ran to completion in silence.
    from ..params.conditioned import check_mapping_fires
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        for m in rate.modifiers:
            if not isinstance(m, Driven):
                continue
            if m.driver == _GENOME_COUNT:
                # a count is numeric: {state: factor} names discrete states a number never equals, and
                # `check_mapping_fires` says exactly that when the states it is given are numbers
                check_mapping_fires(m.mapping, {0}, driver_label=f"{label} (genomes:count)")
            else:
                check_mapping_fires(m.mapping, {"present", "absent"},
                                    driver_label=f"{label} ({m.driver})", exhaustive=True)

    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent, t):
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    copy_counter = 0
    family_counter = 0

    def new_copy(family):
        nonlocal copy_counter
        c = GeneCopy(copy_counter, family)
        copy_counter += 1
        return c

    def new_family():
        nonlocal family_counter
        f = family_counter
        family_counter += 1
        return f

    root = new_node(None, 0.0)
    alive = [root]          # living lineage ids
    gen: list[list] = [[]]  # each lineage's genome (list of GeneCopy), kept in lock-step with `alive`
    species_events: list[SpeciesEvent] = []
    genome_events: list[GeneEdge] = []
    for _ in range(spec.initial_families):  # anonymous families at the origin (t = 0)
        _originate(gen[0], nodes[root], 0.0, genome_events, new_copy, new_family)
    named: dict[str, int] = {}              # a minted id per declared name (the scaled_by("genomes:<name>") handles)
    for name in spec.family_names:
        fid = new_family()
        named[name] = fid
        c = new_copy(fid)
        gen[0].append(c)
        genome_events.append(GeneEdge(0.0, "origination", root, fid, c.id))
    total_copies = len(gen[0])
    genomes_out: dict[int, tuple] = {}

    def driver_value(src, k):
        if src == _GENOME_COUNT:
            return len(gen[k])                                   # a count → a Curve / Scalar
        fid = named[src.split(":", 1)[1]]                 # "genomes:<name>" → presence → a Table
        return "present" if any(c.family == fid for c in gen[k]) else "absent"

    t = 0.0
    ceiling = None if max_lineages is None else max(max_lineages, n_extant or 0)
    while alive:
        nl = len(alive)
        if ceiling is not None and nl > ceiling:
            raise _runaway(ceiling, t)
        drivers = [{s: driver_value(s, k) for s in driver_names} for k in range(nl)]
        wb = [birth_rate.effective(lineages=1, diversity=nl, time=t, drivers=drivers[k]) for k in range(nl)]
        wd = [death_rate.effective(lineages=1, diversity=nl, time=t, drivers=drivers[k]) for k in range(nl)]
        tb, td = sum(wb), sum(wd)
        # the genome's own dynamics are undriven → pooled over the whole live set (per copy / per lineage)
        r_dup = dup.effective(copies=total_copies, lineages=nl, time=t) if total_copies else 0.0
        r_los = los.effective(copies=total_copies, lineages=nl, time=t) if total_copies else 0.0
        r_org = org.effective(copies=total_copies, lineages=nl, time=t)
        total = tb + td + r_dup + r_los + r_org

        next_change = min(birth_rate.next_change(t), death_rate.next_change(t),
                          dup.next_change(t), los.next_change(t), org.next_change(t))
        horizon = next_change if total_time is None else min(next_change, total_time)

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:
                t = t_ev
                if nl == n_extant:
                    break
                r = float(rng.random()) * total
                if r < tb:  # speciation — the genome copies into both daughters (ZOMBI1 re-id)
                    i = _weighted_index(rng, wb, tb)
                    node_id, g = alive[i], gen[i]
                    alive[i] = alive[-1]; alive.pop()
                    gen[i] = gen[-1]; gen.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "speciation"
                    genomes_out[node_id] = tuple(g)
                    total_copies -= len(g)
                    c1, c2 = new_node(node_id, t), new_node(node_id, t)
                    node.children = (c1, c2)
                    for c in (c1, c2):
                        child = []
                        for old in g:
                            nc = new_copy(old.family)
                            child.append(nc)
                            genome_events.append(GeneEdge(t, "speciation", c, old.family, nc.id, parent=old.id))
                        alive.append(c); gen.append(child); total_copies += len(child)
                    species_events.append(SpeciesEvent(t, "speciation", node_id, (c1, c2)))
                elif r < tb + td:  # extinction
                    i = _weighted_index(rng, wd, td)
                    node_id, g = alive[i], gen[i]
                    alive[i] = alive[-1]; alive.pop()
                    gen[i] = gen[-1]; gen.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "extinct"
                    genomes_out[node_id] = tuple(g)
                    total_copies -= len(g)
                    species_events.append(SpeciesEvent(t, "extinction", node_id))
                elif r < tb + td + r_dup:  # duplication (per copy, pooled — the genome's own dynamics)
                    k, j = _pick_copy(rng, gen, total_copies)
                    _duplicate(gen[k], j, nodes[alive[k]], t, genome_events, new_copy)
                    total_copies += 1
                elif r < tb + td + r_dup + r_los:  # loss (per copy, pooled)
                    k, j = _pick_copy(rng, gen, total_copies)
                    _lose_at(gen[k], j, nodes[alive[k]], t, genome_events)
                    total_copies -= 1
                else:  # origination (per lineage, uniform)
                    k = int(rng.integers(nl))
                    _originate(gen[k], nodes[alive[k]], t, genome_events, new_copy, new_family)
                    total_copies += 1
                continue

        if math.isinf(horizon):
            break
        if total_time is not None and horizon == total_time:
            t = total_time
            break
        t = horizon

    for k, node_id in enumerate(alive):  # survivors reach the present
        nodes[node_id].end_time = t
        nodes[node_id].fate = "extant"
        genomes_out[node_id] = tuple(gen[k])
    return Grown(Tree(nodes, root), species_events, genomes=genomes_out,
                 genome_events=genome_events, genome_names=named)


