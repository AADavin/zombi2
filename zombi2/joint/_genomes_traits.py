"""A genome and a discrete trait on a tree the run is handed, each reading the other.

The first joint model whose tree is an input. One Gillespie over the living set, with five
event classes: the genome's duplication, transfer, loss and origination, and the trait's
switch. Transfer works here, unlike the tree-growing models, because on a fixed tree the
set of lineages alive at an instant is already known."""

from __future__ import annotations



from .._runtime.grown import Grown
from ..genomes import FamilyGenome, GeneEdge, GeneFamily
from ..genomes.family import _duplicate, _lose_at, _originate
from ..params.connection import Driven
from ..params.driver import OnTime
from ..params.evaluate import describe
from ..params.scope import PerLineage
from ..traits import Change, DiscreteTrait
from ._runaway import GENOME_COUNT as _GENOME_COUNT

def grow(rng, tree, genome: FamilyGenome, trait: DiscreteTrait, trait_keys,
                         seed) -> Grown:
    """A genome and a discrete trait on a **given** tree, each reading the other as the run goes.

    The first joint model whose tree is an input. The two halves exist already but in different
    shapes: the genome level races its events over every living lineage at once, while the discrete
    trait walks one branch at a time. Neither survives the merge on its own terms — a trait that
    walks a whole branch cannot see a gene lost half-way down it — so the trait moves into the
    genome's race, exactly as it does in `_grow_joint()` against speciation.

    One Gillespie over the living set, then, with five event classes: the genome's
    duplication / transfer / loss / origination, and the trait's switch. Each lineage's rates read
    the other level's state on that lineage — the genome's from ``drivers[k]``, the trait's by
    rebuilding its generator per lineage (`_driven_q`). Both are piecewise-constant between events
    and every event ends the step, so the race is exact and nothing is thinned.

    **Transfer works here**, unlike the tree-growing joint models, and for the reason they refuse it:
    on a tree handed to the run the set of lineages alive at an instant is already known.

    Returns a `~zombi2._runtime.grown.Grown`.
    """
    from ..genomes.family import (_FamilyCounts, _do_transfer, GeneCopy, resolve_families)
    from ..genomes._live import enter, retire, weighted_index
    from ..genomes._transfer import mean_root_to_tip
    from ..traits.discrete import _driven_entries, _driven_q
    from ..params.parameter import as_rate
    from ..params.scope import PerCopy

    # `DiscreteTrait._resolve` settles the generator into one constant matrix, which is exactly what a
    # switch rate reading the genome cannot be. So the alphabet and the root state are taken from the
    # spec here and the generator is left as rate specs, rebuilt per lineage below.
    states = list(trait.states)
    k_states = len(states)
    idx = {s: i for i, s in enumerate(states)}
    if trait.start is None:
        start_i = int(rng.integers(k_states))
    elif trait.start in idx:
        start_i = idx[trait.start]
    else:
        raise ValueError(f"start must be one of states={states} (or None for a uniform draw), "
                         f"got {trait.start!r}")
    # `DiscreteTrait._resolve` checks this, and this engine does not call it (its generator is one
    # constant matrix, which a switch rate reading the genome cannot be), so the check comes along
    at_split = trait.at_speciation
    if at_split is not None and (isinstance(at_split, bool)
                                 or not isinstance(at_split, (int, float))
                                 or not 0.0 <= at_split <= 1.0):
        raise ValueError(f"at_speciation must be a probability in [0, 1] (the shift chance), "
                         f"got {at_split!r}")
    shift = 0.0 if at_split is None else float(at_split)
    entries = _driven_entries(states, trait.switch)          # rate specs, so they can read the genome
    # `FamilyGenome._resolve` refuses a driven rate, and rightly: where the tree is being simulated
    # the genome is what drives it. Here the tree is given and the genome is a target as well, so the
    # rates are resolved on their own terms.
    dup = as_rate(genome.duplication, default_scope=PerCopy, label="duplication")
    los = as_rate(genome.loss, default_scope=PerCopy, label="loss")
    tra = as_rate(genome.transfer, default_scope=PerCopy, label="transfer")
    org = as_rate(genome.origination, default_scope=PerLineage, label="origination")
    for label, rate in (("duplication", dup), ("loss", los), ("transfer", tra),
                        ("origination", org)):
        for m in rate.modifiers:
            if not isinstance(m, (Driven, OnTime)):
                raise ValueError(
                    f"{label} carries {describe(m)}, which this engine does not thread. On a joint "
                    f"run over a given tree a genome rate takes changing_at and scaled_by — the "
                    f"verb that reads the trait simulated beside it.")
    declared, _modules, planted = resolve_families(
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
    st: list[int] = []                       # each lineage's trait state index, parallel to `alive`
    genomes_out: dict[int, tuple] = {}
    node_state: dict[int, int] = {}
    events: list[GeneEdge] = []
    changes: list[Change] = [Change(t, "initial", root.id, None, states[start_i])]
    enter(alive, gen, pos, root.id, [])
    st.append(start_i)
    for _ in range(genome.initial_families):
        _originate(gen[0], root, t, events, new_copy, new_family)
    named: dict[str, int] = {}
    for spec in declared:
        fid = new_family()
        named[spec.name] = fid
        c = new_copy(fid)
        gen[0].append(c)
        events.append(GeneEdge(t, "origination", root.id, fid, c.id))
    total_copies = len(gen[0])
    initial_genome = tuple(gen[0])
    counts = _FamilyCounts(gen)
    depth = mean_root_to_tip(tree)
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)
    si = 0

    def gene_drivers(k):
        """What the trait's switch rate and the genome's own rates read on lineage ``k``."""
        d = {_GENOME_COUNT: len(gen[k])}
        for name, fid in named.items():
            d[f"genomes:{name}"] = "present" if counts.holds(k, fid) else "absent"
        return d

    while si < len(schedule):
        n_alive = len(alive)
        drivers = []
        for k in range(n_alive):
            d = gene_drivers(k)
            for key in trait_keys:
                d[key] = states[st[k]]
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
        # the trait's generator is rebuilt per lineage, because its entries read that lineage's genome
        qs = [_driven_q(entries, k_states, drivers[k], t) for k in range(n_alive)]
        w_sw = [float(-qs[k][st[k], st[k]]) for k in range(n_alive)]
        r_dup, r_los, r_org, r_tra, r_sw = (sum(w) for w in (w_dup, w_los, w_org, w_tra, w_sw))
        total = r_dup + r_los + r_org + r_tra + r_sw

        horizon = min(schedule[si][0], dup.next_change(t), los.next_change(t),
                      org.next_change(t), tra.next_change(t))
        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:
                t = t_ev
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
                elif r < r_dup + r_los + r_org + r_tra:
                    kd = weighted_index(rng, w_tra, r_tra)
                    jd = int(rng.integers(len(gen[kd])))
                    delta, _kr = _do_transfer(rng, tree, alive, gen, counts, kd, jd, t, events,
                                              new_copy, "uniform", False, False, depth, None,
                                              genome.max_family_size, None)
                    total_copies += delta
                else:                                  # the trait switches; the topology is untouched
                    k = weighted_index(rng, w_sw, r_sw)
                    cur = st[k]
                    probs = qs[k][cur].copy()
                    probs[cur] = 0.0
                    probs /= float(-qs[k][cur, cur])
                    new = int(rng.choice(k_states, p=probs))
                    st[k] = new
                    changes.append(Change(t, "on_branch", alive[k], states[cur], states[new]))
                continue

        if horizon == schedule[si][0]:
            t = schedule[si][0]
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                k_out = pos[i]
                g = gen[k_out]
                genomes_out[i] = tuple(g)
                node_state[i] = st[k_out]
                total_copies -= len(g)
                cur = st[k_out]
                retire(alive, gen, pos, k_out)
                st[k_out] = st[-1]; st.pop()           # mirror the swap-remove, or the arrays desync
                inherited = counts.retired(k_out)
                node = tree.nodes[i]
                if node.children:
                    per_daughter = []
                    for c in node.children:
                        child, rows = [], []
                        for old in g:
                            nc = new_copy(old.family)
                            child.append(nc)
                            rows.append(GeneEdge(t, "speciation", c, old.family, nc.id, parent=old.id))
                        per_daughter.append(rows)
                        enter(alive, gen, pos, c, child)
                        counts.entered_like(inherited)
                        total_copies += len(child)
                        d = cur
                        if shift > 0.0 and float(rng.random()) < shift:
                            jj = int(rng.integers(k_states - 1))     # a uniform *other* state
                            d = jj if jj < cur else jj + 1
                            changes.append(Change(t, "on_speciation", c, states[cur], states[d]))
                        st.append(d)
                    for pair in zip(*per_daughter):
                        events.extend(pair)
                si += 1
        else:
            t = horizon
    node_values = {i: states[node_state[i]] for i in tree.nodes}
    return Grown(tree, genomes=genomes_out, genome_events=events, genome_names=named,
                 genome_initial=initial_genome, trait_values=node_values,
                 trait_events=changes)


