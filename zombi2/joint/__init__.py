"""Joint models — the driver co-evolves with what it drives (SPEC §2–4, ``coupling-api.md``).

When the driver **cannot** be grown first — because it is entangled with what it drives as the tree
unfolds — one run must produce both. Two drivers of speciation are grown here, the tree an **output**:

- a **discrete trait** drives speciation (BiSSE / MuSSE), ``P(Species, Traits)`` — birth/death read the
  trait state on each lineage while the trait evolves by its own Mk process on the growing tree;
- **gene content** drives speciation, ``P(Species, Genomes)`` — birth/death read a summary of each
  lineage's live genome (its total gene count, or the presence of a named family) while the genome
  evolves by duplication/loss/origination on the growing tree.

One Gillespie races the event classes over the living lineages at once: **speciation** and
**extinction** (per lineage, driver-read), plus the driver's own events — a **trait switch** (the CTMC
out-rate) or a genome **duplication/loss/origination**. A driver event changes a lineage's state
without touching the topology; a speciation hands the parent's driver state (its trait, its genome) to
both daughters. Because these drivers only change at events, the rate is piecewise-constant between
them and the race is **exact** — no thinning. (A continuously-diffusing driver — QuaSSE — is deferred:
it makes the rate vary continuously, which needs thinning.)

The mechanism is the same ``mod.DrivenBy`` as conditioning; only the ``source`` differs — here a
**live level name** (``"trait"``, ``"genomes:count"``, ``"genomes:<family>"``) rather than a filename.
Driving *both* birth and death recovers full state-dependent diversification (BiSSE's λ and μ)."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..genomes import (
    Event as GenomeEvent,
    GeneCopy,
    GenomesResult,
    UnorderedGenome,
    _duplicate,
    _lose_at,
    _originate,
    _pick_copy,
)
from ..rates.modifiers import DrivenBy, FromParent
from ..rates.rate import as_rate
from ..rates.scope import PerLineage
from ..species import Event, SpeciesResult
from ..tree import Node, Tree
from ..traits import Change, DiscreteTrait, TraitsResult

_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_extant
_GENOME_COUNT = "genomes:count"  # the live gene-content driver source for a count → Curve/Scalar


@dataclass
class JointResult:
    """What :func:`simulate_joint` returns — **both** grown levels of a joint run. ``species`` is the
    grown tree (a :class:`~zombi2.species.SpeciesResult`: ``complete_tree``, ``extant_tree``, the
    speciation/extinction ``events``); the **driver** level that grew with it is either ``trait`` (a
    :class:`~zombi2.traits.TraitsResult`, for a trait→speciation run) or ``genome`` (a
    :class:`~zombi2.genomes.GenomesResult`, for a gene-content→speciation run) — exactly one is set.
    The tree is an output, grown by the driver it carries, so the levels share one ``complete_tree``."""

    species: SpeciesResult
    seed: int | None
    trait: TraitsResult | None = None
    genome: GenomesResult | None = None

    @property
    def complete_tree(self) -> Tree:
        return self.species.complete_tree

    @property
    def extant_tree(self):
        return self.species.extant_tree

    @property
    def n_extant(self) -> int:
        return self.species.n_extant

    @property
    def events(self) -> list:
        """The species events (speciation / extinction). The driver level's own events are
        ``trait.events`` / ``genome.events``."""
        return self.species.events

    def write(self, directory) -> None:
        """Write both levels to ``directory``: the species files (``species_complete.nwk`` /
        ``species_extant.nwk`` / ``species_events.tsv``) and the driver level's — for a trait,
        ``trait_values.tsv`` / ``trait_events.tsv`` / ``trait_tree.nwk``; for a genome,
        ``genome_events.tsv`` / ``profiles.tsv``."""
        self.species.write(directory, outputs=("complete", "extant", "events"))
        if self.trait is not None:
            self.trait.write(directory, outputs=("values", "events", "tree"))
        if self.genome is not None:
            self.genome.write(directory, outputs=("events", "profiles"))


def _weighted_index(rng, weights: list[float], total: float) -> int:
    """Pick a lineage index in proportion to ``weights`` (which sum to ``total``) — the same
    per-lineage pick the species and genome engines use when a rate varies across lineages."""
    r = rng.random() * total
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r < acc:
            return i
    return len(weights) - 1  # floating-point guard: r == total lands on the last lineage


def _grow_joint(rng, birth_rate, death_rate, trait: DiscreteTrait, n_extant, total_time):
    """Grow a forward birth-death tree whose birth/death read a discrete trait that evolves on it.
    Returns ``(tree, species_events, node_values, trait_events)`` — the complete tree, the
    speciation/extinction log, the trait state at every node, and the trait's switch log."""
    states, Q, start_i, shift = trait._resolve(rng)
    k_states = len(states)
    out_rate = [float(-Q[s, s]) for s in range(k_states)]  # the trait's total switch-out rate per state

    # birth/death are driven by the trait; a mapping whose states are none of the trait's would leave
    # every lineage at the default factor — a silently uncoupled run — so refuse it up front
    from ..rates.driver import check_mapping_fires
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        for m in rate.modifiers:
            if isinstance(m, DrivenBy):
                check_mapping_fires(m.mapping, states, source_label=f"{label} (trait)")

    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent, t):
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    root = new_node(None, 0.0)
    alive = [root]      # living lineage ids
    st = [start_i]      # each lineage's trait state index, kept in lock-step with `alive`
    t = 0.0
    species_events: list[Event] = []
    trait_events: list[Change] = []
    end_state: dict[int, int] = {}  # node id → its trait state index when it ended (→ node_values)

    while alive:
        n = len(alive)
        ctx = {"diversity": n, "time": t}
        # per-lineage rates: birth/death read the lineage's trait state (DrivenBy("trait", …)); the
        # trait switch rate is the CTMC out-rate for that state (the trait's own dynamics, undriven).
        wb = [birth_rate.effective(lineages=1, drivers={"trait": states[st[k]]}, **ctx) for k in range(n)]
        wd = [death_rate.effective(lineages=1, drivers={"trait": states[st[k]]}, **ctx) for k in range(n)]
        ws = [out_rate[st[k]] for k in range(n)]
        total_b, total_d, total_s = sum(wb), sum(wd), sum(ws)
        total = total_b + total_d + total_s

        # the trait switch rate is constant between events; only a skyline (OnTime) on birth/death or
        # the total_time limit advances the clock on its own.
        next_change = min(birth_rate.next_change(t), death_rate.next_change(t))
        horizon = next_change if total_time is None else min(next_change, total_time)

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:
                t = t_ev
                if n == n_extant:  # already at the target; stop at this next event's time, unapplied
                    break
                r = float(rng.random()) * total
                if r < total_b:  # speciation
                    i = _weighted_index(rng, wb, total_b)
                    node_id, cur = alive[i], st[i]
                    alive[i] = alive[-1]; alive.pop()          # swap-remove keeps the state array in step
                    st[i] = st[-1]; st.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "speciation"
                    end_state[node_id] = cur
                    c1, c2 = new_node(node_id, t), new_node(node_id, t)
                    node.children = (c1, c2)
                    for c in (c1, c2):  # each daughter inherits the parent's state (+ optional split shift)
                        d = cur
                        if shift > 0.0 and float(rng.random()) < shift:
                            j = int(rng.integers(k_states - 1))  # hop to a uniform *other* state
                            d = j if j < cur else j + 1
                            trait_events.append(Change(t, "on_speciation", c, states[cur], states[d]))
                        alive.append(c); st.append(d)
                    species_events.append(Event(t, "speciation", node_id, (c1, c2)))
                elif r < total_b + total_d:  # extinction
                    i = _weighted_index(rng, wd, total_d)
                    node_id, cur = alive[i], st[i]
                    alive[i] = alive[-1]; alive.pop()
                    st[i] = st[-1]; st.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "extinct"
                    end_state[node_id] = cur
                    species_events.append(Event(t, "extinction", node_id))
                else:  # trait switch — change one lineage's state, no topology change
                    i = _weighted_index(rng, ws, total_s)
                    node_id, cur = alive[i], st[i]
                    probs = Q[cur].copy()
                    probs[cur] = 0.0
                    probs /= out_rate[cur]          # the embedded jump chain: where to, given a jump
                    new = int(rng.choice(k_states, p=probs))
                    st[i] = new
                    trait_events.append(Change(t, "on_branch", node_id, states[cur], states[new]))
                continue

        if math.isinf(horizon):
            break  # nothing scheduled and no skyline change → nothing more can happen
        if total_time is not None and horizon == total_time:
            t = total_time
            break
        t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) birth/death

    for k, node_id in enumerate(alive):  # whoever is still alive reached the present
        nodes[node_id].end_time = t
        nodes[node_id].fate = "extant"
        end_state[node_id] = st[k]

    node_values = {i: states[end_state[i]] for i in nodes}
    return Tree(nodes, root), species_events, node_values, trait_events


def _grow_joint_genome(rng, birth_rate, death_rate, spec: UnorderedGenome, sources, n_extant, total_time):
    """Grow a forward birth-death tree whose birth/death read the genome's **live gene content**, while
    the genome (duplication/loss/origination) evolves on that same growing tree. The species race and
    the genome's own D/L/O race run in one Gillespie over a shared living set. Returns
    ``(tree, species_events, genomes_out, genome_events, family_names)``."""
    dup, los, org = spec._resolve()

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
    species_events: list[Event] = []
    genome_events: list[GenomeEvent] = []
    for _ in range(spec.initial_families):  # anonymous crown families
        _originate(gen[0], nodes[root], 0.0, genome_events, new_copy, new_family)
    family_names: dict[str, int] = {}       # named crown families (the DrivenBy("genomes:<name>") handles)
    for name in spec.families:
        fid = new_family()
        family_names[name] = fid
        c = new_copy(fid)
        gen[0].append(c)
        genome_events.append(GenomeEvent(0.0, "origination", root, fid, c.id))
    total_copies = len(gen[0])
    genomes_out: dict[int, tuple] = {}

    def driver_value(src, k):
        if src == _GENOME_COUNT:
            return len(gen[k])                                   # a count → a Curve / Scalar
        fid = family_names[src.split(":", 1)[1]]                 # "genomes:<name>" → presence → a Table
        return "present" if any(c.family == fid for c in gen[k]) else "absent"

    t = 0.0
    while alive:
        nl = len(alive)
        drivers = [{s: driver_value(s, k) for s in sources} for k in range(nl)]
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
                            genome_events.append(GenomeEvent(t, "speciation", c, old.family, nc.id, parent=old.id))
                        alive.append(c); gen.append(child); total_copies += len(child)
                    species_events.append(Event(t, "speciation", node_id, (c1, c2)))
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
                    species_events.append(Event(t, "extinction", node_id))
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
    return Tree(nodes, root), species_events, genomes_out, genome_events, family_names


def simulate_joint(*, birth, death=0.0, trait=None, genome=None, n_extant=None, total_time=None,
                   seed=None) -> JointResult:
    """Grow a tree **and** the driver that drives its speciation, in one run (SPEC §2–4).

    ``birth`` and ``death`` are rate specs (per lineage). Make either read the driver with
    ``mod.DrivenBy(source, mapping)`` — a **live level name** (not a filename) is what makes this
    *joint* rather than conditioned. Give **exactly one** driver:

    - ``trait = traits.discrete(...)`` — a discrete trait drives speciation (BiSSE / MuSSE), read as
      ``mod.DrivenBy("trait", {"small": 1.0, "large": 2.0})``. Driving both birth and death gives
      state-dependent λ *and* μ.
    - ``genome = genomes.unordered(...)`` — **gene content** drives speciation (``P(Species, Genomes)``),
      read as the total gene count ``mod.DrivenBy("genomes:count", curve)`` or the presence of a named
      family ``mod.DrivenBy("genomes:toxin", {"present": 2.0, "absent": 1.0})`` (declare it with
      ``families=["toxin"]``).

    ::

        joint.simulate_joint(
            birth  = 1.0 * mod.DrivenBy("genomes:toxin", {"present": 3.0, "absent": 1.0}),
            genome = genomes.unordered(origination=0.2, loss=0.1, families=["toxin"]),
            n_extant = 100, seed = 1)

    The driver is an **unexecuted** process spec, grown with the tree. Stop at exactly ``n_extant``
    living lineages (conditioned on survival — a birth-death tree can die out, so it restarts,
    advancing the same generator) **or** at ``total_time`` — give exactly one. Returns a
    :class:`JointResult` carrying the grown tree and the driver level (``.trait`` or ``.genome``).
    Deterministic given ``seed``. Continuous trait→speciation (QuaSSE), clade drift (``FromParent``)
    combined with driving, and gene transfer in a joint run are later slices.
    """
    birth_rate = as_rate(birth, default_scope=PerLineage)
    death_rate = as_rate(death, default_scope=PerLineage)
    if (trait is None) == (genome is None):
        raise TypeError(
            "give exactly one driver: trait=traits.discrete(...) OR genome=genomes.unordered(...)."
        )
    # collect the DrivenBy sources on birth/death (a joint model's diversification must be per lineage)
    sources: list[str] = []
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        if not isinstance(rate.scope, PerLineage):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but a joint diversification rate is "
                f"per lineage — drop the scope wrapper (per lineage is the default)."
            )
        for m in rate.modifiers:
            if isinstance(m, FromParent):
                raise ValueError(
                    f"{label} carries FromParent (clade drift); drift combined with a driven rate is a "
                    f"later slice — use one or the other."
                )
            if isinstance(m, DrivenBy):
                if not isinstance(m.source, str):
                    raise TypeError(
                        f"{label} is driven by a {type(m.source).__name__} object, but a joint model "
                        f"drives from a live level *name* (a string, e.g. \"trait\" / \"genomes:count\"). "
                        f"A grown result object is conditioning — pass it to the target level's run."
                    )
                sources.append(m.source)
    if not sources:
        raise ValueError(
            "a joint model needs the driver to drive something: give birth (or death) a "
            "mod.DrivenBy(...). With neither driven, grow the two levels as independent runs instead."
        )
    # the driver spec must match the sources
    if trait is not None:
        if not isinstance(trait, DiscreteTrait):
            raise TypeError(
                "trait= must be traits.discrete(states=[...], switch=...) — a discrete process spec. "
                "Continuous trait→speciation (QuaSSE) is deferred."
            )
        bad = sorted({s for s in sources if s != "trait"})
        if bad:
            raise ValueError(
                f'with trait=, drive from the live trait — mod.DrivenBy("trait", ...); got source(s) '
                f"{bad}. (A filename source is conditioning, not a joint run.)"
            )
    else:
        if not isinstance(genome, UnorderedGenome):
            raise TypeError("genome= must be genomes.unordered(...) — an unordered-genome process spec.")
        for s in sources:
            if s == _GENOME_COUNT:
                continue
            if s.startswith("genomes:"):
                name = s.split(":", 1)[1]
                if name not in genome.families:
                    raise ValueError(
                        f'DrivenBy("{s}", ...) names family {name!r}, but genomes.unordered was not '
                        f"declared with it — add families=[…, {name!r}]."
                    )
                continue
            raise ValueError(
                f'with genome=, drive from gene content — "genomes:count" or "genomes:<family>"; '
                f"got {s!r}."
            )
    if (n_extant is None) == (total_time is None):
        raise ValueError("give exactly one of n_extant or total_time")
    if n_extant is not None and (isinstance(n_extant, bool) or not isinstance(n_extant, int) or n_extant < 1):
        raise ValueError(f"n_extant must be a positive integer, got {n_extant!r}")
    if total_time is not None and (not isinstance(total_time, (int, float))
                                   or not math.isfinite(total_time) or total_time <= 0):
        raise ValueError(f"total_time must be a positive finite number, got {total_time!r}")

    rng = np.random.default_rng(seed)
    unique_sources = sorted(set(sources))

    def grow_once(target_n, tt) -> tuple[Tree, JointResult]:
        if trait is not None:
            tree, se, nv, te = _grow_joint(rng, birth_rate, death_rate, trait, target_n, tt)
            te.sort(key=lambda c: c.time)
            result = JointResult(SpeciesResult(tree, se, seed, []), seed,
                                 trait=TraitsResult(tree, nv, te, seed, kind="discrete"))
        else:
            tree, se, go, ge, fn = _grow_joint_genome(
                rng, birth_rate, death_rate, genome, unique_sources, target_n, tt)
            result = JointResult(SpeciesResult(tree, se, seed, []), seed,
                                 genome=GenomesResult(tree, go, ge, seed, fn))
        return tree, result

    if total_time is not None:
        return grow_once(None, total_time)[1]

    for _ in range(_MAX_ATTEMPTS):
        tree, result = grow_once(n_extant, None)
        if sum(1 for nd in tree.nodes.values() if nd.fate == "extant") == n_extant:
            return result
    raise RuntimeError(
        f"could not grow a tree to {n_extant} extant lineages in {_MAX_ATTEMPTS} attempts; "
        "birth must comfortably exceed death for large n_extant"
    )


__all__ = ["simulate_joint", "JointResult"]
