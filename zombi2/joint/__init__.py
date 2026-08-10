"""Joint models — the driver co-evolves with what it drives (SPEC §2–4).

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
them and the race is **exact** — no thinning. (A continuously-diffusing driver — QuaSSE — is not
available: it makes the rate vary continuously, which needs thinning.)

The mechanism is the same ``scaled_by`` as conditioning; only the ``driver`` differs — here a
**live level name** (``"trait"``, ``"genomes:count"``, ``"genomes:<family>"``) rather than a filename.
Driving *both* birth and death recovers full state-dependent diversification (BiSSE's λ and μ)."""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass


from .._runtime.draw import weighted_index as _weighted_index
from .._runtime.summary import write_summary
from ..genomes import GeneEdge, GeneCopy, FamilyGenomesResult, FamilyGenome
from ..genomes.family import _duplicate, _lose_at, _originate, _pick_copy  # engine internals
from ..params.mapping import check_not_a_kernel
from ..rng import stream
from ..params.driver import OnTime, OnTotalDiversity
from ..params.evaluate import DRAWN, INHERITED, describe, is_implemented
from ..params.connection import Driven

from ..params.parameter import as_rate
from ..params.scope import PerLineage
from ..species import Event as SpeciesEvent, SpeciesResult
from ..tree import Node, Tree
from ..traits import Change, DiscreteTrait, TraitsResult

#: The rate grammar a joint run supports on ``birth`` / ``death`` (SPEC §5). Declared, like every
#: other level, so the gate below cannot fall behind what the engine threads: the loop passes ``time``
#: and ``diversity`` into every rate and steps its Gillespie at each ``next_change``, so the two
#: covariates are as real here as at the species level, and ``scaled_by`` is what makes the run joint
#: at all. What is missing is missing on purpose — see the rejections in `simulate_joint()`.
IMPLEMENTED_MODIFIERS = (OnTime, OnTotalDiversity, Driven)

_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_extant
_GENOME_COUNT = "genomes:count"  # the live gene-content driver source for a count → Curve/Scalar


@dataclass
class JointResult:
    """What `simulate_joint()` returns — **both** grown levels of a joint run. ``species`` is the
    grown tree (a `SpeciesResult`: ``complete_tree``, ``extant_tree``, the
    speciation/extinction ``events``); the **driver** level that grew with it is either ``trait`` (a
    `TraitsResult`, for a trait→speciation run) or ``genome`` (a
    `FamilyGenomesResult`, for a gene-content→speciation run) — exactly one is set.
    The tree is an output, grown by the driver it carries, so the levels share one ``complete_tree``."""

    species: SpeciesResult
    seed: int | None
    trait: TraitsResult | None = None
    genome: FamilyGenomesResult | None = None

    def __repr__(self) -> str:
        driver = "trait" if self.trait is not None else "genome"
        return (f"JointResult({self.n_extant} extant tips grown with a {driver}, "
                f"seed={self.seed})")

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

    def summary(self) -> dict:
        """What this run produced, as a plain dict — the payload of ``joint_summary.json``.

        A joint run grew two levels at once, so this holds both of their summaries under one roof
        rather than inventing a third vocabulary: ``species`` is the tree that came out, and exactly
        one of ``trait`` / ``genome`` is the driver that shaped it. The tree is an *output* here, which
        is the whole point of the command, so its realised birth and death rates are the numbers worth
        reading — they are what the driver did."""
        out = {"level": "joint", "seed": self.seed,
               "driver": "trait" if self.trait is not None else "genome",
               "species": self.species.summary()}
        if self.trait is not None:
            out["trait"] = self.trait.summary()
        if self.genome is not None:
            out["genome"] = self.genome.summary()
        return out

    def write(self, directory, *, flat: bool = False) -> None:
        """Write both levels to ``directory``: the species files (``species_complete.nwk`` /
        ``species_extant.nwk`` / ``species_events.tsv``) and the driver level's — for a trait,
        ``trait_values.tsv`` / ``trait_events.tsv`` / ``trait_tree.nwk``; for a genome,
        ``genome_events.tsv`` / ``profiles.tsv``. ``flat`` is passed to the driver level, which is
        the only one of the two with a many-files-per-run output."""
        write_summary(pathlib.Path(directory) / "joint_summary.json", self.summary())
        self.species.write(directory, outputs=("complete", "extant", "events"))
        if self.trait is not None:
            self.trait.write(directory, outputs=("values", "events", "tree"))
        if self.genome is not None:
            self.genome.write(directory, outputs=("events", "profiles"), flat=flat)


def _grow_joint(rng, birth_rate, death_rate, trait: DiscreteTrait, n_extant, total_time):
    """Grow a forward birth-death tree whose birth/death read a discrete trait that evolves on it.
    Returns ``(tree, species_events, node_values, trait_events)`` — the complete tree, the
    speciation/extinction log, the trait state at every node, and the trait's switch log."""
    states, Q, start_i, shift = trait._resolve(rng)
    k_states = len(states)
    out_rate = [float(-Q[s, s]) for s in range(k_states)]  # the trait's total switch-out rate per state

    # birth/death are driven by the trait; the trait's declared states are known up front, so the check
    # is exhaustive — every mapping key must be one of them. A key outside the alphabet is a state that
    # can never occur (a typo whose factor would silently never apply), and a mapping matching none of
    # them would be a silently undriven run; both are refused here rather than run as if driven.
    from ..params.conditioned import check_mapping_fires
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        for m in rate.modifiers:
            if isinstance(m, Driven):
                check_mapping_fires(m.mapping, states, driver_label=f"{label} (trait)", exhaustive=True)

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
    species_events: list[SpeciesEvent] = []
    # the initial state at t=0, exactly as the standalone traits engine seeds its log: tree + this row
    # + the switches give the trait on every lineage, which is what lets the log be read back as a
    # driver. Without it a joint run writes a trait_events.tsv no conditioned run can read.
    trait_events: list[Change] = [Change(0.0, "initial", root, None, states[start_i])]
    end_state: dict[int, int] = {}  # node id → its trait state index when it ended (→ node_values)

    while alive:
        n = len(alive)
        ctx = {"diversity": n, "time": t}
        # per-lineage rates: birth/death read the lineage's trait state (scaled_by("trait", …)); the
        # trait switch rate is the CTMC out-rate for that state (the trait's own dynamics, undriven).
        wb = [birth_rate.effective(lineages=1, drivers={"trait": states[st[k]]}, **ctx) for k in range(n)]
        wd = [death_rate.effective(lineages=1, drivers={"trait": states[st[k]]}, **ctx) for k in range(n)]
        ws = [out_rate[st[k]] for k in range(n)]
        total_b, total_d, total_s = sum(wb), sum(wd), sum(ws)
        total = total_b + total_d + total_s

        # the trait switch rate is constant between events; only a skyline (changing_at) on birth/death or
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
                    species_events.append(SpeciesEvent(t, "speciation", node_id, (c1, c2)))
                elif r < total_b + total_d:  # extinction
                    i = _weighted_index(rng, wd, total_d)
                    node_id, cur = alive[i], st[i]
                    alive[i] = alive[-1]; alive.pop()
                    st[i] = st[-1]; st.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "extinct"
                    end_state[node_id] = cur
                    species_events.append(SpeciesEvent(t, "extinction", node_id))
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


def _grow_joint_genome(rng, birth_rate, death_rate, spec: FamilyGenome, driver_names, n_extant, total_time):
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
    while alive:
        nl = len(alive)
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
    return Tree(nodes, root), species_events, genomes_out, genome_events, named


def simulate_joint(*, birth, death=0.0, trait=None, genome=None, n_extant=None, total_time=None,
                   seed=None) -> JointResult:
    """Grow a tree **and** the driver that drives its speciation, in one run (SPEC §2–4).

    ``birth`` and ``death`` are rate specs (per lineage). Make either read the driver with
    ``.scaled_by(driver, mapping)`` — a **live level name** (not a filename) is what makes this
    *joint* rather than conditioned. Give **exactly one** driver:

    - ``trait = traits.discrete(...)`` — a discrete trait drives speciation (BiSSE / MuSSE), read as
      ``.scaled_by("trait", {"small": 1.0, "large": 2.0})``. Driving both birth and death gives
      state-dependent λ *and* μ.
    - ``genome = genomes.family(...)`` — **gene content** drives speciation (``P(Species, Genomes)``),
      read as the total gene count ``.scaled_by("genomes:count", curve)`` or the presence of a named
      family ``.scaled_by("genomes:toxin", {"present": 2.0, "absent": 1.0})`` (declare it with
      ``family_names=["toxin"]``).

    ::

        joint.simulate_joint(
            birth  = PerLineage(1.0).scaled_by("genomes:toxin", {"present": 3.0, "absent": 1.0}),
            genome = genomes.family(origination=0.2, loss=0.1, family_names=["toxin"]),
            n_extant = 100, seed = 1)

    The driver is an **unexecuted** process spec, grown with the tree. Stop at exactly ``n_extant``
    living lineages (conditioned on survival — a birth-death tree can die out, so it restarts,
    advancing the same generator) **or** at ``total_time`` — give exactly one. Returns a
    `JointResult` carrying the grown tree and the driver level (``.trait`` or ``.genome``).
    Deterministic given ``seed``. Continuous trait→speciation (QuaSSE), clade drift (an inherited value)
    combined with driving, and gene transfer in a joint run are not available.
    """
    birth_rate = as_rate(birth, default_scope=PerLineage)
    death_rate = as_rate(death, default_scope=PerLineage)
    if (trait is None) == (genome is None):
        raise TypeError(
            "give exactly one driver: trait=traits.discrete(...) OR genome=genomes.family(...)."
        )
    # collect the Driven driver names on birth/death (a joint model's diversification must be per lineage)
    driver_names: list[str] = []
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        # `Rate.scope` holds the scope **class**, so this is an identity test rather than an
        # isinstance: a scope instance never exists (SPEC §5).
        if rate.scope is not PerLineage:
            assert rate.scope is not None      # `as_rate` above fills in the default scope
            raise ValueError(
                f"{label} has a {rate.scope.__name__} scope, but a joint diversification rate is "
                f"per lineage — write PerLineage(...) (the default, so a bare number is enough)."
            )
        for m in rate.modifiers:
            if m.reads == (DRAWN, "families"):
                # not a missing feature: there is nothing here for it to mean (the species level
                # says the same thing about the same modifier, for the same reason)
                raise ValueError(
                    f"{label} carries a value drawn among families, but a diversification rate has "
                    f"no gene families — varying_among('families', ...) belongs on a genomes rate. "
                    f"To make speciation depend on gene content, drive it: "
                    f"birth = PerLineage(1.0).scaled_by(\"genomes:count\", ...)."
                )
            if m.reads == (INHERITED, "lineages"):
                raise ValueError(
                    f"{label} carries a value inherited among lineages (clade drift); drift and a "
                    f"driven rate are not available together — use one or the other."
                )
            if m.reads == (DRAWN, "lineages"):
                # The species level takes this; the joint engine does not thread it, so accepting it
                # here would run the model without the rate variation the user asked for — and the
                # same `--birth` expression working on `zombi2 species` makes that a trap rather than
                # merely a gap (SPEC §5: reject, never silently ignore).
                raise ValueError(
                    f"{label} carries a value drawn among lineages (independent per-lineage rates); "
                    f"per-lineage rate variation and a driven rate are not available together in a "
                    f"joint run — use one or the other. On its own, "
                    f"varying_among('lineages', ...) works at the species level."
                )
            if not is_implemented(m, IMPLEMENTED_MODIFIERS, "joint"):
                # the backstop: anything this engine does not thread would come back as its default
                # factor of 1.0, which is a run quietly not the model that was asked for (SPEC §5).
                # Declared rather than enumerated here, so a modifier added later cannot slip through.
                raise ValueError(
                    f"{label} carries {describe(m)}, which a joint run does not support. It "
                    f"takes changing_at (skyline), scaled_by(TotalDiversity(cap=...)) "
                    f"(diversity-dependent) and scaled_by (the driver that makes the run joint)."
                )
            if isinstance(m, Driven):
                check_not_a_kernel(m.mapping, label=label)
                if not isinstance(m.driver, str):
                    raise TypeError(
                        f"{label} is driven by a {type(m.driver).__name__} object, but a joint model "
                        f"drives from a live level *name* (a string, e.g. \"trait\" / \"genomes:count\"). "
                        + ("A clade is read off a finished tree, and a joint run grows the tree as it "
                           "goes, so there is no clade to read yet."
                           if type(m.driver).__name__ == "Clade" else
                           "A grown result object is conditioning — pass it to the driven level's run.")
                    )
                driver_names.append(m.driver)
    if not driver_names:
        raise ValueError(
            "a joint model needs the driver to drive something: give birth (or death) a "
            "scaled_by(...). With neither driven, grow the two levels as independent runs instead."
        )
    # the driver spec must match the driver names
    if trait is not None:
        if not isinstance(trait, DiscreteTrait):
            raise TypeError(
                "trait= must be traits.discrete(states=[...], switch=...) — a discrete process spec. "
                "Continuous trait→speciation (QuaSSE) is not available: a continuously varying "
                "rate needs thinning, which this exact race does not do."
            )
        bad = sorted({s for s in driver_names if s != "trait"})
        if bad:
            raise ValueError(
                f'with trait=, drive from the live trait — scaled_by("trait", ...); got driver(s) '
                f"{bad}. (A filename driver is conditioning, not a joint run.)"
            )
    else:
        if not isinstance(genome, FamilyGenome):
            raise TypeError("genome= must be genomes.family(...) — a family-genome process spec.")
        for s in driver_names:
            if s == _GENOME_COUNT:
                continue
            if s.startswith("genomes:"):
                name = s.split(":", 1)[1]
                if name not in genome.family_names:
                    raise ValueError(
                        f'scaled_by("{s}", ...) names family {name!r}, but genomes.family was not '
                        f"declared with it — add family_names=[…, {name!r}]."
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

    rng, seed = stream("joint", seed)       # own stream, and a drawn seed if none was given
    unique_driver_names = sorted(set(driver_names))

    def grow_once(target_n, tt) -> tuple[Tree, JointResult]:
        if trait is not None:
            tree, se, nv, te = _grow_joint(rng, birth_rate, death_rate, trait, target_n, tt)
            te.sort(key=lambda c: c.time)
            result = JointResult(SpeciesResult(tree, se, seed, []), seed,
                                 trait=TraitsResult(tree, nv, te, seed, kind="discrete"))
        else:
            tree, se, go, ge, fn = _grow_joint_genome(
                rng, birth_rate, death_rate, genome, unique_driver_names, target_n, tt)
            result = JointResult(SpeciesResult(tree, se, seed, []), seed,
                                 genome=FamilyGenomesResult(tree, go, ge, seed, fn, {}))
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
