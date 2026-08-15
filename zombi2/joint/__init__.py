"""Joint models — one run simulates two levels at once (SPEC §2–4).

A run is **joint** when neither level can be finished before the other starts, so one run has to
produce both. That is the whole of it, and it says nothing on its own about the species tree: two
levels can drive each other on a tree handed to the run. In the two models built here the species
tree **is** one of the two being simulated, so it comes out of the run rather than going into it:

- a **discrete trait** drives speciation (BiSSE / MuSSE), ``P(Species, Traits)`` — birth/death read the
  trait state on each lineage while the trait evolves by its own Mk process on the growing tree;
- **gene content** drives speciation, ``P(Species, Genomes)`` — birth/death read a summary of each
  lineage's live genome (its total gene count, or the presence of a named family) while the genome
  evolves by duplication/loss/origination on the growing tree;
- a **continuous trait** drives speciation (QuaSSE) — birth/death read a diffusing value on each
  lineage while it diffuses on the growing tree.

One Gillespie races the event classes over the living lineages at once: **speciation** and
**extinction** (per lineage, driver-read), plus the driver's own events — a **trait switch** (the CTMC
out-rate) or a genome **duplication/loss/origination**. A driver event changes a lineage's state
without touching the topology; a speciation hands the parent's driver state (its trait, its genome) to
both daughters. Because a discrete driver only changes at events, the rate is piecewise-constant
between them and the race is **exact** — no thinning. A diffusion is the exception: it moves at every
instant, so that run **slices**, holding the value fixed across a ``step`` the driven rate declares.

The mechanism is the same ``scaled_by`` as conditioning; only the ``driver`` differs — here a
**live level name** (``"trait"``, ``"genomes:count"``, ``"genomes:<family>"``) rather than a filename.
Driving *both* birth and death recovers full state-dependent diversification (BiSSE's λ and μ)."""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass


from .._runtime.draw import weighted_index as _weighted_index
from .._runtime.summary import write_summary
from ..genomes import GeneEdge, GeneCopy, FamilyGenomesResult, FamilyGenome, GeneFamily
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
from ..traits import Change, ContinuousTrait, DiscreteTrait, TraitsResult

#: The rate grammar a joint run supports on ``birth`` / ``death`` (SPEC §5). Declared, like every
#: other level, so the gate below cannot fall behind what the engine threads: the loop passes ``time``
#: and ``diversity`` into every rate and steps its Gillespie at each ``next_change``, so the two
#: covariates are as real here as at the species level, and ``scaled_by`` is what makes the run joint
#: at all. What is missing is missing on purpose — see the rejections in `_simulate_joint()`.
IMPLEMENTED_MODIFIERS = (OnTime, OnTotalDiversity, Driven)

#: `JointResult.write`'s vocabulary. The tokens are the two **levels** and the run's own summary, not
#: their files: a joint run's whole claim is that it writes each level exactly as that level's own
#: command does, so each is written with its own default and there is nothing here to restate.
_WRITE_OUTPUTS = ("summary", "species", "driver")

_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_extant
_GENOME_COUNT = "genomes:count"  # the live gene-content driver source for a count → Curve/Scalar
#: What a trait driver is called when the run holds one and it was given no name. Naming is
#: what lets a run hold two, so the bare form stays for the common case of one.
_ONE_TRAIT = "trait"


@dataclass
class JointResult:
    """What `simulate()` returns — **both** simulated levels of a joint run. ``species`` is the
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

    def write(self, directory, outputs=_WRITE_OUTPUTS, *, flat: bool = False) -> None:
        """Write both levels to ``directory`` (created if needed), each exactly as its own command
        writes it: ``"species"`` → the `SpeciesResult` files (``species_complete.nwk`` /
        ``species_extant.nwk`` / ``species_events.tsv`` / ``species_fates.tsv`` /
        ``species_summary.json``); ``"driver"`` → the level that grew with it, a trait's
        ``trait_values.tsv`` / ``trait_events.tsv`` / ``trait_tree.nwk`` / ``trait_summary.json`` or
        a genome's ``genome_events.tsv`` / ``profiles.tsv`` / ``genomes.tsv`` /
        ``initial_genome.tsv`` / ``gene_trees/`` / ``genome_summary.json``; ``"summary"`` →
        ``joint_summary.json``, the one file that is the joint run's own.

        The tokens are the two **levels**, not their files, because each is written with that level's
        own default — which is what makes a joint run's directory the two runs it stands in for. Pick
        files *within* a level through the level itself: ``result.species.write(d, outputs=…)``,
        ``result.trait.write(d, outputs=…)``. ``flat`` is passed to the driver level, the only one of
        the two with a many-files-per-run output.

        Both levels land in the one directory named here; ``zombi2 joint`` groups them under
        ``species/`` and ``traits/`` / ``genomes/`` instead, and writes the same files."""
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "summary" in outputs:
            write_summary(d / "joint_summary.json", self.summary())
        if "species" in outputs:
            self.species.write(d)
        if "driver" in outputs:
            # two independent tests, not an if/else: a result carrying neither writes neither, rather
            # than reaching for `.write` on None
            if self.trait is not None:
                self.trait.write(d)
            if self.genome is not None:
                self.genome.write(d, flat=flat)


#: What ``max_lineages`` says when it fires. The species engine has had this guard since it grew a
#: tree conditioned on time; a joint run needs it more, not less, because its birth rate can read a
#: driver that its own growth feeds — gene content accumulates, birth rises, more lineages accumulate
#: more gene content — so a rate that looks calm on paper can have no realistic end. It RAISES rather
#: than stopping early, for the reason the species engine gives: a tree cut off at a size is no longer
#: a sample from the process asked for.
def _runaway(ceiling: int, t: float) -> RuntimeError:
    return RuntimeError(
        f"the tree passed {ceiling} standing lineages at time {t:.3g} and is still growing — the "
        f"driven birth rate has no realistic end. Lower the rates, shorten total_time, flatten the "
        f"mapping the driver is read through, or raise max_lineages if the size is what you want "
        f"(max_lineages=None removes the guard).")


def _grow_joint(rng, birth_rate, death_rate, trait: DiscreteTrait, n_extant, total_time,
                max_lineages=None, keys=(_ONE_TRAIT,)):
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

    ceiling = None if max_lineages is None else max(max_lineages, n_extant or 0)
    while alive:
        n = len(alive)
        if ceiling is not None and n > ceiling:
            raise _runaway(ceiling, t)
        ctx = {"diversity": n, "time": t}
        # per-lineage rates: birth/death read the lineage's trait state (scaled_by("trait", …)); the
        # trait switch rate is the CTMC out-rate for that state (the trait's own dynamics, undriven).
        vals = [{key: states[st[k]] for key in keys} for k in range(n)]
        wb = [birth_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
        wd = [death_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
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


#: How many slices a sliced run will walk before it decides nothing is going to happen. A sliced run
#: cannot use the other engines' "nothing scheduled, so stop" test: a rate of zero now says nothing
#: about the rate one slice later, because the driver is still moving. So the walk needs an end of
#: its own, and this is it. It is large enough that no run reaching ``n_extant`` meets it and small
#: enough to answer in seconds when the rates really are dead.
_MAX_SLICES = 2_000_000


def _slice_step(birth_rate, death_rate, trait_keys) -> float:
    """The slice a diffusing driver is held fixed across, taken from the rates that read it.

    ``step`` rides on the connection rather than on the run because it is a property of *this*
    reading: a steep response curve needs a finer step than a flat one, and two participants reading
    the same driver can legitimately disagree. What cannot happen is two rates reading **one live
    trait** at two resolutions — there is one trajectory here, and it either moves at a boundary or
    it does not — so that is refused rather than resolved to the finer of the two.

    There is no default. A step is the size of the approximation being made, and any number this
    function could invent would be a claim about a timescale only the model knows."""
    steps, missing = set(), []
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        for m in rate.modifiers:
            if isinstance(m, Driven) and m.driver in trait_keys:
                (missing.append(label) if m.step is None else steps.add(m.step))
    if missing:
        raise ValueError(
            f"{' and '.join(sorted(set(missing)))} reads a continuously diffusing trait, so it needs "
            f"a step= — the stretch of time the driver is held fixed across. A diffusion moves at "
            f"every instant, so there is no interval where this rate holds still on its own and "
            f"nothing exact to draw against; the run slices instead. Write it on the connection, "
            f'scaled_by("trait", Curve(f), step=0.05), and pick it so the trait moves little within '
            f"one slice — then halve it, rerun the same seed, and see whether the answer moves.")
    if len(steps) > 1:
        raise ValueError(
            f"birth and death read the same live trait at two resolutions, step={sorted(steps)}. "
            f"There is one trajectory in a joint run, and it either moves at a given boundary or it "
            f"does not, so the two readings have to agree on step.")
    return float(next(iter(steps)))


def _grow_joint_continuous(rng, birth_rate, death_rate, trait, step: float, n_extant, total_time,
                           max_lineages=None, keys=(_ONE_TRAIT,)):
    """Grow a forward birth-death tree whose birth/death read a **continuously diffusing** trait —
    QuaSSE. Returns ``(tree, species_events, node_values, trait_events)``.

    Every other joint model here races exactly, because its driver only changes at events and an
    event ends the step. A diffusion changes at every instant, so there is no interval over which
    the birth rate holds still and nothing for a Gillespie step to be drawn against.

    **This one slices.** Time is cut into steps of ``step``; inside a step every lineage's value is
    held where it was, so the rates *are* constant and the race inside the slice is the ordinary
    one. At the boundary each lineage's value is advanced by the exact transition law of its own
    diffusion — ``Normal(0, ∫σ²)`` for Brownian motion, the pull-weighted form under
    Ornstein–Uhlenbeck — over the time since that lineage last moved, which is the slice for an old
    lineage and the remainder of it for one born mid-slice. So the **trait** is exact and only its
    coupling to speciation is approximated: a lineage speciates at the rate its value had at the top
    of the slice rather than the rate it has at that instant.

    The error is first-order in ``step`` and one-sided in a way worth knowing: the value is carried
    forward, never interpolated, because a growing tree has no future to read. The conditioned
    version of the same model — a trait grown first and handed over — interpolates between two known
    ends, so for the same ``step`` it is the more accurate of the two.
    """
    from ..traits.continuous import _accrued_variance
    from ..params.mapping import Table

    rate, start, theta, alpha, jump_sd = trait._resolve()
    is_ou = alpha is not None
    # What to thread the value under. A `Driven` looks its value up by `m.key`, which is the driver
    # name alone only while there is no step — with one it is `(name, step)`, because the same driver
    # read at two resolutions is two trajectories. Threading the bare name here left every factor at
    # 1.0 and the run silently undriven, which is the failure `Driven.factor`'s inert default makes
    # quiet: the tree came out identical whatever curve was written.
    lookup = set()
    for label, r in (("birth", birth_rate), ("death", death_rate)):
        for m in r.modifiers:
            if not isinstance(m, Driven) or m.driver not in keys:
                continue
            lookup.add(m.key)
            # a Curve or a Scalar is defined on the whole line, so there is no alphabet to check it
            # against; a table is the one mapping a diffusing value can never match a key of
            if isinstance(m.mapping, Table):
                raise ValueError(
                    f"{label} reads a continuous trait through a {{state: factor}} table, and a "
                    f"diffusing value is a number that never equals a state name. Map it with a "
                    f"Curve (value → factor) or a Scalar (a log-link): "
                    f'scaled_by("trait", Curve(lambda x: math.exp(0.5 * x)), step=0.05).')

    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent, t):
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    root = new_node(None, 0.0)
    alive = [root]              # living lineage ids
    xs = [float(start)]         # each lineage's trait value, kept in lock-step with `alive`
    since = [0.0]               # …and when it last diffused, so a mid-slice birth is not over-diffused
    t = 0.0
    slice_end = step
    species_events: list[SpeciesEvent] = []
    # a diffusion cannot be rebuilt from its events, so the log carries what the trait level's own log
    # carries: the value at t=0 and each jump at a split. The values themselves ride in node_values.
    trait_events: list[Change] = [Change(0.0, "initial", root, None, float(start))]
    end_value: dict[int, float] = {}

    ceiling = None if max_lineages is None else max(max_lineages, n_extant or 0)
    slices = 0
    while alive:
        n = len(alive)
        if ceiling is not None and n > ceiling:
            raise _runaway(ceiling, t)
        ctx = {"diversity": n, "time": t}
        vals = [{key: xs[k] for key in lookup} for k in range(n)]
        wb = [birth_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
        wd = [death_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
        total_b, total_d = sum(wb), sum(wd)
        total = total_b + total_d

        # the slice boundary is a horizon like any other: the rates change there, so the step stops
        next_change = min(birth_rate.next_change(t), death_rate.next_change(t), slice_end)
        horizon = next_change if total_time is None else min(next_change, total_time)

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:
                t = t_ev
                if n == n_extant:
                    break
                if float(rng.random()) * total < total_b:      # speciation
                    i = _weighted_index(rng, wb, total_b)
                    node_id, x = alive[i], xs[i]
                    alive[i] = alive[-1]; alive.pop()
                    xs[i] = xs[-1]; xs.pop()
                    since[i] = since[-1]; since.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "speciation"
                    end_value[node_id] = x
                    c1, c2 = new_node(node_id, t), new_node(node_id, t)
                    node.children = (c1, c2)
                    for c in (c1, c2):
                        d = x
                        if jump_sd > 0.0:      # the punctuational jump, drawn per daughter
                            d = x + float(rng.normal(0.0, jump_sd))
                            trait_events.append(Change(t, "on_speciation", c, x, d))
                        alive.append(c); xs.append(d); since.append(t)
                    species_events.append(SpeciesEvent(t, "speciation", node_id, (c1, c2)))
                else:                                          # extinction
                    i = _weighted_index(rng, wd, total_d)
                    node_id, x = alive[i], xs[i]
                    alive[i] = alive[-1]; alive.pop()
                    xs[i] = xs[-1]; xs.pop()
                    since[i] = since[-1]; since.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "extinct"
                    end_value[node_id] = x
                    species_events.append(SpeciesEvent(t, "extinction", node_id))
                continue

        if total_time is not None and horizon == total_time:
            t = total_time
            break
        if horizon == slice_end:
            # the slice ends: let every living lineage diffuse over the time since it last moved
            for k in range(len(alive)):
                xs[k] = _advance(rng, rate, xs[k], since[k], slice_end, theta, alpha, is_ou,
                                 _accrued_variance)
                since[k] = slice_end
            t = slice_end
            slice_end += step
            slices += 1
            if slices > _MAX_SLICES:
                raise RuntimeError(
                    f"the run walked {_MAX_SLICES} slices of step={step:g} — to time {t:.6g} — "
                    f"without reaching n_extant={n_extant}. Either the driven birth rate is "
                    f"effectively zero over the values the trait reaches, or step is far finer than "
                    f"the timescale the tree grows on.")
            continue
        t = horizon      # a skyline breakpoint on birth/death: advance and re-read them

    for k, node_id in enumerate(alive):
        nodes[node_id].end_time = t
        # whoever is still alive reached the present, and did so with the last part-slice unrun
        end_value[node_id] = _advance(rng, rate, xs[k], since[k], t, theta, alpha, is_ou,
                                      _accrued_variance)
        nodes[node_id].fate = "extant"

    return Tree(nodes, root), species_events, dict(end_value), trait_events


def _advance(rng, rate, x, t0, t1, theta, alpha, is_ou, accrued):
    """One lineage's diffusion from ``t0`` to ``t1`` — the same transition law `simulate_continuous`
    walks a branch with, over a slice instead of a branch."""
    if t1 <= t0:
        return x
    if is_ou:
        e = math.exp(-alpha * (t1 - t0))
        mean = theta + (x - theta) * e
        var = accrued(rate, t0, t1, pull=alpha)
    else:
        mean, var = x, accrued(rate, t0, t1)
    return mean + (float(rng.normal(0.0, math.sqrt(var))) if var > 0.0 else 0.0)


def _grow_joint_genome(rng, birth_rate, death_rate, spec: FamilyGenome, driver_names, n_extant,
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
    return Tree(nodes, root), species_events, genomes_out, genome_events, named


def _simulate_joint(*, birth, death=0.0, trait=None, genome=None, n_extant=None,
                    total_time=None, seed=None, max_lineages=100_000) -> JointResult:
    """Grow a tree **and** the driver that drives its speciation, in one run (SPEC §2–4).

    ``birth`` and ``death`` are rate specs (per lineage). Make either read the driver with
    ``.scaled_by(driver, mapping)``. Driving speciation has its own function because the driver
    cannot be grown first: it would have to grow on the tree it drives, so the two grow together.
    The driver is a **live level name** (``"trait"``, ``"genomes:count"``); it has to agree with the
    spec given below, and a filename is refused, because a driver read from a file is conditioning.
    Give **exactly one** driver:

    - ``trait = traits.discrete(...)`` — a discrete trait drives speciation (BiSSE / MuSSE), read as
      ``.scaled_by("trait", {"small": 1.0, "large": 2.0})``. Driving both birth and death gives
      state-dependent λ *and* μ.
    - ``trait = traits.continuous(...)`` — a **diffusing** trait drives speciation (QuaSSE), read
      through a `~zombi2.params.Curve` or a `~zombi2.params.Scalar` rather than a table, and with a
      ``step=`` on the connection: a diffusion moves at every instant, so this one run slices.
    - ``genome = genomes.genome(...)`` — **gene content** drives speciation (``P(Species, Genomes)``),
      read as the total gene count ``.scaled_by("genomes:count", curve)`` or the presence of a named
      family ``.scaled_by("genomes:toxin", {"present": 2.0, "absent": 1.0})`` (declare it with
      ``families=[family("toxin")]``).

    The engine behind `simulate()`, which is the way to call it::

        joint.simulate(
            species.birth_death(
                birth = PerLineage(1.0).scaled_by("genomes:toxin", {"present": 3.0, "absent": 1.0}),
                n_extant = 100),
            genomes.genome(origination=0.2, loss=0.1, families=[family("toxin")]), seed = 1)

    ``max_lineages`` (default 100000) stops a run that has no realistic end. A joint birth rate reads
    a driver the run itself grows, so it can feed itself — gene content accumulates, birth rises, and
    more lineages accumulate more gene content — and a rate that looks calm on paper then grows
    without bound. It **raises** rather than truncating, for the reason the species engine gives: a
    tree cut off at a size is no longer a sample from the process asked for. ``max_lineages=None``
    removes the guard.

    The driver is an **unexecuted** process spec, grown with the tree. Stop at exactly ``n_extant``
    living lineages (conditioned on survival — a birth-death tree can die out, so it restarts,
    advancing the same generator) **or** at ``total_time`` — give exactly one. Returns a
    `JointResult` carrying the grown tree and the driver level (``.trait`` or ``.genome``).
    Deterministic given ``seed``. Clade drift (an inherited value) combined with driving, and gene
    transfer in a joint run, are not available.
    """
    birth_rate = as_rate(birth, default_scope=PerLineage)
    death_rate = as_rate(death, default_scope=PerLineage)
    if (trait is None) == (genome is None):
        raise TypeError(
            "give exactly one driver: trait=traits.discrete(...) OR genome=genomes.genome(...)."
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
        if not isinstance(trait, (DiscreteTrait, ContinuousTrait)):
            raise TypeError(
                "trait= must be traits.discrete(states=[...], switch=...) or "
                "traits.continuous(rate=...) — a trait process spec.")
        # a run holds one trait, so the bare "trait" always names it; a name additionally lets a rate
        # say which, which is what a run holding two will need
        trait_keys = {_ONE_TRAIT} | ({f"traits:{trait.name}"} if trait.name else set())
        bad = sorted({s for s in driver_names if s not in trait_keys})
        if bad:
            named = " or ".join(f'"{k}"' for k in sorted(trait_keys))
            raise ValueError(
                f"with a trait participant, drive from that trait — scaled_by({named}, ...); got "
                f"driver(s) {bad}. (A filename driver is conditioning, not a joint run.)"
            )
        if isinstance(trait, ContinuousTrait):
            step = _slice_step(birth_rate, death_rate, trait_keys)
    else:
        if not isinstance(genome, FamilyGenome):
            raise TypeError("genome= must be genomes.genome(...) — a family-genome process spec.")
        if genome.transfer:
            raise ValueError(
                "transfer is not available while the tree is being simulated: a transfer needs the "
                "set of lineages alive at that instant, and on a growing tree that set is still "
                "forming. It works on a tree handed to the run — genomes with traits.")
        for s in driver_names:
            if s == _GENOME_COUNT:
                continue
            if s.startswith("genomes:"):
                name = s.split(":", 1)[1]
                if name not in genome.family_names:
                    raise ValueError(
                        f'scaled_by("{s}", ...) names family {name!r}, but genomes.family was not '
                        f"declared with it — add families=[…, family({name!r})]."
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
        if isinstance(trait, ContinuousTrait):
            tree, se, nv, te = _grow_joint_continuous(rng, birth_rate, death_rate, trait, step,
                                                      target_n, tt, max_lineages,
                                                      tuple(sorted(trait_keys)))
            te.sort(key=lambda c: c.time)
            result = JointResult(SpeciesResult(tree, se, seed, []), seed,
                                 trait=TraitsResult(tree, nv, te, seed))
        elif trait is not None:
            tree, se, nv, te = _grow_joint(rng, birth_rate, death_rate, trait, target_n, tt,
                                           max_lineages, tuple(sorted(trait_keys)))
            te.sort(key=lambda c: c.time)
            result = JointResult(SpeciesResult(tree, se, seed, []), seed,
                                 trait=TraitsResult(tree, nv, te, seed, kind="discrete"))
        else:
            tree, se, go, ge, fn = _grow_joint_genome(
                rng, birth_rate, death_rate, genome, unique_driver_names, target_n, tt,
                max_lineages)
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


def _classify(participants):
    """Sort the things handed to `simulate()` by which level each belongs to.

    A participant is a **process spec** — `~zombi2.species.BirthDeath`,
    `~zombi2.traits.DiscreteTrait`, `~zombi2.genomes.FamilyGenome` — never a finished result. A
    finished result is a driver you already have, which is conditioning, and it belongs in the driven
    level's own run."""
    from ..species import BirthDeath

    kinds = {"species": [], "traits": [], "genomes": []}
    for p in participants:
        if isinstance(p, BirthDeath):
            kinds["species"].append(p)
        elif isinstance(p, (DiscreteTrait, ContinuousTrait)):
            kinds["traits"].append(p)
        elif isinstance(p, FamilyGenome):
            kinds["genomes"].append(p)
        else:
            raise TypeError(
                f"joint.simulate takes process specs — species.birth_death(...), "
                f"traits.discrete(...), traits.continuous(...), genomes.genome(...) — and got "
                f"{p!r}. A finished result is a "
                f"driver you already have, which is conditioning: pass it to the driven level's own "
                f"run instead.")
    return kinds


def _grow_genomes_traits(rng, tree, genome: FamilyGenome, trait: DiscreteTrait, trait_keys,
                         seed) -> "tuple":
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

    Returns ``(genomes_out, genome_events, named, node_values, trait_changes)``.
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
    return genomes_out, events, named, initial_genome, node_values, changes


def _on_a_given_tree(kinds, *, tree, seed) -> JointResult:
    """The joint models whose tree is an **input** — today, a genome and a discrete trait.

    ``tree`` is required here for the reason it is refused when the species level is a participant:
    a run simulates what it is given specs for, and takes everything else."""
    from ..tree import as_tree

    n_traits, n_genomes = len(kinds["traits"]), len(kinds["genomes"])
    if not (n_traits == 1 and n_genomes == 1):
        raise NotImplementedError(
            f"a joint run on a tree you supply is built for one genome and one trait; got "
            f"{n_genomes} genome(s) and {n_traits} trait(s). A level joined to **itself** does not "
            f"come here at all — that is one level and one result, so it stays on that level's own "
            f"function with joint=True.")
    if tree is None:
        raise ValueError(
            "neither participant simulates the species tree, so this run needs one: pass tree=. "
            "Give what you are not simulating.")
    tree = as_tree(tree, level="joint")
    genome, trait = kinds["genomes"][0], kinds["traits"][0]
    if not isinstance(trait, DiscreteTrait):
        raise NotImplementedError(
            "a genome and a CONTINUOUS trait driving each other on a given tree is not built. The "
            "genome's race would have to be sliced against the diffusion, as a continuous trait "
            "driving speciation already is. Use traits.discrete(...) here, or grow one level first "
            "and condition the other on it.")

    # each level must actually read the other, or these are two independent runs wearing one call
    trait_keys = {_ONE_TRAIT} | ({f"traits:{trait.name}"} if trait.name else set())
    reads_trait, reads_genome = False, False
    for rate in (genome.duplication, genome.loss, genome.origination, genome.transfer):
        r = as_rate(rate, default_scope=PerLineage)
        for m in r.modifiers:
            if isinstance(m, Driven) and m.driver in trait_keys:
                reads_trait = True
    from ..traits.discrete import _switch_specs
    for spec in _switch_specs(trait.switch):
        if not isinstance(spec, (int, float)):
            for m in as_rate(spec, default_scope=PerLineage).modifiers:
                if isinstance(m, Driven) and str(m.driver).startswith("genomes:"):
                    reads_genome = True
    if not (reads_trait or reads_genome):
        raise ValueError(
            'neither level reads the other, so this is two independent runs rather than one joint '
            'one. Give a genome rate a scaled_by("trait", ...), or the trait\'s switch a '
            'scaled_by("genomes:<family>", ...) — or run the two levels separately.')

    rng, seed = stream("joint", seed)
    genomes_out, events, named, initial, node_values, changes = _grow_genomes_traits(
        rng, tree, genome, trait, tuple(sorted(trait_keys)), seed)
    changes.sort(key=lambda c: c.time)
    return JointResult(
        # the tree came in rather than out, so its own event log is not this run's to write
        SpeciesResult(tree, [], seed, []), seed,
        trait=TraitsResult(tree, node_values, changes, seed, kind="discrete"),
        genome=FamilyGenomesResult(tree, genomes_out, events, seed, named, {}, initial,
                                   genome.max_family_size))


def simulate(*participants, tree=None, seed=None, max_lineages=100_000) -> JointResult:
    """Simulate two levels **at once**, because neither can be finished before the other starts
    (SPEC §2–4).

    Each participant is a process spec, and you **give what you are not simulating**::

        # the tree is one of the two, so it comes out of the run
        joint.simulate(species.birth_death(birth=faster_if_large, death=0.2, n_extant=100),
                       traits.discrete(name="size", states=["small", "large"], switch=0.1), seed=1)

        joint.simulate(species.birth_death(birth=faster_with_toxin, n_extant=100),
                       genomes.genome(origination=0.2, loss=0.1, families=[family("toxin")]), seed=1)

    A rate reads the other participant by **name**, ``"<level>:<handle>"`` — ``"traits:size"`` for a
    named trait, ``"genomes:toxin"`` for a declared family, ``"genomes:count"`` for a lineage's whole
    gene count. A run holding one unnamed trait also answers to ``"trait"``.

    `~zombi2.traits.continuous` is a participant too, and it is the one driver that does not race
    exactly: a diffusion moves at every instant, so the run holds it fixed across a ``step=`` written
    on the connection and releases it at each boundary.

    The species tree is an output exactly when `~zombi2.species.birth_death` is one of the
    participants; otherwise ``tree`` supplies it. A level driving **itself** does not come here at
    all — that is one level and one result, so it stays on that level's own function with
    ``joint=True``.

    Returns a `JointResult`. Deterministic given ``seed``.
    """
    kinds = _classify(participants)
    n_species, n_traits, n_genomes = (len(kinds[k]) for k in ("species", "traits", "genomes"))
    if n_species == 0:
        return _on_a_given_tree(kinds, tree=tree, seed=seed)
    if tree is not None:
        raise ValueError(
            "the species tree is one of the things this run simulates, so it comes out rather than "
            "going in: drop tree=, or drop species.birth_death(...) and hand the tree over.")
    if n_species > 1:
        raise ValueError("give one species.birth_death(...) — a run grows one tree.")
    if n_traits + n_genomes != 1:
        raise ValueError(
            "give exactly one level for the tree to be simulated with: traits.discrete(...) or "
            f"genomes.genome(...). Got {n_traits} trait(s) and {n_genomes} genome(s).")
    spec = kinds["species"][0]
    driver = kinds["traits"][0] if n_traits else kinds["genomes"][0]
    return _simulate_joint(birth=spec.birth, death=spec.death,
                          n_extant=spec.n_extant, total_time=spec.total_time,
                          seed=seed, max_lineages=max_lineages,
                          **({"trait": driver} if n_traits else {"genome": driver}))


def simulate_joint(**_):
    """Retired. A joint run is written as its **participants** now (SPEC §2–4)::

        joint.simulate(species.birth_death(birth=…, death=…, n_extant=100),
                       traits.discrete(name="size", states=[…], switch=0.1), seed=1)

    The rates and the stop condition move onto `~zombi2.species.birth_death`, which makes the tree
    one of the things being simulated rather than a keyword of the run; and the driver is a
    participant beside it rather than a ``trait=`` / ``genome=`` slot. That is what lets one function
    take every joint model instead of one per pair.
    """
    raise TypeError(
        "simulate_joint is no longer written — a joint run is its participants: "
        "joint.simulate(species.birth_death(birth=…, death=…, n_extant=100), "
        "traits.discrete(name='size', states=[…], switch=0.1), seed=1). The rates and the stop "
        "condition go on species.birth_death, and the driver is a participant beside it.")


__all__ = ["simulate", "JointResult"]
