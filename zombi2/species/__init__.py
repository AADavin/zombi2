"""Species trees — the forward birth-death engine.

Per-lineage birth and death grow a tree forward in time and record every
speciation and extinction. ``birth`` and ``death`` are full **rate specs** — a number,
a scope wrapper, or a product — so ``birth = scope.Global(1.0)`` gives one shared
tree-wide budget (linear growth), ``birth = 1.0 * mod.OnTotalDiversity(cap=100)`` slows the
tree as it fills up, ``birth = 1.0 * mod.OnTime({...})`` runs a skyline (the interval-aware
sampler steps to each breakpoint), and two modifiers make the rate vary from lineage to lineage:
``birth = 1.0 * mod.Inherited(per="lineage", spread=0.2)`` lets it drift down the tree (clade drift, ClaDS) and
``birth = 1.0 * mod.Drawn(per="lineage", spread=0.2)`` gives each lineage an independent draw (*relaxed*
rates). Under either, each lineage threads its own factor and the lineage that speciates or dies is
drawn **weighted** by its effective rate.
"""

from __future__ import annotations

import functools
import math
import pathlib
from dataclasses import dataclass, field


from ..rates.modifiers import (describe, DRAWN, INHERITED, OnTime, OnTotalDiversity,
                               check_one_memory, is_implemented, values_at_birth,
                               values_at_split)
from ..rng import stream
from .._runtime.draw import weighted_index as _weighted_index
from .._runtime.progress import progress_bar
from .._runtime.summary import write_summary
from ..rates.rate import as_rate
from ..rates.scope import Global, PerLineage
from ..tree import Node, Tree, prune

#: The rate grammar this level supports (SPEC §5). Both the engine's gate below and the CLI's help
#: read this, so a modifier can never be advertised without being implemented — or silently ignored.
IMPLEMENTED_SCOPES = (PerLineage, Global)
IMPLEMENTED_MODIFIERS = (OnTime, OnTotalDiversity, (INHERITED, "lineage"), (DRAWN, "lineage"))



@dataclass(frozen=True)
class Event:
    """A recorded event in the true history: a speciation (with its two children) or an extinction."""

    time: float
    kind: str  # "speciation" | "extinction"
    node: int
    children: tuple[int, int] | None = None




_WRITE_OUTPUTS = ("complete", "extant", "events", "fossils", "fates", "summary")  # the write vocabulary the CLI reuses


@dataclass
class SpeciesResult:
    """What ``simulate_species_tree`` returns: the ``complete_tree`` (with the dead) and the derived
    ``extant_tree`` (the observed survivors), the ``events`` log (the recorded true history), the
    ``seed``, and any ``fossils``. (The ``record=`` memory dial lands with the data-heavy levels.)"""

    complete_tree: Tree
    events: list[Event]
    seed: int | None
    #: recovered fossils as ``(lineage_id, time)`` pairs, sorted by time — a side output, present
    #: only when ``fossils`` was set; the fossil's lineage is not removed and is not in the extant tree
    fossils: list[tuple[int, float]] = field(default_factory=list)

    def __repr__(self) -> str:
        fossils = f", {len(self.fossils)} fossils" if self.fossils else ""
        return (f"SpeciesResult({self.n_extant} extant tips, "
                f"{len(self.complete_tree.nodes)} nodes{fossils}, seed={self.seed})")

    @property
    def n_extant(self) -> int:
        """The number of **observed** survivors — the extant tips. Under ``sampling < 1`` this is
        the sampled subset (the rest are ``unsampled``), so it matches the extant tree's tip count."""
        return len(self.complete_tree.extant_leaves())

    @functools.cached_property
    def extant_tree(self) -> Tree | None:
        """The survivors' tree — the complete tree pruned to extant lineages with the
        unifurcations suppressed (dated, bifurcating). ``None`` if nothing survived, which
        ``simulate_species_tree`` refuses to return: a run with no present raises there instead, so a
        result that came from it always has one."""
        return prune(self.complete_tree, keep="extant")

    def summary(self) -> dict:
        """What this run produced, as a plain dict — the payload of ``species_summary.json``.

        Counts, not parameters: the log already says what was asked for. The realised rates are here
        because they are the cheapest check anyone can make on a tree — events divided by the exposure
        that generated them, which is what a declared per-lineage rate means."""
        nodes = self.complete_tree.nodes
        tips = [n for n in nodes.values() if n.children is None]
        extant = self.complete_tree.extant_leaves()
        speciations = sum(1 for e in self.events if e.kind == "speciation")
        extinctions = sum(1 for e in self.events if e.kind == "extinction")
        # total branch length: every node's own branch, which is the exposure a per-lineage rate ran on
        exposure = sum(n.end_time - n.birth_time for n in nodes.values())
        height = max(n.end_time for n in extant) if extant else None
        root = nodes[self.complete_tree.root]
        return {
            "level": "species",
            "seed": self.seed,
            "tips": {"extant": len(extant), "extinct": len(self.complete_tree.extinct_leaves()),
                     "unsampled": len(self.complete_tree.unsampled_leaves()), "total": len(tips)},
            "nodes": len(nodes),
            "events": {"speciation": speciations, "extinction": extinctions},
            "fossils": len(self.fossils),
            "tree": {"height": height, "stem_length": root.end_time - root.birth_time,
                     "total_branch_length": exposure},
            # events per lineage per unit time, as declared rates are counted. A sanity check, not a
            # parameter: it is what the run realised, which a conditioned stop condition can bias.
            "realised_rates": {
                "birth": round(speciations / exposure, 6) if exposure else None,
                "death": round(extinctions / exposure, 6) if exposure else None},
        }

    def write(self, directory, outputs=None) -> None:
        """Write outputs to ``directory``, each file prefixed ``species_``; ``outputs`` selects which
        (default = all applicable): ``"complete"`` → ``species_complete.nwk``, ``"extant"`` →
        ``species_extant.nwk`` (if any survived), ``"events"`` → ``species_events.tsv`` (the
        always-recorded true history, ``time`` · ``kind`` · ``parents`` · ``children``),
        ``"fossils"`` → ``species_fossils.tsv`` (if any recovered),
        ``"fates"`` → ``species_fates.tsv`` (each tip's resolved fate).

        ``species_events.tsv`` names the lineages an event consumed and the lineages it produced, the
        same ``parents`` / ``children`` pair every event file uses: a ``speciation`` row is one parent
        and its two children (``;``-packed), an ``extinction`` row is the dying lineage as the parent
        with no children.

        ``species_fates.tsv`` is the tip-fate table: one ``lineage<TAB>fate`` row per tip, with fate
        one of ``extant`` / ``extinct`` / ``unsampled``. Fate is resolved once, at the end of the run,
        on the same stable ``n<id>`` that keys every other file, so it never renames anything — it is
        a materialised view of information the run already holds. It exists because the ``.nwk`` records
        only branch lengths, from which a reader cannot tell an extinct tip from a survivor that sits at
        the present; this table says so directly, so a downstream level can build the extant set from
        fate rather than guessing from tip depth."""
        if outputs is None:
            outputs = _WRITE_OUTPUTS
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "complete" in outputs:
            (d / "species_complete.nwk").write_text(self.complete_tree.to_newick() + "\n", encoding="utf-8")
        if "extant" in outputs and self.extant_tree is not None:
            (d / "species_extant.nwk").write_text(self.extant_tree.to_newick() + "\n", encoding="utf-8")
        name = self.complete_tree.labels()
        if "events" in outputs:
            # parents / children, not lineage / children: the lineage an event consumed IS its parent
            # (the one that split, or the one that died), so one column pair reads right for both kinds
            rows = ["time\tkind\tparents\tchildren"]
            for e in self.events:
                kids = ";".join(name[c] for c in e.children) if e.children else ""
                rows.append(f"{e.time:.6g}\t{e.kind}\t{name[e.node]}\t{kids}")
            (d / "species_events.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        if "fossils" in outputs and self.fossils:
            # fossils are drawn along every branch, a surviving lineage's as readily as an extinct one's
            rows = ["lineage\ttime"] + [f"{name[i]}\t{t:.6g}" for i, t in self.fossils]
            (d / "species_fossils.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        if "summary" in outputs:
            write_summary(d / "species_summary.json", self.summary())
        if "fates" in outputs:
            # one row per tip (extant / extinct / unsampled); internal nodes are always speciations
            rows = ["lineage\tfate"]
            for n in sorted(self.complete_tree.leaves(), key=lambda nd: nd.id):
                rows.append(f"{name[n.id]}\t{n.fate}")
            (d / "species_fates.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")





_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_extant


def _per_lineage(rate) -> tuple:
    """The modifiers a rate carries **per lineage**, in the order they were written.

    Two modifiers make a rate vary from lineage to lineage, and the engine threads both the same way —
    one factor per living lineage, with the lineage that speciates or dies drawn **weighted** by its
    effective rate rather than uniformly. They differ only in where a daughter's factor comes from,
    which is the whole uncorrelated / autocorrelated split (``SPEC §5``): an inherited value
    inherits it from the parent and nudges it (clade drift, ClaDS), a per-lineage draw
    draws a fresh one with no memory of the parent (*relaxed* rates). A rate carrying neither keeps a
    factor of 1 throughout and its lineage is picked uniformly, exactly as before.

    The list comes from `Rate.carried_modifiers`, which reports what each modifier
    declares it reads rather than testing its class — so a per-lineage modifier this engine has never
    heard of is threaded like the two it has, and **every** one is kept. Taking only the first was
    how a second silently vanished from a run that still reported it."""
    return tuple(m for m, _ in rate.carried_modifiers(unit="lineage"))




def _grow(rng, birth_rate, death_rate, n_extant: int | None, total_time: float | None,
          pulses: list[tuple[float, float]], progress: bool = False,
          max_lineages: int | None = None) -> tuple[Tree, list[Event]]:
    """Grow one forward birth-death tree until it reaches ``n_extant`` living lineages,
    reaches ``total_time``, or dies out. Returns the complete tree and the event log.

    When ``birth`` or ``death`` carries an inherited value or a per-lineage draw
    the rate is *per-lineage*: every lineage threads its own factor, so the lineage that speciates or
    dies is drawn **weighted** by its effective rate rather than uniformly. Under an inherited value a
    daughter's factor is its parent's, nudged at the split (clade drift); under a per-lineage draw it is an
    independent draw with no memory of the parent (relaxed rates). Birth and death vary independently
    of each other. A rate with neither keeps a factor of 1 and picks uniformly, exactly as before.

    ``pulses`` are scheduled mass extinctions as ``(time, survival)`` pairs sorted by time (time runs
    forward from the origin): at each instant every standing lineage is kept with probability
    ``survival`` and otherwise becomes an extinct leaf. They sit at a point on the timeline, so the
    caller passes them only when ``total_time`` is set."""
    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent: int | None, t: float) -> int:
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    birth_drift = _per_lineage(birth_rate)  # the per-lineage modifiers on each rate, possibly none
    death_drift = _per_lineage(death_rate)

    root = new_node(None, 0.0)
    alive = [root]  # a list so picks are reproducible given the seed
    # each lineage's own factors — one per carried modifier — kept in lock-step with `alive` under
    # swap-remove; a rate with no per-lineage modifier holds an empty tuple, whose product is 1.0, so
    # its total is just scope(base) × modifiers over n, picked uniform. The root draws under
    # Drawn(per='lineage') (it is a lineage like any other) and does not under Inherited(per='lineage') (its factor is the
    # ladder's starting point), which is what `values_at_birth` decides.
    root_shared: dict = {}  # the root lineage's cache, so an object on both rates draws once
    inh_b = [values_at_birth(birth_drift, rng, root_shared)]
    inh_d = [values_at_birth(death_drift, rng, root_shared)]
    t = 0.0
    events: list[Event] = []
    pulse_idx = 0  # the next unfired mass extinction in `pulses`

    # a tree grows toward whichever stop condition was given: a tip count, or a time
    bar = progress_bar(n_extant if n_extant is not None else total_time, "species",
                       unit="tip" if n_extant is not None else "time", enabled=progress)
    # A run conditioned on time has no natural ceiling: standing diversity grows like
    # exp((birth - death) * t), so a rate a little too high or a time a little too long is the
    # difference between a thousand lineages and ten million. The guard RAISES rather than stopping
    # early — a tree cut off at a size is no longer a sample from the process asked for, and handing
    # one back as if it were would be worse than not running at all.
    ceiling = None if max_lineages is None else max(max_lineages, n_extant or 0)
    while alive:
        bar.to(len(alive) if n_extant is not None else t)
        n = len(alive)
        if ceiling is not None and n > ceiling:
            bar.close()
            raise RuntimeError(
                f"the tree passed {ceiling} standing lineages at time {t:.3g} and is still growing "
                f"— birth exceeds death by enough that this run has no realistic end. Lower the "
                f"rates, shorten total_time, cap the growth with a modifier "
                f"(birth * OnTotalDiversity(cap=...)), or raise max_lineages if the size is what "
                f"you want (max_lineages=None removes the guard).")
        # standing diversity = the living lineages; OnTotalDiversity/OnTime read `diversity`/`time`
        ctx = {"diversity": n, "time": t}
        # a drifting rate's total is the sum over lineages of each lineage's effective rate —
        # scope(base) × modifiers evaluated per lineage through its inherited factor (lineages=1
        # is one lineage); a non-drifting rate is scope(base) × modifiers once, over all n lineages
        if birth_drift:
            w_b = [birth_rate.effective(lineages=1, carried_factor=math.prod(x), **ctx) for x in inh_b]
            total_birth = sum(w_b)
        else:
            total_birth = birth_rate.effective(lineages=n, **ctx)
        if death_drift:
            w_d = [death_rate.effective(lineages=1, carried_factor=math.prod(x), **ctx) for x in inh_d]
            total_death = sum(w_d)
        else:
            total_death = death_rate.effective(lineages=n, **ctx)
        total = total_birth + total_death
        # the total rate is constant until the next skyline breakpoint, mass extinction, or the total_time
        # limit — advance no further than the earliest of them before re-evaluating
        next_change = min(birth_rate.next_change(t), death_rate.next_change(t))
        next_pulse = pulses[pulse_idx][0] if pulse_idx < len(pulses) else math.inf
        horizon = min(next_change, next_pulse)
        if total_time is not None:
            horizon = min(horizon, total_time)

        if total > 0.0:
            t_event = t + float(rng.exponential(1.0 / total))
            if t_event < horizon:  # an event fires before the rate changes
                t = t_event
                if n == n_extant:
                    # already at the target: stop at the time this next event WOULD fire (do not
                    # apply it), so the present sits *after* the last split and the two newest tips
                    # get a real, non-zero branch length.
                    break
                # birth vs death by their totals; then WHICH lineage — weighted by its effective
                # rate if that rate drifts (so faster lineages are likelier), else uniform
                speciates = rng.random() < total_birth / total
                if speciates:
                    i = _weighted_index(rng, w_b, total_birth) if birth_drift else int(rng.integers(n))
                else:
                    i = _weighted_index(rng, w_d, total_death) if death_drift else int(rng.integers(n))
                node = alive[i]
                parent_b, parent_d = inh_b[i], inh_d[i]
                alive[i] = alive[-1]  # swap-remove keeps picks O(1); the inherited factors move in step
                alive.pop()
                inh_b[i] = inh_b[-1]; inh_b.pop()
                inh_d[i] = inh_d[-1]; inh_d.pop()
                if speciates:
                    nodes[node].end_time = t
                    nodes[node].fate = "speciation"
                    c1, c2 = new_node(node, t), new_node(node, t)
                    nodes[node].children = (c1, c2)
                    alive.extend((c1, c2))
                    # each daughter takes its own factors — the parent's nudged under Inherited(per='lineage'), a
                    # fresh independent draw under Drawn(per='lineage') (an empty tuple when the rate carries
                    # neither, so nothing is drawn and the product stays 1.0). One cache per
                    # daughter, so a modifier written on both rates draws once for that daughter;
                    # the birth pass still runs before the death pass, which keeps the draw order
                    # of a run that shares nothing exactly as it was.
                    d1: dict[int, float] = {}
                    d2: dict[int, float] = {}
                    inh_b.extend((values_at_split(birth_drift, parent_b, rng, d1),
                                  values_at_split(birth_drift, parent_b, rng, d2)))
                    inh_d.extend((values_at_split(death_drift, parent_d, rng, d1),
                                  values_at_split(death_drift, parent_d, rng, d2)))
                    events.append(Event(t, "speciation", node, (c1, c2)))
                else:
                    nodes[node].end_time = t
                    nodes[node].fate = "extinct"
                    events.append(Event(t, "extinction", node))
                continue

        # no stochastic event fired before the horizon
        if math.isinf(horizon):
            break  # nothing scheduled and the rate never changes again → nothing more can happen
        if total_time is not None and horizon == total_time:
            t = total_time
            break
        if pulse_idx < len(pulses) and horizon == next_pulse:
            # a mass extinction: each standing lineage is kept with probability `survival`, the rest
            # become extinct leaves at this instant (their inherited factors leave with them)
            t = next_pulse
            survival = pulses[pulse_idx][1]
            pulse_idx += 1
            kept_a: list[int] = []
            kept_b: list[tuple[float, ...]] = []
            kept_d: list[tuple[float, ...]] = []
            for k, node_id in enumerate(alive):
                if survival >= 1.0 or rng.random() < survival:
                    kept_a.append(node_id)
                    kept_b.append(inh_b[k])
                    kept_d.append(inh_d[k])
                else:
                    nodes[node_id].end_time = t
                    nodes[node_id].fate = "extinct"
                    events.append(Event(t, "extinction", node_id))
            alive[:] = kept_a
            inh_b[:] = kept_b
            inh_d[:] = kept_d
            continue
        t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    bar.close()
    for i in alive:  # whoever is still alive reached the present
        nodes[i].end_time = t
        nodes[i].fate = "extant"

    return Tree(nodes, root), events


def _mass_extinction_pulses(mass_extinctions, total_time: float | None) -> list[tuple[float, float]]:
    """Turn user ``(time, fraction_lost)`` pulses into the engine's ``(time, survival)`` pairs,
    sorted by time. Time runs **forward from the origin** (t=0, the initial lineage), so
    ``(3.0, 0.75)`` = at time 3.0, 75% of
    the standing lineages die (survival 0.25). Empty when none are given. A pulse sits at a point on
    the timeline, so it needs a fixed end — ``total_time`` set — and must fall inside ``(0, total_time)``."""
    if not mass_extinctions:
        return []
    if total_time is None:
        raise ValueError(
            "mass_extinctions need a run with a fixed end — give total_time=..., not n_extant= "
            "(under n_extant= the run can stop before a pulse's time is reached)"
        )
    pulses: list[tuple[float, float]] = []
    for pulse in mass_extinctions:
        time, fraction = pulse
        if (isinstance(time, bool) or not isinstance(time, (int, float))
                or not math.isfinite(time) or not 0.0 < time < total_time):
            raise ValueError(
                f"each mass extinction time must be a number strictly between 0 and total_time ({total_time}), "
                f"got {time!r}"
            )
        if (isinstance(fraction, bool) or not isinstance(fraction, (int, float))
                or not 0.0 <= fraction <= 1.0):
            raise ValueError(f"each mass extinction fraction lost must be in [0, 1], got {fraction!r}")
        pulses.append((time, 1.0 - fraction))
    pulses.sort()
    return pulses


def _apply_sampling(tree: Tree, rho: float, rng) -> None:
    """Incomplete extant sampling: relabel each surviving lineage ``"unsampled"`` with probability
    ``1 - rho`` (in place). The extant tree then prunes to the sampled survivors, while the unsampled
    ones stay in the complete tree, told apart by their fate. ``rho = 1`` observes everyone."""
    if rho >= 1.0:
        return
    for i in sorted(tree.nodes):  # id order + one draw per survivor → reproducible given the seed
        node = tree.nodes[i]
        if node.fate == "extant" and float(rng.random()) >= rho:
            node.fate = "unsampled"


def _recover_fossils(tree: Tree, rate: float, rng) -> list[tuple[int, float]]:
    """Recover fossils along every branch of the complete tree: a branch of length ``L`` yields
    ``Poisson(rate × L)`` fossils, each at a uniform time on the branch. Returns ``(lineage_id,
    time)`` pairs sorted by time. A pure side output — no lineage is removed."""
    if rate <= 0.0:
        return []
    fossils: list[tuple[int, float]] = []
    for i in sorted(tree.nodes):  # id order, then Poisson + uniforms → reproducible given the seed
        node = tree.nodes[i]
        length = node.end_time - node.birth_time
        if length <= 0.0:
            continue
        for _ in range(int(rng.poisson(rate * length))):
            fossils.append((i, node.birth_time + float(rng.random()) * length))
    fossils.sort(key=lambda ft: ft[1])
    return fossils


def simulate_species_tree(birth, death=0.0, *, n_extant=None, total_time=None,
                          mass_extinctions=None, sampling=1.0, fossils=0.0, seed=None,
                          progress=False, max_lineages=100_000) -> SpeciesResult:
    """Grow a forward birth-death tree.

    ``birth`` and ``death`` are rate specs (a number, a ``scope`` wrapper, or a product
    with modifiers); the default scope is **per lineage** (each lineage speciates/dies at
    the base rate, so the tree grows exponentially). Yule = ``death=0``.

    **Rates that vary from lineage to lineage.** ``1.0 * mod.Inherited(per="lineage", spread=σ)`` is *inherited*
    variation — a daughter starts from its parent's rate and is nudged at the split, so fast clades
    stay fast (clade drift; the literature's ClaDS). ``1.0 * mod.Drawn(per="lineage", spread=σ)`` is
    *independent* variation — every lineage draws its own multiplier with no memory of its parent
    (*relaxed* rates), which is the null to compare drift against: the same amount of rate
    heterogeneity, none of it heritable, so the tree-shape signature that heritability leaves
    (lopsidedness — fast clades hoarding the tips) is absent. Both draws are mean-corrected, so
    widening ``spread`` spreads lineages out without moving the average one off the base rate; both
    make the lineage that speciates or dies drawn weighted by its own rate; and both must be counted
    per lineage. They answer the same question and a rate carrying both is refused.

    Stop at exactly ``n_extant`` living lineages, **or** at ``total_time`` — give exactly
    one. ``n_extant`` is **conditioned on survival**: a birth-death tree can die out, so we
    restart (advancing the same generator) until one reaches ``n_extant``. Deterministic given
    ``seed``.

    **Where ``n_extant`` puts the present.** The run stops the first moment ``n_extant`` lineages
    are alive together, then draws one more waiting time and places the present where that next
    event *would* have fired, without applying it — so the two newest tips get a real branch length
    rather than a zero-length one. Under pure birth this is exactly the general sampling approach
    (Hartmann, Wong & Stadler 2010): the tree is the process observed at an instant drawn uniformly
    over the time it spends holding n lineages. **With extinction it is not.** It is a *first
    hitting* rule — the run stops the first time it touches n, so an interval at n lineages reached
    by falling back from n+1 is never sampled — and the trees are correspondingly shallower than the
    birth-death process conditioned on n tips. The gap grows with turnover and shrinks with n: it is
    within noise at ``death=0``, around a tenth of the tree height at n=10 with ``death/birth=0.4``,
    and roughly a third to a half at n=10 with ``death/birth=0.8``; by n=50 at moderate turnover it
    is back in the noise. If you publish trees grown this way, say which rule made them — a rate
    estimator applied to them will otherwise look broken for reasons that are not its fault.

    ``total_time`` is not conditioned on survival: it can die out, and then it **raises** rather than
    handing back a tree with no present. Looping over seeds and skipping the failures is survival
    conditioning by another name, and changes the distribution of everything downstream.

    ``mass_extinctions`` is a list of ``(time, fraction_lost)`` pulses — e.g. ``[(3.0, 0.75)]`` culls
    75% of the lineages alive at time 3.0 (time runs forward from the origin, t=0). It is a point-in-time
    intervention on the process (not a rate) placed on the timeline, so it needs a fixed end:
    give ``total_time`` (not ``n_extant``), with each time strictly inside ``(0, total_time)``.

    ``sampling`` (ρ, default 1.0) is incomplete extant sampling: each survivor is observed with
    probability ρ, the rest relabelled ``unsampled``. It prunes the **extant tree** to the sampled
    survivors (the unsampled ones remain only in the complete tree). ``n_extant`` still stops at that
    many *survivors*; sampling then thins what you observe, so ``result.n_extant`` can be smaller.

    ``fossils`` is a recovery rate along the branches: each branch of length ``L`` yields
    ``Poisson(fossils × L)`` fossils, returned as ``result.fossils`` = ``(lineage, time)`` pairs. A
    **side output** — the fossil's lineage is not removed and does not enter the extant tree.
    """
    birth_rate = as_rate(birth, default_scope=PerLineage)
    death_rate = as_rate(death, default_scope=PerLineage)
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        # a modifier this engine does not thread would return its default factor of 1.0 — a run that
        # is quietly not the model asked for — so reject it (SPEC §5, the genome engine's discipline)
        if not isinstance(rate.scope, IMPLEMENTED_SCOPES):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the species engine counts "
                f"lineages — use PerLineage(...) (the default, so a bare number is enough) or "
                f"Global(...) for one shared budget."
            )
        for m in rate.modifiers:
            if m.reads == (DRAWN, "family"):
                # not a missing feature: there is nothing here for it to mean
                raise ValueError(
                    f"{label} carries Drawn(per='family'), but a species tree has no gene families — Drawn(per='family') "
                    f"belongs on a genomes rate. For per-lineage heterogeneity here use Drawn(per='lineage') "
                    f"(independent) or Inherited(per='lineage') (inherited)."
                )
            if not is_implemented(m, IMPLEMENTED_MODIFIERS, "species"):
                raise ValueError(
                    f"{label} carries {describe(m)}, which the species engine does not "
                    f"support. It takes OnTime (skyline), OnTotalDiversity (diversity-dependent), "
                    f"Inherited(per='lineage') (inherited rate drift, ClaDS) and Drawn(per='lineage') (independent "
                    f"per-lineage rates). A birth or death that reads an evolved value cannot be "
                    f"conditioned — the tree and its driver grow together — so it is a joint run: "
                    f"joint.simulate_joint(birth=1.0 * mod.Driven('trait', {{...}}), ...)."
                )
        # SPEC §5: one memory structure per axis. Drawn(per='lineage') has none and Inherited(per='lineage') has a continuous
        # one, so a rate carrying both asks for a lineage's factor to be independent of its parent's
        # and inherited from it at once — there is no model there to implement, so say so rather than
        # silently letting whichever comes first win.
        check_one_memory(_per_lineage(rate), label=label, unit="lineage")
        if (per_lineage := _per_lineage(rate)) and not isinstance(rate.scope, PerLineage):
            raise ValueError(
                f"{label} carries {describe(per_lineage[0])} (a factor per lineage) but its scope "
                f"is {type(rate.scope).__name__}; a rate that varies by lineage must be counted per "
                f"lineage — drop the scope wrapper (per lineage is the default) or use PerLineage(...)"
            )
    if (n_extant is None) == (total_time is None):
        raise ValueError("give exactly one of n_extant or total_time")
    if n_extant is not None and (isinstance(n_extant, bool) or not isinstance(n_extant, int) or n_extant < 1):
        raise ValueError(f"n_extant must be a positive integer, got {n_extant!r}")
    if total_time is not None and (not isinstance(total_time, (int, float)) or not math.isfinite(total_time) or total_time <= 0):
        raise ValueError(f"total_time must be a positive finite number, got {total_time!r}")
    if isinstance(fossils, bool) or not isinstance(fossils, (int, float)) or not math.isfinite(fossils) or fossils < 0:
        raise ValueError(f"fossils must be a non-negative finite rate, got {fossils!r}")
    if isinstance(sampling, bool) or not isinstance(sampling, (int, float)) or not 0.0 < sampling <= 1.0:
        raise ValueError(f"sampling must be a fraction in (0, 1], got {sampling!r}")
    pulses = _mass_extinction_pulses(mass_extinctions, total_time)  # [] unless mass_extinctions given (needs total_time)

    rng, seed = stream("species", seed)     # own stream, and a drawn seed if none was given

    def _finish(tree: Tree, events: list[Event]) -> SpeciesResult:
        # observe (sampling relabels survivors) then recover fossils along the grown branches
        alive = sum(1 for nd in tree.nodes.values() if nd.fate == "extant")
        _apply_sampling(tree, sampling, rng)
        # Sampling can take none of them, and then the run has no present — the same dead end the
        # extinction guard above refuses, reached by the other road. Refusing here too keeps the two
        # "nothing observed" outcomes consistent: neither hands back a result whose extant tree is
        # None for a downstream level to trip over. Only reachable with sampling < 1, since both
        # callers guarantee a survivor before this point.
        if not any(nd.fate == "extant" for nd in tree.nodes.values()):
            raise RuntimeError(
                f"sampling={sampling:g} observed none of the {alive} survivor"
                f"{'' if alive == 1 else 's'}, so the run has no present to grow a genome, sequence "
                f"or trait along. This is the sampling process, not a bad parameter — it has "
                f"probability {(1 - sampling) ** alive:.3g} here — so raise sampling, ask for more "
                f"survivors, or draw another seed.")
        return SpeciesResult(tree, events, seed, _recover_fossils(tree, fossils, rng))

    if total_time is not None:
        tree, events = _grow(rng, birth_rate, death_rate, None, total_time, pulses, progress,
                             max_lineages)
        # A time-conditioned run is not conditioned on survival, so with death ≥ birth it can reach
        # total_time with nothing alive. An empty tree is not a sample anyone can use — the extant
        # tree is None and every downstream level would otherwise mistake the last-dying tip for a
        # survivor — so refuse it here rather than hand back a tree with no present.
        if not any(nd.fate == "extant" for nd in tree.nodes.values()):
            raise RuntimeError(
                f"the run went extinct before total_time={total_time:g}: no lineage is alive at the "
                f"present, so there is nothing to grow a genome, sequence or trait along. With death "
                f"close to or above birth, total extinction is likely — lower death, shorten "
                f"total_time, or use n_extant=... (which is conditioned on survival).")
        return _finish(tree, events)

    for _ in range(_MAX_ATTEMPTS):
        tree, events = _grow(rng, birth_rate, death_rate, n_extant, None, [], progress,
                             max_lineages)
        if sum(1 for nd in tree.nodes.values() if nd.fate == "extant") == n_extant:  # survivors (pre-sampling)
            return _finish(tree, events)
    raise RuntimeError(
        f"could not grow a tree to {n_extant} extant lineages in {_MAX_ATTEMPTS} attempts; "
        "birth must comfortably exceed death for large n_extant"
    )


__all__ = ["simulate_species_tree", "SpeciesResult", "Event"]
