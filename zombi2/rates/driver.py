"""A conditioned `DrivenBy`'s file-backing (SPEC §2).

When ``DrivenBy``'s ``source`` is a **filename**, the relation is *conditioned*: the driver was grown
first and written to a file, and two ordinary runs in order do the rest
(``loss = 0.25 * mod.DrivenBy("habitat.tsv", {...})``). This module — living beside ``DrivenBy`` in
``rates`` because it is that modifier's file end — turns the written driver into the per-lineage lookup
the target engine queries as it walks the (already-grown) tree. (Conditioning needs no engine of its
own: it *folds into the target level's* run; only genuinely-joint models get a dedicated engine,
``zombi2.joint``.)

The driver file is the trait **event log** (``trait_events.tsv``, written by
`zombi2.traits.TraitsResult.write()` with ``outputs=("events",)``): an ``initial`` row giving the
state at t=0, then every switch — ``time · kind · lineage · from · to``. The driver ran on the same
complete tree the target now runs on, so replaying the log **against that tree** rebuilds each
lineage's branch as constant stretches (a discrete driver switches *mid-branch*, so this is the exact
stochastic character map, not one value per branch). `DriverTrajectory` then answers both
*what is the driver on this lineage now?* (`value()`) and *when does it next
change?* (`next_change()`, so the target's Gillespie steps at each switch).

The join key is the **species node id**: ``node n7`` in the log is lineage 7 in the target run.
"""

from __future__ import annotations

import warnings

import bisect
import math
import pathlib

from .modifiers import DrivenBy
from ..tree import node_from_label, node_label


class DriverTrajectory:
    """A driver's value along every lineage, as a piecewise-constant function of time — the
    per-lineage lookup a conditioned `DrivenBy` reads.

    Built from segments ``{node_id: [(start_time, state), …]}`` (each lineage's branch cut into
    constant stretches, sorted by start). The engine calls `value()` to get a lineage's driver
    state at the current instant and `next_change()` to learn when it next switches (a horizon
    breakpoint, so the Gillespie re-evaluates the driven rate exactly at each switch)."""

    def __init__(self, segments: dict[int, list[tuple[float, object]]]) -> None:
        self._starts: dict[int, list[float]] = {}
        self._states: dict[int, list[object]] = {}
        for node_id, segs in segments.items():
            ordered = sorted(segs)  # by start time
            self._starts[node_id] = [s for s, _ in ordered]
            self._states[node_id] = [v for _, v in ordered]

    def states(self) -> set:
        """Every state the driver actually takes, anywhere on the tree — what a discrete mapping's
        keys are checked against, so a mapping that names none of them can be caught."""
        return {s for states in self._states.values() for s in states}

    def value(self, node_id: int, time: float) -> object:
        """The driver's state on lineage ``node_id`` at ``time`` — the segment whose start is the
        latest at or before ``time`` (right-continuous: at a switch instant the new state applies)."""
        starts = self._starts.get(node_id)
        if starts is None:
            raise KeyError(
                f"the driver file has no lineage {node_label(node_id)}; the driver must be grown "
                f"on the SAME "
                f"complete tree the target runs on (node ids must match)."
            )
        i = bisect.bisect_right(starts, time) - 1
        if i < 0:  # a query before the branch's first segment: clamp to the first (branch-start) state
            i = 0
        return self._states[node_id][i]

    def next_change(self, node_id: int, time: float) -> float:
        """The next time strictly after ``time`` at which lineage ``node_id``'s driver switches, else
        ``inf`` (it stays constant for the rest of the branch). Feeds the target Gillespie's horizon."""
        starts = self._starts.get(node_id)
        if starts is None:
            raise KeyError(f"the driver file has no lineage {node_label(node_id)} "
                           "(node ids must match the target tree).")
        i = bisect.bisect_right(starts, time)
        return starts[i] if i < len(starts) else math.inf


def load_driver(path, tree, *, step: float | None = None) -> DriverTrajectory:
    """Read a written driver and rebuild it against ``tree``. Two files can be a driver, and which
    one you need depends on which kind of trait was grown.

    A **discrete** trait's driver is its event log, ``trait_events.tsv``
    (``time · kind · lineage · from · to``: an ``initial`` row then every switch). The log alone is
    not enough — a switch says *when* the state changed, not what each branch started in — so the tree
    supplies branch birth/end times and the topology, and the reconstruction walks parent-before-child:
    the root begins in the ``initial`` row's state, every other lineage in its own ``on_speciation``
    state if it has one else its parent's ending state, and ``on_branch`` rows cut the branch into
    constant stretches. That is the exact stochastic character map.

    A **continuous** trait's driver is its value table, ``trait_values.tsv``
    (``node · kind · trait``). A diffusion has no switches to log — its event file holds only the
    ``initial`` row — so the values at the nodes are what carries it, and the path between them is
    interpolated at a resolution of ``step`` (`interpolated_segments`). Pointing a conditioned rate at
    a continuous trait's *event* log used to be accepted and to yield a driver frozen at the root
    value for the whole tree; it now raises and names the file to use instead.

    The join key is the species node id either way: ``n7`` in the file is lineage 7 here."""
    try:
        text = pathlib.Path(path).read_text(encoding="utf-8")
    except FileNotFoundError:
        raise FileNotFoundError(
            f"DrivenBy driver file not found: {str(path)!r}. A conditioned rate points at this file, but "
            f"it is not there — check the path (it is relative to where you run zombi2), or grow the "
            f"driver first (run the level that writes it) so there is something to condition on."
        ) from None
    rows = [line for line in text.splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"driver file {str(path)!r} is empty")
    header = rows[0].split("\t")
    expected = ["time", "kind", "lineage", "from", "to"]
    if header == ["node", "kind", "trait"]:
        return _load_values_driver(path, rows, tree, step)
    if header != expected:
        raise ValueError(
            f"driver file {str(path)!r} must be a trait event log with header {expected}, or a trait "
            f"value table with header ['node', 'kind', 'trait']; got {header}. Write one with "
            f"TraitsResult.write(dir, outputs=('events',)) or outputs=('values',)."
        )
    initial_state = None
    clado: dict[int, object] = {}                       # lineage -> its on-speciation start state
    switches: dict[int, list[tuple[float, object]]] = {}   # lineage -> [(time, to_state), …]
    for line in rows[1:]:
        parts = line.split("\t")
        if len(parts) != 5:
            raise ValueError(f"driver file {str(path)!r} row is not 5 columns: {line!r}")
        time_s, kind, node_s, _from, to = parts
        node_id = node_from_label(node_s)
        if kind == "initial":
            initial_state = to
        elif kind == "on_speciation":
            clado[node_id] = to
        else:
            switches.setdefault(node_id, []).append((float(time_s), to))
    if initial_state is None:
        raise ValueError(
            f"driver file {str(path)!r} has no 'initial' row, so the state at t=0 is unknown and the "
            "trajectory cannot be reconstructed. Re-write it with a current ZOMBI2 (the t=0 row used "
            "to be spelled 'root')."
        )
    if not switches and not clado:
        # A log with an initial row and nothing else is either a discrete trait that happened never to
        # switch — fine, a constant driver is the truth — or a CONTINUOUS trait, whose diffusion is not
        # in this file at all. Replaying the latter silently produced a driver frozen at the root value
        # for the whole tree, so a run that looked conditioned was the undriven model with a constant
        # factor. The sibling value table tells the two apart, and is the file a diffusion belongs in.
        values_file = pathlib.Path(path).with_name("trait_values.tsv")
        if values_file.exists() and _looks_continuous(values_file):
            raise ValueError(
                f"driver file {str(path)!r} is a CONTINUOUS trait's event log, which records no "
                f"switches — a diffusion has none — so it carries only the value at t=0 and would "
                f"drive every lineage at that one constant value.\n"
                f"Point the conditioned rate at {str(values_file)!r} instead: a continuous trait's "
                f"driver is its value table, and the path between nodes is interpolated from it."
            )
    return DriverTrajectory(_replay(tree, initial_state, clado, switches))


def _looks_continuous(values_file: pathlib.Path) -> bool:
    """Whether a ``trait_values.tsv`` holds numbers (a diffusion) rather than state labels.

    Read rather than guessed from the event log, because a discrete trait's states may themselves be
    spelled ``0`` and ``1`` and a never-switching discrete log is indistinguishable from a continuous
    one on its own."""
    try:
        rows = [r for r in values_file.read_text(encoding="utf-8").splitlines() if r.strip()]
    except OSError:
        return False
    if len(rows) < 2 or rows[0].split("\t") != ["node", "kind", "trait"]:
        return False
    for row in rows[1:]:
        parts = row.split("\t")
        if len(parts) != 3:
            return False
        try:
            float(parts[2])
        except ValueError:
            return False
    return True


def _load_values_driver(path, rows, tree, step) -> DriverTrajectory:
    """A continuous trait's ``trait_values.tsv`` → a `DriverTrajectory`, interpolated at ``step``."""
    values: dict[int, float] = {}
    for line in rows[1:]:
        parts = line.split("\t")
        if len(parts) != 3:
            raise ValueError(f"driver file {str(path)!r} row is not 3 columns: {line!r}")
        node_s, _kind, value_s = parts
        try:
            values[node_from_label(node_s)] = float(value_s)
        except ValueError:
            raise ValueError(
                f"driver file {str(path)!r} is a value table whose entries are not numbers "
                f"({value_s!r} on {node_s}). A DISCRETE trait's driver is its event log "
                f"(trait_events.tsv), which carries the switch times this table has lost; only a "
                f"CONTINUOUS trait is conditioned on its values."
            ) from None
    missing = [node_label(i) for i in tree.nodes if i not in values]
    if missing:
        raise ValueError(
            f"driver file {str(path)!r} has no value for {len(missing)} lineage(s) of the tree this "
            f"run is on (e.g. {', '.join(missing[:5])}). The driver must have been grown on the SAME "
            f"complete tree — including the lineages that went extinct."
        )
    return DriverTrajectory(interpolated_segments(tree, values, step))


def _replay(tree, initial_state, clado, switches) -> dict[int, list[tuple[float, object]]]:
    """Rebuild each lineage's constant stretches ``{node: [(start_time, state), …]}`` from the tree and
    the parsed log. Parent before child, so a lineage can read its parent's ending state."""
    segments: dict[int, list[tuple[float, object]]] = {}
    end_state: dict[int, object] = {}
    stack = [tree.root]
    while stack:                                        # pre-order: a parent is popped before its kids
        i = stack.pop()
        node = tree.nodes[i]
        if node.parent is None:
            start = initial_state
        elif i in clado:
            start = clado[i]
        else:
            start = end_state[node.parent]
        segs, t, state = [], node.birth_time, start
        for when, to in sorted(switches.get(i, ())):
            segs.append((t, state))
            t, state = when, to
        segs.append((t, state))
        segments[i] = segs
        end_state[i] = state
        if node.children is not None:
            stack.extend(node.children)
    return segments


#: default per-branch resolution for a continuous driver — how many constant stretches each branch is
#: cut into. Higher is a finer approximation of the continuous path (and more Gillespie breakpoints).
#: Default continuous-driver resolution, as a fraction of the tree's height: a stretch lasts 1% of
#: the run, wherever on the tree it sits. A fraction rather than an absolute duration because a tree
#: may be measured in expected substitutions or in millions of years, and a default in "time units"
#: would be meaninglessly fine on one and uselessly coarse on the other.
CONTINUOUS_DRIVER_FRACTION = 0.01


def driver_from_result(result, *, step: float | None = None) -> DriverTrajectory:
    """Build a `DriverTrajectory` **directly from a grown trait result** — the same per-lineage lookup
    `load_driver()` builds from a file, but skipping the file round-trip. This is how a conditioned
    ``DrivenBy(trait, …)`` reads a trait grown in the same Python session: still conditioning (the
    driver was grown first and is held fixed), just handed over in memory.

    A **discrete** trait (`traits.simulate_discrete`) has a stochastic character map, so each branch is
    cut into its exact constant segments. A **continuous** trait (`traits.simulate_continuous`) has no
    such map, so it is handled by `driver_from_continuous_result()` — a piecewise-constant
    approximation whose stretches last at most ``step`` time units."""
    tree = getattr(result, "complete_tree", None)
    history = getattr(result, "history", None)
    if tree is None:
        raise ValueError(
            "a conditioned driver object must be a grown trait result (from traits.simulate_discrete or "
            f"simulate_continuous), carrying its complete tree; got {type(result).__name__}.")
    if history is None:                                       # no character map -> continuous (or threshold)
        if getattr(result, "node_values", None) is not None:
            return driver_from_continuous_result(result, step=step)
        raise ValueError(
            "a conditioned driver object must be a DISCRETE trait result (with a stochastic character "
            "map) or a CONTINUOUS one (with per-node values); got "
            f"{type(result).__name__} with neither.")
    segments: dict[int, list[tuple[float, object]]] = {}
    for i, node in tree.nodes.items():
        t = node.birth_time
        segs: list[tuple[float, object]] = []
        for state, dur in history[i]:
            segs.append((t, state))
            t += dur
        segments[i] = segs
    return DriverTrajectory(segments)


def tree_height(tree) -> float:
    """Origin to present — what a ``step=None`` continuous driver takes its default resolution from."""
    root_birth = tree.nodes[tree.root].birth_time
    return max(n.end_time for n in tree.nodes.values()) - root_birth


def default_step(tree) -> float:
    """The default driver resolution for ``tree``: `CONTINUOUS_DRIVER_FRACTION` of its height."""
    height = tree_height(tree)
    if not (height > 0.0) or not math.isfinite(height):
        raise ValueError(f"cannot pick a default driver step for a tree of height {height!r}; "
                         f"pass step= explicitly")
    return height * CONTINUOUS_DRIVER_FRACTION


def interpolated_segments(tree, values: dict, step: float | None = None) -> dict:
    """Cut every branch into stretches of **at most ``step`` time units** and give each the trait
    linearly interpolated between the branch's endpoints, read at the stretch's midpoint.

    The cut is per unit of TIME, not per branch. Cutting each branch into a fixed number of pieces
    makes the approximation as coarse as the branch is long: a branch ten times another gets stretches
    ten times cruder, so the error is uneven across the tree and worst exactly where the driver has
    had most time to move. A fixed time step gives every stretch the same duration wherever it sits,
    so the resolution means the same thing everywhere and refining it is one number.

    A branch shorter than ``step`` gets a single stretch at its midpoint value. The number of
    stretches on a branch of length ``L`` is ``ceil(L / step)``, and they divide ``L`` evenly."""
    if step is None:
        step = default_step(tree)
    step = float(step)
    if not (step > 0.0) or not math.isfinite(step):
        raise ValueError(f"a driver step is a duration and must be finite and positive, got {step!r}")

    segments: dict[int, list[tuple[float, object]]] = {}
    for i, node in tree.nodes.items():
        end_v = float(values[i])
        start_v = float(values[node.parent]) if node.parent is not None else end_v  # root: no earlier value
        t0 = node.birth_time
        dt = node.end_time - node.birth_time
        if dt <= 0:
            segments[i] = [(t0, end_v)]
            continue
        n = max(1, math.ceil(dt / step))
        segments[i] = [(t0 + k * dt / n, start_v + (end_v - start_v) * (k + 0.5) / n)
                       for k in range(n)]
    return segments


def driver_from_continuous_result(result, *, step: float | None = None) -> DriverTrajectory:
    """Build a `DriverTrajectory` from a **continuous** trait result (`traits.simulate_continuous`).

    A continuous trait has no discrete switches, so there is no exact stochastic map to replay.
    Instead the path is approximated by a piecewise-constant one — the same shape a discrete trait's
    map has, which is what lets one engine consume both — with each stretch lasting at most ``step``
    time units (`interpolated_segments`). The engine then steps its Gillespie to every stretch
    boundary, so within a stretch the rate really is constant and the exponential draw is exact
    there; nothing is thinned or rejected.

    ``step`` is the resolution, in the tree's own time units, and it is the knob that trades accuracy
    for speed: halving it doubles the breakpoints and so the work. ``None`` takes
    `CONTINUOUS_DRIVER_FRACTION` of the tree's height.

    Two things this approximation is not, worth knowing before reading a driven run. It is the
    **straight line** between a branch's endpoint values, so it is the mean of the Brownian bridge
    between them with the excursions dropped — a real path wanders either side of that line and this
    driver does not. And under a non-linear response curve those dropped excursions do not average
    out, so a smaller ``step`` is not only more precise, it removes a bias.

    The driver's values are **floats**, so its `DrivenBy` needs a continuous mapping (a
    `~zombi2.rates.mapping.Curve` or `~zombi2.rates.mapping.Scalar`), not a discrete
    ``{state: factor}`` `~zombi2.rates.mapping.Table` (which would match no float and never fire)."""
    tree = result.complete_tree
    values = result.node_values
    if not values:
        raise ValueError("continuous driver result has no node values to interpolate")
    sample = next(iter(values.values()))
    if isinstance(sample, dict):
        raise ValueError(
            "a continuous driver must be a SINGLE-trait result; got multi-trait node values "
            f"({sorted(sample)}). Grow one trait as the driver, or select a component before conditioning.")
    return DriverTrajectory(interpolated_segments(tree, values, step))


def check_mapping_fires(mapping, available_states, *, source_label: str, exhaustive: bool = False) -> None:
    """Raise if a **discrete** (`Table`) mapping's states do not line up
    with the states the driver can take. Such a mismatch leaves lineages at the table's default factor —
    a rate that is never touched — so the run drifts from the model the log records. It is almost always
    a typo or a stale / mismatched driver, so it is refused. Continuous mappings (Curve / Scalar) apply
    to every value and have nothing to mismatch.

    ``available_states`` is the set the mapping is checked against, and ``exhaustive`` says what that set
    *is*:

    - ``exhaustive=False`` (the default) — ``available_states`` are the states the driver was **observed**
      to take (e.g. replayed from a written trait file). At least **one** named state must occur, but a
      mapping may still list a state this particular realisation never reached (a legitimate partial
      mapping), so only an *empty* overlap is an error. A key that matches nothing still **warns**:
      it is far more often a typo than a deliberate partial mapping, and the failure it causes is the
      worst kind — the run completes, reports that it was driven, and applies the driver to only some
      of the lineages the user meant. Silence there is what lets a wrong result be published.
    - ``exhaustive=True`` — ``available_states`` is the driver's **complete declared alphabet**, known up
      front (e.g. a joint trait's declared states). Then **every** named state must be one of them: a key
      outside the alphabet is a state that can never occur, so its factor could never apply — an
      unambiguous typo, refused even when other keys do match."""
    from .mapping import Table

    if not isinstance(mapping, Table):
        return
    if available_states and all(isinstance(s, (int, float)) and not isinstance(s, bool)
                                for s in available_states):
        raise ValueError(
            f"DrivenBy on {source_label}: the driver is CONTINUOUS (its values are numbers), so its "
            "mapping must be a Curve (value -> factor) or a Scalar (a log-link), not a {state: factor} "
            "table — a table names discrete states, which a continuous value never equals.")
    named = set(mapping.per_state)
    have = {str(s) for s in available_states}
    if exhaustive:
        stray = named - have
        if stray:
            raise ValueError(
                f"DrivenBy on {source_label}: the mapping names state(s) {sorted(stray)} that are not "
                f"among the driver's states {sorted(have)} — a factor for a state that can never occur, "
                f"so it would silently never apply. Check for a typo in the state names.")
        return
    if not (named & have):
        raise ValueError(
            f"DrivenBy on {source_label}: the mapping's states {sorted(named)} match none of the "
            f"driver's states {sorted(have)}, so the mapping would silently do nothing — every "
            f"lineage falls to the default factor and the rate is never driven. Check for a typo in "
            f"the state names, or a stale or mismatched driver file.")
    stray = named - have
    if stray:
        # Some keys matched, so the rate IS driven and the run is a legitimate model — which is why
        # this warns rather than raises. But a key matching nothing is far more often a typo than a
        # deliberate partial mapping, and the shape of that failure is the dangerous one: the run
        # completes, the summary says it was driven, and the factor the user cared about was applied
        # to nobody. Saying so costs one line; not saying it is how a wrong result gets published.
        warnings.warn(
            f"DrivenBy on {source_label}: the mapping names state(s) {sorted(stray)} that the "
            f"driver never takes (it takes {sorted(have)}), so those factors were never applied. "
            f"The states it did match are still driving the rate. Check for a typo — this is a "
            f"warning rather than an error only because a mapping may legitimately name a state "
            f"that this particular run did not reach.",
            stacklevel=2)


def driven_mods(rate) -> list:
    """The `DrivenBy` modifiers a rate carries, or ``[]`` when it carries none. A non-empty list means
    the rate reads another level on each lineage, so the engine must thread a ``drivers`` value and
    step where the driver switches.

    It lives here rather than in one level's package because every level that reads a driver asks the
    same question first, and the answer is a fact about the rate grammar (SPEC §5), not about any one
    engine. (The trait and genome engines still each carry a private copy from before this existed;
    folding them onto this one is a tidy-up, not a behaviour change.)"""
    return [m for m in rate.modifiers if isinstance(m, DrivenBy)]


def names_a_live_level(source: object) -> bool:
    """Whether a ``DrivenBy`` ``source`` names a **level growing beside the run** rather than a
    finished driver.

    SPEC §5: "a finished result makes the run conditioned, and the name of a level growing beside it
    makes the run joint". One modifier, one spelling, and the *source* is what tells the two apart —
    so this is the predicate that reads the source, not a judgement about what the target level then
    does with it. The live names are the ones `zombi2.joint` accepts: ``"trait"``, ``"genomes:count"``
    and ``"genomes:<family>"``.

    A level that cannot be joined with the driver at all (Traits–Sequences, SPEC §3) uses this to say
    so in the modelling terms, instead of letting the string fall through to `load_driver()` and come
    back as a missing file called ``'trait'``."""
    return isinstance(source, str) and (source == "trait" or source.startswith("genomes:"))


def refuse_wrong_direction(source, level: str | None) -> None:
    """Raise when ``level``'s engine may not read ``source`` at all — a direction SPEC §3 rules out,
    not something merely unimplemented.

    A driver declares this by answering ``refuses(level)`` with the reason, or ``None``. It lives on
    the driver for the same reason `resolve_driver()`'s protocol does: this module serves every level
    and must not import from any of them.

    ``level`` is the engine name each level already answers to (``"genomes.family"``, ``"sequences"``,
    … — `zombi2.rates.modifiers.Modifier.implemented_for` lists them); ``None`` checks nothing."""
    if level is None:
        return
    refuses = getattr(source, "refuses", None)
    if refuses is None:
        return
    why = refuses(level)
    if why:
        raise ValueError(why)


def resolve_driver(source, tree, *, step: float | None = None,
                   level: str | None = None) -> DriverTrajectory:
    """Resolve a conditioned ``DrivenBy`` ``source`` into a `DriverTrajectory` — a **filename**
    (str) via `load_driver()` (replayed against ``tree``, the target run's own species tree), an
    object that answers for itself through ``as_driver_trajectory(tree, step=…)`` (a genome run's
    ``presence("name")``, a sequence run's ``gc()``), or an **in-memory** trait result via
    `driver_from_result()` (which carries its own tree).
    Both are conditioning (the driver grown first); the object form just spares you the ``write``/read
    step in a single session.

    ``step`` is the continuous-driver resolution (see `interpolated_segments`); it is ignored by a
    discrete driver, whose stretches are exact. ``level`` names the engine doing the reading, which is
    what lets a driver refuse a level that sits above it (`refuse_wrong_direction`)."""
    if isinstance(source, str):
        return load_driver(source, tree, step=step)
    if hasattr(source, "as_driver_trajectory"):
        # a level that knows how to answer "what state was lineage L in at time t?" for itself —
        # `genomes.presence("tox")` is the first. The protocol is one method rather than an isinstance
        # branch per level so this module stays free of imports from the levels it serves.
        refuse_wrong_direction(source, level)
        return source.as_driver_trajectory(tree, step=step)
    return driver_from_result(source, step=step)


__all__ = ["DriverTrajectory", "load_driver", "driver_from_result", "resolve_driver",
           "refuse_wrong_direction", "check_mapping_fires", "driven_mods", "names_a_live_level"]
