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

import bisect
import math
import pathlib

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


def load_driver(path, tree) -> DriverTrajectory:
    """Read a trait **event log** (``trait_events.tsv``: ``time · kind · lineage · from · to``, an
    ``initial`` row then the switches) and **replay it against ``tree``** into a `DriverTrajectory`.

    The log alone is not enough — a switch says *when* the state changed, not what each branch started
    in — so the tree supplies branch birth/end times and the topology, and the reconstruction walks
    parent-before-child: the root begins in the ``initial`` row's state, every other lineage in its own
    ``on_speciation`` state if it has one else its parent's ending state, and ``on_branch`` rows cut
    the branch into constant stretches. This is the same tree the target level runs on, so ``node n7``
    in the log is lineage 7 here. (``tree`` is the run's own species tree, always in hand where a
    conditioned rate is resolved.)"""
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
    if header != expected:
        raise ValueError(
            f"driver file {str(path)!r} must be a trait event log with header {expected}, got "
            f"{header} — write it with TraitsResult.write(dir, outputs=('events',)). (The old "
            "node·start·end·state driver table was retired: the event log is the driver now.)"
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
    return DriverTrajectory(_replay(tree, initial_state, clado, switches))


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
CONTINUOUS_DRIVER_STEPS = 8


def driver_from_result(result, *, steps: int = CONTINUOUS_DRIVER_STEPS) -> DriverTrajectory:
    """Build a `DriverTrajectory` **directly from a grown trait result** — the same per-lineage lookup
    `load_driver()` builds from a file, but skipping the file round-trip. This is how a conditioned
    ``DrivenBy(trait, …)`` reads a trait grown in the same Python session: still conditioning (the
    driver was grown first and is held fixed), just handed over in memory.

    A **discrete** trait (`traits.simulate_discrete`) has a stochastic character map, so each branch is
    cut into its exact constant segments. A **continuous** trait (`traits.simulate_continuous`) has no
    such map, so it is handled by `driver_from_continuous_result()` — a piecewise-constant
    approximation (``steps`` per branch)."""
    tree = getattr(result, "complete_tree", None)
    history = getattr(result, "history", None)
    if tree is None:
        raise ValueError(
            "a conditioned driver object must be a grown trait result (from traits.simulate_discrete or "
            f"simulate_continuous), carrying its complete tree; got {type(result).__name__}.")
    if history is None:                                       # no character map -> continuous (or threshold)
        if getattr(result, "node_values", None) is not None:
            return driver_from_continuous_result(result, steps=steps)
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


def driver_from_continuous_result(result, *, steps: int = CONTINUOUS_DRIVER_STEPS) -> DriverTrajectory:
    """Build a `DriverTrajectory` from a **continuous** trait result (`traits.simulate_continuous`).

    A continuous trait has no discrete switches, so there is no exact stochastic map. Instead each
    branch is cut into ``steps`` equal stretches whose value is the trait **linearly interpolated** from
    the branch's start (its parent's node value) to its end (this node's value), sampled at each
    stretch's midpoint — a piecewise-constant approximation of the continuous path that the same engine
    consumes, converging as ``steps`` grows. It adds ``steps`` breakpoints per branch, so a
    continuously-driven run steps its Gillespie more often than a discrete one.

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
    steps = max(1, int(steps))
    segments: dict[int, list[tuple[float, object]]] = {}
    for i, node in tree.nodes.items():
        end_v = float(values[i])
        start_v = float(values[node.parent]) if node.parent is not None else end_v  # root: no earlier value
        t0 = node.birth_time
        dt = node.end_time - node.birth_time
        if dt <= 0:
            segments[i] = [(t0, end_v)]
            continue
        segments[i] = [(t0 + k * dt / steps, start_v + (end_v - start_v) * (k + 0.5) / steps)
                       for k in range(steps)]
    return DriverTrajectory(segments)


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
      mapping), so only an *empty* overlap is an error.
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


def resolve_driver(source, tree) -> DriverTrajectory:
    """Resolve a conditioned ``DrivenBy`` ``source`` into a `DriverTrajectory` — a **filename**
    (str) via `load_driver()` (replayed against ``tree``, the target run's own species tree), or an
    **in-memory** discrete trait result via `driver_from_result()` (which carries its own tree).
    Both are conditioning (the driver grown first); the object form just spares you the ``write``/read
    step in a single session."""
    if isinstance(source, str):
        return load_driver(source, tree)
    return driver_from_result(source)


__all__ = ["DriverTrajectory", "load_driver", "driver_from_result", "resolve_driver",
           "check_mapping_fires"]
