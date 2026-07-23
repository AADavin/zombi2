"""A conditioned :class:`~zombi2.rates.modifiers.DrivenBy`'s file-backing (SPEC §2, ``coupling-api.md``).

When ``DrivenBy``'s ``source`` is a **filename**, the coupling is *conditioned*: the driver was grown
first and written to a file, and two ordinary runs in order do the rest
(``loss = 0.25 * mod.DrivenBy("habitat.tsv", {...})``). This module — living beside ``DrivenBy`` in
``rates`` because it is that modifier's file end — turns the written driver into the per-lineage lookup
the target engine queries as it walks the (already-grown) tree. (Conditioning needs no engine of its
own: it *folds into the target level's* run; only genuinely-joint models get a dedicated engine,
``zombi2.joint``.)

The driver file is the trait **event log** (``trait_events.tsv``, written by
:meth:`zombi2.traits.TraitsResult.write` with ``outputs=("events",)``): a ``root`` row giving the
initial state, then every switch — ``time · kind · lineage · from · to``. The driver ran on the same
complete tree the target now runs on, so replaying the log **against that tree** rebuilds each
lineage's branch as constant stretches (a discrete driver switches *mid-branch*, so this is the exact
stochastic character map, not one value per branch). :class:`DriverTrajectory` then answers both
*what is the driver on this lineage now?* (:meth:`~DriverTrajectory.value`) and *when does it next
change?* (:meth:`~DriverTrajectory.next_change`, so the target's Gillespie steps at each switch).

The join key is the **species node id**: ``node n7`` in the log is lineage 7 in the target run.
"""

from __future__ import annotations

import bisect
import math
import pathlib


class DriverTrajectory:
    """A driver's value along every lineage, as a piecewise-constant function of time — the
    per-lineage lookup a conditioned :class:`~zombi2.rates.modifiers.DrivenBy` reads.

    Built from segments ``{node_id: [(start_time, state), …]}`` (each lineage's branch cut into
    constant stretches, sorted by start). The engine calls :meth:`value` to get a lineage's driver
    state at the current instant and :meth:`next_change` to learn when it next switches (a horizon
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
                f"the driver file has no lineage n{node_id}; the driver must be grown on the SAME "
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
            raise KeyError(f"the driver file has no lineage n{node_id} (node ids must match the target tree).")
        i = bisect.bisect_right(starts, time)
        return starts[i] if i < len(starts) else math.inf


def load_driver(path, tree) -> DriverTrajectory:
    """Read a trait **event log** (``trait_events.tsv``: ``time · kind · lineage · from · to``, a
    ``root`` row then the switches) and **replay it against ``tree``** into a :class:`DriverTrajectory`.

    The log alone is not enough — a switch says *when* the state changed, not what each branch started
    in — so the tree supplies branch birth/end times and the topology, and the reconstruction walks
    parent-before-child: the root begins in the ``root`` row's state, every other lineage in its own
    ``on_speciation`` state if it has one else its parent's ending state, and ``on_branch`` rows cut
    the branch into constant stretches. This is the same tree the target level runs on, so ``node n7``
    in the log is lineage 7 here. (``tree`` is the run's own species tree, always in hand where a
    conditioned rate is resolved.)"""
    text = pathlib.Path(path).read_text()
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
    root_state = None
    clado: dict[int, object] = {}                       # lineage -> its on-speciation start state
    switches: dict[int, list[tuple[float, object]]] = {}   # lineage -> [(time, to_state), …]
    for line in rows[1:]:
        parts = line.split("\t")
        if len(parts) != 5:
            raise ValueError(f"driver file {str(path)!r} row is not 5 columns: {line!r}")
        time_s, kind, node_s, _from, to = parts
        node_id = int(node_s[1:]) if node_s.startswith("n") else int(node_s)
        if kind == "root":
            root_state = to
        elif kind == "on_speciation":
            clado[node_id] = to
        else:
            switches.setdefault(node_id, []).append((float(time_s), to))
    if root_state is None:
        raise ValueError(
            f"driver file {str(path)!r} has no 'root' row, so the initial state is unknown and the "
            "trajectory cannot be reconstructed. Re-write it with a current ZOMBI2."
        )
    return DriverTrajectory(_replay(tree, root_state, clado, switches))


def _replay(tree, root_state, clado, switches) -> dict[int, list[tuple[float, object]]]:
    """Rebuild each lineage's constant stretches ``{node: [(start_time, state), …]}`` from the tree and
    the parsed log. Parent before child, so a lineage can read its parent's ending state."""
    segments: dict[int, list[tuple[float, object]]] = {}
    end_state: dict[int, object] = {}
    stack = [tree.root]
    while stack:                                        # pre-order: a parent is popped before its kids
        i = stack.pop()
        node = tree.nodes[i]
        if node.parent is None:
            start = root_state
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


def driver_from_result(result) -> DriverTrajectory:
    """Build a :class:`DriverTrajectory` **directly from a discrete trait result** — the same
    per-lineage lookup :func:`load_driver` builds from a file, but skipping the file round-trip. This
    is how a conditioned ``DrivenBy(habitat, …)`` reads a trait grown in the same Python session: still
    conditioning (the driver was grown first and is held fixed), just handed over in memory rather than
    written out. Needs a **discrete** trait (its stochastic character map cuts each branch into the
    constant segments a driver needs); a continuous / threshold trait has no such map."""
    history = getattr(result, "history", None)
    tree = getattr(result, "complete_tree", None)
    if history is None or tree is None:
        raise ValueError(
            "a conditioned driver object must be a DISCRETE trait result (from traits.simulate_discrete), "
            "whose stochastic character map cuts each branch into constant segments; got "
            f"{type(result).__name__} with no such map. Driving with a continuous trait is a later slice."
        )
    segments: dict[int, list[tuple[float, object]]] = {}
    for i, node in tree.nodes.items():
        t = node.birth_time
        segs: list[tuple[float, object]] = []
        for state, dur in history[i]:
            segs.append((t, state))
            t += dur
        segments[i] = segs
    return DriverTrajectory(segments)


def check_mapping_fires(mapping, available_states, *, source_label: str) -> None:
    """Raise if a **discrete** (:class:`~zombi2.rates.mapping.Table`) mapping names none of the states
    the driver can actually take. Such a mapping leaves every lineage at the table's default factor —
    a rate that is never touched — so the run is the fully *uncoupled* model while the log records it as
    driven. That is almost always a typo or a stale / mismatched driver file, so it is refused.

    At least **one** named state must occur; a mapping may still list a state this particular
    realisation never reached (a legitimate partial mapping), so only an *empty* overlap is an error.
    Continuous mappings (Curve / Scalar) apply to every value and have nothing to mismatch."""
    from .mapping import Table

    if not isinstance(mapping, Table):
        return
    named = set(mapping.per_state)
    have = {str(s) for s in available_states}
    if not (named & have):
        raise ValueError(
            f"DrivenBy on {source_label}: the mapping's states {sorted(named)} match none of the "
            f"driver's states {sorted(have)}, so the coupling would silently do nothing — every "
            f"lineage falls to the default factor and the rate is never driven. Check for a typo in "
            f"the state names, or a stale or mismatched driver file.")


def resolve_driver(source, tree) -> DriverTrajectory:
    """Resolve a conditioned ``DrivenBy`` ``source`` into a :class:`DriverTrajectory` — a **filename**
    (str) via :func:`load_driver` (replayed against ``tree``, the target run's own species tree), or an
    **in-memory** discrete trait result via :func:`driver_from_result` (which carries its own tree).
    Both are conditioning (the driver grown first); the object form just spares you the ``write``/read
    step in a single session."""
    if isinstance(source, str):
        return load_driver(source, tree)
    return driver_from_result(source)


__all__ = ["DriverTrajectory", "load_driver", "driver_from_result", "resolve_driver",
           "check_mapping_fires"]
