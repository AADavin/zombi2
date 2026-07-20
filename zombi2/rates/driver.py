"""A conditioned :class:`~zombi2.rates.modifiers.DrivenBy`'s file-backing (SPEC §2, ``coupling-api.md``).

When ``DrivenBy``'s ``source`` is a **filename**, the coupling is *conditioned*: the driver was grown
first and written to a file, and two ordinary runs in order do the rest
(``loss = 0.25 * mod.DrivenBy("habitat.tsv", {...})``). This module — living beside ``DrivenBy`` in
``rates`` because it is that modifier's file end — turns the written driver into the per-lineage lookup
the target engine queries as it walks the (already-grown) tree. (Conditioning needs no engine of its
own: it *folds into the target level's* run; only genuinely-joint models get a dedicated engine,
``zombi2.joint``.)

The driver file is a **segment table** (``trait_driver.tsv``, written by
:meth:`zombi2.traits.TraitsResult.write` with ``outputs=("driver",)``): one row per constant
stretch of a lineage's branch, ``node · start · end · state``. A discrete driver switches
*mid-branch*, so the table is the exact stochastic character map, not one value per branch —
:class:`DriverTrajectory` answers both *what is the driver on this lineage now?*
(:meth:`~DriverTrajectory.value`) and *when does it next change?*
(:meth:`~DriverTrajectory.next_change`, so the target's Gillespie steps at each switch).

The join key is the **species node id** — the driver ran on the same complete tree the target now
runs on, so ``node n7`` in the file is lineage 7 in the target run.
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


def load_driver(path) -> DriverTrajectory:
    """Read a driver **segment table** (``node · start · end · state``, tab-separated, one header
    row) into a :class:`DriverTrajectory`. Node ids are written ``n<id>`` (matching every other
    ZOMBI2 tip label) and read back to the integer lineage id. States are kept as written (a
    discrete label such as ``aquatic``). This is the file :class:`~zombi2.rates.modifiers.DrivenBy`
    names as its ``source`` when the driver was grown first (conditioning)."""
    text = pathlib.Path(path).read_text()
    rows = [line for line in text.splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"driver file {str(path)!r} is empty")
    header = rows[0].split("\t")
    expected = ["node", "start", "end", "state"]
    if header != expected:
        raise ValueError(
            f"driver file {str(path)!r} must have header {expected}, got {header} — write it with "
            f"TraitsResult.write(dir, outputs=('driver',))."
        )
    segments: dict[int, list[tuple[float, object]]] = {}
    for line in rows[1:]:
        parts = line.split("\t")
        if len(parts) != 4:
            raise ValueError(f"driver file {str(path)!r} row is not 4 columns: {line!r}")
        node_s, start_s, _end_s, state = parts
        node_id = int(node_s[1:]) if node_s.startswith("n") else int(node_s)
        segments.setdefault(node_id, []).append((float(start_s), state))
    return DriverTrajectory(segments)


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


def resolve_driver(source) -> DriverTrajectory:
    """Resolve a conditioned ``DrivenBy`` ``source`` into a :class:`DriverTrajectory` — a **filename**
    (str) via :func:`load_driver`, or an **in-memory** discrete trait result via
    :func:`driver_from_result`. Both are conditioning (the driver grown first); the object form just
    spares you the ``write``/read step in a single session."""
    if isinstance(source, str):
        return load_driver(source)
    return driver_from_result(source)


__all__ = ["DriverTrajectory", "load_driver", "driver_from_result", "resolve_driver"]
