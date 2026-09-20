"""A conditioned `Driven`'s file-backing (SPEC §2).

When ``Driven``'s ``driver`` is a **filename**, the relation is *conditioned*: the driver was grown
first and written to a file, and two ordinary runs in order do the rest
(``loss = PerCopy(0.25).scaled_by("habitat.tsv", {...})``). This module — living beside `Driven` in
``params`` because it is that modifier's file end — turns the written driver into the per-lineage lookup
the driven engine queries as it walks the (already-grown) tree. (Conditioning needs no engine of its
own: it *folds into the driven level's* run; only genuinely-joint models get a dedicated engine,
``zombi2.joint``.)

Two files can be a driver. The usual one is the trait **event log** (``trait_events.tsv``, written by
`zombi2.traits.TraitsResult.write()` with ``outputs=("events",)``); a continuous trait's value table
(``trait_values.tsv``) is the other, and `load_driver()` dispatches on the header. The event log
holds an ``initial`` row giving the state at t=0, then every switch — ``time · kind · lineage · from · to``. The driver ran on the same
complete tree the driven run now walks, so replaying the log **against that tree** rebuilds each
lineage's branch as constant stretches (a discrete driver switches *mid-branch*, so this is the exact
stochastic character map, not one value per branch). `DriverTrajectory` then answers both
*what is the driver on this lineage now?* (`value()`) and *when does it next
change?* (`next_change()`, so the driven Gillespie steps at each switch).

The join key is the **species node id**: ``node n7`` in the log is lineage 7 in the driven run.
"""

from __future__ import annotations

import warnings

import bisect
import math
import pathlib
import zlib
from typing import cast

from .connection import Driven
from ..rng import branch_stream
from ..tree import node_from_label, node_label


class DriverTrajectory:
    """A driver's value along every lineage, as a piecewise-constant function of time — the
    per-lineage lookup a conditioned `Driven` reads.

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
                f"complete tree the driven run walks (node ids must match)."
            )
        i = bisect.bisect_right(starts, time) - 1
        if i < 0:  # a query before the branch's first segment: clamp to the first (branch-start) state
            i = 0
        return self._states[node_id][i]

    def next_change(self, node_id: int, time: float) -> float:
        """The next time strictly after ``time`` at which lineage ``node_id``'s driver switches, else
        ``inf`` (it stays constant for the rest of the branch). Feeds the driven Gillespie's horizon."""
        starts = self._starts.get(node_id)
        if starts is None:
            raise KeyError(f"the driver file has no lineage {node_label(node_id)} "
                           "(node ids must match the driven run's tree).")
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
            f"driver file not found: {str(path)!r}. A conditioned rate points at this file, but "
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


def _load_path_file(values_file: pathlib.Path):
    """The written path beside a ``trait_values.tsv``, as a `~zombi2.traits.path.WrittenPath`, the
    branch start values it carries, and a seed — or ``(None, None, None)`` when there is no such file.

    ``trait_path.tsv`` is ``node · time · trait · variance``, one row per point, and every branch's
    two endpoints are among its rows. Reading it is what makes a run driven from files take the same
    path a run driven in memory takes; without it the only thing a file carries is the value at each
    node, and the reader has to fall back to the straight line between them."""
    from ..traits.path import WrittenPath

    path_file = values_file.with_name("trait_path.tsv")
    try:
        text = path_file.read_text(encoding="utf-8")
    except OSError:
        return None, None, None
    rows = [r for r in text.splitlines() if r.strip()]
    if not rows or rows[0].split("\t") != ["node", "time", "trait", "variance"]:
        raise ValueError(
            f"{str(path_file)!r} must be a trait path table with header "
            f"['node', 'time', 'trait', 'variance']; got {rows[0].split(chr(9)) if rows else []}. "
            f"Write one with TraitsResult.write(dir) — a continuous trait writes it by default.")
    points: dict[int, list[tuple[float, float, float]]] = {}
    for line in rows[1:]:
        parts = line.split("\t")
        if len(parts) != 4:
            raise ValueError(f"{str(path_file)!r} row is not 4 columns: {line!r}")
        node_s, time_s, value_s, var_s = parts
        try:
            point = (float(time_s), float(value_s), float(var_s))
        except ValueError:
            raise ValueError(
                f"{str(path_file)!r} row has a non-numeric time, value or variance: {line!r}"
            ) from None
        points.setdefault(node_from_label(node_s), []).append(point)
    for segs in points.values():
        segs.sort()
    starts = {i: segs[0][1] for i, segs in points.items() if segs}
    return WrittenPath(points), starts, zlib.crc32(text.encode("utf-8"))


def _load_values_driver(path, rows, tree, step) -> DriverTrajectory:
    """A continuous trait's ``trait_values.tsv`` → a `DriverTrajectory`, read at ``step``.

    The values table gives the trait at each node. What it happened to do **between** them is in
    ``trait_path.tsv`` beside it, written by default since #454, and this reads that file when it is
    there — the same sideways look `load_driver` already does for the values table itself. Without
    it the reader can only take the straight line between node values, which drops the path's
    excursions and biases every non-linear mapping, so it says so rather than doing it quietly."""
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
    law, starts, path_seed = _load_path_file(pathlib.Path(path))
    if law is None:
        warnings.warn(
            f"no trait_path.tsv beside {str(path)!r}, so this continuous driver is read as the "
            f"STRAIGHT LINE between node values. That line is the mean of the trait's path with "
            f"the excursions dropped, and under a non-linear mapping they do not average out, so a "
            f"smaller step does not remove the bias. Re-write the driver with a current ZOMBI2 "
            f"(a continuous trait writes trait_path.tsv by default) to read its real path.",
            RuntimeWarning, stacklevel=3)
    return DriverTrajectory(interpolated_segments(tree, values, step, law=law,
                                                  seed=path_seed,
                                                  starts=starts))


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
        if node.children:
            stack.extend(node.children)
    return segments


#: Default continuous-driver resolution, as a fraction of the tree's height: a stretch lasts 1% of
#: the run, wherever on the tree it sits. A fraction rather than an absolute duration because a tree
#: may be measured in expected substitutions or in millions of years, and a default in "time units"
#: would be meaninglessly fine on one and uselessly coarse on the other.
CONTINUOUS_DRIVER_FRACTION = 0.01


def driver_from_result(result, *, step: float | None = None, path: bool = True) -> DriverTrajectory:
    """Build a `DriverTrajectory` **directly from a grown trait result** — the same per-lineage lookup
    `load_driver()` builds from a file, but skipping the file round-trip. This is how a conditioned
    ``scaled_by(trait, …)`` reads a trait grown in the same Python session: still conditioning (the
    driver was grown first and is held fixed), just handed over in memory.

    A **discrete** trait (`traits.simulate_discrete`) has a stochastic character map, so each branch is
    cut into its exact constant segments. A **continuous** trait (`traits.simulate_continuous`) has no
    such map, so it is handled by `driver_from_continuous_result()` — a piecewise-constant reading
    whose stretches last at most ``step`` time units. ``path`` is passed on there; a discrete
    trait's map is exact and ignores it."""
    tree = getattr(result, "complete_tree", None)
    history = getattr(result, "history", None)
    if tree is None:
        raise ValueError(
            "a conditioned driver object must be a grown trait result (from traits.simulate_discrete or "
            f"simulate_continuous), carrying its complete tree; got {type(result).__name__}.")
    if history is None:                                       # no character map -> continuous (or threshold)
        if getattr(result, "node_values", None) is not None:
            return driver_from_continuous_result(result, step=step, path=path)
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


def branch_starts(result) -> dict:
    """Each branch's **left endpoint** — the trait when the branch begins — read off a continuous
    trait result's event log.

    Not the same number as the parent's node value. A run with ``at_speciation`` adds a jump at the
    split before the branch diffuses, so the branch starts from the POST-jump value; the log records
    that jump as an ``on_speciation`` row whose ``to_state`` is it. Reading the parent's value
    instead made every branch of such a run start from the wrong place, which pulled the whole
    within-branch reading toward the parent. The root's left endpoint is the ``initial`` row, the
    trait at t=0, which is not in ``node_values`` at all (``node_values[root]`` is the value at the
    first split).

    A branch with neither row is not in the result, and the caller falls back to the parent's node
    value."""
    starts: dict[int, float] = {}
    for e in getattr(result, "events", None) or []:
        if e.kind in ("initial", "on_speciation"):
            try:
                starts[int(e.lineage)] = float(cast("float", e.to_state))
            except (TypeError, ValueError):
                return {}               # a discrete trait's log: its states are not numbers
    return starts


def path_points(tree, values: dict, step: float | None = None, *,
                law=None, seed: int | None = None, starts: dict | None = None) -> dict:
    """Every branch's path as points in time: ``{node: [(time, value, variance), …]}``.

    A branch is cut into stretches of **at most ``step`` time units** and the trait is read at each
    stretch's midpoint. The list a branch gets runs left endpoint → midpoints → node value, so it
    stands on its own: the branch's two ends are in it, and nothing outside is needed to read
    anywhere on it. ``variance`` on a point is what the trait accrued since the point before it,
    which is what lets a later reader refine between two of them
    (`zombi2.traits.path.WrittenPath`).

    ``law`` decides what a midpoint value **is**, and there are two answers:

    - With a `~zombi2.traits.path.BridgeLaw` (the usual case for a continuous trait grown in this
      session), the trait's real path is **drawn**, conditioned on both of the branch's endpoint
      values. The excursions either side of the straight line are then in the reading, which is what
      a non-linear mapping needs: dropping them biases the answer, and a smaller ``step`` does not
      remove that bias.
    - With ``law=None``, the value is the straight line between the endpoints — the mean of that
      bridge with the excursions dropped, and every variance zero.

    ``seed`` is the run's seed, and each branch's path is drawn from its own sub-stream off it
    (`zombi2.rng.branch_stream`), so a branch's path does not depend on the traversal order, on how
    many other branches there are, or on which reader asked first.

    ``starts`` gives each branch's left endpoint where it is known (`branch_starts`); a branch
    missing from it falls back to the parent's node value."""
    if step is None:
        step = default_step(tree)
    step = float(step)
    if not (step > 0.0) or not math.isfinite(step):
        raise ValueError(f"a driver step is a duration and must be finite and positive, got {step!r}")
    if law is not None and seed is None:
        law = None                      # nothing to key a reproducible path on; take the line

    points: dict[int, list[tuple[float, float, float]]] = {}
    for i, node in tree.nodes.items():
        end_v = float(values[i])
        if starts is not None and i in starts:
            start_v = float(starts[i])
        elif node.parent is not None:
            start_v = float(values[node.parent])
        else:
            start_v = end_v             # root with no recorded initial value: nothing earlier to read
        t0, t1 = node.birth_time, node.end_time
        dt = t1 - t0
        if dt <= 0:
            points[i] = [(t0, end_v, 0.0)]
            continue
        n = max(1, math.ceil(dt / step))
        mids = [t0 + (k + 0.5) * dt / n for k in range(n)]
        if law is None:
            drawn = [(start_v + (end_v - start_v) * (k + 0.5) / n, 0.0) for k in range(n)]
            tail = 0.0
        else:
            assert seed is not None
            drawn = law.sample(i, t0, t1, start_v, end_v, mids,
                               branch_stream("trait_paths", seed, i))
            tail = law.variance(i, mids[-1], t1)
        row = [(t0, start_v, 0.0)]
        row += [(mids[k], drawn[k][0], drawn[k][1]) for k in range(n)]
        row.append((t1, end_v, tail))
        points[i] = row
    return points


def interpolated_segments(tree, values: dict, step: float | None = None, *,
                          law=None, seed: int | None = None, starts: dict | None = None) -> dict:
    """Cut every branch into stretches of **at most ``step`` time units** and give each the trait
    read at the stretch's midpoint — the piecewise-constant shape a `DriverTrajectory` is built from.

    The cut is per unit of TIME, not per branch. Cutting each branch into a fixed number of pieces
    makes the reading as coarse as the branch is long: a branch ten times another gets stretches ten
    times cruder, so the error is uneven across the tree and worst exactly where the driver has had
    most time to move. A fixed time step gives every stretch the same duration wherever it sits, so
    the resolution means the same thing everywhere and refining it is one number.

    A branch shorter than ``step`` gets a single stretch at its midpoint value. The number of
    stretches on a branch of length ``L`` is ``ceil(L / step)``, and they divide ``L`` evenly.

    ``law``, ``seed`` and ``starts`` are `path_points`': the midpoint values are its drawn path, or
    the straight line between the branch's endpoints when there is no law. A **file-backed** driver
    passes a `~zombi2.traits.path.WrittenPath` here, so it reads the path the run wrote; a written
    run with no ``trait_path.tsv`` beside it has no law and takes the line."""
    points = path_points(tree, values, step, law=law, seed=seed, starts=starts)
    segments: dict[int, list[tuple[float, object]]] = {}
    for i, node in tree.nodes.items():
        row = points[i]
        if len(row) == 1:                              # a branch with no length
            segments[i] = [(row[0][0], row[0][1])]
            continue
        t0, t1 = node.birth_time, node.end_time
        n = len(row) - 2                               # the midpoints, between the two endpoints
        dt = t1 - t0
        segments[i] = [(t0 + k * dt / n, row[k + 1][1]) for k in range(n)]
    return segments


def driver_from_continuous_result(result, *, step: float | None = None,
                                  path: bool = True) -> DriverTrajectory:
    """Build a `DriverTrajectory` from a **continuous** trait result (`traits.simulate_continuous`).

    A continuous trait has no discrete switches, so there is no stochastic map to replay. Instead
    the trait is read as a piecewise-constant function of time — the same shape a discrete trait's
    map has, which is what lets one engine consume both — with each stretch lasting at most ``step``
    time units (`interpolated_segments`). The engine then steps its Gillespie to every stretch
    boundary, so within a stretch the rate really is constant and the exponential draw is exact
    there; nothing is thinned or rejected.

    ``step`` is the resolution, in the tree's own time units, and it is the knob that trades accuracy
    for speed: halving it doubles the breakpoints and so the work. ``None`` takes
    `CONTINUOUS_DRIVER_FRACTION` of the tree's height.

    ``path`` says what the trait is **between** the nodes, and it is on by default:

    - ``True`` draws the trait's real path within each branch, conditioned on both of the branch's
      endpoint values (`zombi2.traits.path`). It is drawn here, after the run, off its own random
      stream keyed by the branch, so the trait's node values are the same numbers either way.
    - ``False`` takes the **straight line** between the endpoints, which is the mean of that bridge
      with the excursions dropped. A real path wanders either side of the line, and under a
      non-linear mapping those excursions do not average out, so the line is biased and a smaller
      ``step`` does not remove the bias. This was the only reading before #454.

    ``True`` falls back to the line where the result carries no bridge law — a trait grown with
    ``regimes=`` or with a driven optimum, whose θ is piecewise constant along the branch. A
    file-backed driver (``trait_values.tsv``) always takes the line: the file carries the node
    values, not the law that produced them.

    The driver's values are **floats**, so its `Driven` needs a continuous mapping (a
    `~zombi2.params.mapping.Curve` or `~zombi2.params.mapping.Scalar`), not a discrete
    ``{state: factor}`` `~zombi2.params.mapping.Table` (which would match no float and never fire)."""
    tree = result.complete_tree
    values = result.node_values
    if not values:
        raise ValueError("continuous driver result has no node values to interpolate")
    sample = next(iter(values.values()))
    if isinstance(sample, dict):
        raise ValueError(
            "a continuous driver must be a SINGLE-trait result; got multi-trait node values "
            f"({sorted(sample)}). Grow one trait as the driver, or select a component before conditioning.")
    law = getattr(result, "path_law", None) if path else None
    return DriverTrajectory(interpolated_segments(tree, values, step, law=law,
                                                  seed=getattr(result, "seed", None),
                                                  starts=branch_starts(result)))


def check_mapping_fires(mapping, available_states, *, driver_label: str, exhaustive: bool = False) -> None:
    """Raise if a **discrete** (`Table`) mapping's states do not line up
    with the states the driver can take. Such a mismatch leaves lineages at the table's default factor —
    a number that is never touched — a rate's, an extent's or a ``transfer_to`` weight's — so the run
    drifts from the model the log records. It is almost always
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
            f"the driver on {driver_label}: the driver is CONTINUOUS (its values are numbers), so its "
            "mapping must be a Curve (value -> factor) or a Scalar (a log-link), not a {state: factor} "
            "table — a table names discrete states, which a continuous value never equals.")
    named = set(mapping.per_state)
    have = {str(s) for s in available_states}
    if exhaustive:
        stray = named - have
        if stray:
            raise ValueError(
                f"the driver on {driver_label}: the mapping names state(s) {sorted(stray)} that are not "
                f"among the driver's states {sorted(have)} — a factor for a state that can never occur, "
                f"so it would silently never apply. Check for a typo in the state names.")
        return
    if not (named & have):
        raise ValueError(
            f"the driver on {driver_label}: the mapping's states {sorted(named)} match none of the "
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
            f"the driver on {driver_label}: the mapping names state(s) {sorted(stray)} that the "
            f"driver never takes (it takes {sorted(have)}), so those factors were never applied. "
            f"The states it did match are still driving the rate. Check for a typo — this is a "
            f"warning rather than an error only because a mapping may legitimately name a state "
            f"that this particular run did not reach.",
            stacklevel=2)


def driven_mods(rate) -> list:
    """The `Driven` modifiers a rate carries, or ``[]`` when it carries none. A non-empty list means
    the rate reads an evolved value on each lineage, so the engine must thread a ``drivers`` value and
    step where the driver switches.

    It lives here rather than in one level's package because every level that reads a driver asks the
    same question first, and the answer is a fact about the rate grammar (SPEC §5), not about any one
    engine. (The trait and genome engines still each carry a private copy from before this existed;
    folding them onto this one is a tidy-up, not a behaviour change.)"""
    return [m for m in rate.modifiers if isinstance(m, Driven)]


def names_a_live_level(driver: object) -> bool:
    """Whether a ``Driven`` ``driver`` names a **level growing beside the run** rather than a
    finished driver.

    SPEC §5: "a finished result makes the run conditioned, and the name of a level growing beside it
    makes the run joint". One modifier, one spelling, and the *driver* is what tells the two apart —
    so this is the predicate that reads it, not a judgement about what the driven level then
    does with it. The live names are the ones the joint runs read — `zombi2.joint.simulate` across
    levels, a level's own ``joint=True`` form within one: ``"trait"``, ``"traits:<name>"``,
    ``"genomes:count"``, ``"genomes:<family>"`` and ``"sequences:<name>"``.

    An engine that resolves its drivers up front uses this to answer a live name in the modelling
    terms — the run asked for is joint, and it belongs on the joint entry point — instead of letting
    the string fall through to `load_driver()` and come back as a missing file called ``'trait'``."""
    return isinstance(driver, str) and (
        driver == "trait" or driver.startswith(("genomes:", "traits:", "sequences:")))


def refuse_wrong_direction(driver, level: str | None) -> None:
    """Raise when ``level``'s engine may not read ``driver`` at all — a direction SPEC §3 rules out,
    not something merely unimplemented.

    A driver declares this by answering ``refuses(level)`` with the reason, or ``None``. It lives on
    the driver for the same reason `resolve_driver()`'s protocol does: this module serves every level
    and must not import from any of them.

    ``level`` is the engine name each level already answers to (``"genomes.family"``, ``"sequences"``,
    … — `zombi2.params.evaluate.Modifier.implemented_for` lists them); ``None`` checks nothing."""
    if level is None:
        return
    refuses = getattr(driver, "refuses", None)
    if refuses is None:
        return
    why = refuses(level)
    if why:
        raise ValueError(why)


def resolve_driver(driver, tree, *, step: float | None = None,
                   level: str | None = None, path: bool = True) -> DriverTrajectory:
    """Resolve a conditioned `Driven`'s ``driver`` into a `DriverTrajectory` — a **filename**
    (str) via `load_driver()` (replayed against ``tree``, the driven run's own species tree), an
    object that answers for itself through ``as_driver_trajectory(tree, step=…)`` (a genome run's
    ``presence("name")``, a sequence run's ``gc()``), or an **in-memory** trait result via
    `driver_from_result()` (which carries its own tree).
    All three are conditioning (the driver grown first); the object forms just spare you the
    ``write``/read step in a single session.

    ``step`` is the continuous-driver resolution (see `interpolated_segments`); it is ignored by a
    discrete driver, whose stretches are exact. ``path`` says whether a continuous trait handed over
    in memory is read along its drawn path (the default) or along the straight line between its node
    values (`driver_from_continuous_result`); the other two forms ignore it. ``level`` names the
    engine doing the reading, which is what lets a driver refuse a level that sits above it
    (`refuse_wrong_direction`).

    A **live level name** is refused here, whichever engine is asking: it is the joint spelling of a
    driver (`names_a_live_level`), read as the run goes by a joint engine, and there is nothing grown
    yet for this to resolve. The refusal lives at this one choke point so that no engine that resolves
    up front can let the name fall through to `load_driver()` as a filename."""
    if names_a_live_level(driver):
        raise ValueError(
            f"the driver {driver!r} names a level growing beside the run — the joint spelling of a "
            "driver (SPEC §5) — and this run resolves its drivers before it starts, which takes a "
            "finished one: a result grown first, or the file it wrote. Two levels that drive each "
            "other are one run: joint.simulate(...) across levels, or the level's own function with "
            "joint=True within one.")
    if isinstance(driver, str):
        # A file-backed driver keeps the straight line. ``trait_values.tsv`` carries the node values
        # but not the law that produced them — the σ², the pull, the optimum, the schedule — so
        # there is nothing here to draw a bridge from, and writing the law into the file is an open
        # question rather than part of #454.
        return load_driver(driver, tree, step=step)
    if hasattr(driver, "as_driver_trajectory"):
        # a level that knows how to answer "what state was lineage L in at time t?" for itself —
        # `genomes.presence("tox")` is the first. The protocol is one method rather than an isinstance
        # branch per level so this module stays free of imports from the levels it serves.
        refuse_wrong_direction(driver, level)
        return driver.as_driver_trajectory(tree, step=step)
    return driver_from_result(driver, step=step, path=path)


__all__ = ["DriverTrajectory", "load_driver", "driver_from_result", "branch_starts",
           "path_points", "interpolated_segments", "resolve_driver",
           "refuse_wrong_direction", "check_mapping_fires", "driven_mods", "names_a_live_level"]
