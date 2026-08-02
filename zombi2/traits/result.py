"""Traits — the result object and its serialisation: Change, TraitsResult, and the writers."""

from __future__ import annotations

import collections
import pathlib
from dataclasses import dataclass, field
from functools import cached_property
from typing import cast


from .._runtime.summary import _stats, write_summary
from ..genomes.events import _name
from ..tree import Node, Tree

_WRITE_OUTPUTS = ("values", "events", "tree", "summary")  # write vocabulary; "events" = the trait event log

@dataclass(frozen=True)
class Change:
    """A realized trait change — one entry of the event log, the trait twin of the genome level's
    `Event`. On lineage ``lineage`` at ``time`` (origin-forward, the species-tree
    clock) the state went from ``from_state`` to ``to_state``. ``kind`` is ``"on_branch"`` — a switch
    *along* a branch (an Mk transition) — ``"on_speciation"`` — a jump *at* a speciation node (from
    ``at_speciation``; for a continuous trait ``from_state`` / ``to_state`` are the pre- and post-jump
    values) — or ``"initial"``, one synthetic entry at t=0 giving the **initial state** the run
    started in (``from_state`` ``None``, ``time`` the root's ``birth_time``). That row is what lets the
    log stand on its own: the tree plus the initial state plus the switches determines the trait on
    every lineage at every instant, so no separate driver file is needed."""

    time: float
    kind: str
    lineage: int
    from_state: object
    to_state: object



@dataclass
class TraitsResult:
    """What ``simulate_continuous`` / ``simulate_discrete`` returns: the ``complete_tree`` it ran on,
    ``node_values`` at **every** node (the value at each node — extant, extinct, and internal alike; a
    float for a continuous trait, a state label for a discrete / threshold one, a per-trait dict for
    correlated traits), the ``events`` log, the ``seed``, and the ``kind`` (``"continuous"`` /
    ``"discrete"`` / ``"threshold"``). The observed dataset is the extant tips, ``.values``.

    ``events`` is the timestamped event log — the **same shape as the genome level's** and the source
    of truth for a discrete (Mk) trait, from which ``history`` (the per-branch stochastic character
    map) is derived. A continuous trait has no along-branch events, so its log holds only the
    on-speciation jumps (empty without ``at_speciation``) while ``node_values`` carries the diffusion;
    a threshold trait's crossings are un-timed, so its log is empty and it has no map.
    """

    complete_tree: Tree
    node_values: dict[int, object]
    events: list[Change] = field(default_factory=list)
    seed: int | None = None
    kind: str = "continuous"

    def __repr__(self) -> str:
        return (f"TraitsResult(a {self.kind} trait over {len(self.values)} extant tips, "
                f"{len(self.node_values)} nodes, {len(self.events)} events, seed={self.seed})")

    @property
    def values(self) -> dict[str, object]:
        """The observed trait dataset — the value at each **extant** tip (the comparative-data
        vector), keyed by the **tip name** the tree writes: ``n5``, or ``e5`` for a lineage that died.

        Keyed by name because the only thing anyone does with this is join it to the tree, and the
        tree names its tips. It used to be keyed by the bare node id (``5``), which the written
        ``trait_values.tsv`` and every Newick label do not use — so a comparative dataset built in
        Python shared **no keys at all** with the tree beside it, and nothing said so. The two files
        a ``write()`` produces always did match; it was the in-memory pair that did not.

        `values_by_id` is the old view, for code that joins on node ids. Internal and extinct nodes
        keep their exact ancestral / lineage values in ``node_values``, which stays id-keyed: it is
        the run's own record, not a dataset to export."""
        name = self.complete_tree.labels()
        return {name[n.id]: self.node_values[n.id] for n in self.complete_tree.extant_leaves()}

    @property
    def values_by_id(self) -> dict[int, object]:
        """`values`, keyed by node id rather than tip name — the shape `values` had before it was
        keyed to match the tree. For joining against ``node_values``, ``complete_tree.nodes`` or
        anything else that works in ids."""
        return {n.id: self.node_values[n.id] for n in self.complete_tree.extant_leaves()}

    @cached_property
    def history(self) -> dict[int, list] | None:
        """The per-branch **stochastic character map** — ``{node: [(state, duration), …]}`` whose
        durations sum to the branch length — **derived from the event log** (a discrete / Mk trait
        only). ``None`` for a continuous trait (a diffusion has no map) and for a threshold trait
        (its liability crossings are un-timed)."""
        if self.kind != "discrete":
            return None
        return _history_from_events(self.complete_tree, self.node_values, self.events)

    def summary(self) -> dict:
        """What this run produced, as a plain dict — the payload of ``trait_summary.json``.

        A trait's shape depends on its kind, so the summary does too. A **discrete** trait — and a
        **threshold** one, which reads a discrete state off a continuous liability — is described
        by its switches: how many, and how the tips ended up distributed over the states — which is the
        thing you look at first, because a run whose tips are all in one state has told you nothing.
        A **continuous** one is described by where the values got to, since there are no along-branch
        events to count; its log holds only the on-speciation jumps, and that count is here so an
        empty one is visibly empty rather than ambiguous."""
        values = list(self.values.values())
        switches = sum(1 for e in self.events if e.kind == "on_branch")
        jumps = sum(1 for e in self.events if e.kind == "on_speciation")
        out: dict[str, object] = {
            "level": "traits",
            "seed": self.seed,
            "kind": self.kind,
            "tips": len(values),
            "nodes": len(self.node_values),
            "events": {"on_branch": switches, "on_speciation": jumps},
        }
        # continuous is the only NUMERIC kind. A threshold trait reads a discrete state off a
        # continuous liability, so its values are states like a discrete trait's — taking a mean of
        # them is what the first version of this tried to do.
        if self.kind != "continuous":
            counts = collections.Counter(str(v) for v in values)
            out["states"] = dict(sorted(counts.items()))
            out["states_at_tips"] = len(counts)
            # a run whose tips all share one state is degenerate, and the number that says so
            out["most_common_share"] = (round(counts.most_common(1)[0][1] / len(values), 6)
                                        if values else None)
        else:
            numeric = [float(cast(float, v)) for v in values]
            out["values"] = _stats(numeric)
            # the root NODE, i.e. after diffusing along the stem — not the value the run started
            # from, which is `start` and belongs to no node
            out["value_at_root_node"] = float(
                cast(float, self.node_values[self.complete_tree.root]))
        return out

    def write(self, directory, outputs=("values",)) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed): ``"values"`` →
        ``trait_values.tsv`` (the ``node<TAB>kind<TAB>trait`` table over **every** node — tips, extinct
        lineages and internal nodes; ``kind`` is the tip's fate — ``extant`` / ``extinct`` (/ ``unsampled``
        under incomplete sampling) — or ``ancestor`` for an internal node, so the extant tips filter out
        with ``kind == "extant"``); ``"events"`` →
        ``trait_events.tsv``, the event log (``time · kind · lineage · from · to``) — one ``initial``
        row at t=0 giving the initial state, then every switch in time order; ``"tree"`` →
        ``trait_tree.nwk``, the complete tree as Newick with **every** node annotated ``[&trait=…]``
        (a *trait tree*, carrying the exact ancestral values; opens in FigTree / iTOL).

        ``trait_events.tsv`` is also the **conditioning file**: a genome / sequence run drives a rate
        with ``mod.DrivenBy("trait_events.tsv", …)``, replaying it against the shared tree. A
        **discrete** trait's log reconstructs its state on every lineage exactly (that is what the
        ``initial`` row and the switch times are for); a continuous trait's diffusion cannot be rebuilt
        from events, so it carries only the ``initial`` row and any on-speciation jumps."""
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        names = self.complete_tree.labels()   # e<id> for a lineage that died; n<id> for the rest
        if "values" in outputs:
            # every node — extant tips, extinct lineages (e<id>) and internal nodes (n<id>) alike, each
            # with its exact value and a `kind` column (the tip's fate, or `ancestor`) so a comparative
            # vector of exactly the extant tips is one `kind == "extant"` filter
            kinds = {i: _node_kind(n) for i, n in self.complete_tree.nodes.items()}
            (d / "trait_values.tsv").write_text(_values_tsv(self.node_values, names, kinds),
                                                encoding="utf-8")
        if "events" in outputs:
            (d / "trait_events.tsv").write_text(_events_tsv(self.events, names), encoding="utf-8")
        if "summary" in outputs:
            write_summary(d / "trait_summary.json", self.summary())
        if "tree" in outputs:
            (d / "trait_tree.nwk").write_text(
                _trait_newick(self.complete_tree, self.node_values) + "\n", encoding="utf-8")



def _fmt(v) -> str:
    """A trait value for TSV: a continuous float compactly, a discrete state label as-is."""
    return f"{v:.6g}" if isinstance(v, float) else str(v)



def _trait_annotation(v) -> str:
    """The BEAST/FigTree ``[&…]`` comment for a node value — ``[&trait=…]`` for a single trait,
    one ``key=value`` per trait for a correlated (dict) value."""
    if isinstance(v, dict):
        return "[&" + ",".join(f"{k}={_fmt(x)}" for k, x in v.items()) + "]"
    return f"[&trait={_fmt(v)}]"



def _trait_newick(tree: "Tree", node_values: dict) -> str:
    """The complete tree as Newick with **every** node annotated with its trait value (a *trait
    tree*). Mirrors `zombi2.species.Tree.to_newick()` — branch length ``end_time − birth_time``,
    leaves and internals named ``n<id>`` (``e<id>`` for a lineage that died), the root carrying its
    stem — and adds the ``[&trait=…]``
    comment at each node, so the exact ancestral states ride along the tree."""
    name = tree.labels()

    def emit(i: int) -> str:
        node = tree.nodes[i]
        bl = node.end_time - node.birth_time
        tag = f"{name[i]}{_trait_annotation(node_values[i])}"
        if node.children is None:
            return f"{tag}:{bl:.7g}"
        return f"({','.join(emit(c) for c in node.children)}){tag}:{bl:.7g}"

    root = tree.nodes[tree.root]
    tag = f"{name[tree.root]}{_trait_annotation(node_values[tree.root])}"
    stem = root.end_time - root.birth_time
    if root.children is None:
        return f"{tag}:{stem:.7g};"
    return f"({','.join(emit(c) for c in root.children)}){tag}:{stem:.7g};"



def _node_kind(node: "Node") -> str:
    """A node's ``kind`` for the values table: an internal node is an ``ancestor``; a tip is labelled
    by its fate — ``extant``, ``extinct``, or ``unsampled`` (a survivor not observed under incomplete
    sampling). A present-day tip whose fate was never resolved (the bare ``"alive"`` default) reads as
    ``extant``. So ``kind == "extant"`` isolates exactly the observed tips a comparative method wants."""
    if not node.is_leaf:
        return "ancestor"
    return node.fate if node.fate in ("extinct", "unsampled") else "extant"


def _values_tsv(values: dict[int, object], names: dict | None = None,
                kinds: dict | None = None) -> str:
    """Node values as a ``node<TAB>kind<TAB>…`` table, one row per node in id order (``n<id>``, or
    ``e<id>`` for a lineage that died, to match the Newick). ``kind`` (from ``kinds``, keyed by node id)
    is a tip's fate — ``extant`` / ``extinct`` / ``unsampled`` — or ``ancestor`` for an internal node,
    so the extant tips filter out of the all-nodes table with ``kind == "extant"``. A single trait gives
    one ``trait`` column; correlated traits give one per trait."""
    kinds = kinds or {}
    kind = lambda i: kinds.get(i, "ancestor")
    first = next(iter(values.values())) if values else None
    if isinstance(first, dict):                                   # correlated / multi-trait
        per_trait = cast("dict[int, dict[str, object]]", values)
        cols = list(first)
        rows = ["node\tkind\t" + "\t".join(str(c) for c in cols)]
        for i in sorted(per_trait):
            rows.append(f"{_name(names, i)}\t{kind(i)}\t"
                        + "\t".join(_fmt(per_trait[i][c]) for c in cols))
        return "\n".join(rows) + "\n"
    rows = ["node\tkind\ttrait"]
    for i in sorted(values):
        rows.append(f"{_name(names, i)}\t{kind(i)}\t{_fmt(values[i])}")
    return "\n".join(rows) + "\n"



def _events_tsv(changes: list[Change], names: dict | None = None) -> str:
    """The event log as ``time<TAB>kind<TAB>lineage<TAB>from<TAB>to`` (``kind`` = initial / on_branch /
    on_speciation), the ``initial`` row first, then the switches in time order — the trait twin of
    ``genome_events.tsv``, and the conditioning file a driven run replays.

    ``lineage · from · to`` is the shape a *state change* wants, not the ``parents`` / ``children`` of
    a birth: nothing is born or dies here, one lineage's trait simply moves from one value to another.

    Times are written at **full float precision** (``repr``), not rounded: a driven run steps its
    Gillespie exactly at each switch, so a rounded time would make the file-driven run diverge from the
    in-memory one. The ``initial`` row's ``from`` is empty.

    **Correlated traits widen the table** rather than repeating a row per trait, exactly as
    ``trait_values.tsv`` does: the ``from`` / ``to`` pair becomes one pair per trait,
    ``from:<trait>`` · ``to:<trait>``. A correlated jump moves every trait at once, so it is one
    event, and one row per trait would suggest they were several."""
    multi = bool(changes) and isinstance(changes[0].to_state, dict)
    if not multi:
        rows = ["time\tkind\tlineage\tfrom\tto"]
        for c in changes:
            frm = "" if c.from_state is None else _fmt(c.from_state)  # the initial row leads from nothing
            rows.append(f"{c.time!r}\t{c.kind}\t{_name(names, c.lineage)}\t{frm}\t{_fmt(c.to_state)}")
        return "\n".join(rows) + "\n"
    cols = list(cast("dict[str, object]", changes[0].to_state))
    rows = ["time\tkind\tlineage\t"
            + "\t".join(f"from:{c}" for c in cols) + "\t" + "\t".join(f"to:{c}" for c in cols)]
    for c in changes:
        to_state = cast("dict[str, object]", c.to_state)      # every change of a multi-trait run
        from_state = cast("dict[str, object] | None", c.from_state)
        cells = ["" for _ in cols] if from_state is None else [_fmt(from_state[t]) for t in cols]
        rows.append(f"{c.time!r}\t{c.kind}\t{_name(names, c.lineage)}\t"
                    + "\t".join(cells) + "\t" + "\t".join(_fmt(to_state[t]) for t in cols))
    return "\n".join(rows) + "\n"



def _history_from_events(tree: "Tree", node_values: dict, events: list) -> dict:
    """Reconstruct the per-branch stochastic character map ``{node: [(state, duration), …]}`` from the
    event log — the inverse of how the Gillespie writes the log. On each branch the on-branch events
    (sorted) cut it into segments; an on-speciation event sets the branch's *start* state, and
    ``node_values[node]`` is its *end* state (the constant value when a branch has no events)."""
    ana: dict[int, list] = {i: [] for i in tree.nodes}
    clado_to: dict[int, object] = {}
    for e in events:
        if e.kind == "initial":
            continue                               # the t=0 marker; node_values covers it here
        if e.kind == "on_speciation":
            clado_to[e.lineage] = e.to_state
        else:
            ana[e.lineage].append(e)
    history: dict[int, list] = {}
    for i in tree.nodes:
        node = tree.nodes[i]
        evs = sorted(ana[i], key=lambda e: e.time)
        state = evs[0].from_state if evs else clado_to.get(i, node_values[i])  # branch-start state
        segs, t = [], node.birth_time
        for e in evs:
            segs.append((state, e.time - t))
            state, t = e.to_state, e.time
        segs.append((state, node.end_time - t))  # final segment to the node's end
        history[i] = segs
    return history


