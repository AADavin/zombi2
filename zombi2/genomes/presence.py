"""Gene presence as a conditioning driver — is this family in this lineage, right now?

A driver answers one question: *what state was lineage L in at time t?* Traits answer it
from their own event log. A **gene family** can answer it too — present or absent — and that is the
other direction of the same relation: a trait can already make a genome rate faster, and this is what
lets a genome make a trait's rate faster.

The trajectory is read off the family's **gene tree**. Each `GeneNode` records the species branch it
lived on and when it ended, and its parent records when it began, so a gene is an interval on a branch
and the family is present there exactly while at least one interval covers it.

One reader therefore serves **all three resolutions**: each recovers gene trees, so what differs
between them is already behind us here. What a resolution still decides is how a family is *named*
(`_declared_names`) and what takes it away — a whole-copy loss at family and ordered, a deletion
covering the gene at nucleotide.

That is not a way of avoiding the event log — it is the log after the step that turns it into
genealogies. Presence needs to know when each copy *began* and *ended*, and neither form of the log
says so directly: ``result.events`` is one row per gene-tree edge (two per duplication, two per
transfer), and the written ``genome_events.tsv`` is one row per event with ``parents`` / ``children``
gene ids. Either way, getting from ids to intervals means rebuilding the genealogy, which is exactly
`gene_trees_from_edges`. Doing it here would be a second implementation of a tested one.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..params.conditioned import DriverTrajectory

#: what the driver's states are called. A mapping written against this driver names these, so they
#: are the vocabulary and not an implementation detail: ``{"present": 4.0, "absent": 1.0}``.
PRESENT, ABSENT = "present", "absent"


@dataclass(frozen=True)
class GenePresence:
    """Whether a named gene family is in each lineage's genome, over time.

    Built by ``result.presence("name")`` and read with ``scaled_by`` like a grown trait::

        g = simulate_genomes_family(sp, families=[family("tox")], loss=0.3, seed=1)
        simulate_discrete(sp.complete_tree, states=["harmless", "pathogenic"], start="harmless",
                          switch=PerLineage(0.1).scaled_by(g.presence("tox"),
                                                           {"present": 5.0, "absent": 1.0}), seed=2)

    This is **conditioning**: the genome is grown first and held fixed while the trait reads it. For
    the two growing together — a gene whose presence shapes the tree the genome is evolving on —
    that is `zombi2.joint`, and it drives speciation rather than a trait.
    """

    result: object
    name: str

    def _segments(self, tree):
        _same_tree(self.result, tree)
        family = _named(self.result, self.name)
        return _segments(tree, {family: _spans(self.result, family)},
                         lambda covering: PRESENT if covering else ABSENT)

    def as_driver_trajectory(self, tree, *, step: float | None = None) -> DriverTrajectory:
        """The per-lineage present/absent trajectory, for `zombi2.params.conditioned.resolve_driver`.
        ``step`` is ignored: the gene tree says exactly when a copy was gained or lost, so there is no
        resolution left to choose."""
        return DriverTrajectory(self._segments(tree))

    def history(self, tree) -> dict[int, list]:
        """``{node id: [(state, duration), …]}`` — the per-branch map, durations summing to the
        branch length.

        The same shape `zombi2.traits.TraitsResult.history` has, so anything that draws a trait's
        history down a tree draws this one too. `as_driver_trajectory` is what a rate reads; this is
        what a figure reads."""
        return _history(self._segments(tree), tree)

    def __repr__(self) -> str:
        return f"presence({self.name!r})"


@dataclass(frozen=True)
class ModuleCompletion:
    """How much of a **module** — a named group of gene families — a lineage carries, over time.

    A number in ``[0, 1]``: the fraction of the module's families with at least one copy. Built by
    ``result.completion("name")`` and read with a `Curve`, the way any continuous driver is::

        g = simulate_genomes_family(sp, seed=1, families=[
            family(n, module="flagellum") for n in ("flgA", "flgB", "flgC")])
        simulate_discrete(sp.complete_tree, states=["sessile", "motile"], start="sessile",
                          switch=PerLineage(0.05).scaled_by(g.completion("flagellum"),
                                                            Curve(lambda f: 0.2 + 20.0 * f ** 4)), seed=2)

    **A fraction rather than a yes/no, deliberately.** Under independent loss the chance that *every*
    family of a module survives falls off geometrically with its size — measured on a 200-tip tree, a
    module of three was complete at 189 tips and one of six at none of them — so a complete/incomplete
    driver would be a constant for anything but the smallest modules, and the run would say nothing.
    The fraction is always informative, and a threshold is expressible where every other response
    shape already lives: in the `Curve`. ``lambda f: 8.0 if f > 0.8 else 1.0`` is "eight times faster
    once four fifths of it is there".
    """

    result: object
    name: str

    def _segments(self, tree):
        _same_tree(self.result, tree)
        members = self.result.modules[self.name]
        ids = {_named(self.result, m): m for m in members}
        per_family = {fid: _spans(self.result, fid) for fid in ids}
        n = len(ids)
        return _segments(tree, per_family, lambda covering: len(covering) / n)

    def as_driver_trajectory(self, tree, *, step: float | None = None) -> DriverTrajectory:
        """The per-lineage completion trajectory, for `zombi2.params.conditioned.resolve_driver`. ``step``
        is ignored: the fraction is a number, but it moves one family at a time and the gene trees
        give every step exactly."""
        return DriverTrajectory(self._segments(tree))

    def history(self, tree) -> dict[int, list]:
        """``{node id: [(fraction, duration), …]}`` — the per-branch map, as `GenePresence.history`
        gives it, with a number in ``[0, 1]`` in place of a state."""
        return _history(self._segments(tree), tree)

    def __repr__(self) -> str:
        return f"completion({self.name!r})"


def _spans(result, family: int) -> dict[int, list[tuple[float, float]]]:
    """``{branch: [(began, ended), …]}`` — every gene of ``family``, as an interval on the branch it
    lived on. A `GeneNode` records when it *ended*; it began when its parent ended, and the root when
    the family originated."""
    out: dict[int, list[tuple[float, float]]] = {}
    gt = result.gene_trees.get(family)
    if gt is None:
        return out
    stack = [(gt.complete, gt.origination)]
    while stack:
        node, began = stack.pop()
        out.setdefault(node.species, []).append((began, node.time))
        for child in node.children:
            stack.append((child, node.time))
    return out


def _segments(tree, per_family, value_of) -> dict[int, list[tuple[float, object]]]:
    """Cut every branch wherever any of the families appears or disappears on it, and label each
    stretch with ``value_of(covering)`` — the set of families with a copy alive there.

    **Every** lineage of ``tree`` gets segments, because a driver has to answer for whatever branch
    the engine is on; one that never held any of the families simply reads the empty-set value for
    its whole life rather than being missing.
    """
    segments: dict[int, list[tuple[float, object]]] = {}
    for node_id, nd in tree.nodes.items():
        start, end = nd.birth_time, nd.end_time
        here = {fam: spans.get(node_id, ()) for fam, spans in per_family.items()}
        cuts = sorted({start} | {t for spans in here.values() for span in spans for t in span
                       if start < t < end})
        segs: list[tuple[float, object]] = []
        for t in cuts:
            covering = frozenset(f for f, spans in here.items() if any(b <= t < e for b, e in spans))
            state = value_of(covering)
            if not segs or segs[-1][1] != state:        # only where it actually changes
                segs.append((t, state))
        segments[node_id] = segs
    return segments


def _same_tree(result, tree) -> None:
    if getattr(result, "complete_tree", None) is not tree and \
            set(result.complete_tree.nodes) != set(tree.nodes):
        raise ValueError(
            "the genome run and the run reading it are on different species trees, so their lineage "
            "ids do not refer to the same branches. Grow both on the same complete tree.")


def _declared_names(result) -> dict[str, int]:
    """``{name: family id}`` — the run's declared families, under whichever handle its resolution
    keeps them: ``family_names`` at family and ordered, ``gene_names`` at nucleotide, where a declared
    family *is* one gene and its name came off a GFF. Asked once here, not at each call site."""
    names = getattr(result, "family_names", None)
    return getattr(result, "gene_names", {}) if names is None else names


def _named(result, name: str) -> int:
    names = _declared_names(result)
    if name not in names:
        raise KeyError(
            f"no named family {name!r} in this genome run; declared families are {sorted(names)}. A "
            f"family must be declared to be named here — `families=[family(...)]` at the family and ordered "
            f"resolutions, a GFF `ID` / `Name` at the nucleotide one.")
    return names[name]


def _history(segments, tree) -> dict[int, list[tuple[object, float]]]:
    """``{branch: [(state, duration), …]}``, durations summing to the branch length."""
    out: dict[int, list[tuple[object, float]]] = {}
    for node_id, segs in segments.items():
        bounds = [t for t, _ in segs] + [tree.nodes[node_id].end_time]
        out[node_id] = [(state, bounds[i + 1] - bounds[i]) for i, (_, state) in enumerate(segs)]
    return out
