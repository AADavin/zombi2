"""Gene presence as a conditioning driver — is this family in this lineage, right now?

A `DrivenBy` driver answers one question: *what state was lineage L in at time t?* Traits answer it
from their own event log. A **gene family** can answer it too — present or absent — and that is the
other direction of the same relation: a trait can already make a genome rate faster, and this is what
lets a genome make a trait's rate faster.

The trajectory is read off the family's **gene tree**. Each `GeneNode` records the species branch it
lived on and when it ended, and its parent records when it began, so a gene is an interval on a branch
and the family is present there exactly while at least one interval covers it.

That is not a way of avoiding the event log — it is the log after the step that turns it into
genealogies. Presence needs to know when each copy *began* and *ended*, and neither form of the log
says so directly: ``result.events`` is one row per gene-tree edge (two per duplication, two per
transfer), and the written ``genome_events.tsv`` is one row per event with ``parents`` / ``children``
gene ids. Either way, getting from ids to intervals means rebuilding the genealogy, which is exactly
`gene_trees_from_events`. Doing it here would be a second implementation of a tested one.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..rates.driver import DriverTrajectory

#: what the driver's states are called. A mapping written against this driver names these, so they
#: are the vocabulary and not an implementation detail: ``{"present": 4.0, "absent": 1.0}``.
PRESENT, ABSENT = "present", "absent"


@dataclass(frozen=True)
class GenePresence:
    """Whether a named gene family is in each lineage's genome, over time.

    Built by ``result.presence("name")`` and handed to `DrivenBy` like a grown trait::

        g = simulate_genomes_family(sp, family_names=["tox"], loss=0.3, seed=1)
        simulate_discrete(sp.complete_tree, states=["harmless", "pathogenic"], start="harmless",
                          switch=0.1 * mod.DrivenBy(g.presence("tox"),
                                                    {"present": 5.0, "absent": 1.0}), seed=2)

    This is **conditioning**: the genome is grown first and held fixed while the trait reads it. For
    the two growing together — a gene whose presence shapes the tree the genome is evolving on —
    that is `zombi2.joint`, and it drives speciation rather than a trait.
    """

    result: object
    name: str

    def as_driver_trajectory(self, tree) -> DriverTrajectory:
        """The per-lineage present/absent trajectory, for `zombi2.rates.driver.resolve_driver`."""
        return DriverTrajectory(_presence_segments(self.result, self.name, tree))

    def __repr__(self) -> str:
        return f"presence({self.name!r})"


def _presence_segments(result, name: str, tree) -> dict[int, list[tuple[float, object]]]:
    """``{node id: [(start time, state), …]}`` — each lineage's branch cut where the family's copy
    count crosses zero.

    ``tree`` is the target run's species tree and decides which lineages get a segment: **every**
    one does, because a driver has to answer for any branch the engine asks about, and a lineage that
    never held the family is simply absent for its whole life rather than missing.
    """
    if name not in result.family_names:
        raise KeyError(
            f"no named family {name!r} in this genome run; declared families are "
            f"{sorted(result.family_names)}. A family has to be declared with `family_names=` to be "
            f"named here — an anonymous one has an id but nothing stable to call it by.")
    family = result.family_names[name]
    if getattr(result, "complete_tree", None) is not tree and \
            set(result.complete_tree.nodes) != set(tree.nodes):
        raise ValueError(
            "the genome run and the run reading it are on different species trees, so their lineage "
            "ids do not refer to the same branches. Grow both on the same complete tree.")

    # every gene of the family, as (branch, when it began, when it ended). A node records when it
    # *ended*; it began when its parent ended, and the root when the family originated.
    spans: dict[int, list[tuple[float, float]]] = {}
    gt = result.gene_trees.get(family)
    if gt is not None:
        stack = [(gt.complete, gt.origination)]
        while stack:
            node, began = stack.pop()
            spans.setdefault(node.species, []).append((began, node.time))
            for child in node.children:
                stack.append((child, node.time))

    segments: dict[int, list[tuple[float, object]]] = {}
    for node_id, nd in tree.nodes.items():
        start, end = nd.birth_time, nd.end_time
        here = spans.get(node_id, ())
        # the branch is cut at every instant a gene of this family appears or disappears on it; the
        # state between two cuts is whether any gene covers that stretch
        cuts = sorted({start} | {t for span in here for t in span if start < t < end})
        segs: list[tuple[float, object]] = []
        for t in cuts:
            state = PRESENT if any(b <= t < e for b, e in here) else ABSENT
            if not segs or segs[-1][1] != state:        # only where it actually changes
                segs.append((t, state))
        segments[node_id] = segs
    return segments
