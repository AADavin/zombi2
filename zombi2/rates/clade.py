"""`Clade` — a value read off the **tree itself**, so a rate can depend on where a lineage sits.

Everything else a rate reads was grown by another level: a trait, a genome's gene content, a
sequence's composition. A clade is different — it is a fact about the tree the run is already walking,
so there is nothing to grow first and nothing to condition on::

    loss = 0.2 * ScaledBy(Clade({"fast": ["n12", "n27"]}), {"fast": 3.0})

Membership is **constant along a branch**, which is what makes this cheap: a driver that switches
mid-branch forces the Gillespie to stop at each switch, and this one never switches. It is painted
once per run, and from the engine's point of view it is an ordinary conditioned driver — the value is
resolved before the run and looked up per lineage, exactly like a trait read from a file.

**Naming a clade.** Either a list of tips, in which case the clade is the subtree below their most
recent common ancestor, or a single node id (an ``int``, or an ``"n<id>"`` label), in which case it is
that node's whole subtree. Clades must be disjoint; a lineage in none of them is in the implicit group
``"rest"``, which a mapping may name like any other state.

`Clades` — the plural — is the `transfer_to` rule that weights a
**pair**, the donor's clade against the recipient's. This is the singular value: one lineage, one
label, readable by any rate. The two share `resolve_groups`, so a clade means the same
thing whichever way it is read.
"""

from __future__ import annotations


class Clade:
    """Which named clade a lineage belongs to, as a categorical value on that lineage.

    A **conditioned** driver in the ordinary sense — its value is known before the run — but with no
    file and no earlier run behind it, because the tree is already an input. It therefore works at
    every level that reads a driver at all, and on a growing tree it does not: a clade is only defined
    once the tree exists, which is why a joint run refuses it.
    """

    def __init__(self, groups: dict) -> None:
        if not isinstance(groups, dict) or not groups:
            raise ValueError(
                "Clade needs a non-empty {label: clade} dict, where a clade is a list of tips (the "
                "subtree below their MRCA) or a single node id — e.g. Clade({'A': ['n1', 'n2']})")
        for label in groups:
            if not isinstance(label, str) or not label.strip():
                raise ValueError(f"clade labels must be non-empty strings, got {label!r}")
        self.groups = groups

    def as_driver_trajectory(self, tree, *, step: float | None = None):
        """The per-lineage lookup a driven rate reads — one stretch per lineage, starting at its
        birth, because membership never changes along a branch.

        ``step`` is the resolution a *continuous* driver is read at and is meaningless here: a clade
        label is categorical and its stretches are already exact, so nothing is approximated and
        nothing is gained by cutting the branch finer.
        """
        from ..genomes._transfer import resolve_groups
        from .driver import DriverTrajectory

        painted = resolve_groups(tree, self.groups)
        return DriverTrajectory({i: [(tree.nodes[i].birth_time, painted[i])] for i in tree.nodes})

    def written_form(self) -> str:
        """A clade is built from literals — labels, node ids, tip names — so unlike every other
        driver it can be written into a run's log and pasted back. `DrivenBy`
        asks for this when recording the rate."""
        return repr(self)

    def __repr__(self) -> str:
        return f"Clade({self.groups!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Clade) and other.groups == self.groups

    def __hash__(self) -> int:
        return hash((Clade, tuple(sorted(self.groups))))


__all__ = ["Clade"]
