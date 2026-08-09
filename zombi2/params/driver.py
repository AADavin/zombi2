"""The **drivers** — everything a parameter can read (SPEC §5).

A parameter that is not a constant reads a driver, and a driver is two facts: what it is attached to
(its **unit** — the run, a lineage, a gene family) and how its number is made. Those are independent,
and keeping them apart is what stops the grammar needing a new class per model::

    Random("families", LogNormal(0.0, 0.5))          # made by a draw,      among families
    Random("lineages", Drift(LogNormal(0.0, 0.2)))   # made by inheritance, among lineages
    Time()                                           # measured from the run
    Clade({"fast": [...]})                           # read off the tree itself

`connection` is the other half — what a parameter *does* with the number.

**The unit is an argument, not a class.** A draw among families and a draw among lineages are one
model at two attachments, so ``Random("families", …)`` and ``Random("lineages", …)`` are one class
and two cells of a grid rather than two names to remember. The cells are models the field already
has names for — rate heterogeneity across gene families, the relaxed clock, ClaDS — and the prose
uses those.

What the grid buys is the next cell. A per-chromosome draw needs no new class and no invented name:
it is ``Random("chromosomes", …)``, which constructs, and a cell no engine carries refuses at that
level's gate naming itself rather than reading as a typo.

`Drawn` and `Inherited`, the two objects `Random` builds, are in `law` beside the law that chooses
between them. `Time` is here because it is a driver that is **not** a modifier: it is read through a
verb, where a draw is already a factor and needs none.
"""

from __future__ import annotations

from .modifiers import UNITS, Drawn, Inherited




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
        from .conditioned import DriverTrajectory

        painted = resolve_groups(tree, self.groups)
        return DriverTrajectory({i: [(tree.nodes[i].birth_time, painted[i])] for i in tree.nodes})

    def written_form(self) -> str:
        """A clade is built from literals — labels, node ids, tip names — so unlike every other
        driver it can be written into a run's log and pasted back. `Driven.written_call` asks for
        this when recording the rate."""
        return repr(self)

    def __repr__(self) -> str:
        return f"Clade({self.groups!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Clade) and other.groups == self.groups

    def __hash__(self) -> int:
        return hash((Clade, tuple(sorted(self.groups))))





class Measured:
    """Base for a value the run already knows and can be asked for at any moment — no draw, no
    inheritance, nothing carried per unit. Abstract: use `Time`.

    A measured value is **not** a `Modifier`, and that is the
    distinction the two halves of the grammar rest on. A drawn value is already a dimensionless
    factor, so it multiplies a base on its own. A measured one is a time, a count, an age — a number
    in its own units — so it reaches a rate only through a verb and a mapping.
    """

    unit: str = "run"


class Time(Measured):
    """The run's clock, as a value::

        birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})    # a skyline, in multiples of 0.5
        birth = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})   # the rates themselves

    Attached to the whole run, so every parameter may read it — the run is the coarsest unit, and a
    parameter's units always include it.

    The clock is written with ``changing_at`` when the schedule holds factors, and it is the one
    driver `set_by` also takes, when the schedule holds the numbers themselves. Those two are the
    only spellings: ``scaled_by(Time(), …)`` is refused and names the shortcut, so there is one way
    to write each reading rather than two.
    """

    unit: str = "run"

    def __repr__(self) -> str:
        return "Time()"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Time)

    def __hash__(self) -> int:
        return hash(Time)

#: The drivers a parameter may be **written** with. `Measured` is absent because it is
#: abstract, and `Drawn` / `Inherited` because a law names them, not a driver.
WRITABLE = ("Clade", "Time")

__all__ = ["UNITS", "WRITABLE", "Clade", "Drawn", "Inherited", "Measured", "Time"]
