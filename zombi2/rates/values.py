"""Values — the *what is read* half of a driven parameter (SPEC §5).

A parameter that is not a constant reads a **value**, and a value is two facts: what it is attached
to (its **unit** — the run, a lineage, a gene family) and how its number is **made**. Those two are
independent, and keeping them apart is what stops the grammar needing a new class per model::

    Drawn(per="family", spread=0.5)       # made by a draw,      attached to a family
    Drawn(per="lineage", spread=0.3)      # made by a draw,      attached to a lineage
    Inherited(per="lineage", spread=0.2)  # made by inheritance, attached to a lineage
    Time()                                # measured from the run
    Clade({"fast": [...]})                # read off the tree itself

The `verbs` module is the other half — what a parameter *does* with the number.

`Drawn` and `Inherited` live in `modifiers` beside the classes an engine
dispatches on, because that is what they are; they are re-exported here so the grammar can be
imported as one surface. `Time` is here because it is a value that is **not** a modifier: it is read
through a verb (``ScaledBy(Time(), …)``), where a draw is already a factor and needs none.

**The named cells are the same six models the ``mod.`` names have always given you**, written as a
grid rather than as a list: ``Drawn(per="family", …)`` *is* ``ByFamily(…)``, and both build the
identical object, so a run is unchanged and a script written either way behaves the same. The old
names stay; they are the names of cells, kept because each names a model people cite.

The reason to prefer this spelling is what happens next. A per-chromosome draw needs no new class and
no invented name — it is ``Drawn(per="chromosome", …)``, a cell of a grid you can already see — and a
cell no engine carries refuses at that level's gate, naming itself (*"drawn per chromosome, which the
species engine does not support"*), rather than reading as a typo.
"""

from __future__ import annotations

from .clade import Clade
from .modifiers import UNITS, Drawn, Inherited


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

        birth = 1.0 * ScaledBy(Time(), {0: 1.0, 3: 0.3})     # a skyline

    Attached to the whole run, so every parameter may read it — the run is the coarsest unit, and a
    parameter's units always include it.
    """

    unit: str = "run"

    def __repr__(self) -> str:
        return "Time()"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Time)

    def __hash__(self) -> int:
        return hash(Time)


__all__ = ["UNITS", "Clade", "Drawn", "Inherited", "Measured", "Time"]
