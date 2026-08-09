"""Values — the *what is read* half of a driven parameter (SPEC §5).

A parameter that is not a constant reads a **value**, and a value is two facts: what it is attached
to (its **unit** — the run, a lineage, a gene family) and how its number is **made**. Those two are
independent, and keeping them apart is what stops the grammar needing a new class per model::

    Random("families", LogNormal(0.0, 0.5))          # made by a draw,      among families
    Random("lineages", LogNormal(0.0, 0.3))          # made by a draw,      among lineages
    Random("lineages", Drift(LogNormal(0.0, 0.2)))   # made by inheritance, among lineages
    Time()                                # measured from the run
    Clade({"fast": [...]})                # read off the tree itself

The `verbs` module is the other half — what a parameter *does* with the number.

`Drawn` and `Inherited`, the two objects `Random` builds, live in `modifiers` beside the classes an
engine dispatches on, because that is what they are; they are re-exported here so the grammar can be
imported as one surface. `Time` is here because it is a value that is **not** a modifier: it is read
through a verb (``changing_at({0: 1.0, 3: 0.3})``, ``set_by(Time(), …)``), where a draw is already a
factor and needs none.

**The unit is an argument, not a class.** A draw among families and a draw among lineages are one
model at two attachments, so ``Random("families", …)`` and ``Random("lineages", …)`` are one class
and two cells of a grid rather than two names to remember. The cells are models the field already
has names for — rate heterogeneity across gene families, the relaxed clock, ClaDS — and the prose
uses those.

What the grid buys is the next cell. A per-chromosome draw needs no new class and no invented name:
it is ``Random("chromosomes", …)``, which constructs, and a cell no engine carries refuses at that
level's gate naming itself — *"drawn among chromosomes, which the species engine does not
support"* — rather than reading as a typo.
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


#: The values a rate may be **written** with. `Measured` is absent because it is abstract, and
#: `Drawn` / `Inherited` because they are modifiers and `modifiers.WRITABLE` already lists them.
WRITABLE = ("Clade", "Time")

__all__ = ["UNITS", "WRITABLE", "Clade", "Drawn", "Inherited", "Measured", "Time"]
