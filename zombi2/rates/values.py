"""Values — the *what is read* half of a driven parameter (SPEC §5).

A parameter that is not a constant reads a **value**, and a value is two facts: what it is attached
to (its **unit** — the run, a lineage, a gene family) and how its number is **made**. Those two are
independent, and keeping them apart is what stops the grammar needing a new class per model::

    Drawn(per="family", spread=0.5)      # made by a draw,      attached to a family
    Drawn(per="lineage", spread=0.3)     # made by a draw,      attached to a lineage
    Inherited(per="lineage", spread=0.2) # made by inheritance, attached to a lineage
    Time()                               # measured from the run

The `verbs` module is the other half — what a parameter *does* with the number.

**These are the same six models the ``mod.`` names have always given you**, written as a grid rather
than as a list: ``Drawn(per="family", …)`` is ``ByFamily(…)``, ``Inherited(per="lineage", …)`` is
``FromParent(…)``, and both spellings build the identical object, so a run is unchanged and a script
written either way behaves the same. The old names stay; they are the names of cells.

The reason to prefer this spelling is what happens next. A per-chromosome draw needs no new class and
no new name here — it is ``Drawn(per="chromosome", …)``, a cell of a grid you can already see — where
under the old scheme it would need a ``ByChromosome`` invented, documented and remembered. A cell no
engine implements refuses at the level's gate, naming itself, rather than being silently ignored.
"""

from __future__ import annotations

from .modifiers import ByFamily, ByLineage, FromParent, Modifier

#: The units a value can be attached to, in the order they nest. A value's unit decides what may read
#: it: a parameter may read a value when the parameter's own units include the value's (SPEC §5), so a
#: trait on a lineage can drive gene loss, and a family's tempo cannot drive speciation.
UNITS = ("run", "lineage", "chromosome", "family", "copy", "site")

#: Which unit each *drawn* cell is built from today. A unit that is absent is not a different kind of
#: model — it is a cell nobody has needed yet, and adding one is a line here plus the engine threading
#: it, rather than a new modifier class.
_DRAWN = {"family": ByFamily, "lineage": ByLineage}

#: The same, for *inherited* values. Inheritance follows a genealogy, and the only genealogy an engine
#: walks today is the species tree, so "lineage" is the only cell.
_INHERITED = {"lineage": FromParent}


def _unit(per: str, built: dict, made: str) -> type:
    if not isinstance(per, str):
        raise TypeError(f"per= names the unit a value is attached to, as a string; got {per!r}")
    if per in built:
        return built[per]
    if per in UNITS:
        raise ValueError(
            f"a {made} value per {per!r} is not implemented. It is a cell of the grid rather than a "
            f"different kind of model, so it needs the engine to carry one number per {per}; today "
            f"{made} values exist per {' and per '.join(sorted(built))}.")
    raise ValueError(f"unknown unit {per!r}; a value is attached to one of {list(UNITS)}")


def Drawn(*, per: str, spread: float, dist: str = "lognormal") -> Modifier:
    """A value **drawn once per unit** and then fixed for that unit's life — i.i.d. heterogeneity,
    with no memory of anything else::

        loss = 0.25 * Drawn(per="family", spread=0.5)     # each family its own tempo
        substitution = 1.0 * Drawn(per="lineage", spread=0.3)   # the uncorrelated relaxed clock

    ``spread`` is the width and ``dist`` the shape — ``"lognormal"`` (σ on the log scale) or
    ``"gamma"`` (σ as the coefficient of variation). Both are mean-corrected, so widening ``spread``
    spreads the units out without moving the average one off the base.

    **Writing one object or two decides what varies together.** One `Drawn` read by
    several rates is one draw for that unit, so a family that loses fast also duplicates fast; two
    separately built ones are independent even with the same spread (SPEC §5).
    """
    return _unit(per, _DRAWN, "drawn")(spread=spread, dist=dist)


def Inherited(*, per: str, spread: float, bins: int | None = None) -> Modifier:
    """A value **inherited from a parent and perturbed** at each split — continuous memory down the
    tree, the autocorrelated twin of `Drawn`::

        birth = 1.0 * Inherited(per="lineage", spread=0.2)       # rate drift down a clade (ClaDS)
        substitution = 1.0 * Inherited(per="lineage", spread=0.3)  # the autocorrelated clock

    The per-split factor is lognormal and mean-corrected, so the rate does not creep upward down the
    tree. ``bins`` discretises the drift onto a ladder, where a daughter moves to a neighbouring rung
    or stays — the rate-category form of the same model.

    A value is either drawn afresh or inherited, never both on one unit: they are two answers to the
    same question, and a rate carrying both raises (SPEC §5).
    """
    return _unit(per, _INHERITED, "inherited")(spread=spread, bins=bins)


class Measured:
    """Base for a value the run already knows and can be asked for at any moment — no draw, no
    inheritance, nothing carried. Abstract: use `Time`."""

    unit: str = "run"


class Time(Measured):
    """The run's clock, as a value::

        birth = 1.0 * ScaledBy(Time(), {0: 1.0, 3: 0.3})     # a skyline

    Attached to the whole run, so anything can read it. Read through a mapping like any other value,
    which is what distinguishes it from `Drawn`: a draw is already a factor, while a
    time is a number that a mapping has to turn into one."""

    unit: str = "run"

    def __repr__(self) -> str:
        return "Time()"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Time)

    def __hash__(self) -> int:
        return hash(Time)


__all__ = ["UNITS", "Drawn", "Inherited", "Measured", "Time"]
