"""Scopes — the "per what?" of a rate (SPEC §5).

Every rate is written **from its scope**, so *per what?* is answered on the page rather than by a
default nobody types::

    birth        = PerLineage(0.5)   # each lineage speciates at 0.5 -> total = 0.5 × (lineages alive)
    birth        = Global(0.5)       # one shared budget for the whole tree -> total = 0.5 (constant)
    loss         = PerCopy(0.25)     # each gene copy is lost at 0.25 -> total = 0.25 × (copies present)
    substitution = PerSite(0.01)

Calling one **is** the rate: `PerCopy(0.25)` evaluates to a `Rate`, and the verbs chain onto it
(``PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))``). There is no intermediate object.

A scope carries no number of its own. It is a marker for the unit the base is counted per, and
`Rate.scope` holds the **class** rather than an instance — because while a scope carried a base the
same number lived in two places (``Rate.base`` and ``Rate.scope.base``) and nothing made them agree.

``PerCopy()`` — with no number — is what a `Rate.set_by` is written from: replacing *how fast* says
nothing about *per what*, so the scope still stands while the base does not.

The word *"per"* is reserved for these. A driver never starts with "per", and the unit a value
varies among is written with ``varying_among`` (SPEC §5).

There is deliberately **no** ``PerGenome``: one genome lives in one lineage, so "per genome" is
``PerLineage``.

A bare number (``birth = 1.0``) is coerced by each level to its natural default scope — species
birth/death and gene origination per lineage, duplication/transfer/loss per copy, substitution per
site. The scopes here are the explicit override, and the only spelling where more than one is legal.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, ClassVar

if TYPE_CHECKING:                      # `Rate` imports this module, so the real import is deferred
    from .parameter import Rate


class Scope:
    """The unit a rate's base is counted per, as a marker class.

    Abstract: use one of `Global`, `PerLineage`, `PerCopy`, `PerSite`, `PerChromosome`. Calling one
    builds a `Rate`; a ``Scope`` instance never exists, which is why `total_of` is a classmethod
    and `Rate.scope` holds the class.
    """

    #: the count keyword this scope multiplies by; ``None`` = a constant total (`Global`).
    unit: ClassVar[str | None] = None

    # The return type is deliberately not a `Scope`: `PerCopy(0.25)` *is* the rate, and there is no
    # second object for a level to coerce later. mypy flags a `__new__` returning something outside
    # the class, which is exactly the unusual thing being done here on purpose.
    def __new__(cls, base: float | None = None) -> "Rate":   # type: ignore[misc]
        from .parameter import Rate

        if cls is Scope:
            raise TypeError(
                "Scope is the abstract marker — write the one that says per what: PerLineage, "
                "PerCopy, PerSite, PerChromosome, or Global.")
        return Rate(base, cls, ())

    @classmethod
    def total_of(cls, base: float, **counts: Any) -> float:
        """The total rate for ``base``, given the current counts.

        ``counts`` supplies the units in scope right now (``lineages``, ``copies``, ``sites``,
        ``chromosomes``); each scope reads only the one it needs and ignores the rest. `Global`
        reads none.

        The base is an argument rather than a field because a `SetBy` supplies one: a scope answers
        *per what?*, and that question is unchanged when a driver replaces the number — "0.5 per
        copy in cave lineages" still multiplies by the copies present.
        """
        if cls.unit is None:
            return base
        try:
            return base * counts[cls.unit]
        except KeyError:
            raise KeyError(
                f"{cls.__name__} needs a {cls.unit!r} count; got {sorted(counts)}"
            ) from None


class Global(Scope):
    """One shared budget for the whole system: the total does not scale with anything.

    ``Global`` (capitalised — ``global`` is a Python keyword) makes a process run at a
    constant total rate: linear growth, not exponential.
    """

    unit: ClassVar[str | None] = None


class PerLineage(Scope):
    """Per lineage — the total scales with the number of lineages present.

    The default for species birth/death and gene origination. Within a single genome
    there is one lineage, so this reads as a constant per-genome budget; across the
    species tree it is ``base × (lineages alive)`` (exponential diversification).
    """

    unit: ClassVar[str | None] = "lineages"


class PerCopy(Scope):
    """Per gene copy — the total scales with family/genome size (duplication, transfer, loss).

    A large family therefore turns over faster: ``base × (copies present)``.
    """

    unit: ClassVar[str | None] = "copies"


class PerSite(Scope):
    """Per sequence site — the total scales with the number of sites (substitutions)."""

    unit: ClassVar[str | None] = "sites"


class PerChromosome(Scope):
    """Per chromosome — the total scales with the number of chromosomes (fission/fusion/loss)."""

    unit: ClassVar[str | None] = "chromosomes"


__all__ = ["Scope", "Global", "PerLineage", "PerCopy", "PerSite", "PerChromosome"]
