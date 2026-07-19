"""Scope wrappers — the "per what?" of a rate (SPEC §5).

Every rate in ZOMBI2 is ``scope(base) × modifiers``. A *scope wrapper* answers
*per what?*: it tags a base rate with the unit it applies to, so the engine can
turn that per-unit base into a **total** rate by multiplying by however many of
the unit are present right now::

    birth = scope.PerLineage(1.0)   # each lineage speciates at 1.0 -> total = 1.0 × (lineages alive)
    birth = scope.Global(1.0)       # one shared budget for the whole tree -> total = 1.0 (constant)
    loss  = scope.PerCopy(0.25)     # each gene copy is lost at 0.25 -> total = 0.25 × (copies present)

Wrappers **wrap** a base; they do **not** multiply — multiplying is what
modifiers (``zombi2.modifiers``) do. The word *"per"* is reserved for these
wrappers; a modifier never starts with "per".

There is deliberately **no** ``PerGenome``: one genome lives in one lineage, so
"per genome" is ``PerLineage``.

A bare number (``birth = 1.0``) is coerced by each level to its natural default
wrapper — species birth/death and gene origination default to per lineage,
duplication/transfer/loss to per copy, substitution to per site. The wrappers
here are the explicit override.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar


@dataclass(frozen=True)
class Scope:
    """A base rate tagged with the unit it applies to.

    Abstract: use one of :class:`Global`, :class:`PerLineage`, :class:`PerCopy`,
    :class:`PerSite`, :class:`PerChromosome`.
    """

    base: float

    #: the :meth:`total` keyword this wrapper multiplies by; ``None`` = a constant total.
    unit: ClassVar[str | None] = None

    def __post_init__(self) -> None:
        if isinstance(self.base, bool) or not isinstance(self.base, (int, float)):
            raise TypeError(f"a rate base must be a real number, got {self.base!r}")
        if not math.isfinite(self.base) or self.base < 0:
            raise ValueError(f"a rate base must be finite and non-negative, got {self.base!r}")

    def total(self, **counts: float) -> float:
        """The total rate given the current counts.

        ``counts`` supplies the units in scope right now (``lineages``, ``copies``,
        ``sites``, ``chromosomes``); each wrapper reads only the one it needs and
        ignores the rest. :class:`Global` reads none.
        """
        if self.unit is None:
            return self.base
        try:
            return self.base * counts[self.unit]
        except KeyError:
            raise KeyError(
                f"{type(self).__name__} needs a {self.unit!r} count; got {sorted(counts)}"
            ) from None

    def __mul__(self, other: object):
        # composing a scope with a modifier builds a Rate (internal plumbing, see zombi2.rate)
        from .modifiers import Modifier
        from .rate import Rate

        if isinstance(other, Modifier):
            return Rate(self.base, self, (other,))
        return NotImplemented


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
