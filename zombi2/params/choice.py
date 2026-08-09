"""The **choice** — which candidate receives, rather than how often or how much (SPEC §5).

A rate says how often an event fires and an extent says how much it takes. A choice says *who gets
it*, and `transfer_to` is the only one. Its numbers are per-candidate **weights**, compared against
each other and normalised across the contemporaneous lineages, so they change neither how fast nor
how many transfers happen — only which lineage receives.

Four rules, and the two written here as classes are the **topological** ones: they read a fact about
the species tree rather than a value another level evolved.

- ``"uniform"`` — every contemporaneous lineage equally likely;
- `Distance` — closer relatives likelier, in units of tree depth;
- `Clades` — weight by the **pair** (donor's named clade, recipient's named clade);
- ``Recipients().weighted_by(driver, mapping)`` — weight by another level, which is the rest of the
  grammar.

They live beside the rate grammar and not with the transfer engine because they are things a user
**writes**, and everything a user writes belongs to one vocabulary: the same expression has to read
in Python, on the command line and in a ``--params`` file (SPEC §5, *one written form*). While they
sat at the genome level the parser could not see them, so `Distance(decay=…)` was Python-only and a
non-default decay could not be typed into a flag at all. Turning the *tree* into group membership
still belongs to the engine — that needs a tree, which a written rule does not.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import connection as verbs
from .mapping import Between
from .connection import Driven

#: The choice rules a `transfer_to` may be **written** with, for the parser's whitelist. The two
#: named strings (``"uniform"`` / ``"distance"``) are not here: they are values, not calls.
WRITABLE = ("Recipients", "Distance", "Clades")

__all__ = ["Choice", "Recipients", "Distance", "Clades", "WRITABLE"]


@dataclass(frozen=True, repr=False)     # repr=False: `__repr__` below is the written form
class Choice:
    """Which candidate receives — the parameter `Recipients()` opens and `weighted_by` fills in.

    A choice has no base: only the ratios between candidates are read, so there is nothing to write
    in front of the first verb. ``Recipients()`` on its own is the uniform rule, every
    contemporaneous lineage equally likely, which is what a choice with no weights says.
    """

    weights: tuple[Driven, ...] = ()

    def weighted_by(self, driver: object, mapping: object = None, *,
                    step: float | None = None) -> "Choice":
        """Weight the candidates by ``driver`` — see `verbs.weighted_by`.

        Weights **multiply and are then normalised across the candidates**, so chaining two is
        meaningful: prefer close relatives *and* run a highway between two distant clades. Whether a
        given engine reads more than one is that level's declaration.
        """
        return Choice(self.weights + (verbs.weighted_by(driver, mapping, step=step),))

    def scaled_by(self, driver: object, mapping: object = None, *,
                  step: float | None = None) -> "Choice":
        """Refused. A choice has no base to scale."""
        from .rate import RateCompositionError
        raise RateCompositionError(
            "a choice has no base to scale: its numbers are weights, compared against each other "
            "and normalised across the candidates. The verb is weighted_by — the same driver and "
            "the same mapping, read as a weight.")

    def set_by(self, driver: object, mapping: object = None, *,
               step: float | None = None) -> "Choice":
        """Refused. A choice has no base to replace."""
        from .rate import RateCompositionError
        raise RateCompositionError(
            "a choice has no base to replace: its numbers are weights, compared against each other "
            "and normalised across the candidates. The verb is weighted_by.")

    def __repr__(self) -> str:
        """The expression that produced this rule, written with no base in front — because a choice
        has none, and the flag refuses one there."""
        return "Recipients()" + "".join(f".{w.written_call()}" for w in self.weights)


def Recipients() -> Choice:
    """Open a ``transfer_to`` rule: the candidates that could receive a transfer::

        transfer_to = Recipients().weighted_by(competence, {"competent": 3.0, "normal": 1.0})

    On its own it is the uniform rule. It is a function rather than a class for the same reason a
    scope constructor returns a rate: what you get back is the parameter, ready for a verb.
    """
    return Choice()


@dataclass(frozen=True)
class Distance:
    """A ``transfer_to`` weighting by relatedness: a recipient at patristic distance ``d`` from the
    donor gets weight ``exp(-decay × d / depth)``, where ``depth`` is the tree's mean root-to-tip
    time — so ``decay`` is **scale-free** (in units of tree depth), meaning the same across trees of
    different absolute timescales. ``transfer_to="distance"`` is ``Distance(decay=1.0)``."""

    decay: float = 1.0

    def __post_init__(self) -> None:
        if isinstance(self.decay, bool) or not isinstance(self.decay, (int, float)) \
                or not math.isfinite(self.decay) or self.decay < 0:
            raise ValueError(f"Distance decay must be a finite non-negative number, got {self.decay!r}")


@dataclass
class Clades:
    """A ``transfer_to`` weighting by **named clades** — the topological, *donor-conditioned* sibling of
    `Distance`. Each group is a clade of the species tree, and a
    `Between` kernel weights a candidate recipient by the **pair** (donor's
    clade, recipient's clade), so a transfer can be steered to run *between* two clades rather than
    within them — which the per-recipient weight of a `Driven` cannot
    express::

        transfer_to = Clades({"A": ["n12", "n27"], "B": 40},
                             Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0))

    A clade is named either by a **set of tips** (a list — the clade is the subtree below their MRCA) or
    by a single **node id** (an int, or an ``"n<id>"`` label — the clade is that node's whole subtree).
    Groups must be disjoint; a lineage in none of them is in the implicit group ``"rest"``, usable as a
    kernel key. Membership is read from the **tree** (a clade is a fact about the tree, not another
    level), so this is a topological rule like ``"distance"``, resolved once per run — **not** a
    ``weighted_by`` reading another level, and needing no driver file."""

    groups: dict
    between: object

    def __post_init__(self) -> None:
        if not isinstance(self.groups, dict) or not self.groups:
            raise ValueError(
                "Clades needs a non-empty {label: clade} dict, where a clade is a list of tips (its "
                "MRCA's subtree) or a single node id — e.g. Clades({'A': ['n1', 'n2'], 'B': 40}, ...)")
        for label in self.groups:
            if not isinstance(label, str) or not label.strip():
                raise ValueError(f"clade labels must be non-empty strings, got {label!r}")
            if label == "rest":
                raise ValueError(
                    "'rest' is reserved for lineages in no named clade — name your clade something else")
        if isinstance(self.between, dict):
            self.between = Between(self.between)
        if not isinstance(self.between, Between):
            raise ValueError(
                "Clades takes a Between kernel (or a plain {(from, to): weight} dict) as its second "
                f"argument — the per-pair recipient weights — got {self.between!r}")
        unknown = self.between.groups() - (set(self.groups) | {"rest"})
        if unknown:
            raise ValueError(
                f"the Between kernel names groups {sorted(unknown)} that are not defined clades; "
                f"defined clades are {sorted(self.groups)} (plus the implicit 'rest')")

    def __repr__(self) -> str:
        """The written form: the constructor call the API and the parser both take, not the
        ``Clades(groups=…, between=…)`` a dataclass gives. A run's log records this, and a record
        the flag will not read back is not a record."""
        return f"Clades({self.groups!r}, {self.between!r})"
