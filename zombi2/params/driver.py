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

from dataclasses import dataclass

import math
from typing import Any, ClassVar, Mapping

from .evaluate import (CHANGING_AT, MEASURED, SCALED_BY, SET_BY, UNITS, VARYING_UNITS,
                      Modifier, _WRITTEN_AS)
from .law import Drawn, Drift, Inherited
from .retired import check_no_retired_keywords




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

    def resolve(self, tree) -> dict[str, list[int]]:
        """Which lineages each named group covers on ``tree`` — ``{label: [node ids]}``, in the order
        the groups were written, with ``"rest"`` last when some lineage is in no named clade::

            Clade({"fast": ["n27", "n51"]}).resolve(tree)
            # {'fast': [24, 27, 28, 51, 52], 'rest': [0, 1, 2, ...]}

        A read-back, for checking a clade **before** trusting a run that reads it. Three things it
        makes visible, none of them derivable from the tip names alone: the **MRCA's own branch is
        inside the clade** (``n24`` above, which nobody named), a clade holds the extinct and
        internal lineages of its subtree as well as its tips, and a lineage in no named clade is in
        ``"rest"``.

        It cannot disagree with the run, because it paints with `resolve_groups` — the one function
        the engine paints membership with, for this driver and for the `Clades` transfer rule alike.
        ``zombi2 tools tree TREE --clades`` is the other half: it lists the clades a tree offers to
        name (Appendix C)."""
        from ..genomes._transfer import resolve_groups
        from ..tree import as_tree

        painted = resolve_groups(as_tree(tree, level="clade"), self.groups)
        covers: dict[str, list[int]] = {label: [] for label in self.groups}
        for i in sorted(painted):
            covers.setdefault(painted[i], []).append(i)
        return covers

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

class OnTime(Modifier):
    """The rate changes in time — a skyline / episodic schedule. Written ``changing_at``::

        birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})    # 1.0 on [0, 3), then 0.3 on

    ``schedule`` maps each interval's start time to a relative factor, dimensionless: on a base of
    ``2.0`` the schedule scales it. Before the earliest breakpoint the earliest factor applies
    (define the schedule from time 0 to avoid surprise).

    The **same object** carries the other half of `Time`, ``set_by(Time(), ...)``, where the
    schedule holds the rates themselves rather than multiples of a base::

        birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})    # 30% of what it was
        birth = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})   # the rates themselves

    Those two are the same model, because a schedule on a base of 1.0 *is* the rate — which is why
    `Rate.set_by` builds this rather than a `SetBy` and there is nothing for an engine to learn.
    What differs is the sentence the reader typed, so `verb` records which, exactly as a `Driven`
    records its own: a run's log has to say back what was written, not an equivalent.
    """

    reads: ClassVar[tuple[str, str] | None] = (MEASURED, "run")

    #: `CHANGING_AT`, or `SET_BY` when the schedule holds the rates rather than factors.
    verb: str = CHANGING_AT

    def __init__(self, schedule: Mapping[float, float], *, verb: str | None = None) -> None:
        steps = tuple(sorted((float(t), float(f)) for t, f in schedule.items()))
        if not steps:
            raise ValueError(
                "a time schedule cannot be empty, e.g. .changing_at({0: 1.0, 3: 0.3})")
        for t, f in steps:
            if not math.isfinite(t):
                raise ValueError(f"schedule times must be finite, got {t!r}")
            if not math.isfinite(f) or f < 0:
                raise ValueError(f"schedule values must be finite and non-negative, got {f!r}")
        self._steps = steps
        if verb is not None:
            self.verb = verb

    # `factor` narrows the base signature: this modifier cannot answer without its key, and
    # giving it a default would make a level that forgot to thread it return a plausible
    # wrong number in silence. `IMPLEMENTED_MODIFIERS` is what guarantees the key arrives — a level
    # that does not thread it rejects the modifier outright rather than reaching here.
    def factor(self, *, time: float, **_: Any) -> float:  # type: ignore[override]
        f = self._steps[0][1]  # before the first breakpoint, the earliest factor applies
        for t, fac in self._steps:
            if t <= time:
                f = fac
            else:
                break
        return f

    def next_change(self, time: float) -> float:
        for t, _ in self._steps:  # steps are sorted; the first breakpoint strictly after `time`
            if t > time:
                return t
        return math.inf

    def written_call(self) -> str:
        # `repr(float)` rather than `:g`, which rounds to six significant figures: this text is what
        # a run's log records and what a reader pastes back into a flag, so a rate that prints
        # rounded is a log naming a model that is not the one that ran.
        inner = ", ".join(f"{float(t)!r}: {float(f)!r}" for t, f in self._steps)
        if self.verb == SET_BY:
            return f"{SET_BY}(Time(), {{{inner}}})"
        return f"{CHANGING_AT}({{{inner}}})"

    def __eq__(self, other: object) -> bool:
        """By the schedule alone. `verb` is outside this on purpose: the two spellings run the same
        model, and two rates that behave identically should compare equal."""
        return isinstance(other, OnTime) and other._steps == self._steps

    def __hash__(self) -> int:
        return hash((OnTime, self._steps))


@dataclass(frozen=True)
class TotalDiversity:
    """The lineages standing right now, as a **driver** a rate can be scaled by::

        birth = PerLineage(1.0).scaled_by(TotalDiversity(cap=100))

    The carrying capacity is written on the driver rather than in a mapping beside it, and that is
    a limit rather than a design: the factor an engine reads is the linear fall to a cap, and a
    general curve of diversity would have to be integrated rather than read at a point, exactly as
    a smooth function of time would (SPEC §5). So there is one shape, and it takes its one number
    here. `scaled_by` refuses a mapping alongside, rather than accepting one and ignoring it.
    """

    cap: float | None = None

    def __post_init__(self) -> None:
        if self.cap is None:
            raise ValueError(
                "TotalDiversity needs its carrying capacity — TotalDiversity(cap=100), the "
                "diversity at which the rate reaches zero. The linear fall to a cap is the one "
                "shape an engine reads, so there is no mapping to write instead.")
        if isinstance(self.cap, bool) or not isinstance(self.cap, (int, float)):
            raise TypeError(f"TotalDiversity cap must be a real number, got {self.cap!r}")
        if not math.isfinite(self.cap) or self.cap <= 0:
            raise ValueError(f"TotalDiversity cap must be finite and positive, got {self.cap!r}")

    def __repr__(self) -> str:
        assert self.cap is not None            # __post_init__ refuses a driver without its cap
        return f"TotalDiversity(cap={float(self.cap)!r})"


@dataclass(frozen=True, repr=False)            # repr=False: the written form, from `Modifier`
class OnTotalDiversity(Modifier):
    """The rate slows as standing diversity grows — diversity-dependence. Built by
    ``scaled_by(TotalDiversity(cap=100))``.

    The factor falls linearly from 1 toward 0 as diversity rises to ``cap`` (a carrying
    capacity), and stays 0 beyond it: a cap of 100 halves the rate at 50 lineages and stops it
    at 100.
    """

    reads: ClassVar[tuple[str, str] | None] = (MEASURED, "run")

    cap: float

    def __post_init__(self) -> None:
        if isinstance(self.cap, bool) or not isinstance(self.cap, (int, float)):
            raise TypeError(f"TotalDiversity cap must be a real number, got {self.cap!r}")
        if not math.isfinite(self.cap) or self.cap <= 0:
            raise ValueError(f"TotalDiversity cap must be finite and positive, got {self.cap!r}")

    # `factor` narrows the base signature: this modifier cannot answer without its key, and
    # giving it a default would make a level that forgot to thread it return a plausible
    # wrong number in silence. `IMPLEMENTED_MODIFIERS` is what guarantees the key arrives — a level
    # that does not thread it rejects the modifier outright rather than reaching here.
    def factor(self, *, diversity: float, **_: Any) -> float:  # type: ignore[override]
        return max(0.0, 1.0 - diversity / self.cap)

    def written_call(self) -> str:
        return f"{SCALED_BY}(TotalDiversity(cap={float(self.cap)!r}))"


def Random(unit: str | None = None, law: object = None, **retired: Any) -> Modifier:
    """A value drawn for each unit of that kind — the one driver that is not measured anywhere.

    ``unit`` is plural (``'lineages'``, ``'families'``, ``'copies'``, ``'sites'``,
    ``'chromosomes'``). The **law** says what happens to the value afterwards, which is a separate
    question from what it starts as::

        Random('families', LogNormal(0.0, 0.5))          # drawn once, held for that family's life
        Random('lineages', Drift(LogNormal(0.0, 0.3)))   # the parent's, perturbed at each split
        Random('lineages', Drift(LogNormal(0.0, 0.3), bins=8))    # the rate-category clock

    A bare distribution is deliberate, not an oversight: it follows the convention the grammar
    already uses everywhere — a bare dict is a table, a bare function a curve — where the plain case
    is written plainly and anything else is named.

    Usually written through the verb, ``rate.varying_among('families', law)``, which builds one of
    these and attaches it. Building it **by name** is how two rates share one draw::

        family_speed = Random('families', LogNormal(0.0, 0.5))
        duplication  = PerCopy(0.20).varying_among(family_speed)
        loss         = PerCopy(0.10).varying_among(family_speed)   # exactly half, in every family

    because the engine caches a unit's draw by object identity. Two separately built ``Random``
    objects are two draws even with identical arguments: the question is whether you wrote one
    thing or two.

    ``**retired`` catches ``spread=`` and ``per=``, the two keywords this replaced, so Python
    answers them with the same sentence a flag does rather than with "unexpected keyword argument".
    The unit has a default for that reason alone: ``Random(per='family', …)`` writes it into a
    keyword, and a required positional would make Python complain about the missing argument before
    anything here could say what ``per=`` became.
    """
    check_no_retired_keywords(retired, where="Random")
    if unit is None:
        raise TypeError(
            f"a Random needs the plural unit its value varies among — Random('families', "
            f"LogNormal(0.0, 0.5)); one of {list(VARYING_UNITS)}.")
    if isinstance(law, Drift):
        return Inherited(unit, law.dist, law.bins)
    return Drawn(unit, law)


_WRITTEN_AS.update({OnTime: CHANGING_AT,
                    OnTotalDiversity: f"{SCALED_BY}(TotalDiversity(...))"})


#: The drivers a parameter may be **written** with. `Measured` is absent because it is
#: abstract, and `Drawn` / `Inherited` because a law names them, not a driver.
WRITABLE = ("Clade", "Time", "Random", "TotalDiversity")

__all__ = ["UNITS", "WRITABLE", "Clade", "Measured", "OnTime", "OnTotalDiversity", "Random",
           "Time", "TotalDiversity"]
#: `Drawn` and `Inherited` are imported above because `Random` builds them, and are deliberately
#: absent here: they belong to `law`, which is where the law that chooses between them lives.
#: One concept, one public home — re-exporting them was a habit carried over from `values.py`.
