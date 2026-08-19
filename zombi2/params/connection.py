"""The **link** — what joins a driver to the parameter that reads it (SPEC §5, §7).

The manual calls the whole written dependency a *connection* (driver, link, target); this module
is the link, the verb-and-mapping part. The filename predates that vocabulary.

A parameter is written from its scope and a driver says what it reads. What sits between them is
this: a verb saying what the number *does*, and the object that verb builds.

    scaled_by    multiplies the base                      -> Driven
    set_by       replaces it, in the parameter's own units -> SetBy
    weighted_by  compares the candidates of a choice       -> Driven, recorded as a weight

Two drivers keep a verb of their own because they are written constantly — ``varying_among`` for a
`Random` and ``changing_at`` for the clock — and `scaled_by` refuses both by name, so there is
exactly one spelling for each.

Which verb is legal is a fact about the **parameter**, not about the driver: a rate can be scaled or
replaced, an extent only scaled, a choice only weighted. That is why the verbs are methods over in
`parameter`, and why what they build lives here rather than there — one file for the joining, one for
the things being joined.
"""


from __future__ import annotations

import math
from typing import Any, ClassVar, Mapping

from .driver import OnTime, OnTotalDiversity, Random, TotalDiversity
from .evaluate import CHANGING_AT, DRIVEN, Modifier, SCALED_BY, SET_BY, VARYING_AMONG, WEIGHTED_BY, _WRITTEN_AS, _driver_form, describe
from .law import Drawn, Drift, Inherited
from .retired import check_no_retired_keywords
from .driver import Measured, Time

#: The verb names, and so the attribute names `zombi2.params.parse` will follow in an expression.
#: The written form is a call on an attribute, so this is the parser's whitelist for the attribute
#: half — the other half being `_NAMES`, the things it may call by name.
#:
#: Built from the constants defined in `zombi2.params.evaluate`, beside the ``verb`` field that
#: records which verb was written, and re-exported here so a reader outside this package compares
#: against ``verbs.SCALED_BY`` rather than against a string it spelled out itself.
VERBS = (SCALED_BY, SET_BY, WEIGHTED_BY, VARYING_AMONG, CHANGING_AT)

__all__ = ["VERBS", "SCALED_BY", "SET_BY", "WEIGHTED_BY", "VARYING_AMONG", "CHANGING_AT",
           "written_with", "scaled_by", "set_by", "weighted_by", "varying_among", "changing_at"]


def written_with(m: object, verb: str) -> bool:
    """Whether ``m`` records ``verb`` as the verb that wrote it.

    The same object serves several verbs — a `Driven` is what `scaled_by` and `weighted_by` both
    build — so which one was typed is a fact only the object remembers, and it is what tells a
    mismatched verb from a right one. Asking through here rather than by reaching for ``.verb``
    keeps the reading in one place and works for a modifier that records none.
    """
    return getattr(m, "verb", None) == verb


def _refuse_a_factor(driver: object, verb: str) -> None:
    """Something that is already a dimensionless factor, handed to a verb as though it were a
    driver. A `Random` has its own verb; anything else here is what a verb *builds*, so in neither
    case is there a value for a mapping to turn into a number."""
    if isinstance(driver, (Drawn, Inherited, Drift)):
        raise TypeError(
            f"a Random is already a factor and has its own verb: write "
            f".varying_among('families', LogNormal(0.0, 0.5)) rather than {verb}(...). Verbs are "
            f"for drivers a mapping has to turn into a number — a trait, a clock, a count.")
    if isinstance(driver, Modifier):
        raise TypeError(
            f"{describe(driver)} is what a verb builds, not something a verb reads: it is already a "
            f"dimensionless factor, so {verb}(...) has nothing to map. Write the verb that produces "
            f"it on the parameter itself.")


def _refuse_time(driver: object, verb: str) -> None:
    """`Time` handed to a general verb. Its two readings are `changing_at` (factors) and
    ``set_by(Time(), ...)`` (the rates themselves), and those are the only two spellings."""
    if isinstance(driver, Time):
        raise ValueError(
            f"the run's clock has its own verb: write .changing_at({{0: 1.0, 3: 0.3}}) rather than "
            f"{verb}(Time(), ...). If you meant the schedule to hold the rates themselves rather "
            f"than multiples of a base, that is PerLineage().set_by(Time(), {{0: 0.5, 3: 0.15}}).")


def _refuse_a_whole_rule(driver: object) -> None:
    """`Distance` or `Clades` handed to `weighted_by` as though it were a driver.

    Both are complete ``transfer_to`` rules — a rule read off the tree, mapping and all — rather
    than a value a mapping turns into a weight, and each is written on its own. Passing one as a
    driver builds a `Driven` whose lookup key no engine ever fills, and which renders as
    ``weighted_by(<Distance>, ...)``: a run's log recording a rule nobody can paste back, and one
    the round-trip check cannot catch because the text never parses. Splitting `Distance` into a
    driver and a mapping of its own is designed and not built, so this refuses rather than accepting
    something that would quietly weight nothing."""
    from .choice import Clades, Distance

    if isinstance(driver, Distance):
        raise ValueError(
            "Distance() is a whole transfer_to rule, not a driver a mapping reads: write "
            "transfer_to=Distance(decay=1.0), which weights a recipient by exp(-decay × d / depth). "
            "Choosing the shape of that fall yourself — weighted_by(Distance(), lambda d: 1 / (1 + "
            "d ** 2)) — is not implemented; today the decay is the one knob.")
    if isinstance(driver, Clades):
        raise ValueError(
            "Clades({...}, Between({...})) is a whole transfer_to rule, not a driver: it already "
            "carries the kernel that weights the (donor, recipient) clade pair, and it is written "
            "on its own — transfer_to=Clades({'A': ['n12'], 'B': ['n40']}, "
            "Between({('A', 'B'): 10.0})). "
            "weighted_by reads a value another level evolved, and clade membership is a fact about "
            "the tree.")


def _schedule(mapping: object, verb: str) -> dict:
    """A time driver's mapping. A schedule — ``{0: 1.0, 3: 0.3}``, the value from each breakpoint on
    — is the only shape an engine reads, because a rate that is piecewise-constant in time can be
    stepped to exactly. A smooth curve of time is a real model and is refused rather than
    approximated: it makes the rate vary continuously, which needs the engine to integrate its
    hazard rather than sample it."""
    if isinstance(mapping, dict):
        return mapping
    raise ValueError(
        f"{verb} takes a schedule — {{0: 1.0, 3: 0.3}}, the value from each breakpoint on — because "
        f"a rate that changes in steps can be stepped to exactly. A smooth function of time is not "
        f"implemented: it makes the rate vary continuously between events, which needs the engine "
        f"to integrate the rate rather than read it at a point.")


def scaled_by(driver: object, mapping: object = None, *, step: float | None = None) -> Modifier:
    """Multiply the parameter's base by a factor read from ``driver``.

    The factor is dimensionless, and almost every parameter takes one::

        loss  = PerCopy(0.25).scaled_by(habitat, {"cave": 4.0, "surface": 1.0})   # a grown trait
        birth = PerLineage(1.0).scaled_by(TotalDiversity(cap=100))                # the standing LTT

    ``mapping`` turns the driver's value into that factor, and its shape follows the driver's
    **type**: a categorical driver takes a table (a dict), a numerical one a curve (a callable) or a
    ``Scalar`` log-link.

    ``step`` is the resolution a **continuous** driver is read at, in the tree's own time units. A
    categorical driver switches at moments the engine can step to exactly and ignores it.
    """
    _refuse_time(driver, "scaled_by")
    _refuse_a_factor(driver, "scaled_by")
    if isinstance(driver, TotalDiversity):
        if mapping is not None:
            raise ValueError(
                "TotalDiversity carries its own shape — the linear fall to a cap — so there is no "
                "mapping to write beside it: scaled_by(TotalDiversity(cap=100)). A general curve "
                "of standing diversity is not implemented (SPEC §5).")
        assert driver.cap is not None            # TotalDiversity refuses a driver without its cap
        return OnTotalDiversity(driver.cap)
    if isinstance(driver, Measured):
        raise ValueError(
            f"scaled_by({type(driver).__name__}(), ...) is not implemented — that driver exists in "
            f"the grammar but no engine supplies it yet.")
    if mapping is None:
        raise ValueError(
            "scaled_by(driver, mapping) needs a mapping: a dict for a categorical driver, a "
            "callable for a numerical one.")
    return Driven(driver, mapping, step, verb=SCALED_BY)


def set_by(driver: object, mapping: object = None, *, step: float | None = None) -> Modifier:
    """Replace the parameter's base with a number read from ``driver``, in the parameter's own
    units::

        loss  = PerCopy().set_by("habitat.tsv", {"aquatic": 1.0, "terrestrial": 0.25})
        birth = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})

    Written with no base in front, because the driver supplies the whole number; the scope still
    stands, because replacing *how fast* says nothing about *per what*.

    The clock is the one driver that can replace as well as scale, and it needs no new machinery to:
    a schedule on a base of 1.0 *is* the rate, so this builds the same `OnTime` `changing_at` does
    and records which verb wrote it (`OnTime.verb`), so a run's log says back what was typed.
    """
    if isinstance(driver, Time):
        return OnTime(_schedule(mapping, "set_by(Time(), ...)"), verb=SET_BY)
    _refuse_a_factor(driver, "set_by")
    if isinstance(driver, (TotalDiversity, Measured)):
        raise ValueError(
            f"set_by({type(driver).__name__}(), ...) is not implemented — no engine can be handed "
            f"a base from that driver. Scale a base you write yourself instead: "
            f"scaled_by({type(driver).__name__}(...)).")
    if mapping is None:
        raise ValueError(
            "set_by(driver, mapping) needs a mapping: a dict for a categorical driver, a callable "
            "for a numerical one. Its numbers are the rate itself, not factors.")
    return SetBy(driver, mapping, step, verb=SET_BY)


def weighted_by(driver: object, mapping: object = None, *, step: float | None = None) -> Driven:
    """Weight the candidates of a **choice** — an argument that decides *who*, not how fast.

    ``transfer_to``, the recipient of a horizontal transfer, is the only choice today. A choice has
    no base, because only the ratios between candidates are read, and a weight of zero means that
    candidate cannot be chosen::

        transfer_to = Recipients().weighted_by(competence, {"competent": 3.0, "normal": 1.0})

    A weight may read **both ends** — the donor's group and the recipient's — through a ``Between``
    kernel, which is the mapping for a driver that sits on a pair rather than on one lineage.
    """
    _refuse_time(driver, "weighted_by")
    _refuse_a_factor(driver, "weighted_by")
    _refuse_a_whole_rule(driver)
    if isinstance(driver, (TotalDiversity, Measured)):
        raise ValueError(
            f"weighted_by({type(driver).__name__}(), ...) is not implemented — a choice weights "
            f"each candidate by something that candidate has, and that driver is a property of the "
            f"run rather than of a lineage, so every candidate would weigh the same.")
    if mapping is None:
        raise ValueError(
            "weighted_by(driver, mapping) needs a mapping: a dict of per-candidate weights, a "
            "callable, or a Between kernel to weight the (donor, recipient) pair.")
    return Driven(driver, mapping, step, verb=WEIGHTED_BY)


def varying_among(among: object = None, law: object = None, **retired: object) -> Modifier:
    """Let the parameter vary at random among the units of one kind (SPEC §5)::

        loss = PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))
        rate = PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))

    ``among`` is the plural unit name and ``law`` says what happens to the drawn value afterwards —
    a bare distribution for a value drawn and held, a `Drift` for one carried down the tree and
    perturbed at each split.

    It also takes a **named** `Random`, with no second argument, which is how two rates share one
    draw: the engine caches a unit's value by object identity, so one object read twice is one
    number and two objects are two.

    ``**retired`` catches ``per=`` and ``spread=``, the keywords this verb replaced, so Python
    answers them with the same sentence a flag does. Every parameter's ``varying_among`` passes its
    own through to here, so the answer cannot be good at one level and absent at another. The unit
    has a default for that reason alone: ``varying_among(per='families', …)`` writes the unit into a
    keyword, and a required positional would make Python complain about the missing argument before
    anything here could say what ``per=`` became.
    """
    check_no_retired_keywords(retired, where="varying_among")
    if among is None:
        raise TypeError(
            "varying_among takes the plural unit to vary among and the law it follows — "
            "varying_among('families', LogNormal(0.0, 0.5)) — or a Random built by name, "
            "varying_among(family_speed), which is how two rates share one draw.")
    if isinstance(among, (Drawn, Inherited)):
        if law is not None:
            raise TypeError(
                "a named Random already carries its law, so varying_among takes it alone: "
                "varying_among(family_speed). Sharing one object is what makes two rates share one "
                "draw — build a second one and there is nothing to share.")
        return among
    if isinstance(among, Modifier):
        raise TypeError(
            f"varying_among takes a Random — a unit and a law — and {describe(among)} is not one.")
    if not isinstance(among, str):
        raise TypeError(
            f"varying_among takes the plural unit to vary among and the law it follows — "
            f"varying_among('families', LogNormal(0.0, 0.5)) — or a Random built by name. "
            f"Got {among!r}.")
    return Random(among, law)


def changing_at(schedule: object) -> OnTime:
    """Let the parameter change in time — a skyline, the run's clock read as a schedule of factors::

        birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})   # 1.0, then 30% of it from time 3

    The numbers are multiples of the base. For the other reading — the schedule holding the rates
    themselves — write ``set_by(Time(), ...)``, which builds the same thing on a base of 1.0.
    """
    return OnTime(_schedule(schedule, "changing_at"))


class Driven(Modifier):
    """The factor is read from **another evolved value** — the one mechanism behind both conditioning and
    joining (SPEC §2).

    **You do not write this class; you write a verb.** Which verb says what the number does, and is
    decided by what you are attaching it to: `scaled_by` multiplies a rate or an extent,
    `weighted_by` compares the candidates of a choice, `set_by` replaces a base. All three build
    this::

        loss        = PerCopy(0.25).scaled_by("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})
        birth       = PerLineage(1.0).scaled_by("trait", {"small": 1.0, "large": 2.0})  # joint
        transfer_to = Recipients().weighted_by("competence.tsv", {"competent": 3.0})

    It is Chapter 8's definition made literal: *a rate that reads a value which varies from lineage
    to lineage, rather than a fixed number*. It reads the driver's value on each lineage and the
    mapping turns it into a number.

    ``driver`` says where the driven value comes from, and that single choice splits *conditioned*
    from *joint* — the chapter's spine, *can the driver be grown first?*:

    - a **filename** (``"habitat.tsv"``), a **grown driver result** (a ``TraitsResult``, discrete or
      continuous), or a genome result's ``presence(...)`` / ``completion(...)`` —
      the driver was grown first and handed over (**conditioned**): two ordinary runs. The result
      object is the file's in-memory shortcut — same conditioning, no ``write``/read step;
    - a **level name** (``"trait"``, ``"genomes:count"``) — the driver co-evolves in one run
      (**joint**): neither level can be grown first.

    ``mapping`` says how the driver's value becomes the factor — a `Table`
    (a dict, for a discrete driver), a `Curve` (a callable, continuous),
    a `Scalar` (a log-link coefficient), or a `Between` (a weight per donor/recipient pair, which
    only ``transfer_to`` takes); a raw dict / callable / number is coerced (`as_mapping()`).

    What a ``Driven`` can be attached to comes in **three kinds**, and only the first is a rate:
    *how often* an event fires (a rate, e.g. ``loss``), *how much* it takes (an extent, e.g.
    ``loss_extent``, at the ordered and nucleotide resolutions), and a **choice** of who receives it
    (``transfer_to``, a weight per candidate rather than a multiplier). It always maps a value to a
    number; it never drives a *value*, such as an OU optimum.

    Like a carried modifier, a ``Driven`` reads
    a value the **engine** threads per lineage — here a ``drivers`` mapping ``{key: value}`` — and
    is otherwise dumb: it just maps the value to a factor. The engine owns *where* the value comes from
    (a file it loaded, or the live level growing beside the tree) and *when* it changes (a discrete
    driver switches mid-branch, so the engine steps its Gillespie at each switch); a rate reaching an
    engine that has not threaded its ``driver`` gets a factor of 1.0 (inert).
    """

    #: A driven value is *not* carried: the engine resolves it per lineage into ``drivers`` and this
    #: modifier maps it, whether it was recorded beforehand (conditioned) or is growing alongside
    #: (joint). Which of those it is depends on the ``driver`` argument, not on the class, so the
    #: kind here is the pair's shared name rather than one of them.
    reads: ClassVar[tuple[str, str] | None] = (DRIVEN, "lineages")

    #: The verb that wrote this one, and so how it prints. A driven value is written with the verb
    #: that says what its number does — and only the writer knows which, because it follows from
    #: what the value is attached to, not from anything the value itself holds. Recording it is what
    #: lets a run's log say back exactly what was typed.
    verb: str = SCALED_BY

    def __init__(self, driver: object, mapping: object, step: float | None = None, *,
                 verb: str | None = None) -> None:
        from .mapping import as_mapping

        if isinstance(driver, str):
            if not driver.strip():
                raise ValueError("a driver must be a non-empty string (a filename or level name)")
            base: object = driver                    # a string driver is its own context key
        else:
            base = id(driver)                        # an in-memory driver result (conditioning): key by identity
        if step is not None:
            step = float(step)
            if not (step > 0.0) or step == float("inf"):
                raise ValueError(
                    f"step is the resolution a CONTINUOUS driver is read at, in the tree's own "
                    f"time units, so it must be finite and positive; got {step!r}.")
        # the step is part of the key: the same driver read at two resolutions is two trajectories, and
        # keying on the driver alone would silently resolve it once and share the first one
        self.key: object = base if step is None else (base, step)
        self.driver = driver
        self.step = step
        self.mapping = as_mapping(mapping)
        if verb is not None:
            self.verb = verb

    def factor(self, *, drivers: Mapping | None = None, time: float | None = None,
               **_: Any) -> float:
        """The mapped multiplier for this lineage's driver value — the engine threads the value under
        ``drivers[key]`` (``key`` is the driver string, or the identity of an in-memory driver). No
        ``drivers`` (or this driver absent) ⇒ 1.0, so an unthreaded rate is inert (the engine is
        responsible for supplying the value where a driven rate is supported).

        ``time`` rides in from the same context every engine already passes, and is used only where
        the mapping has a `Schedule` entry — this driver state, but only after t."""
        if drivers is None:
            return 1.0
        value = drivers.get(self.key)
        if value is None:
            return 1.0
        return self.mapping.multiplier(value, time=time)

    def next_change(self, time: float) -> float:
        """A scheduled mapping entry changes on its own, so its breakpoints have to reach the
        engine's horizon or the Gillespie steps straight past them. Everything else answers ``inf``,
        which is `Modifier.next_change`'s default and what every mapping but `Table` returns."""
        nc = getattr(self.mapping, "next_change", None)
        return math.inf if nc is None else nc(time)

    def written_call(self) -> str:
        """``step`` is written whenever it is set, because a written form that omits an argument
        records a different model. Leaving it out meant a driven rate with a step rendered as one
        without, reparsed as one without, and compared equal to one without — so a run's log said
        something the run had not done, and every round-trip check agreed."""
        step = f", step={self.step!r}" if self.step is not None else ""
        return f"{self.verb}({_driver_form(self.driver)}, {self.mapping!r}{step})"

    def __eq__(self, other: object) -> bool:
        # By the **driver**, not by `key`. `key` is a runtime lookup handle: it is `(path, step)`
        # for a file and `id()` for a driver that is an object — so comparing keys made two rates
        # reading the same `Clade` unequal, and a clade is the one driver written from literals
        # precisely so that it round-trips. Two equal clades describe one partition of one tree, and
        # nothing is drawn, so there is no sense in which they could be two different drivers. (That
        # is what separates this from `Drawn`, where writing one object or two is the model.)
        #
        # ``step`` is compared explicitly because it used to ride along inside `key`: two drivers
        # read at different resolutions are different models, and dropping it here would have made
        # them equal.
        return (isinstance(other, Driven) and other.driver == self.driver
                and other.mapping == self.mapping and other.step == self.step)

    def __hash__(self) -> int:
        # by the driver alone: a mapping is a dict or a callable and need not be hashable, so this is
        # coarser than __eq__ rather than inconsistent with it, which is all a hash owes. A driver
        # that is itself unhashable falls back to the class, keeping a Rate carrying one hashable.
        try:
            return hash((Driven, self.driver, self.step))
        except TypeError:
            return hash(Driven)


class SetBy(Driven):
    """**Replace** the parameter's base with a value read from a driver, rather than multiplying it.
    What ``set_by`` builds::

        loss = PerCopy().set_by(habitat, {"cave": 1.0, "surface": 0.25})   # the rate itself

    Written with **no base in front**, because there is none to write: the driver supplies the whole
    number, in the parameter's own units rather than as a dimensionless factor. That is what the
    literature usually means — "the loss rate is 1.0 in caves", not "four times a background nobody
    stated" — and spelling an absolute statement as a multiple of an invented background is the kind
    of quiet mismatch this grammar exists to avoid.

    **The scope still applies.** ``set_by`` replaces the base, not the *per what?*: a per-copy rate
    set to 1.0 is still 1.0 per copy, so it is multiplied by the copies present exactly as a written
    base would be. Only the number changes, which is why the scope is still written in front:
    ``PerCopy()``.

    It is a `Driven`, so every engine that resolves drivers resolves this one too — the trajectory,
    the mid-branch switches, the mapping checks are all the same machinery. What differs is one line
    in `Rate.effective`, which asks a ``SetBy`` for the base and every
    other modifier for a factor. The two compose: a replaced base may still be scaled.

        loss = PerCopy().set_by(habitat, {...}).scaled_by(size, Scalar(0.5))

    A rate may carry **one** ``SetBy``, written first. Two would be two answers to the same
    question, and neither order of application is more right than the other, so it raises rather
    than picking.
    """

    verb: str = SET_BY

    #: What a level's gate reads instead of an isinstance check (`Modifier.replaces_base`).
    #: Replacing a base is a capability an engine has or has not, and only three declare it — asking
    #: the object rather than its class is what lets this live beside the verb that writes it.
    replaces_base: ClassVar[bool] = True


_WRITTEN_AS.update({Driven: SCALED_BY, SetBy: SET_BY})
