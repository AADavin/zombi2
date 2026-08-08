"""The mapping of a `Driven` — what turns the driver's value into a number (SPEC §5).

The ``ScaledBy`` modifier reads the driver's value on a lineage; the **mapping** turns that value
into the number the modifier contributes — a dimensionless multiplier on a rate or an extent, a
normalised weight on ``transfer_to``. There are four shapes:

- `Table`  — a **discrete** driver → a dict of factors: ``{"aquatic": 3.0, "terrestrial": 1.0}``.
- `Curve`  — a **continuous** driver → a function: ``lambda x: math.exp(0.5 * x)``.
- `Scalar` — a single log-link coefficient: ``multiplier = exp(strength · value)``.
- `Between` — a weight per **(donor group, recipient group)** pair, which only ``transfer_to``
  takes: it reads the driver at both ends, so a rate or an extent refuses it
  (`check_not_a_kernel()`).

You rarely name the first three — pass a raw ``dict`` / callable / number as ``mapping=`` and
`as_mapping()` coerces it (a dict → ``Table``, a callable → ``Curve``, a number →
``Scalar``), exactly as `as_rate()` coerces a rate spec.

**Jump** (a burst fired *at an event*, e.g. a pulse of gene change at each split) is not a mapping:
it changes a state at a moment rather than scaling a number, so it does not live here and is not
reachable through ``ScaledBy`` (SPEC §4).
"""

from __future__ import annotations

import math
from typing import cast

_MAX_EXPONENT = 40.0  # clamp the log-link argument so a large driver value cannot overflow exp()


class Mapping:
    """Base for a driver-value → factor mapping. Abstract — use `Table`,
    `Curve`, or `Scalar` (or pass a raw dict / callable / number, which
    `as_mapping()` coerces). A mapping returns a **dimensionless, non-negative** factor — a
    multiplier on a rate or an extent, a weight on ``transfer_to``."""

    def multiplier(self, value: object) -> float:
        raise NotImplementedError


class Table(Mapping):
    """A **discrete** driver → a lookup of factors, one per driver state::

        Table({"aquatic": 3.0, "terrestrial": 1.0})   # 3× the rate in aquatic lineages

    ``default`` (1.0) is the factor for any state not named — so an unlisted state leaves the
    rate unchanged. This is the primary ``ScaledBy`` mapping (MuSSE-style per-state rates).

    States are matched by their **string form** — ``Table({0: 3.0, 1: 1.0})`` and ``Table({"0":
    3.0, "1": 1.0})`` behave identically, and both match a driver whose value is ``0`` or ``"0"``.
    A conditioned driver arrives from a text file (always a string), and a live joint driver arrives
    as its native label; string-matching makes the two agree, so an ``int``-labelled trait does not
    silently miss its mapping."""

    def __init__(self, per_state, default: float = 1.0) -> None:
        if not isinstance(per_state, dict) or not per_state:
            raise ValueError(f"Table needs a non-empty {{state: factor}} dict, got {per_state!r}")
        table = {}
        for state, factor in per_state.items():
            key = str(state)  # states matched by string form (a driver file is text); see the class docstring
            if key in table:
                raise ValueError(
                    f"Table states collide as strings: {state!r} and an earlier key both map to {key!r}"
                )
            table[key] = _check_factor(factor, f"Table factor for {state!r}")
        self.per_state = table
        self.default = _check_factor(default, "Table default")

    def multiplier(self, value: object) -> float:
        return self.per_state.get(str(value), self.default)

    def __repr__(self) -> str:
        # `repr(float)`, not `:g` — see `OnTime.__repr__`: six significant figures in a run's log
        # is a record of a different model.
        inner = ", ".join(f"{s!r}: {float(f)!r}" for s, f in self.per_state.items())
        tail = "" if self.default == 1.0 else f", default={float(self.default)!r}"
        return f"Table({{{inner}}}{tail})"

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, Table) and other.per_state == self.per_state
                and other.default == self.default)


class Curve(Mapping):
    """A **continuous** driver → an arbitrary function of the value, optionally capped::

        Curve(lambda x: math.exp(0.5 * x))          # exponential response
        Curve(lambda x: 1.0 + x, bound=5.0)          # linear, capped at 5×

    ``bound`` (a ceiling on the factor) is what an exact Gillespie thinner needs when the
    driver is unbounded; omit it for a naturally-bounded ``fn``. The function must return a
    finite, non-negative number for every driver value it sees (a rate cannot go negative)."""

    def __init__(self, fn, bound: float | None = None) -> None:
        if not callable(fn):
            raise TypeError(f"Curve needs a callable value→factor function, got {fn!r}")
        if bound is not None:
            if isinstance(bound, bool) or not isinstance(bound, (int, float)) \
                    or not math.isfinite(bound) or bound < 0:
                raise ValueError(f"Curve bound must be a finite non-negative number, got {bound!r}")
            bound = float(bound)
        self.fn = fn
        self.bound = bound

    def multiplier(self, value: object) -> float:
        f = self.fn(_numeric(value, "Curve"))
        if isinstance(f, bool) or not isinstance(f, (int, float)) or not math.isfinite(f) or f < 0:
            raise ValueError(
                f"the Curve returned {f!r} for driver value {value!r}; a mapping's factor must be a "
                f"finite non-negative number"
            )
        f = float(f)
        return min(f, self.bound) if self.bound is not None else f

    def __repr__(self) -> str:
        tail = "" if self.bound is None else f", bound={self.bound:g}"
        return f"Curve({getattr(self.fn, '__name__', 'fn')}{tail})"


class Scalar(Mapping):
    """A single log-link coefficient — ``multiplier = exp(strength · value)``::

        Scalar(0.0)    # null: factor 1 for every value
        Scalar(0.7)    # a binary 0/1 driver gives factor 1 (off) or exp(0.7) ≈ 2.0 (on)

    The natural response when the driver is already a 0/1 indicator or a single continuous
    covariate: one knob, ``strength`` (0 ⇒ the driver does not change the rate). The exponent is
    clamped so a large value cannot overflow."""

    def __init__(self, strength: float) -> None:
        if isinstance(strength, bool) or not isinstance(strength, (int, float)) \
                or not math.isfinite(strength):
            raise ValueError(f"Scalar strength must be a finite number, got {strength!r}")
        self.strength = float(strength)

    def multiplier(self, value: object) -> float:
        x = self.strength * _numeric(value, "Scalar")
        x = max(-_MAX_EXPONENT, min(_MAX_EXPONENT, x))  # guard exp() against overflow
        return math.exp(x)

    def __repr__(self) -> str:
        return f"Scalar(strength={float(self.strength)!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Scalar) and other.strength == self.strength


class Between:
    """A weight over ordered **(donor-group, recipient-group)** pairs — the 2-D kernel of the transfer
    **choice** of who receives (SPEC §5), the donor-conditioned sibling of `Table`::

        Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)   # A↔B only, nothing else receives
        Between({("A", "B"): 3.0})                                 # A→B 3× baseline, every other pair 1×

    A `Table` weights a candidate recipient by *that candidate's* state alone; a ``Between``
    weights it by the **pair** — the donor's group and the recipient's — which is what lets a transfer
    be steered to run *between* two groups rather than within them. It is therefore **not** a
    `Mapping` (a ``Mapping.multiplier`` reads one value): its `weight()` reads two, and the
    engine passes both. It is used in ``transfer_to`` — on its own as the kernel of a
    `Clades` rule (groups from the tree), or inside a
    `ScaledBy` (groups from a trait). It is **not** a rate multiplier:
    a rate has no donor to condition on, so a ``Between`` on a rate is refused.

    Keys are ``(from_group, to_group)`` pairs matched by **string form**, exactly like ``Table``'s
    states, so an integer-labelled group still finds its entry. ``default`` (1.0) is the weight for any
    pair not named — ``default=0.0`` gives the "only the flows I name can happen" idiom, reusing the
    rule that a **weight of 0 means the donor cannot send to that recipient group**;
    when every candidate weighs 0 the transfer has nowhere to land and does not fire."""

    def __init__(self, per_pair, default: float = 1.0) -> None:
        if not isinstance(per_pair, dict) or not per_pair:
            raise ValueError(
                f"Between needs a non-empty {{(from_group, to_group): weight}} dict, got {per_pair!r}")
        table: dict[tuple[str, str], float] = {}
        for pair, weight in per_pair.items():
            if not (isinstance(pair, tuple) and len(pair) == 2):
                raise ValueError(
                    f"Between keys are (from_group, to_group) pairs, got {pair!r} — write "
                    f"Between({{('A', 'B'): 1.0}}), the donor group first, the recipient group second")
            key = (str(pair[0]), str(pair[1]))  # groups matched by string form, like Table's states
            if key in table:
                raise ValueError(
                    f"Between pairs collide as strings: {pair!r} and an earlier key both map to {key!r}")
            table[key] = _check_factor(weight, f"Between weight for {pair!r}")
        self.per_pair = table
        self.default = _check_factor(default, "Between default")

    def weight(self, from_group: object, to_group: object) -> float:
        """The weight for a transfer from a ``from_group`` donor to a ``to_group`` recipient — the
        named pair's weight, or `default` if the pair is unnamed."""
        return self.per_pair.get((str(from_group), str(to_group)), self.default)

    def groups(self) -> set:
        """Every group named on either side of a pair — what a fires-check tests against the groups
        that actually occur, so a kernel naming only absent groups (a typo) can be caught."""
        return {g for pair in self.per_pair for g in pair}

    def __repr__(self) -> str:
        inner = ", ".join(f"({a!r}, {b!r}): {float(w)!r}" for (a, b), w in self.per_pair.items())
        tail = "" if self.default == 1.0 else f", default={float(self.default)!r}"
        return f"Between({{{inner}}}{tail})"

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, Between) and other.per_pair == self.per_pair
                and other.default == self.default)


def check_kernel_fires(kernel: Between, available_groups, *, driver_label: str) -> None:
    """Raise if a `Between` names **no pair whose two groups both occur** among
    ``available_groups`` — the recipient-weight twin of
    `check_mapping_fires()`. Such a kernel weights every candidate at its
    ``default``, so the recipient choice is secretly *uniform* while the run records it as steered —
    almost always a typo in a group name or a stale driver. A kernel may still name a pair this
    realisation never realises (a legitimate partial kernel), so only an *empty* overlap is refused."""
    have = {str(g) for g in available_groups}
    if not any(a in have and b in have for a, b in kernel.per_pair):
        raise ValueError(
            f"Between on {driver_label}: the kernel's groups {sorted(kernel.groups())} include no pair "
            f"whose two groups both occur in {sorted(have)}, so the weighting would silently do nothing "
            f"— every candidate falls to the default weight and the recipient is drawn uniformly. Check "
            f"for a typo in the group names, or a stale or mismatched driver.")


def _check_factor(x: object, where: str) -> float:
    """Coerce ``x`` to a finite, non-negative float (a rate multiplier) or raise naming ``where``."""
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0:
        raise ValueError(f"{where} must be a finite non-negative number, got {x!r}")
    return float(x)


def _numeric(value: object, cls: str) -> float:
    """A driver value as a float, for a **continuous** mapping (``Curve`` / ``Scalar``). Raises a
    clear error if the driver is a discrete label — the usual sign a discrete driver was given a
    continuous mapping (use a ``Table`` / dict for a discrete driver)."""
    try:
        return float(value)  # type: ignore[arg-type]  # the except below is the check
    except (TypeError, ValueError):
        raise ValueError(
            f"{cls} is a continuous-driver mapping but got the discrete driver value {value!r}; use a "
            f"Table (a dict mapping) for a discrete driver such as a habitat state."
        ) from None


def as_mapping(spec: object) -> Mapping:
    """Coerce a ``Driven`` mapping spec into a `Mapping`.

    Accepts an already-built mapping (returned unchanged), a ``dict`` (→ `Table`), a
    callable (→ `Curve`), or a number (→ `Scalar`). Mirrors
    `as_rate()` / `as_distribution()`.
    """
    if isinstance(spec, (Mapping, Between)):
        # a Between is a choice's kernel, not a rate multiplier; carried through here so
        # ScaledBy(..., Between(...)) works, and refused on a rate or an extent by the engine — which is
        # why the declared return type is the one every *rate* caller may rely on.
        return cast(Mapping, spec)
    if isinstance(spec, dict):
        return Table(spec)
    if isinstance(spec, bool):
        raise TypeError(f"a Driven mapping cannot be a bool, got {spec!r}")
    if isinstance(spec, (int, float)):
        return Scalar(float(spec))
    if callable(spec):
        return Curve(spec)
    raise TypeError(
        f"a Driven mapping must be a dict (Table), a callable (Curve), a number (Scalar), a "
        f"Table/Curve/Scalar, or a Between (a transfer_to kernel), got {spec!r}"
    )


def check_not_a_kernel(mapping, *, label: str) -> None:
    """Raise if a **rate** (or an extent) is driven through a `Between` kernel.

    A ``Between`` weights a recipient by the ``(donor, recipient)`` group pair, so it answers *who
    receives* and belongs in ``transfer_to``. A rate is read on one lineage and has no donor to
    condition on, so a kernel there has nothing to be a pair with: `Between` deliberately implements
    no ``multiplier``, and an engine that does not check first dies part-way through a run with
    ``AttributeError: 'Between' object has no attribute 'multiplier'`` — a traceback from inside the
    engine, naming neither the rate nor the mistake.

    Every engine that accepts ``ScaledBy`` on a rate or an extent calls this, so the message is the
    same one wherever the kernel was put."""
    if isinstance(mapping, Between):
        raise ValueError(
            f"{label} carries ScaledBy(…, Between(…)); a Between kernel is donor-conditioned — it "
            f"weights a recipient by the (donor, recipient) group pair — so it belongs in transfer_to "
            f"(who RECEIVES) and never in a rate or an extent, which are read on one lineage and have "
            f"no donor to condition on. Drive this with a Table (a plain dict) or a Curve, and put the "
            f"kernel in transfer_to=Weights(driver, Between({{...}})).")


__all__ = ["Mapping", "Table", "Curve", "Scalar", "Between", "check_kernel_fires",
           "check_not_a_kernel", "as_mapping"]
