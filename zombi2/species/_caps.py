"""Per-class capability descriptors for species-tree models.

Replaces the scattered ``isinstance`` / ``getattr`` dispatch in ``sim.py`` / ``forward.py``:
each model class declares a :class:`SpeciesCaps` as a class attribute, and the simulators read
it via :func:`species_caps`. An unregistered model — or an unknown growth engine — is then a
loud error at dispatch time, not a silently wrong tree (the failure mode the old
``isinstance(model, ...)`` ladders allowed). See ``docs/design/model-architecture.md``.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass
from typing import Any, Callable


class GrowthEngine(enum.Enum):
    """Which forward growth loop simulates a model."""

    THINNING = "thinning"      # time-varying rates -> the thinning loop ``_grow``
    GILLESPIE = "gillespie"    # rates constant between events -> ``_grow_gillespie``


@dataclass(frozen=True)
class Feature:
    """A forward-only model feature and the test for whether it is *active* on an instance.

    ``is_active`` knows the feature's neutral value **and** its shape (e.g. fossilization is a
    scalar on :class:`~zombi2.BirthDeath` but a per-epoch list on
    :class:`~zombi2.EpisodicBirthDeath`), so adding a feature is one declaration rather than a
    new scalar/list branch in the backward-mode guard.
    """

    name: str
    is_active: Callable[[Any], bool]


def _fossilization_active(v) -> bool:
    if v is None:
        return False
    return (sum(v) if isinstance(v, list) else v) > 0


FOSSILIZATION = Feature("fossilization", _fossilization_active)
REMOVAL = Feature("removal", lambda v: v is not None and v != 1.0)
MASS_EXTINCTIONS = Feature("mass_extinctions", bool)


@dataclass(frozen=True)
class SpeciesCaps:
    """What a species-tree model can do — declared once, on the class."""

    growth: GrowthEngine
    supports_backward: bool = False              # has a closed-form reconstructed CDF (sim.py)
    supports_ghosts: bool = False                # add_ghost_lineages implemented (ghosts.py)
    supports_n_tips: bool = True                 # forward ``n_tips`` stop mode allowed (forward.py)
    incomplete_sampling_backward: bool = False   # ρ<1 allowed in backward mode (sim.py)
    forward_only_features: tuple[Feature, ...] = ()   # rejected in backward mode


def species_caps(model) -> SpeciesCaps:
    """The :class:`SpeciesCaps` in force for ``model``.

    An **instance** may override its class caps by setting ``self._caps`` — used when a model's
    engine depends on its configuration (e.g. ``BirthDeath(per="shared")`` is forward-only Gillespie
    where the per-lineage default is backward-capable). Otherwise walks the MRO so a subclass (e.g.
    ``Yule(BirthDeath)``) inherits its parent's caps. Raises :class:`TypeError` for a model with no
    declared caps — the loud error that replaces a forgotten ``isinstance`` branch silently routing a
    model into the wrong engine.
    """
    inst = model.__dict__.get("_caps")
    if isinstance(inst, SpeciesCaps):
        return inst
    for klass in type(model).__mro__:
        caps = klass.__dict__.get("_caps")
        if isinstance(caps, SpeciesCaps):
            return caps
    raise TypeError(
        f"{type(model).__name__} has no SpeciesCaps; declare one as a class attribute "
        "(see zombi2/species/_caps.py)"
    )


def active_forward_features(model, caps: SpeciesCaps) -> list[str]:
    """Names of ``caps.forward_only_features`` actually switched on for ``model`` — the uniform
    replacement for the scalar/list/``getattr`` probing in the backward-mode guard."""
    return [f.name for f in caps.forward_only_features
            if f.is_active(getattr(model, f.name, None))]
