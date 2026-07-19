"""Bridge: compile a grammar coupling onto the rate engine's emission seam.

The declarative grammar (:mod:`zombi2.coevolve.grammar`) says *what* a coupling is; this module says
*how* a **rate**-target coupling runs — as a :class:`CouplingModifier`, a
:class:`zombi2.genomes.rates.Modifier` whose :meth:`~CouplingModifier.factor` reads a
:class:`~zombi2.coevolve.grammar.DriverSignal`'s value along the tree and applies the coupling's
:class:`~zombi2.coevolve.grammar.Response`. Composed via :class:`~zombi2.genomes.rates.ModifiedRates`
over any base :class:`~zombi2.genomes.rates.RateModel`, it multiplies **only** its target event's
rate; every other event is left at ``1.0``, so a null coupling (``response = 0``) leaves the base
weights byte-identical.

This is the one piece of the grammar that touches the rate engine — the pure language stays in
``grammar.py``. It subclasses the *stable* ``Modifier`` / ``ModifiedRates`` API (untouched by the
``per=`` opportunity-knob rename, which only renames the base-rate classes). See
``docs/design/coevolve-grammar.md`` §4.3.
"""

from __future__ import annotations

from zombi2.coevolve.grammar import DriverSignal, Response
from zombi2.genomes.events import EventType
from zombi2.genomes.rates import Modifier

#: Grammar genomes target-variable → the :class:`EventType` it bends. Only the unambiguous ones —
#: ``"gain"`` is context-dependent (transfer vs origination), so a caller passes an explicit
#: ``EventType`` for it rather than relying on this map.
EVENT_FOR_VARIABLE: dict[str, EventType] = {
    "loss": EventType.LOSS,
    "duplication": EventType.DUPLICATION,
    "transfer": EventType.TRANSFER,
    "origination": EventType.ORIGINATION,
}


class CouplingModifier(Modifier):
    """A grammar coupling on a genome **rate** target, compiled onto the emission seam.

    Parameters
    ----------
    signal:
        The driver's realised history (a :class:`~zombi2.coevolve.grammar.DriverSignal` — e.g. a
        ``TraitTrajectory``); its ``value(branch, time)`` is the driver value the response reads.
    response:
        The grammar :class:`~zombi2.coevolve.grammar.Response` (``rate_multiplier``).
    event:
        The :class:`EventType` this coupling bends. The factor applies only to this event; every
        other event returns ``1.0`` (no effect).
    weights:
        Optional per-family scaling of the driver value, the way the trait→gene panel does: family
        ``f`` sees ``weights[f] · driver_value`` (so ``loss = base · exp(-strength · w_f · s)``). A
        family absent from ``weights`` (weight ``0``) is unaffected. When ``weights`` is given the
        modifier keys on family (:attr:`keys_on_family`), so :class:`ModifiedRates` expands aggregate
        ``family=None`` base weights into per-family ones for the factor to key on.
    """

    def __init__(self, signal: DriverSignal, response: Response, event: EventType, *, weights=None):
        self.signal = signal
        self.response = response
        self.event = event
        self.weights = dict(weights) if weights is not None else None
        # Per-family weights need per-family base weights → key on family (see ModifiedRates).
        self.keys_on_family = weights is not None

    def factor(self, event, family, branch, time) -> float:
        if event is not self.event:
            return 1.0
        driver_value = self.signal.value(branch, time)
        if self.weights is not None:
            w = self.weights.get(family, 0.0)
            if w == 0.0:
                return 1.0
            driver_value = w * driver_value
        return self.response.rate_multiplier(driver_value)

    def refresh_times(self, t0, t1):
        # The coupled rate changes exactly when the driver does — forward the signal's change points
        # so the simulator refreshes the branch there (no blunt time_dependent full refresh).
        return self.signal.refresh_times(t0, t1)
