"""Rate models — the coupling seam.

A rate model consumes the **whole genome** and emits a list of weighted candidate
events, each tagged with the family it acts on (or ``None`` to act on a uniformly chosen
copy). This one interface handles every planned variation as a subclass:

* uniform (v1)      — every family shares D/T/L; emit ``family=None`` entries scaled by
  genome size (target copy chosen uniformly).
* per-family sampled — each family draws its own D/T/L at first sighting; emit one entry
  per family so the target family is chosen *weighted by its own rate*.
* genome-wise (future) — size-independent totals; emit ``family=None`` constant entries.
* Potts / coupled (future) — read ``genome.presence_vector(order)``.

Emitting per-family entries is what makes per-family rates *correct*: a plain uniform
copy draw would pick the target family without regard to its rate.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import namedtuple

from .distributions import as_distribution
from .events import EventType, TargetParams

#: A weighted candidate event. ``family`` is the family id to act on, or ``None`` to act
#: on a uniformly chosen gene copy (used for uniform/genome-wise rates and origination).
EventWeight = namedtuple("EventWeight", ["event", "family", "rate"])


class RateModel(ABC):
    """Abstract rate model: turns a genome into weighted candidate events."""

    @abstractmethod
    def event_weights(self, genome, branch: str, time: float) -> list[EventWeight]:
        ...

    def target_params(self, event: EventType, genome, branch: str, time: float) -> TargetParams:
        """Parameters handed to :meth:`Genome.draw_target`. v1: the trivial default."""
        return TargetParams()

    def bind_rng(self, rng) -> None:
        """Called once at the start of a simulation. Default: no-op.

        Stateful models (e.g. per-family sampled rates) use this to seed their sampling
        and reset any per-run cache.
        """


class UniformRates(RateModel):
    """Every gene family shares the same per-copy D/T/L rates (v1 default).

    D/T/L totals scale with genome size; the target copy is chosen uniformly, so the
    acting family is implicitly ~ copy number. Origination is a per-branch rate.
    """

    def __init__(self, duplication: float = 0.0, transfer: float = 0.0,
                 loss: float = 0.0, origination: float = 0.0):
        for name, value in (("duplication", duplication), ("transfer", transfer),
                            ("loss", loss), ("origination", origination)):
            if value < 0:
                raise ValueError(f"{name} rate must be >= 0, got {value}")
        self.duplication = float(duplication)
        self.transfer = float(transfer)
        self.loss = float(loss)
        self.origination = float(origination)

    def event_weights(self, genome, branch, time):
        n = genome.size()
        out: list[EventWeight] = []
        if n > 0:
            if self.duplication > 0:
                out.append(EventWeight(EventType.DUPLICATION, None, self.duplication * n))
            if self.transfer > 0:
                out.append(EventWeight(EventType.TRANSFER, None, self.transfer * n))
            if self.loss > 0:
                out.append(EventWeight(EventType.LOSS, None, self.loss * n))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        return out


class FamilySampledRates(RateModel):
    """Each gene family draws its OWN D/T/L rates from distributions (ZOMBI-1 style).

    A family's rates are sampled once, the first time it is seen, and kept for the life
    of the family (inherited by its descendants). ``duplication``/``transfer``/``loss``
    accept a built-in :class:`~zombi2.distributions.Distribution`, a plain float, a
    scipy.stats frozen distribution, or a callable ``rng -> float``. ``origination`` is
    a single per-branch rate (the birth rate of *new* families).
    """

    def __init__(self, duplication=0.0, transfer=0.0, loss=0.0, origination: float = 0.0):
        self._dup = as_distribution(duplication)
        self._trans = as_distribution(transfer)
        self._loss = as_distribution(loss)
        if origination < 0:
            raise ValueError(f"origination rate must be >= 0, got {origination}")
        self.origination = float(origination)
        self._rng = None
        self._family_rates: dict[str, tuple[float, float, float]] = {}

    def bind_rng(self, rng) -> None:
        self._rng = rng
        self._family_rates = {}

    def rates_for(self, family: str) -> tuple[float, float, float]:
        """The (dup, transfer, loss) rates for a family, sampled and cached on first use."""
        cached = self._family_rates.get(family)
        if cached is None:
            rng = self._rng
            cached = (
                max(0.0, self._dup.sample(rng)),
                max(0.0, self._trans.sample(rng)),
                max(0.0, self._loss.sample(rng)),
            )
            self._family_rates[family] = cached
        return cached

    def event_weights(self, genome, branch, time):
        out: list[EventWeight] = []
        for family in genome.families():
            cn = genome.copy_number(family)
            if cn == 0:
                continue
            d, t, l = self.rates_for(family)
            if d > 0:
                out.append(EventWeight(EventType.DUPLICATION, family, d * cn))
            if t > 0:
                out.append(EventWeight(EventType.TRANSFER, family, t * cn))
            if l > 0:
                out.append(EventWeight(EventType.LOSS, family, l * cn))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        return out
