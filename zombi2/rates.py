"""Rate models — the coupling seam, plus duplication growth control.

A rate model consumes the **whole genome** and emits a list of weighted candidate
events, each tagged with the family it acts on (or ``None`` to act on a uniformly chosen
copy). This one interface handles every planned variation as a subclass:

* uniform (v1)      — every family shares D/T/L; emit ``family=None`` entries scaled by
  genome size (target copy chosen uniformly).
* per-family sampled — each family draws its own D/T/L at first sighting; emit one entry
  per family so the target family is chosen *weighted by its own rate*.
* genome-wise (future) — size-independent totals; emit ``family=None`` constant entries.
* Potts / coupled (future) — read ``genome.presence_vector(order)``.

**Growth control.** A family's copy number is a birth-death process; with duplication >
loss it grows like ``e^{(d-l)t}``. Two mechanisms:

* ``carrying_capacity`` (K, per rate model, optional) — logistic density dependence: the
  per-copy duplication rate is scaled by ``max(0, 1 - n/K)`` (a soft, mechanistic model).
* ``max_family_size`` (a *simulation* parameter, delivered here via :meth:`bind`) — a hard
  ceiling on family size. It suppresses duplication at the ceiling; transfers are capped
  separately by the simulator (an over-cap transfer becomes a replacement). The hard cap
  is the one that can be set as a fraction of the number of species (see the driver).

Any growth control makes duplication per-family (its rate depends on that family's size),
so the duplication term flips to one entry per family; loss/transfer/origination are
unaffected.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import namedtuple

from .distributions import as_distribution
from .events import EventType, TargetParams

#: A weighted candidate event. ``family`` is the family id to act on, or ``None`` to act
#: on a uniformly chosen gene copy (used for uniform/genome-wise rates and origination).
EventWeight = namedtuple("EventWeight", ["event", "family", "rate"])


def duplication_factor(n: int, carrying_capacity: float | None, max_copies: int | None) -> float:
    """Multiplier on the per-copy duplication rate for a family with ``n`` copies."""
    if max_copies is not None and n >= max_copies:
        return 0.0
    if carrying_capacity is not None:
        return max(0.0, 1.0 - n / carrying_capacity)
    return 1.0


class RateModel(ABC):
    """Abstract rate model: turns a genome into weighted candidate events."""

    _rng = None
    _max_family_size: int | None = None

    @abstractmethod
    def event_weights(self, genome, branch: str, time: float) -> list[EventWeight]:
        ...

    def target_params(self, event: EventType, genome, branch: str, time: float) -> TargetParams:
        return TargetParams()

    def bind(self, rng, max_family_size: int | None = None) -> None:
        """Called once at the start of a simulation with the RNG and the resolved
        hard family-size cap (or ``None``). Subclasses may extend; call ``super().bind``."""
        self._rng = rng
        self._max_family_size = max_family_size


class UniformRates(RateModel):
    """Every gene family shares the same per-copy D/T/L rates (v1 default)."""

    def __init__(self, duplication: float = 0.0, transfer: float = 0.0,
                 loss: float = 0.0, origination: float = 0.0,
                 *, inversion: float = 0.0, transposition: float = 0.0,
                 carrying_capacity: float | None = None):
        rates = (("duplication", duplication), ("transfer", transfer), ("loss", loss),
                 ("origination", origination), ("inversion", inversion),
                 ("transposition", transposition))
        for name, value in rates:
            if value < 0:
                raise ValueError(f"{name} rate must be >= 0, got {value}")
        self.duplication = float(duplication)
        self.transfer = float(transfer)
        self.loss = float(loss)
        self.origination = float(origination)
        # rearrangements: only fired by genomes that support them (e.g. OrderedGenome)
        self.inversion = float(inversion)
        self.transposition = float(transposition)
        self.carrying_capacity = carrying_capacity

    def _regulated(self) -> bool:
        return self.carrying_capacity is not None or self._max_family_size is not None

    def event_weights(self, genome, branch, time):
        n = genome.size()
        out: list[EventWeight] = []

        if self.duplication > 0 and n > 0:
            if self._regulated():  # per-family duplication (rate depends on family size)
                for family in genome.families():
                    cn = genome.copy_number(family)
                    f = duplication_factor(cn, self.carrying_capacity, self._max_family_size)
                    if f > 0:
                        out.append(EventWeight(EventType.DUPLICATION, family, self.duplication * cn * f))
            else:  # aggregate fast path
                out.append(EventWeight(EventType.DUPLICATION, None, self.duplication * n))

        if n > 0:
            if self.transfer > 0:
                out.append(EventWeight(EventType.TRANSFER, None, self.transfer * n))
            if self.loss > 0:
                out.append(EventWeight(EventType.LOSS, None, self.loss * n))
            if self.inversion > 0:
                out.append(EventWeight(EventType.INVERSION, None, self.inversion * n))
            if self.transposition > 0:
                out.append(EventWeight(EventType.TRANSPOSITION, None, self.transposition * n))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        return out


class GenomeWiseRates(RateModel):
    """Genome-wise rates: each event type fires at a **constant per-genome rate**,
    independent of how many gene copies the genome holds.

    Contrast :class:`UniformRates`, where the total duplication/transfer/loss rate scales
    with genome size (per-copy rates). Here the totals are fixed, so when an event fires a
    target copy is chosen uniformly. A useful consequence: family sizes grow *linearly*
    rather than exponentially, so genome-wise models are intrinsically far less prone to
    runaway growth. Origination is per branch.
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
        out: list[EventWeight] = []
        if genome.size() > 0:
            if self.duplication > 0:
                out.append(EventWeight(EventType.DUPLICATION, None, self.duplication))
            if self.transfer > 0:
                out.append(EventWeight(EventType.TRANSFER, None, self.transfer))
            if self.loss > 0:
                out.append(EventWeight(EventType.LOSS, None, self.loss))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        return out


class FamilySampledRates(RateModel):
    """Each gene family draws its OWN D/T/L rates from distributions (ZOMBI-1 style)."""

    def __init__(self, duplication=0.0, transfer=0.0, loss=0.0, origination: float = 0.0,
                 *, carrying_capacity: float | None = None):
        self._dup = as_distribution(duplication)
        self._trans = as_distribution(transfer)
        self._loss = as_distribution(loss)
        if origination < 0:
            raise ValueError(f"origination rate must be >= 0, got {origination}")
        self.origination = float(origination)
        self.carrying_capacity = carrying_capacity
        self._family_rates: dict[str, tuple[float, float, float]] = {}

    def bind(self, rng, max_family_size: int | None = None) -> None:
        super().bind(rng, max_family_size)
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
                f = duplication_factor(cn, self.carrying_capacity, self._max_family_size)
                if f > 0:
                    out.append(EventWeight(EventType.DUPLICATION, family, d * cn * f))
            if t > 0:
                out.append(EventWeight(EventType.TRANSFER, family, t * cn))
            if l > 0:
                out.append(EventWeight(EventType.LOSS, family, l * cn))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        return out
