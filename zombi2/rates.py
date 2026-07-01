"""Rate models — the coupling seam.

The single most important design decision in ZOMBI2: a rate model consumes the **whole
genome** and returns per-event propensities. This one signature absorbs every planned
rate variation as a subclass overriding :meth:`RateModel.propensities`:

* gene-wise (v1)      -> ``rate * genome.size()``
* genome-wise         -> ``rate`` (size-independent)
* length-dependent    -> ``rate * genome.total_length()``
* Potts / coupled     -> reads ``genome.presence_vector(order)``

Nothing else in the codebase changes when a new rate model is added.
"""

from __future__ import annotations

from dataclasses import dataclass

from .events import EventType, TargetParams


@dataclass(frozen=True)
class EventRates:
    """Base per-event rates. Interpretation depends on the rate model.

    In the v1 gene-wise model, D/T/L are *per gene copy* and O is *per branch*.
    """

    duplication: float
    transfer: float
    loss: float
    origination: float

    def validate(self) -> None:
        for name, value in vars(self).items():
            if value < 0:
                raise ValueError(f"{name} rate must be >= 0, got {value}")


class RateModel:
    """v1 default and base class: constant, gene-wise rates.

    Subclass and override :meth:`propensities` for genome-wise, length-dependent or
    coupled (Potts) rates — the ``(genome, branch, time)`` signature already gives you
    everything those variants need. (Alias :data:`ConstantGeneWiseRates`.)
    """

    def __init__(self, rates: EventRates, *, per_family: bool = False):
        rates.validate()
        self.rates = rates
        self.per_family = per_family  # reserved for future per-family rate draws

    def propensities(self, genome, branch: str, time: float) -> dict[EventType, float]:
        """Instantaneous rates for each event type on ``genome``."""
        n = genome.size()
        r = self.rates
        return {
            EventType.DUPLICATION: r.duplication * n,
            EventType.TRANSFER: r.transfer * n,
            EventType.LOSS: r.loss * n,
            EventType.ORIGINATION: r.origination,  # per-branch, size-independent
        }

    def target_params(self, event: EventType, genome, branch: str, time: float) -> TargetParams:
        """Parameters handed to :meth:`Genome.draw_target`. v1: the trivial default."""
        return TargetParams()


#: The spec name for the v1 default rate model.
ConstantGeneWiseRates = RateModel
