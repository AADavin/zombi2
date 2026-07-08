"""Rate models — the coupling seam, plus duplication growth control.

A rate model consumes the **whole genome** and emits a list of weighted candidate
events, each tagged with the family it acts on (or ``None`` to act on a uniformly chosen
copy). This one interface handles every planned variation as a subclass:

* uniform (v1)      — every family shares D/T/L; emit ``family=None`` entries scaled by
  genome size (target copy chosen uniformly).
* per-family sampled — each family draws its own D/T/L at first sighting; emit one entry
  per family so the target family is chosen *weighted by its own rate*.
* genome-wise (future) — size-independent totals; emit ``family=None`` constant entries.

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

import math
from abc import ABC, abstractmethod
from collections import namedtuple

from zombi2.distributions import as_distribution
from zombi2.genomes.events import EventType, TargetParams

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


def _check_carrying_capacity(carrying_capacity: float | None) -> float | None:
    """Validate a soft carrying capacity ``K``: ``None`` (disabled) or a finite positive
    value. A non-positive ``K`` is nonsensical — ``K=0`` divides by zero mid-simulation and
    a negative ``K`` inverts the density term into runaway amplification — so reject it up
    front (mirroring the species-side ``DiversityDependent`` guard)."""
    if carrying_capacity is None:
        return None
    k = float(carrying_capacity)
    if not math.isfinite(k) or k <= 0:
        raise ValueError(
            f"carrying_capacity must be a finite positive number, got {carrying_capacity}"
        )
    return k


class RateModel(ABC):
    """Abstract rate model: turns a genome into weighted candidate events."""

    _rng = None
    _max_family_size: int | None = None

    #: The simulator caches each branch's weights and refreshes only branches whose genome
    #: changed. Set this True in a subclass whose weights vary *continuously* with ``time``
    #: within a branch interval, to force a full refresh every event (correctness over speed).
    time_dependent = False

    @abstractmethod
    def event_weights(self, genome, branch: str, time: float) -> list[EventWeight]:
        ...

    def target_params(self, event: EventType, genome, branch: str, time: float) -> TargetParams:
        return TargetParams()

    def refresh_times(self, t0: float, t1: float) -> list[tuple[float, str]]:
        """Times in ``(t0, t1)`` at which a branch's weights change on their own — i.e. *not*
        as a consequence of a gene event — and which branch changes. The simulator pauses the
        Gillespie loop at each, refreshes that branch's cached weights, and continues, so a
        rate model whose weights follow an externally-imposed schedule (e.g. a trait value that
        drifts/jumps along the branch, see :mod:`zombi2.coevolve.trait_coupling`) stays exact without
        the blunt ``time_dependent`` full-refresh-every-event. Default: none — weights change
        only at events, so no extra refresh points are needed."""
        return []

    def establishment_probability(self, selection, recipient_genome, time: float) -> float:
        """Probability that a horizontally transferred ``selection`` *establishes* in
        ``recipient_genome`` — rolled once per transfer, *after* a real donor is found, so gain
        stays donor-limited. Default ``1.0``: every transfer establishes, the caller skips the
        RNG draw, and this default never inspects ``selection`` (whose ``genes`` are only
        populated at extraction for some genome types) — so existing models keep byte-identical
        streams. A subclass may override this to bias establishment by the recipient genome's
        current content, putting part of any coupling on the *gain* channel."""
        return 1.0

    def bind(self, rng, max_family_size: int | None = None, tree=None) -> None:
        """Called once at the start of a simulation with the RNG, the resolved hard
        family-size cap (or ``None``), and the species ``tree`` (some models, e.g.
        autocorrelated branch rates, need the topology). Subclasses may extend; call
        ``super().bind``."""
        self._rng = rng
        self._max_family_size = max_family_size


class SharedRates(RateModel):
    """Every gene family shares the same per-copy D/T/L rates (v1 default)."""

    def __init__(self, duplication: float = 0.0, transfer: float = 0.0,
                 loss: float = 0.0, origination: float = 0.0,
                 *, inversion: float = 0.0, transposition: float = 0.0,
                 insertion: float = 0.0, deletion: float = 0.0,
                 carrying_capacity: float | None = None):
        rates = (("duplication", duplication), ("transfer", transfer), ("loss", loss),
                 ("origination", origination), ("inversion", inversion),
                 ("transposition", transposition), ("insertion", insertion),
                 ("deletion", deletion))
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
        # intergenic indels: only fired by the nucleotide genome (per-nucleotide rate)
        self.insertion = float(insertion)
        self.deletion = float(deletion)
        self.carrying_capacity = _check_carrying_capacity(carrying_capacity)

    def _regulated(self) -> bool:
        # Only the *soft* carrying capacity needs per-family duplication weights. A hard
        # max_family_size is enforced cheaply by the simulator (it skips an over-cap
        # duplication), so we keep the O(1) aggregate duplication term.
        return self.carrying_capacity is not None

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
            if self.insertion > 0:
                out.append(EventWeight(EventType.INSERTION, None, self.insertion * n))
            if self.deletion > 0:
                out.append(EventWeight(EventType.DELETION, None, self.deletion * n))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        return out


class PerGenomeRates(RateModel):
    """Genome-wise rates: each event type fires at a **constant per-genome rate**,
    independent of how many gene copies the genome holds.

    Contrast :class:`SharedRates`, where the total duplication/transfer/loss rate scales
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
        self.carrying_capacity = _check_carrying_capacity(carrying_capacity)
        self._family_rates: dict[str, tuple[float, float, float]] = {}

    def bind(self, rng, max_family_size: int | None = None, tree=None) -> None:
        super().bind(rng, max_family_size, tree)
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


class BranchRates(RateModel):
    """Make rates vary per species-tree branch by scaling a base rate model.

    Wraps any base rate model (``SharedRates``, ``FamilySampledRates``, ...) and
    multiplies its duplication/transfer/loss weights on each branch by a per-branch
    factor (a single scalar scaling D/T/L together; origination is left unscaled). This
    composes with the base model, so branch heterogeneity and family/uniform rates
    combine into a two-factor model.

    Provide exactly one factor source:

    * ``autocorr_sigma`` — a **relaxed clock**: the factor evolves down the tree as
      ``factor(child) = factor(parent) * exp(N(0, sigma * sqrt(branch_length)))``, so
      closely related lineages have similar rates. Needs the tree (passed via ``bind``).
      ``sigma = 0`` reduces to the base model.
    * ``per_branch`` — a distribution (built-in / scipy / callable) drawn **i.i.d. per
      branch**, independently.
    * ``factors`` — an explicit ``{branch_name: factor}`` map (branches not listed use
      ``root_rate``).

    ``root_rate`` is the factor of the root branch (default 1.0).
    """

    def __init__(self, base: RateModel, *, autocorr_sigma: float | None = None,
                 per_branch=None, factors: dict | None = None, root_rate: float = 1.0):
        sources = [autocorr_sigma is not None, per_branch is not None, factors is not None]
        if sum(sources) != 1:
            raise ValueError("specify exactly one of autocorr_sigma, per_branch, or factors")
        if autocorr_sigma is not None and autocorr_sigma < 0:
            raise ValueError("autocorr_sigma must be >= 0")
        self.base = base
        self.autocorr_sigma = autocorr_sigma
        self.per_branch = as_distribution(per_branch) if per_branch is not None else None
        self.explicit = dict(factors) if factors is not None else None
        self.root_rate = float(root_rate)
        self._factor: dict[str, float] = {}

    _SCALED = (EventType.DUPLICATION, EventType.TRANSFER, EventType.LOSS)

    def bind(self, rng, max_family_size: int | None = None, tree=None) -> None:
        super().bind(rng, max_family_size, tree)
        self.base.bind(rng, max_family_size, tree)
        self._factor = {}
        if self.explicit is not None:
            self._factor = dict(self.explicit)
        elif self.autocorr_sigma is not None:
            if tree is None:
                raise ValueError("autocorrelated branch rates need the species tree")
            for node in tree.nodes_preorder():
                if node.parent is None:
                    self._factor[node.name] = self.root_rate
                else:
                    scale = self.autocorr_sigma * math.sqrt(max(node.branch_length(), 0.0))
                    drift = math.exp(rng.normal(0.0, scale)) if scale > 0 else 1.0
                    self._factor[node.name] = self._factor[node.parent.name] * drift
        # per-branch i.i.d.: sampled lazily in _branch_factor

    def _branch_factor(self, branch: str) -> float:
        f = self._factor.get(branch)
        if f is None:
            f = self.per_branch.sample(self._rng) if self.per_branch is not None else self.root_rate
            self._factor[branch] = f
        return f

    def event_weights(self, genome, branch, time):
        factor = self._branch_factor(branch)
        out = []
        for ew in self.base.event_weights(genome, branch, time):
            if ew.event in self._SCALED and factor != 1.0:
                out.append(EventWeight(ew.event, ew.family, ew.rate * factor))
            else:
                out.append(ew)
        return out

    def target_params(self, event, genome, branch, time):
        return self.base.target_params(event, genome, branch, time)
