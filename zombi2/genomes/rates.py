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
import warnings
from abc import ABC, abstractmethod
from collections import namedtuple

from zombi2.distributions import as_distribution
from zombi2.genomes.events import EventType, TargetParams

#: A weighted candidate event. ``family`` is the family id to act on, or ``None`` to act
#: on a uniformly chosen gene copy (used for uniform/genome-wise rates and origination).
EventWeight = namedtuple("EventWeight", ["event", "family", "rate"])

#: A per-event **opportunity override**. Attach a specific opportunity ``unit`` — ``"copy"`` /
#: ``"lineage"`` / ``"shared"`` — to a single event's ``rate``, overriding the model-level ``per``.
#: e.g. ``Rates(duplication=Per("shared", 0.5), loss=0.3)`` is a shared duplication clock with per-copy
#: loss — a *self-limiting* family (births capped tree-wide, deaths growing with copy number). A bare
#: number uses the model-level ``per``. See ``docs/design/opportunity-knob.md``.
Per = namedtuple("Per", ["unit", "rate"])


def _opportunity(value, default_unit):
    """Split a rate spec into ``(rate, unit)``: a bare number uses ``default_unit``; a :class:`Per`
    overrides it with its own unit."""
    if isinstance(value, Per):
        return float(value.rate), str(value.unit)
    return float(value), default_unit


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

    def shared_event_weights(self, family) -> list[tuple]:
        """For a ``per="shared"`` model: the **tree-wide** base rate of each per-family event —
        ``[(EventType, rate), ...]`` — for ``family``. The simulator runs one such constant-rate
        clock per family present anywhere on the tree (independent of copy number), localising a
        fire to a copy chosen uniformly across the whole family. Default: none (not a shared model),
        so the simulator's shared pool stays empty and every other model is untouched."""
        return []

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


class Rates(RateModel):
    """Gene-family event rates, with a selectable **opportunity** (``per``) — the unit the clock rides.

    - ``per="copy"`` (default) — one clock **per gene copy**: the total D/T/L rate scales with the
      genome's copy number, so families grow **exponentially**. This is the built-in model the Rust
      engine implements.
    - ``per="lineage"`` — one clock **per genome**: the total D/T/L rate is a constant regardless of
      how many copies the genome holds (the target copy is chosen uniformly), so families grow
      **linearly**. Only the D/T/L/origination rates apply; rearrangement/chromosome events are a
      per-copy notion and are rejected here.

    ``per`` is the opportunity axis; it is orthogonal to *per-family* heterogeneity (which lives in the
    base :class:`FamilySampledRates` or a :class:`FamilyModifier`) — see ``docs/guide/rates.md`` and
    ``docs/design/opportunity-knob.md``. (``per="shared"`` — one clock for a whole family across the
    tree — is planned; not yet implemented for genomes.)

    An optional per-copy ``conversion`` rate adds **intra-genome gene conversion**: one copy of a
    family overwrites ("converts") another copy of the *same* family in the same genome —
    non-reciprocal and copy-number-neutral, the intra-genome analogue of transfer. It fires only on
    families holding two or more copies (a family with ``n`` copies converts at total rate
    ``conversion · n``), pulling within-family coalescences toward the present (concerted evolution).
    Pair it with a :class:`~zombi2.genomes.conversion.ConversionModel` (via
    ``simulate_genomes(..., conversions=...)``) to bias the donor; without one, conversion is
    unbiased. Conversion shapes gene *trees* rather than copy-number profiles, so a model with
    ``conversion > 0`` runs on the Python engine (never the Rust counts-only path).
    """

    def __init__(self, duplication: float = 0.0, transfer: float = 0.0,
                 loss: float = 0.0, origination: float = 0.0,
                 *, per: str = "copy",
                 inversion: float = 0.0, transposition: float = 0.0,
                 translocation: float = 0.0,
                 insertion: float = 0.0, deletion: float = 0.0,
                 conversion: float = 0.0,
                 chromosome_origination: float = 0.0, chromosome_loss: float = 0.0,
                 fission: float = 0.0, fusion: float = 0.0,
                 carrying_capacity: float | None = None):
        if per not in ("copy", "lineage", "shared"):
            raise ValueError(f"per must be 'copy', 'lineage', or 'shared', got {per!r}")
        self.per = per
        # per-event opportunity: duplication/transfer/loss may each be a Per(unit, rate) that
        # overrides the model-level `per` (e.g. duplication=Per("shared", ...) with per-copy loss).
        duplication, dup_unit = _opportunity(duplication, per)
        transfer, transfer_unit = _opportunity(transfer, per)
        loss, loss_unit = _opportunity(loss, per)
        self._units = {EventType.DUPLICATION: dup_unit, EventType.TRANSFER: transfer_unit,
                       EventType.LOSS: loss_unit}
        for u in self._units.values():
            if u not in ("copy", "lineage", "shared"):
                raise ValueError(f"opportunity unit must be 'copy', 'lineage', or 'shared', got {u!r}")
        rates = (("duplication", duplication), ("transfer", transfer), ("loss", loss),
                 ("origination", origination), ("inversion", inversion),
                 ("transposition", transposition), ("translocation", translocation),
                 ("insertion", insertion),
                 ("deletion", deletion), ("conversion", conversion),
                 ("chromosome_origination", chromosome_origination),
                 ("chromosome_loss", chromosome_loss), ("fission", fission), ("fusion", fusion))
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
        # translocation: move an arc to another chromosome (nucleotide model; needs >= 2 chromosomes)
        self.translocation = float(translocation)
        # intergenic indels: only fired by the nucleotide genome (per-nucleotide rate)
        self.insertion = float(insertion)
        self.deletion = float(deletion)
        # intra-genome gene conversion: only fired by multiset genomes with >= 2 copies of a family
        self.conversion = float(conversion)
        # chromosome-tier events: only fired by genomes that support them (OrderedGenome). Off by
        # default (rate 0), so a run that does not set them is byte-identical to before.
        self.chromosome_origination = float(chromosome_origination)
        self.chromosome_loss = float(chromosome_loss)
        self.fission = float(fission)
        self.fusion = float(fusion)
        self.carrying_capacity = _check_carrying_capacity(carrying_capacity)
        if not self._all_copy:
            per_copy_only = {
                "inversion": self.inversion, "transposition": self.transposition,
                "translocation": self.translocation, "insertion": self.insertion,
                "deletion": self.deletion, "conversion": self.conversion,
                "chromosome_origination": self.chromosome_origination,
                "chromosome_loss": self.chromosome_loss, "fission": self.fission,
                "fusion": self.fusion,
            }
            set_now = [k for k, v in per_copy_only.items() if v]
            if set_now or self.carrying_capacity is not None:
                bad = set_now + (["carrying_capacity"] if self.carrying_capacity is not None else [])
                raise ValueError(
                    "rearrangements / carrying_capacity require every event per-copy (no "
                    f"per='lineage'/'shared' and no Per(...) overrides); remove {', '.join(bad)}")
        if self.transfer > 0 and self._units[EventType.TRANSFER] == "shared":
            raise ValueError("a shared transfer clock is not yet supported; use per-copy/lineage "
                             "transfer (duplication and loss can be shared)")

    @property
    def has_shared(self) -> bool:
        """True if any per-family event rides the ``shared`` opportunity (a tree-wide clock the
        simulator drives through its shared pool)."""
        return any(u == "shared" for u in self._units.values())

    @property
    def _all_copy(self) -> bool:
        """True if every per-family event is per-copy (the default) — the plain model the Rust engine
        implements and the only one that carries rearrangements / carrying_capacity."""
        return all(u == "copy" for u in self._units.values())

    def _mixed_weights(self, genome):
        """Per-branch weights for a lineage/shared/mixed model: per-copy events scale by copy count,
        per-lineage events are constant, shared events go to the simulator's pool (not emitted here)."""
        out: list[EventWeight] = []
        n = genome.size()
        for event, rate in ((EventType.DUPLICATION, self.duplication),
                             (EventType.TRANSFER, self.transfer), (EventType.LOSS, self.loss)):
            if rate > 0 and n > 0:
                unit = self._units[event]
                if unit == "copy":
                    out.append(EventWeight(event, None, rate * n))
                elif unit == "lineage":
                    out.append(EventWeight(event, None, rate))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        return out

    def _regulated(self) -> bool:
        # Only the *soft* carrying capacity needs per-family duplication weights. A hard
        # max_family_size is enforced cheaply by the simulator (it skips an over-cap
        # duplication), so we keep the O(1) aggregate duplication term.
        return self.carrying_capacity is not None

    def event_weights(self, genome, branch, time):
        if not self._all_copy:
            # any lineage/shared/mixed per-event opportunity → the simple per-event router (shared
            # events go to the simulator's pool, not here). The full per-copy body below (rearrangements,
            # regulation, conversion, chromosomes) runs only when every event is per-copy.
            return self._mixed_weights(genome)

        n = genome.size()   # every per-family event is per-copy (the plain, Rust-eligible model)
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
        if self.conversion > 0 and n > 0:
            # intra-genome gene conversion: a same-family copy overwrites another. It needs a donor
            # AND a recipient, so only families with >= 2 copies qualify; a family with cn copies
            # converts at total rate conversion * cn. Family-specific (never the family=None fast
            # path), so the simulator is told which family to convert.
            for family in genome.families():
                cn = genome.copy_number(family)
                if cn >= 2:
                    out.append(EventWeight(EventType.CONVERSION, family, self.conversion * cn))
        if self.origination > 0:
            out.append(EventWeight(EventType.ORIGINATION, None, self.origination))
        # chromosome-tier events (OrderedGenome opts in via supported_events; every branch is gated
        # on rate > 0, so a model that leaves these at 0 adds nothing and stays byte-identical).
        n_chrom = len(getattr(genome, "chromosomes", ()))
        if self.chromosome_origination > 0:
            out.append(EventWeight(EventType.CHROMOSOME_ORIGINATION, None, self.chromosome_origination))
        if self.fission > 0 and n > 0:  # a chromosome needs genes to be worth splitting
            out.append(EventWeight(EventType.FISSION, None, self.fission * n_chrom))
        if n_chrom >= 2:  # loss / fusion / translocation need a second chromosome
            if self.chromosome_loss > 0:
                out.append(EventWeight(EventType.CHROMOSOME_LOSS, None, self.chromosome_loss * n_chrom))
            if self.fusion > 0:
                out.append(EventWeight(EventType.FUSION, None, self.fusion * n_chrom))
            if self.translocation > 0 and n > 0:  # per-nucleotide: move an arc to another chromosome
                out.append(EventWeight(EventType.TRANSLOCATION, None, self.translocation * n))
        return out

    def shared_event_weights(self, family):
        out = []
        for event, rate in ((EventType.DUPLICATION, self.duplication), (EventType.LOSS, self.loss)):
            if rate > 0 and self._units[event] == "shared":
                out.append((event, rate))
        return out


class PerCopyRates(Rates):
    """Deprecated preset for ``Rates(per="copy")`` — one clock per gene copy (exponential families,
    the built-in Rust model). The opportunity is now a knob on :class:`Rates`, so this named class is
    redundant: it still works but warns and has left ``zombi2.__all__``; removed in 0.4.0."""

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "PerCopyRates is deprecated; use Rates(..., per='copy') (the default). "
            "The old name still works but is removed in 0.4.0.",
            DeprecationWarning, stacklevel=2,
        )
        super().__init__(*args, per="copy", **kwargs)


class PerLineageRates(Rates):
    """Deprecated preset for ``Rates(per="lineage")`` — one clock per genome (constant D/T/L totals,
    a uniformly chosen target, so families grow *linearly* not exponentially). The opportunity is now
    a knob on :class:`Rates`: this named class is redundant, warns, and has left ``zombi2.__all__``;
    removed in 0.4.0."""

    def __init__(self, duplication: float = 0.0, transfer: float = 0.0,
                 loss: float = 0.0, origination: float = 0.0):
        warnings.warn(
            "PerLineageRates is deprecated; use Rates(..., per='lineage'). "
            "The old name still works but is removed in 0.4.0.",
            DeprecationWarning, stacklevel=2,
        )
        super().__init__(duplication, transfer, loss, origination, per="lineage")


#: Backwards-compatible alias for the (now deprecated) per-copy preset. ``SharedRates`` was the
#: original name for ``PerCopyRates``; both now point at the same deprecated class object.
SharedRates = PerCopyRates


class FamilySampledRates(RateModel):
    """Each gene family has its OWN duplication/transfer/loss rates (ZOMBI-1 style).

    By default every family **draws** its ``(d, t, l)`` from the given distributions at first
    sighting. Pass ``rates`` to instead **fix** the rates of specific families by name: a
    ``{family_id: (dup, transfer, loss)}`` map (e.g. ``{"1": (3, 2, 1), "2": (4, 0, 1)}``). Listed
    families use their tabulated triple; families **not** listed fall back to drawing from the
    distributions — so with the default rates of ``0`` (``Fixed(0)``) the unlisted families are
    simply inert, and the tabulated ones are exactly as specified. This is the hand-specified
    per-family model reachable from the CLI via ``--family-rates FILE`` (see
    :func:`~zombi2.genomes.read_rates.read_family_rates`).
    """

    def __init__(self, duplication=0.0, transfer=0.0, loss=0.0, origination: float = 0.0,
                 *, per: str = "copy", rates: dict | None = None,
                 carrying_capacity: float | None = None):
        if per not in ("copy", "lineage"):
            raise ValueError(f"per must be 'copy' or 'lineage', got {per!r}")
        self.per = per
        self._dup = as_distribution(duplication)
        self._trans = as_distribution(transfer)
        self._loss = as_distribution(loss)
        if origination < 0:
            raise ValueError(f"origination rate must be >= 0, got {origination}")
        self.origination = float(origination)
        self.carrying_capacity = _check_carrying_capacity(carrying_capacity)
        if self.per == "lineage" and self.carrying_capacity is not None:
            raise ValueError("per='lineage' does not support carrying_capacity (a per-copy notion)")
        self._fixed: dict[str, tuple[float, float, float]] = {}
        for fam, triple in (rates or {}).items():
            d, t, l = triple
            if min(d, t, l) < 0:
                raise ValueError(f"family {fam!r}: rates must be >= 0, got {triple}")
            self._fixed[str(fam)] = (float(d), float(t), float(l))
        self._family_rates: dict[str, tuple[float, float, float]] = {}

    def bind(self, rng, max_family_size: int | None = None, tree=None) -> None:
        super().bind(rng, max_family_size, tree)
        self._family_rates = {}

    def rates_for(self, family: str) -> tuple[float, float, float]:
        """The (dup, transfer, loss) rates for a family: the tabulated triple if the family is in
        the ``rates`` map, otherwise sampled from the distributions and cached on first use."""
        fixed = self._fixed.get(family)
        if fixed is not None:
            return fixed
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
        lineage = self.per == "lineage"
        for family in genome.families():
            cn = genome.copy_number(family)
            if cn == 0:
                continue
            d, t, l = self.rates_for(family)
            if lineage:
                # one clock per genome for this family: constant totals (not × copy number),
                # a uniformly chosen copy within the family
                if d > 0:
                    out.append(EventWeight(EventType.DUPLICATION, family, d))
                if t > 0:
                    out.append(EventWeight(EventType.TRANSFER, family, t))
                if l > 0:
                    out.append(EventWeight(EventType.LOSS, family, l))
                continue
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


_EVENT_BY_NAME = {
    "duplication": EventType.DUPLICATION,
    "transfer": EventType.TRANSFER,
    "loss": EventType.LOSS,
}


def _resolve_events(events) -> frozenset:
    """Normalise a set of scalable event kinds (names or ``EventType``) to a frozenset. ``None``
    means all of duplication/transfer/loss (the default)."""
    if events is None:
        return frozenset((EventType.DUPLICATION, EventType.TRANSFER, EventType.LOSS))
    out = set()
    for e in events:
        if isinstance(e, EventType):
            out.add(e)
        else:
            key = str(e).lower()
            if key not in _EVENT_BY_NAME:
                raise ValueError(f"unknown event {e!r}; use duplication / transfer / loss")
            out.add(_EVENT_BY_NAME[key])
    if not out:
        raise ValueError("events must name at least one of duplication / transfer / loss")
    return frozenset(out)


class Modifier(ABC):
    """A context-keyed multiplier on an event's rate — the composable building block behind
    per-branch, per-family and lineage-pair heterogeneity (see the *Modifiers* section of the rates
    primer, ``docs/guide/rates.md``).

    A modifier reads whatever context it cares about — the ``event`` kind, the ``family`` it acts on
    (or ``None``), the ``branch``, the ``time`` — and returns a multiplier on the base rate; ``1.0``
    means no effect. Modifiers compose **multiplicatively** and default to 1, so stacking them is
    order-independent. This is the *emission-seam* modifier (it scales how often an event fires); its
    *recipient-seam* counterpart, transfer donor->recipient bias, lives on
    :class:`~zombi2.genomes.transfers.TransferModel`. See ``docs/design/rate-modifiers.md``.
    """

    time_dependent = False

    #: True if this modifier keys on the *family* an event acts on. When any modifier in a
    #: :class:`ModifiedRates` stack sets this, the composer expands a base's aggregate ``family=None``
    #: weights into one weight per family (family f's share of a uniform-target total rate ``r`` is
    #: ``r * n_f / n``), so the family factor has a family to attach to.
    keys_on_family = False

    def bind(self, rng, tree=None) -> None:
        """One-time setup at the start of a simulation (e.g. autocorrelated factors need the tree)."""

    @abstractmethod
    def factor(self, event: EventType, family, branch: str, time: float) -> float:
        ...

    def refresh_times(self, t0: float, t1: float) -> list[tuple[float, str]]:
        """Times in ``(t0, t1)`` at which this modifier's factor changes on its own
        (see :meth:`RateModel.refresh_times`). Default: none."""
        return []


class ModifiedRates(RateModel):
    """A base rate model with a stack of emission-seam :class:`Modifier` s applied to its weights.

    Each candidate event's rate is multiplied by the product of the modifiers' factors for that
    ``(event, family, branch, time)``. An empty stack is exactly ``base``; ``target_params``,
    ``establishment_probability`` and ``refresh_times`` compose base and modifiers. This is the
    emission-seam composer of the rate-modifier design (``docs/design/rate-modifiers.md``); the plain
    :class:`PerCopyRates` fast path is untouched — any ``ModifiedRates`` runs on the Python engine.
    """

    def __init__(self, base: RateModel, modifiers):
        # Modifiers scale the *per-branch* weights; a per="shared" base routes duplication/loss through
        # a tree-wide pool that a ModifiedRates would not see (its shared_event_weights would be the
        # empty default), silently switching the shared clock off — so reject the combination.
        if getattr(base, "per", None) == "shared":
            raise ValueError("modifiers (LineageRates / ModifiedRates / FamilyModifier) do not yet "
                             "compose with a per='shared' base; the shared clock has no per-branch "
                             "weights for a modifier to scale")
        self.base = base
        self.modifiers = list(modifiers)
        # A per-family modifier needs per-family weights, so expand aggregate family=None weights.
        self._expand_families = any(getattr(m, "keys_on_family", False) for m in self.modifiers)

    @property
    def time_dependent(self) -> bool:
        return self.base.time_dependent or any(m.time_dependent for m in self.modifiers)

    def bind(self, rng, max_family_size: int | None = None, tree=None) -> None:
        super().bind(rng, max_family_size, tree)
        self.base.bind(rng, max_family_size, tree)
        for m in self.modifiers:
            m.bind(rng, tree)

    def event_weights(self, genome, branch, time):
        out = []
        for ew in self.base.event_weights(genome, branch, time):
            for w in (self._per_family(ew, genome) if self._expand_families else (ew,)):
                factor = 1.0
                for m in self.modifiers:
                    factor *= m.factor(w.event, w.family, branch, time)
                out.append(w if factor == 1.0 else EventWeight(w.event, w.family, w.rate * factor))
        return out

    @staticmethod
    def _per_family(ew, genome):
        """Expand an aggregate ``family=None`` weight (a total rate with a uniformly chosen target
        copy) into one weight per family, ``rate * n_f / n`` — distributionally identical but with a
        family a :class:`FamilyModifier` can key on. Origination (which mints a *new* family) and
        weights that already name a family pass through unchanged."""
        if ew.family is not None or ew.event is EventType.ORIGINATION:
            return (ew,)
        n = genome.size()
        if n <= 0:
            return (ew,)
        return tuple(EventWeight(ew.event, f, ew.rate * genome.copy_number(f) / n)
                     for f in genome.families())

    def target_params(self, event, genome, branch, time):
        return self.base.target_params(event, genome, branch, time)

    def establishment_probability(self, selection, recipient_genome, time):
        return self.base.establishment_probability(selection, recipient_genome, time)

    def refresh_times(self, t0, t1):
        times = list(self.base.refresh_times(t0, t1))
        for m in self.modifiers:
            times.extend(m.refresh_times(t0, t1))
        return times


class LineageModifier(Modifier):
    """Per-lineage factor — the emission-seam modifier behind :class:`LineageRates` (a relaxed clock
    is exactly this: a rate multiplier that varies from lineage to lineage down the tree).

    Provide exactly one source: ``autocorr_sigma`` (relaxed clock), ``per_branch`` (i.i.d. per
    lineage), or ``factors`` (an explicit ``{lineage: factor}`` map). ``events`` selects which event
    kinds it scales (default duplication/transfer/loss); ``root_rate`` is the root lineage's factor
    and the fallback for lineages an explicit map omits.
    """

    def __init__(self, *, autocorr_sigma: float | None = None, per_branch=None,
                 factors: dict | None = None, root_rate: float = 1.0, events=None):
        sources = [autocorr_sigma is not None, per_branch is not None, factors is not None]
        if sum(sources) != 1:
            raise ValueError("specify exactly one of autocorr_sigma, per_branch, or factors")
        if autocorr_sigma is not None and autocorr_sigma < 0:
            raise ValueError("autocorr_sigma must be >= 0")
        self.autocorr_sigma = autocorr_sigma
        self.per_branch = as_distribution(per_branch) if per_branch is not None else None
        self.explicit = dict(factors) if factors is not None else None
        self.root_rate = float(root_rate)
        self._scaled = _resolve_events(events)
        self._factor: dict[str, float] = {}
        self._rng = None

    def bind(self, rng, tree=None) -> None:
        self._rng = rng
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

    def factor(self, event, family, branch, time):
        if event not in self._scaled:
            return 1.0
        return self._branch_factor(branch)


class LineageRates(ModifiedRates):
    """Make rates vary from lineage to lineage by scaling a base rate model.

    Wraps any base rate model (``PerCopyRates``, ``FamilySampledRates``, ...) and multiplies its
    weights on each lineage by a per-lineage factor. By default the factor scales
    duplication/transfer/loss together (origination is left unscaled); pass ``events`` to restrict it
    to specific event kinds — e.g. ``events=("transfer",)`` makes a branch more (or less) prone to
    **donating** a transfer without touching its duplication/loss. This is the transfer-*emission*
    dial; its counterpart, transfer *receptivity* (how likely a branch is to **receive**), lives on
    :class:`~zombi2.genomes.transfers.TransferModel`. It composes with the base model, so branch
    heterogeneity and family/uniform rates combine into a two-factor model.

    Provide exactly one factor source:

    * ``autocorr_sigma`` — a **relaxed clock**: the factor evolves down the tree as
      ``factor(child) = factor(parent) * exp(N(0, sigma * sqrt(branch_length)))``, so
      closely related lineages have similar rates. Needs the tree (passed via ``bind``).
      ``sigma = 0`` reduces to the base model.
    * ``per_branch`` — a distribution (built-in / scipy / callable) drawn **i.i.d. per
      branch**, independently.
    * ``factors`` — an explicit ``{branch_name: factor}`` map (branches not listed use
      ``root_rate``).

    ``root_rate`` (default 1.0) is the root-branch factor in ``autocorr_sigma`` mode and the fallback
    factor for branches absent from an explicit ``factors`` map. In ``per_branch`` mode every branch
    is drawn i.i.d. (``root_rate`` is unused) — harmless, since the root's zero-length branch carries
    no events. ``events`` selects which event kinds the factor scales (default: duplication, transfer,
    loss).
    """

    def __init__(self, base: RateModel, *, autocorr_sigma: float | None = None,
                 per_branch=None, factors: dict | None = None, root_rate: float = 1.0,
                 events=None):
        super().__init__(base, [LineageModifier(
            autocorr_sigma=autocorr_sigma, per_branch=per_branch, factors=factors,
            root_rate=root_rate, events=events)])

    @property
    def _factor(self) -> dict:
        """The per-branch factor map — kept for backward compatibility; it now lives on the
        underlying :class:`LineageModifier`."""
        return self.modifiers[0]._factor


class FamilyModifier(Modifier):
    """Per-gene-family multiplier — per-family rate heterogeneity as a composable **overlay**.

    Contrast :class:`FamilySampledRates`, which bakes per-family rates into a base model; a
    ``FamilyModifier`` instead multiplies any base's per-family rate by a family factor, so it stacks
    with other modifiers. Composed on :class:`PerLineageRates` it expresses **per-genome × per-family**
    rates (the combination the old ``--rate-model`` enum could not name); composed with a
    :class:`LineageModifier` it gives **family × branch** heterogeneity.

    Provide exactly one source: ``factors`` (an explicit ``{family: factor}`` map; families absent
    from it use ``root_factor``) or ``per_family`` (a distribution drawn i.i.d. per family at first
    sight and cached). ``events`` selects which event kinds it scales (default duplication/transfer/
    loss). Because it keys on the family, a :class:`ModifiedRates` carrying one runs on the Python
    engine (the aggregate weights are expanded per family).
    """

    keys_on_family = True

    def __init__(self, *, factors: dict | None = None, per_family=None,
                 root_factor: float = 1.0, events=None):
        if (factors is None) == (per_family is None):
            raise ValueError("specify exactly one of factors or per_family")
        self.explicit = ({str(k): float(v) for k, v in factors.items()}
                         if factors is not None else None)
        self.per_family = as_distribution(per_family) if per_family is not None else None
        self.root_factor = float(root_factor)
        self._scaled = _resolve_events(events)
        self._factor: dict[str, float] = {}
        self._rng = None

    def bind(self, rng, tree=None) -> None:
        self._rng = rng
        self._factor = dict(self.explicit) if self.explicit is not None else {}

    def _family_factor(self, family: str) -> float:
        f = self._factor.get(family)
        if f is None:
            f = (self.per_family.sample(self._rng)
                 if self.per_family is not None else self.root_factor)
            self._factor[family] = f
        return f

    def factor(self, event, family, branch, time):
        if family is None or event not in self._scaled:
            return 1.0
        return self._family_factor(str(family))


#: Backwards-compatible aliases — the rate vocabulary standardises on "lineage" (see
#: ``docs/design/rate-vocabulary.md``): "per genome" was always "per lineage" (one genome per
#: lineage), and the per-branch heterogeneity modifier is a per-lineage one. The old names remain
#: as the **same class objects**, so existing code and any ``type(...)`` checks keep working.
PerGenomeRates = PerLineageRates
BranchModifier = LineageModifier
BranchRates = LineageRates
