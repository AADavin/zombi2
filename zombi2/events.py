"""Shared vocabulary for ZOMBI2 events.

Everything in this module is deliberately representation-agnostic: it is the common
language spoken by the genome, the rate model, the simulator loop and the I/O layer.
The types are intentionally shaped for the *future* (segment-based events, transfers
that carry an insertion point, per-gene log rows) even though v1 only ever uses the
trivial, single-gene forms.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # avoid a runtime import cycle events <-> genome
    from .genome import Gene


class EventType(Enum):
    """The kinds of events that can appear in an event log.

    The *stochastic* subset (:data:`STOCHASTIC_EVENTS`) is what a rate model can fire.
    ``SPECIATION`` and ``LEAF`` are structural markers written by the driver, not
    sampled. New event kinds (e.g. ``INVERSION``, ``TRANSPOSITION``) are added by
    appending members here — the simulator loop iterates over the events a genome
    declares it supports, so no call site changes.
    """

    ORIGINATION = "O"
    DUPLICATION = "D"
    TRANSFER = "T"
    LOSS = "L"
    # --- structural markers (logged, never sampled) ---
    SPECIATION = "S"
    LEAF = "F"
    # --- future: INVERSION = "I", TRANSPOSITION = "P" ---


#: The events a rate model may fire (order-free v1 set).
STOCHASTIC_EVENTS: tuple[EventType, ...] = (
    EventType.ORIGINATION,
    EventType.DUPLICATION,
    EventType.TRANSFER,
    EventType.LOSS,
)


@dataclass(frozen=True)
class Region:
    """Optional positional metadata for a :class:`Selection`.

    Unused in v1 (selections carry ``region=None``). The ordered-genome extension
    fills this in with a contiguous run; kept here so that adding order is additive.
    """

    chromosome: int
    start: int
    length: int
    strand: int = 1


@dataclass(frozen=True)
class Selection:
    """What an event acts on.

    ``genes`` is always populated (v1: exactly one gene). ``region`` carries positional
    information only for ordered genomes. Never replace this with a ``SingleGene`` type:
    keeping it segment-shaped is what lets segment events slot in without a signature
    change.
    """

    genes: tuple["Gene", ...]
    region: Region | None = None

    @property
    def family(self) -> str:
        return self.genes[0].family


@dataclass(frozen=True)
class TargetParams:
    """Parameters handed from the rate model to :meth:`Genome.draw_target`.

    v1 draws a single gene and ignores every field. ``extension`` (segment length) and
    ``extent_mode`` ("gene_count" vs "nucleotide") are read by future extensions.
    """

    extension: float | None = None
    extent_mode: str = "gene_count"


@dataclass
class GeneOp:
    """One gene's involvement in an event — a single row of an :class:`EventRecord`.

    ``role`` describes what happened to this gene, e.g. ``"origin"``, ``"parent"``,
    ``"copy"``, ``"lost"``, ``"donor_kept"``, ``"transferred"``. Optional fields
    (``orientation``, ``length``) stay ``None`` in v1.
    """

    gid: str
    family: str
    role: str
    orientation: int | None = None
    length: int | None = None


class InsertionPoint(Enum):
    """Where a transferred segment lands in the recipient.

    An order-free multiset has no meaningful position, so v1 always uses ``ANYWHERE``.
    Ordered genomes return richer insertion descriptions from
    :meth:`Genome.choose_insertion_point`.
    """

    ANYWHERE = "anywhere"


@dataclass
class TransferSegment:
    """A portable, self-contained copy of the genes moved by a transfer.

    Carries freshly minted gene ids so the recipient can adopt them directly. ``genes``
    is the copy that lands in the recipient; the donor keeps its own copy when
    ``extract_segment`` is called with ``keep_copy=True``.
    """

    family: str
    genes: tuple["Gene", ...]
    replacement: bool = False


@dataclass
class EventRecord:
    """A single logged event.

    In v1 every record concerns exactly one gene family (all its :class:`GeneOp` rows
    share a family), so per-family logs are just a filter. ``genes`` is a *list* of rows
    — length 1 in v1, longer once segment events exist — which is why the log schema
    never has to change.
    """

    event: EventType
    branch: str
    time: float
    genes: list[GeneOp]
    donor: str | None = None
    recipient: str | None = None
    insertion: object | None = None

    @property
    def family(self) -> str:
        return self.genes[0].family


@dataclass
class EventLog:
    """A chronological collection of :class:`EventRecord`s with per-family views."""

    records: list[EventRecord] = field(default_factory=list)

    def add(self, record: EventRecord) -> None:
        self.records.append(record)

    def by_family(self) -> dict[str, list[EventRecord]]:
        out: dict[str, list[EventRecord]] = {}
        for r in self.records:
            out.setdefault(r.family, []).append(r)
        return out

    def __iter__(self):
        return iter(self.records)

    def __len__(self) -> int:
        return len(self.records)
