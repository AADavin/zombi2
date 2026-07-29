"""The gene-genealogy event log — the shared source of truth every resolution writes.

An `Event` records one moment in a gene family's history; the per-family gene trees are
*derived* from a run's events (see `gene_trees`), identically whether the genome was an
multiset of families or an ordered set of chromosomes. Position and orientation are **not** here —
they live in the genome snapshots and the rearrangement log — because an event is about gene
*identity and descent*, which is resolution-blind. So this module is imported by both the family
core and the ordered engine, and neither owns it.
"""

from __future__ import annotations

from dataclasses import dataclass

# the node spelling belongs to the tree, not to this level — re-exported here because every
# genome module already reaches for it alongside the gene labels below
from ..tree import node_from_label, node_label  # noqa: F401


@dataclass(frozen=True)
class Event:
    """A recorded genome event — the true history every per-family gene tree is derived from. Gene
    ids are **per segment** (the ZOMBI1 model): every event ends a gene and starts fresh ids for its
    descendants, so an id belongs to exactly one species branch and every node's genome has its own
    ids. ``lineage`` is the species-tree node the event fired on; ``time`` is when (crown-forward,
    the species tree's clock). By kind:

    - ``"origination"`` — ``copy`` is a founding gene of a fresh family (``parent`` ``None``): a root.
    - ``"duplication"`` — the gene ``parent`` ends; ``copy`` is one of its **two** descendants (the
      continuation and the new copy — two rows, same ``parent``), both on ``lineage``.
    - ``"transfer"`` — the donor gene ``parent`` ends; ``copy`` is one of its two descendants: the
      continuation on the donor ``lineage``, or the transferred copy on the ``recipient`` lineage (a
      horizontal edge). Two rows, same ``parent``, and **both name ``donor``** — the branch the
      material left. Without it the arriving row said only where the copy landed (twice over:
      ``lineage`` and ``recipient`` are the same branch there), so reading who donated to whom meant
      pairing the two rows on ``(time, parent)``. A transfer is an edge; each of its rows names both
      ends.
    - ``"speciation"`` — the gene ``parent`` ends at a split; ``copy`` is its descendant in daughter
      species ``lineage`` (one row per daughter — two, same ``parent``).
    - ``"loss"`` — the gene ``copy`` ends with no descendant (``parent`` ``None``).
    """

    time: float
    kind: str  # "origination" | "duplication" | "loss" | "transfer" | "speciation"
    lineage: int  # the species-tree node id where it fired (for a transfer: the donor lineage)
    family: int
    copy: int  # the copy born (origination / duplication / transfer) or removed (loss)
    parent: int | None = None  # duplication & transfer: the source copy (which survives)
    recipient: int | None = None  # transfer only: the species lineage the new copy is born on
    #: transfer only: the species lineage the material left. Set on **both** rows, so either names
    #: the whole edge; on the donor's own row it repeats ``lineage``, which is the price of that.
    donor: int | None = None

    @property
    def event(self) -> int:
        """Which event this row belongs to: **the gene copy whose fate the event is**.

        A row is a gene-tree *edge*, so an event that leaves two descendants writes two rows. Those
        rows agree here and nowhere else that is safe to group on — ``time`` is a float, and pairing
        on it invites both false joins and float-equality bugs. A copy ends exactly once, so within a
        kind this is unique: ``sort -u`` on it counts events rather than edges, which counting rows
        does not. For the kinds that end a copy (duplication, transfer, speciation) that copy is
        ``parent``; a loss ends ``copy`` itself, and an origination begins one with no parent, so
        those name ``copy``. Both are already one row per event — the column exists so a single rule
        works for every kind and no reader has to know which of the two to group on."""
        return self.copy if self.parent is None else self.parent


_COLS = ("time", "kind", "lineage", "family", "copy", "parent", "recipient", "donor", "event")

#: Columns holding a species-tree node, written as ``n<id>``.
_NODE_COLS = frozenset({"lineage", "recipient", "donor"})

#: Columns holding a gene copy, written as ``g<id>``. This table names the species in its own column,
#: so the copy cell carries the copy alone; a file with no such column names it
#: ``n<species>_g<copy>`` instead (`copy_label()`). Either way the ``g<id>`` half is the same
#: token, so a copy joins across every file without translation. ``parent`` is a gene copy too (the
#: source copy a duplication/transfer descends from), so it is g-labelled; the species columns are not.
_GENE_COLS = frozenset({"copy", "parent", "event"})


def gene_label(copy_id: int | None) -> str:
    """A gene copy as every ZOMBI2 **table** writes it: ``g<id>``. Empty for ``None``.

    A table names the copy's species in a column of its own, so this is the whole cell. Where there is
    no such column — a Newick leaf, a FASTA record — the copy is written ``n<species>_g<copy>``
    (`copy_label()`), which embeds this unchanged, so a copy still joins across files without
    translation. The ``g`` also keeps a copy id from being read as a bare number in a Newick leaf,
    where that is ambiguous with a support value or a branch length."""
    return "" if copy_id is None else f"g{copy_id}"


def copy_label(species: int | None, copy_id: int | None, fate: str | None = None) -> str:
    """A gene copy **where it sits**: ``n<species>_g<copy>``, the token every file naming a copy of a
    known species uses — the gene-tree and phylogram leaves, the alignment FASTA records and the
    homology tables.

    The copy id alone is unique, so the species is redundant for joining; it is here because a
    sequence or a tip is read by people and by tools that never see the rest of the run, and ``g2179``
    alone cannot say which genome it came from. That forced everyone benchmarking orthology to join
    the alignments back to ``genomes.tsv`` themselves. The two halves are each still their own label,
    so splitting on the single ``_`` recovers them; ``_`` rather than ``|`` because a FASTA record
    name goes on to be parsed by aligners and tree builders that treat ``|`` as a field separator."""
    return f"{node_label(species, fate)}_{gene_label(copy_id)}"


def gene_from_label(cell: str) -> int:
    """The inverse of `gene_label()`. A bare integer is accepted too, so a table written before
    the copy columns carried their ``g`` still replays."""
    return int(cell[1:] if cell[:1] == "g" else cell)


#: the event-log header line — the column names, tab-joined. Shared so a streamed shard and the
#: whole-log writer put the same header on the same columns.
EVENTS_HEADER = "\t".join(_COLS)


def _name(names: "dict[int, str] | None", node_id: int | None) -> str:
    """One node's written name: from the run's ``names`` map when there is one, else plain ``n<id>``.

    The map is the only thing that knows a lineage died, so a writer that can name a dead branch takes
    it; one that provably cannot (a profile's columns, an alignment's records — extant tips both) does
    not need to, and gets the same answer without it."""
    return node_label(node_id) if names is None or node_id is None else names[node_id]


def _cell(e: Event, col: str, names: dict[int, str] | None) -> str:
    v = getattr(e, col)
    if v is None:
        return ""
    if col in _NODE_COLS:
        return _name(names, v)
    if col in _GENE_COLS:
        return gene_label(v)
    return str(v)


def event_rows(events: list[Event], names: dict[int, str] | None = None) -> list[str]:
    """The event rows **without** the header — one tab-joined line per event. The one row format, so a
    streamed per-worker shard and `events_tsv()` cannot drift; the shard writes only rows and the
    finalize prepends `EVENTS_HEADER` once.

    ``names`` is the run's node names (`Tree.labels()`) — needed because a lineage that went
    extinct is written ``e<id>``, and an event log names *every* branch, the dead ones included.
    Omitting it names them all ``n<id>``, which is right only where no dead branch can appear."""
    return ["\t".join(_cell(e, c, names) for c in _COLS) for e in events]


def events_tsv(events: list[Event], names: dict[int, str] | None = None) -> str:
    """The event log as TSV — one row per event; empty cells for the fields a kind does not use.
    ``names`` as in `event_rows()`."""
    return "\n".join([EVENTS_HEADER, *event_rows(events, names)]) + "\n"


def events_from_tsv(text: str) -> list[Event]:
    """Parse the TSV `events_tsv()` writes back into a ``list[Event]`` — the deserializer twin, so
    a written ``genome_events.tsv`` can be replayed (a downstream level's gene trees are derived from
    the log by `gene_trees_from_events()`). ``time`` is a float, the
    id columns are ints, and the optional ``parent`` / ``recipient`` are ints or ``None`` (empty)."""
    lines = text.splitlines()
    if not lines:
        raise ValueError("empty genome event log — is the file empty?")
    header = lines[0].split("\t")
    if tuple(header[:len(_COLS)]) == _COLS:
        # the ordered resolution writes the same genealogy with each event's position beside it, and
        # its rearrangements in the same table. The extra columns are ignored here and the
        # rearrangement rows skipped: this reader is about identity and descent, which they do not
        # touch. (A genome level that needs the positions reads them itself.)
        return _parse(lines, header)
    if tuple(header) != _COLS:
        # `block_events.tsv` — a nucleotide run's own record, keyed by ancestral interval. Every
        # resolution writes its genealogy here in one format, so that file is the only table left
        # that looks like this one and is not; naming it saves reading the columns to find out.
        hint = ("; this looks like a nucleotide run's block_events.tsv, whose rows are ancestral "
                "intervals rather than gene-tree edges. Read one with "
                "zombi2.genomes.nucleotide.read_nucleotide_genomes. That run's genome_events.tsv is "
                "the genealogy, in this format, and is what this reader wants"
                if "source" in header and "family" not in header else "")
        if not hint and set(header) < set(_COLS):
            # a genuinely short header: every column it has is one of ours, so it is a log this
            # version no longer reads rather than a different table
            hint = (f"; it is missing {', '.join(c for c in _COLS if c not in header)}. Re-run "
                    f"'zombi2 genomes' to write a log this version reads")
        raise ValueError(f"unexpected genome-event columns {header}; expected {list(_COLS)}{hint}")
    return _parse(lines, header)


#: the kinds this log records about gene identity and descent. A wider table may also carry the
#: ancestry-**neutral** rearrangements; they end no gene lineage, so they are not events here.
_GENEALOGY = frozenset({"origination", "duplication", "loss", "transfer", "speciation"})


def _parse(lines: list[str], header: list[str]) -> list[Event]:
    """Read the rows by column **name**, so a table carrying more than the canonical columns parses
    unchanged and only the genealogy rows come back. ``event`` is derived from ``parent``/``copy``,
    so it is written but never read back: there is one definition of it, on `Event`."""
    at = {c: i for i, c in enumerate(header)}
    events: list[Event] = []
    for lineno, raw in enumerate(lines[1:], 2):
        if not raw:                                     # tolerate a trailing blank line
            continue
        cells = raw.split("\t")
        if len(cells) != len(header):
            raise ValueError(f"genome event log line {lineno}: expected {len(header)} columns, "
                             f"got {len(cells)}")
        if cells[at["kind"]] not in _GENEALOGY:
            continue
        get = lambda c: cells[at[c]]                    # noqa: E731
        events.append(Event(
            time=float(get("time")), kind=get("kind"), lineage=node_from_label(get("lineage")),
            family=int(get("family")), copy=gene_from_label(get("copy")),
            parent=gene_from_label(get("parent")) if get("parent") else None,
            recipient=node_from_label(get("recipient")) if get("recipient") else None,
            donor=node_from_label(get("donor")) if get("donor") else None))
    return events
