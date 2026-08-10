"""The gene-genealogy event log — the shared source of truth every resolution writes.

An `GeneEdge` records one moment in a gene family's history; the per-family gene trees are
*derived* from a run's events (see `gene_trees`), identically whether the genome was a multiset of
families, an ordered set of chromosomes or a coordinate space of blocks. Position and orientation
are **not** here — they live in the genome snapshots and the rearrangement log — because an event is
about gene *identity and descent*, which is resolution-blind. So this module is imported by all
three genome engines, and none of them owns it.

**The file is one row per event**, five columns wide — ``time kind family parents children`` — at the
family and nucleotide resolutions; the ordered one appends the event's coordinates to those five. An
event that ends one gene and starts two writes one row with one parent and two children, not a row
per descendant. The participants are written ``n<species>_g<copy>``, each carrying the branch it
lived on inside the token, which is what lets the ``lineage`` / ``recipient`` / ``donor`` columns go:
a transfer's row says who donated and who received by naming the two copies where they sit.

**In memory an `GeneEdge` is still one gene-tree edge.** That is what the engines emit and what
`gene_trees_from_edges()` reads — a gene tree is a parent→children graph, so the edge is the natural
object. The aggregation happens here, at write time, and `edges_from_tsv()` undoes it exactly, so
the file and the objects can each have the shape that suits them without either being a translation
of the other everywhere else.
"""

from __future__ import annotations

import collections
from dataclasses import dataclass

# the node spelling belongs to the tree, not to this level — re-exported here because every
# genome module already reaches for it alongside the gene labels below
from ..tree import node_from_label, node_label  # noqa: F401


@dataclass(frozen=True)
class GeneEdge:
    """A recorded genome event — the true history every per-family gene tree is derived from. Gene
    ids are **per segment** (the ZOMBI1 model): every event ends a gene and starts fresh ids for its
    descendants, so an id belongs to exactly one species branch and every node's genome has its own
    ids. ``lineage`` is the species-tree node the event fired on; ``time`` is when (origin-forward,
    the species tree's clock). By kind:

    - ``"origination"`` — ``copy`` is a founding gene of a fresh family (``parent`` ``None``): a root.
    - ``"duplication"`` — the gene ``parent`` ends; ``copy`` is one of its **two** descendants (the
      continuation and the new copy — two of these, same ``parent``), both on ``lineage``.
    - ``"transfer"`` — the donor gene ``parent`` ends; ``copy`` is one of its two descendants: the
      continuation on the donor ``lineage``, or the transferred copy on the ``recipient`` lineage (a
      horizontal edge). Two of these, same ``parent``, and **both name ``donor``** — the branch the
      material left. Without it the arriving edge said only where the copy landed (twice over:
      ``lineage`` and ``recipient`` are the same branch there), so reading who donated to whom meant
      pairing the two on ``(time, parent)``.
    - ``"speciation"`` — the gene ``parent`` ends at a split; ``copy`` is its descendant in daughter
      species ``lineage`` (one per daughter — two, same ``parent``).
    - ``"loss"`` — the gene ``copy`` ends with no descendant (``parent`` ``None``).

    These are the objects, one per gene-tree **edge**, which is the shape a gene tree is read in. The
    written log gathers them into one row per event (see the module docstring) and splits ``transfer``
    by ``replaced``: ``transfer_additive`` when the arriving copy joined the recipient's genome,
    ``transfer_replacing`` when it overwrote a resident there.
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
    #: replacing transfer only: the resident copy the arriving one overwrote, on **both** rows (as
    #: ``donor`` is). That copy dies, so it also has its own ``loss`` event — this names it, which is
    #: what lets the writer fold the two into one ``transfer_replacing`` row and the reader put them
    #: back. The alternative was to recognise the pair by their shared timestamp, and pairing rows on
    #: a float is exactly the bug this field exists to avoid.
    replaced: int | None = None

    @property
    def event(self) -> int:
        """Which event this row belongs to: **the gene copy whose fate the event is**.

        A row is a gene-tree *edge*, so an event that leaves two descendants is two of them. Those
        rows agree here and nowhere else that is safe to group on — ``time`` is a float, and pairing
        on it invites both false joins and float-equality bugs. A copy ends exactly once, so within a
        kind this is unique: it is what `event_rows()` groups the edges by to write one row per event.
        For the kinds that end a copy (duplication, transfer, speciation) that copy is ``parent``; a
        loss ends ``copy`` itself, and an origination begins one with no parent, so those name
        ``copy``. Both are already one event apiece — the property exists so a single rule works for
        every kind and no reader has to know which of the two to group on."""
        return self.copy if self.parent is None else self.parent


@dataclass(frozen=True)
class Event:
    """One genome event — **the same thing one row of ``genome_events.tsv`` is**.

    A `GeneEdge` is one edge of a gene tree, so a duplication is two of them (the copy that continues
    and the copy that is new) and a transfer likewise. An `Event` is the event itself: one object,
    with the copies it ended in ``parents`` and the copies it began in ``children``. Counting
    duplications in Python and counting them in the file give the same number, and a filter written
    against one works on the other — ``kind`` here is the file's vocabulary, so a transfer is
    ``transfer_additive`` or ``transfer_replacing`` rather than the edge form's bare ``transfer``.

    ``parents`` and ``children`` are gene ids in the order the row writes them: for a transfer, the
    copy left on the donor's branch leads. An origination has no parents and a loss no children.
    """

    time: float
    kind: str
    family: int
    parents: tuple[int, ...]
    children: tuple[int, ...]


_COLS = ("time", "kind", "family", "parents", "children")

#: TRANSITIONAL — the columns of the old log, one row per gene-tree **edge**, with the branch of each
#: participant in a column of its own. Nothing writes them any more — every resolution writes
#: `_COLS` — so `edges_from_tsv` reads them only to replay a log written by an older version.
_EDGE_COLS = ("time", "kind", "lineage", "family", "copy", "parent", "recipient", "donor", "event")

#: what separates the copies inside the ``parents`` / ``children`` cells. A cell is a *list*, so it
#: needs a separator that is neither the column separator nor part of a copy token — and ``;`` reads
#: as a list to a person and splits in one call in every language a reader might use.
_PACK = ";"

#: the kinds the log is written with. ``transfer`` is split by whether the arriving copy replaced a
#: resident, because the two are different events to anyone counting: an additive transfer grows the
#: recipient's genome, a replacing one does not, and only the second kills a gene.
_ADDITIVE, _REPLACING = "transfer_additive", "transfer_replacing"
_WRITTEN_KINDS = frozenset({"origination", "duplication", "loss", "speciation",
                            _ADDITIVE, _REPLACING})


def gene_label(copy_id: int | None) -> str:
    """A gene copy as every ZOMBI2 **table** writes it: ``g<id>``. Empty for ``None``.

    Where a table names the copy's species in a column of its own this is the whole cell. Where there
    is none — a Newick leaf, a FASTA record, this log's ``parents`` / ``children`` — the copy is
    written ``n<species>_g<copy>`` (`copy_label()`), which embeds this unchanged, so a copy still
    joins across files without translation. The ``g`` also keeps a copy id from being read as a bare
    number in a Newick leaf, where that is ambiguous with a support value or a branch length."""
    return "" if copy_id is None else f"g{copy_id}"


def _copy_cell(node_name: str, copy_id: int | None) -> str:
    """`copy_label()`'s spelling, from a node name already written. The one place the two halves are
    joined, so the log (whose names come from the run's ``labels()`` map) and the trees (whose names
    come from a fate) cannot spell a copy differently."""
    return f"{node_name}_{gene_label(copy_id)}"


def copy_label(species: int | None, copy_id: int | None, fate: str | None = None) -> str:
    """A gene copy **where it sits**: ``n<species>_g<copy>``, the token every file naming a copy of a
    known species uses — the gene-tree and phylogram leaves, the alignment FASTA records, the homology
    tables, and this log's ``parents`` and ``children``.

    The copy id alone is unique, so the species is redundant for joining; it is here because a
    sequence or a tip is read by people and by tools that never see the rest of the run, and ``g2179``
    alone cannot say which genome it came from. That forced everyone benchmarking orthology to join
    the alignments back to ``genomes.tsv`` themselves. The two halves are each still their own label,
    so splitting on the single ``_`` recovers them; ``_`` rather than ``|`` because a FASTA record
    name goes on to be parsed by aligners and tree builders that treat ``|`` as a field separator."""
    return _copy_cell(node_label(species, fate), copy_id)


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


def _branches(events: list[GeneEdge]) -> dict[int, int]:
    """``{copy id: the species branch it lived on}``, read off the log itself.

    A copy is born exactly once, on the branch it spends its whole life on, and its birth is some
    event's ``copy`` on that event's ``lineage`` — so the log already carries this and no engine has
    to pass it in. It is what lets a participant be written ``n2_g34``: the row a copy *dies* in
    (a speciation's parent, most of all) names the daughter branches, not the branch the parent
    lived on."""
    return {e.copy: e.lineage for e in events}


def event_counts(edges: list[GeneEdge], origin_time: float) -> dict[str, int]:
    """``{kind: how many events}`` — one number per *event*, which is what a row of
    ``genome_events.tsv`` is, and what every resolution's ``genome_summary.json`` reports.

    Counted from `GeneEdge`, which is one per gene-tree **edge**: a duplication, a transfer and a
    speciation each end one gene and start two, so counting edges inflates them exactly 2×. A gene
    ends at exactly one event, so distinct parents *are* the events. Origination and loss begin or end
    a single lineage and already have one edge apiece.

    ``loss`` therefore counts every gene that died, which under ``replacement`` is **more than the
    log's ``loss`` rows**: a copy displaced by an arriving transfer has no row of its own — it is the
    second parent of that ``transfer_replacing`` row, because its death and the transfer are one
    event. That gap is exactly what the migration guide warns a returning ZOMBI v1 user about, and
    this is the corrected number it tells them to trust.

    ``origination`` is split from ``initial``: the starting genome is logged as origination at the
    root's own start time, so a bare count is de-novo arrivals plus ``initial_families`` — a number
    nobody asked for. ``origin_time`` is the root's ``birth_time``, the cut between the two.

    One function because all three resolutions must agree here by construction. Only the family
    resolution used to report any of it, which left the two resolutions with the *larger* undercount
    (64% at ordered, measured) with no corrected figure to consult."""
    per_kind: dict[str, set] = collections.defaultdict(set)
    singles: collections.Counter = collections.Counter()
    for e in edges:
        if e.parent is None:                  # origination and loss: one row apiece already
            singles[e.kind] += 1
        else:
            per_kind[e.kind].add(e.parent)    # two rows, one parent, one event
    counts = {k: len(v) for k, v in sorted(per_kind.items())}
    counts.update(sorted(singles.items()))
    initial = sum(1 for e in edges if e.kind == "origination" and e.time <= origin_time)
    return {"initial": initial,
            "origination": counts.get("origination", 0) - initial,
            **{k: counts.get(k, 0) for k in ("duplication", "transfer", "loss", "speciation")}}


def events_from_edges(edges: list[GeneEdge]) -> list[Event]:
    """The edges gathered into ``(time, kind, family, parents, children)``, one per event, in the
    order the events were recorded.

    The key is `GeneEdge.event` — the copy whose fate the event is — paired with the kind, because a
    copy that an origination *begins* is the same id a later duplication *ends*. Grouping on ids
    rather than on ``time`` is the point: the two edges of one event carry the same copy id, and
    pairing them on a float is a bug waiting for the run that produces two events in the same
    microsecond.

    A replacing transfer's displaced copy is folded in here: it is named by ``replaced`` on the
    transfer, so its own ``loss`` edge is skipped and becomes the transfer's second parent. That is
    the whole of what the ``transfer_replacing`` kind means — one event, one row, two parents (the
    donor's copy and the copy it overwrote), and no phantom loss a reader would have to tell from a
    real one by matching timestamps."""
    folded = {e.replaced for e in edges if e.replaced is not None}
    rows: dict[tuple[str, int], tuple[float, str, int, list[int], list[int]]] = {}
    for e in edges:
        if e.kind == "loss" and e.copy in folded:
            continue                            # its death IS the replacing transfer (see above)
        key = (e.kind, e.event)
        row = rows.get(key)
        if row is None:
            kind, parents = e.kind, []
            if e.kind == "loss":
                parents.append(e.copy)          # a loss ends the copy it names: that is its parent
            elif e.parent is not None:
                parents.append(e.parent)
            if e.kind == "transfer":
                kind = _REPLACING if e.replaced is not None else _ADDITIVE
                if e.replaced is not None:
                    parents.append(e.replaced)
            rows[key] = row = (e.time, kind, e.family, parents, [])
        if e.kind == "loss":
            continue                            # nothing descends from it
        # a transfer reads donor → recipient: the copy left on the donor's branch leads, whichever
        # order the engine happened to emit the two edges in. `recipient` is set on the arriving row
        # alone, which is what tells the two sides apart.
        if e.kind == "transfer" and e.recipient is None:
            row[4].insert(0, e.copy)
        else:
            row[4].append(e.copy)
    return [Event(time, kind, family, tuple(parents), tuple(children))
            for time, kind, family, parents, children in rows.values()]


def event_rows(events: list[GeneEdge], names: dict[int, str] | None = None) -> list[str]:
    """The event rows **without** the header — one tab-joined line per event. The one row format, so a
    streamed per-worker shard and `events_tsv()` cannot drift; the shard writes only rows and the
    finalize prepends `EVENTS_HEADER` once.

    ``names`` is the run's node names (`Tree.labels()`) — needed because a lineage that went
    extinct is written ``e<id>``, and an event log names *every* branch, the dead ones included.
    Omitting it names them all ``n<id>``, which is right only where no dead branch can appear."""
    where = _branches(events)

    def cell(copies: tuple[int, ...]) -> str:
        return _PACK.join(_copy_cell(_name(names, where[c]), c) for c in copies)

    return [f"{e.time}\t{e.kind}\t{e.family}\t{cell(e.parents)}\t{cell(e.children)}"
            for e in events_from_edges(events)]


def events_tsv(events: list[GeneEdge], names: dict[int, str] | None = None) -> str:
    """The event log as TSV — one row per event; an empty cell where a kind has no parents (an
    origination) or no children (a loss). ``names`` as in `event_rows()`."""
    return "\n".join([EVENTS_HEADER, *event_rows(events, names)]) + "\n"


def edges_from_tsv(text: str) -> list[GeneEdge]:
    """Parse the TSV `events_tsv()` writes back into a ``list[GeneEdge]`` — the deserializer twin, so
    a written ``genome_events.tsv`` can be replayed (a downstream level's gene trees are derived from
    the log by `gene_trees_from_edges()`). Each row is expanded back into one `GeneEdge` per gene-tree
    edge, which is what every reader downstream of here expects."""
    lines = text.splitlines()
    if not lines:
        raise ValueError("empty genome event log — is the file empty?")
    header = lines[0].split("\t")
    if tuple(header[:len(_COLS)]) == _COLS:
        # a resolution may write more beside the genealogy — its events' positions, its
        # rearrangements in the same table. The extra columns are ignored here and rows of a kind
        # this log does not know are skipped: this reader is about identity and descent, which they
        # do not touch. (A genome level that needs the positions reads them itself.)
        return _parse(lines, header)
    if tuple(header[:len(_EDGE_COLS)]) == _EDGE_COLS:
        return _parse_edges(lines, header)      # the pre-aggregation table (see `_EDGE_COLS`)
    # `block_events.tsv` — a nucleotide run's own record, keyed by ancestral interval. Every
    # resolution writes its genealogy in one format, so that file is the only table left that looks
    # like this one and is not; naming it saves reading the columns to find out.
    hint = ("; this looks like a nucleotide run's block_events.tsv, whose rows are ancestral "
            "intervals rather than genome events. Read one with "
            "zombi2.genomes.nucleotide.read_nucleotide_genomes. That run's genome_events.tsv is "
            "the genealogy, in this format, and is what this reader wants"
            if "source" in header and "family" not in header else "")
    for cols in () if hint else (_COLS, _EDGE_COLS):
        # a genuinely short header: every column it has is one of a log's, so it is a log this
        # version no longer reads rather than a different table. Named against the current columns
        # first, since that is what a truncated file of this schema is missing.
        if set(header) < set(cols):
            hint = (f"; it is missing {', '.join(c for c in cols if c not in header)}. Re-run "
                    f"'zombi2 genomes' to write a log this version reads")
            break
    raise ValueError(f"unexpected genome-event columns {header}; expected {list(_COLS)}{hint}")


#: the kinds the pre-aggregation log recorded about gene identity and descent. A wider table may also
#: carry the ancestry-**neutral** rearrangements; they end no gene lineage, so they are not events here.
_GENEALOGY = frozenset({"origination", "duplication", "loss", "transfer", "speciation"})


def _copies(cell: str, lineno: int, col: str) -> list[tuple[int, int]]:
    """One ``parents`` / ``children`` cell as ``[(species, copy), …]`` — the branch each participant
    lived on and its id, which is everything the columns this log dropped used to say."""
    if not cell:
        return []
    out = []
    for token in cell.split(_PACK):
        species, sep, gene = token.rpartition("_")
        if not sep:
            raise ValueError(f"genome event log line {lineno}: {col} entry {token!r} is not a copy "
                             f"of the form n<species>_g<copy>")
        out.append((node_from_label(species), gene_from_label(gene)))
    return out


def _need(what: list, n: int, lineno: int, kind: str, col: str) -> None:
    if len(what) != n:
        raise ValueError(f"genome event log line {lineno}: a {kind} has {n} {col}, got {len(what)}")


def _parse(lines: list[str], header: list[str]) -> list[GeneEdge]:
    """Read the rows by column **name**, so a table carrying more than the canonical columns parses
    unchanged and only the genealogy rows come back, and expand each into its gene-tree edges — the
    exact inverse of `events_from_edges()`, down to the order the edges come back in."""
    at = {c: i for i, c in enumerate(header)}
    events: list[GeneEdge] = []
    for lineno, raw in enumerate(lines[1:], 2):
        if not raw:                                     # tolerate a trailing blank line
            continue
        cells = raw.split("\t")
        if len(cells) != len(header):
            raise ValueError(f"genome event log line {lineno}: expected {len(header)} columns, "
                             f"got {len(cells)}")
        kind = cells[at["kind"]]
        if kind not in _WRITTEN_KINDS:
            continue
        t, fam = float(cells[at["time"]]), int(cells[at["family"]])
        parents = _copies(cells[at["parents"]], lineno, "parents")
        children = _copies(cells[at["children"]], lineno, "children")
        if kind == "origination":
            _need(parents, 0, lineno, kind, "parents")
            _need(children, 1, lineno, kind, "children")
            (lineage, copy), = children
            events.append(GeneEdge(t, kind, lineage, fam, copy))
        elif kind == "loss":
            _need(parents, 1, lineno, kind, "parents")
            _need(children, 0, lineno, kind, "children")
            (lineage, copy), = parents
            events.append(GeneEdge(t, kind, lineage, fam, copy))
        elif kind in ("duplication", "speciation"):
            _need(parents, 1, lineno, kind, "parents")
            _need(children, 2, lineno, kind, "children")
            (_, source), = parents
            events.extend(GeneEdge(t, kind, lineage, fam, copy, parent=source)
                          for lineage, copy in children)
        else:                                           # a transfer, additive or replacing
            _need(parents, 2 if kind == _REPLACING else 1, lineno, kind, "parents")
            _need(children, 2, lineno, kind, "children")
            (donor, cont), (recipient, arrived) = children   # donor's copy leads (events_from_edges)
            (_, source), *rest = parents
            replaced = rest[0][1] if rest else None
            if replaced is not None:
                # the copy the arriving one overwrote: it dies here, and its death is an edge of the
                # gene tree like any other, so it comes back as the `loss` the file no longer spends
                # a row on. Written first, which is where the engines record it.
                events.append(GeneEdge(t, "loss", recipient, fam, replaced))
            events.append(GeneEdge(t, "transfer", donor, fam, cont, parent=source, donor=donor,
                                replaced=replaced))
            events.append(GeneEdge(t, "transfer", recipient, fam, arrived, parent=source,
                                recipient=recipient, donor=donor, replaced=replaced))
    return events


def _parse_edges(lines: list[str], header: list[str]) -> list[GeneEdge]:
    """The pre-aggregation table (`_EDGE_COLS`): every row is already one edge, so this is a
    field-by-field read. Kept so a log written before the aggregation still replays."""
    at = {c: i for i, c in enumerate(header)}
    events: list[GeneEdge] = []
    for lineno, raw in enumerate(lines[1:], 2):
        if not raw:
            continue
        cells = raw.split("\t")
        if len(cells) != len(header):
            raise ValueError(f"genome event log line {lineno}: expected {len(header)} columns, "
                             f"got {len(cells)}")
        if cells[at["kind"]] not in _GENEALOGY:
            continue
        get = lambda c: cells[at[c]]                    # noqa: E731
        events.append(GeneEdge(
            time=float(get("time")), kind=get("kind"), lineage=node_from_label(get("lineage")),
            family=int(get("family")), copy=gene_from_label(get("copy")),
            parent=gene_from_label(get("parent")) if get("parent") else None,
            recipient=node_from_label(get("recipient")) if get("recipient") else None,
            donor=node_from_label(get("donor")) if get("donor") else None))
    return events
