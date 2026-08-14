"""Chromosomes — the two logs every genome resolution with chromosomes writes.

Two different things happen to a chromosome, and they are two files because they are two kinds of
statement.

A `ChromosomeEvent` is an **edge of the chromosome genealogy**: a chromosome id is re-minted at every
event that reshapes it, so a chromosome lineage begins, splits, merges or ends, and the edge list *is*
the **chromosome network** — the genealogy between the species tree and the gene trees.

A **rearrangement** — an inversion, a transposition, a translocation — moves a segment around inside a
genome. It begins and ends no lineage (the chromosomes it touches keep their ids, the genes keep
theirs), so it is an edge of nothing and belongs in a log of its own: ``rearrangement_events.tsv``,
whose rows are *segments* and whose columns are coordinates.

Both are resolution-agnostic — the ids are chromosome ids whether a chromosome holds ordered gene
tokens or nucleotide blocks — so both live here, one home, imported by the ordered and the nucleotide
engines. Each engine keeps its own ``Inversion`` / ``Transposition`` / ``Translocation`` records,
because their coordinates are in different units (genes, base pairs); what they share is the table,
which is written here so the two resolutions cannot spell one differently.
"""

from __future__ import annotations

from dataclasses import dataclass

from .events import _name, node_from_label, node_label  # noqa: F401


@dataclass(frozen=True)
class ChromosomeEvent:
    """One edge of the **chromosome genealogy** — a chromosome lineage's birth, split, merge, or
    death, fired on species branch ``lineage`` at ``time``. ``parents`` → ``children`` are chromosome
    ids and the arity names the event: ``"initial"`` (``()`` → one child: a replicon the run *started*
    with) and ``"origination"`` (``()`` → one child: a de-novo replicon), both **roots**;
    ``"speciation"`` and ``"fission"`` (one parent → two children, a **bifurcation**), ``"fusion"``
    (two parents → one child, the **reticulation** — in-degree 2, what makes this a network and not a
    tree), ``"loss"`` (one parent → ``()``, a **leaf**). The edge list is the network's ground truth; a
    graph serialisation is derived from it (never eNewick — a multi-rooted, reticulating graph is not a
    tree).

    ``lineage`` is the branch the event *fired* on, which for a speciation is the branch that split —
    its children live on the daughters. The written file says where each chromosome lived by naming it
    ``n<species>_c<chromosome>`` (`chromosome_label()`), so it carries no ``lineage`` column at all."""

    time: float
    kind: str  # "initial" | "origination" | "speciation" | "fission" | "fusion" | "loss"
    lineage: int
    parents: tuple[int, ...]
    children: tuple[int, ...]


#: ``chromosome_events.tsv``. Every participant carries its own branch in its token, so this table
#: is the same four columns the species and gene logs are: when, what, what ended, what began.
_COLS = ("time", "kind", "parents", "children")

#: what separates the ids inside a ``parents`` / ``children`` cell — the same separator the gene log
#: packs its copies with, because a cell means the same thing in both.
_PACK = ";"


def chromosome_label(species: int | None, chromosome_id: int, names=None) -> str:
    """A chromosome **where it sits**: ``n<species>_c<chromosome>`` — the token every table naming a
    chromosome lineage uses. A chromosome id is re-minted at each event, so it is unique on its own and
    the species is redundant for joining; it is here for the same reason a gene copy carries one (see
    `copy_label`): a row is read by people, and ``c7`` alone cannot say whose chromosome it was."""
    return f"{_name(names, species)}_c{chromosome_id}"


def chromosome_from_label(cell: str) -> tuple[int, int]:
    """The inverse of `chromosome_label()`: ``n0_c1`` → ``(0, 1)``, the branch and the chromosome."""
    species, sep, chrom = cell.rpartition("_")
    if not sep or chrom[:1] != "c":
        raise ValueError(f"{cell!r} is not a chromosome of the form n<species>_c<chromosome>")
    return node_from_label(species), int(chrom[1:])


def _branches(chromosome_events, tree) -> dict[int, int]:
    """``{chromosome id: the species branch it lived on}``, read off the log itself.

    A chromosome id is minted once, as the child of exactly one edge, and lives its whole life on one
    branch — so the log already carries this. The one edge whose children are *not* on the branch it
    fired on is a **speciation**: the daughters' chromosomes live on the daughters, in the tree's own
    child order (which is the order the engines mint them in), so that is where the tree is needed."""
    where: dict[int, int] = {}
    for e in chromosome_events:
        if e.kind == "speciation":
            for daughter, child in zip(tree.nodes[e.lineage].children or (), e.children):
                where[child] = daughter
        else:
            for child in e.children:
                where[child] = e.lineage
    return where


def chromosome_events_tsv(chromosome_events: list[ChromosomeEvent], tree, names=None) -> str:
    """The chromosome network as TSV — one row per edge, the ids of a multi-ended side joined by
    ``;`` (a fusion has two parents, a fission two children). Both resolutions write this file, so
    the writer lives with the record it writes.

    ``tree`` is the species tree the run happened in, needed only to put a speciation's children on
    the daughters they were minted for; ``names`` is its ``labels()`` map, so a chromosome on a
    lineage that died is written ``e<id>`` like everywhere else."""
    where = _branches(chromosome_events, tree)

    def cell(ids) -> str:
        return _PACK.join(chromosome_label(where[c], c, names) for c in ids)

    rows = [f"{e.time}\t{e.kind}\t{cell(e.parents)}\t{cell(e.children)}" for e in chromosome_events]
    return "\n".join(["\t".join(_COLS), *rows]) + "\n"


#: ``rearrangement_events.tsv``. A rearrangement acts on a **segment**, and a segment is not an entity
#: with an id — nothing is born and nothing dies — so unlike every other log here the branch cannot
#: ride inside a token and stays a column of its own. ``start`` (not ``position``) because a
#: ``length`` travels with it: the two name a run together. ``dest_chromosome`` is set only by a
#: translocation (a transposition lands on the chromosome it left), ``dest_position`` only where the
#: record has one, and ``flipped`` only by the two that move a block, which may land inverted.
REARRANGEMENT_COLS = ("time", "kind", "lineage", "chromosome", "start", "length",
                      "dest_chromosome", "dest_position", "flipped", "cuts")
#: How a list of ancestral breakpoints is written into one cell: ``source:position``, semicolon
#: separated, empty where a record has none. The nucleotide resolution fills it — the partition has
#: to tell an indel breakpoint from an ordinary one, and these coordinates are the only ones in the
#: same frame as the block log. The ordered resolution has no ancestral coordinates to give and
#: leaves it blank.
CUTS_SEP, CUTS_PAIR = ";", ":"


def cuts_cell(cuts) -> str:
    """``((0, 40), (0, 60))`` → ``"0:40;0:60"``; empty for none."""
    return CUTS_SEP.join(f"{src}{CUTS_PAIR}{at}" for (src, at) in cuts)


def cuts_from_cell(cell: str) -> tuple:
    """The inverse of `cuts_cell()`."""
    if not cell:
        return ()
    return tuple(tuple(int(x) for x in pair.split(CUTS_PAIR)) for pair in cell.split(CUTS_SEP))


def _rearrangement_cells(r) -> tuple:
    """One record as ``(kind, chromosome, start, length, dest_chromosome, dest_position, flipped)``.

    Keyed on the record's own class name, which the two resolutions share — they hold the same three
    kinds in the same words, differing only in the unit their coordinates are in (genes / base pairs),
    and the unit is the engine's business, not the table's. A nucleotide translocation records no
    landing position (its blocks keep their source coordinates and the arc is placed by the engine),
    so that one cell is empty there — the one field the two resolutions do not both carry."""
    kind = type(r).__name__.lower()
    if kind == "inversion":
        return (kind, r.chromosome, r.start, r.length, None, None, None)
    if kind == "transposition":
        return (kind, r.chromosome, r.start, r.length, None, r.dest, int(r.flipped))
    if kind == "translocation":
        return (kind, r.source, r.start, r.length, r.dest,
                getattr(r, "dest_position", None), int(r.flipped))
    raise AssertionError(f"unhandled rearrangement {type(r).__name__}")


def rearrangement_events_tsv(rearrangements, names=None) -> str:
    """The rearrangement log as TSV — one row per segment moved, in the order it happened (the engines
    append as they fire, and a run is one forward pass). ``names`` is the run's ``labels()`` map, so a
    rearrangement on a lineage that died is written ``e<id>``.

    Empty cells where a kind has no such field; the columns are the same for every row, which is the
    point of the file existing. These rows used to sit in the genealogy table with its nine columns
    blank, which made that table mostly empty and this record hard to read out of it."""
    rows = []
    for r in rearrangements:
        kind, *rest = _rearrangement_cells(r)
        rows.append("\t".join([str(r.time), kind, _name(names, r.lineage),
                               *("" if c is None else str(c) for c in rest),
                               cuts_cell(getattr(r, "cuts", ()))]))
    return "\n".join(["\t".join(REARRANGEMENT_COLS), *rows]) + "\n"


__all__ = ["ChromosomeEvent", "chromosome_events_tsv", "chromosome_label",
           "chromosome_from_label", "rearrangement_events_tsv", "REARRANGEMENT_COLS",
           "cuts_cell", "cuts_from_cell"]
