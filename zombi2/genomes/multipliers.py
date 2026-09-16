"""The drawn rate multipliers of a run, as tables: ``family_multipliers.tsv`` and
``lineage_multipliers.tsv``.

A rate varies along one of two axes, and they are separate questions that multiply.

``varying_among("families", law)`` draws a number for each **family** when the family is born, and
the family keeps it for its whole life: the family's rate is the run's rate multiplied by it. A
ribosomal protein family and a phage tail family then turn over at different speeds in the same
genome.

``varying_among("lineages", law)`` draws a number for each **species branch**, before the run
starts, and it multiplies the rate of every family on that branch: one branch sets the rate of a
whole genome. This is the axis the sequences level's clock rides (SPEC §7), read here on the
genome's own event rates. `draw_lineage_multipliers` is the one entry point.

Either way the run records the numbers and writes them with its other files, so a method that fits a
rate per family, or per branch, can be compared against the rates it was simulated with.

Two levels draw a per-family number, so both write ``family_multipliers.tsv``, in the same format
and under the same name: the genomes level, whose columns are its event rates, and the sequences
level, whose one column is the substitution rate (`zombi2.sequences.multipliers`). Each writes it
into its own run directory, because each drew its own numbers. Which columns a file has is what says
which level wrote it.
"""

from __future__ import annotations

import math

from ..params.evaluate import values_at_birth, values_at_split

#: the parameters a family resolution run draws a multiplier for, in the order the columns list them
FAMILY_TARGETS = ("duplication", "transfer", "loss")

#: the same at the ordered resolution, which also draws one for each rearrangement
ORDERED_TARGETS = (*FAMILY_TARGETS, "inversion", "transposition", "translocation")

#: the sequences level, which has one rate to draw among families: the substitution rate
SEQUENCE_TARGETS = ("substitution",)


def multipliers_of(tables, targets, own=None) -> dict[int, dict[str, "float | None"]]:
    """``{family: {target: multiplier}}``, one entry per family, from an engine's tables
    ``{target: {family: multiplier}}``.

    ``own`` is ``{target: family ids}``: the families whose own rate replaces the run's rate for that
    target, so no drawn number reaches them there, and their entry is ``None``."""
    own = own or {}
    families = sorted({f for target in targets for f in tables.get(target, {})})
    return {f: {target: (None if f in own.get(target, ()) else float(tables[target][f]))
                for target in targets}
            for f in families}


def _cell(value: "float | None") -> str:
    """A multiplier as text: the ``repr`` of the float, so it reads back as the same number, and an
    empty cell where the family's own rate took its place."""
    return "" if value is None else repr(float(value))


def multipliers_header(targets) -> str:
    """The header of ``family_multipliers.tsv``."""
    return "\t".join(("family", *targets))


def multipliers_row(family: int, row, targets) -> str:
    """One family's row of ``family_multipliers.tsv``."""
    return "\t".join((str(family), *(_cell(row[target]) for target in targets)))


def multipliers_tsv(multipliers, targets) -> str:
    """``family_multipliers.tsv``: the header, then one row per family in family order. A run whose
    rates do not vary among families writes the header alone."""
    rows = [multipliers_header(targets)]
    rows += [multipliers_row(f, multipliers[f], targets) for f in sorted(multipliers)]
    return "\n".join(rows) + "\n"


def multipliers_from_tsv(text: str) -> dict[int, dict[str, "float | None"]]:
    """The multipliers `multipliers_tsv` wrote, read back."""
    lines = [line for line in text.splitlines() if line.strip()]
    header = lines[0].split("\t") if lines else []
    targets = tuple(header[1:])
    known = (FAMILY_TARGETS, ORDERED_TARGETS, SEQUENCE_TARGETS)
    if not header or header[0] != "family" or targets not in known:
        raise ValueError("family_multipliers.tsv must start with one of the headers "
                         + ", ".join(repr(multipliers_header(t)) for t in known))
    out: dict[int, dict[str, "float | None"]] = {}
    for line in lines[1:]:
        family, *cells = line.split("\t")
        out[int(family)] = {target: (float(cell) if cell else None)
                            for target, cell in zip(targets, cells)}
    return out


# --- the other axis: a multiplier per species branch, shared by every family ----------------------

#: the parameters a genome run draws a per-lineage multiplier for, in the order the columns list
#: them. Origination is here and absent from `FAMILY_TARGETS` for the reason the engines give: a
#: family does not exist when origination is read, so it can carry no per-family draw — but the
#: lineage does, so it can carry a per-lineage one.
LINEAGE_TARGETS = ("duplication", "transfer", "loss", "origination")

#: the same at the ordered resolution, which counts the rearrangements and the chromosome events too.
#: Wider than `ORDERED_TARGETS`, its per-family counterpart, because a multiplier per branch scales
#: whatever the event acts on, while a per-family one has to reach the genes a segment covers and so
#: cannot apply to an event on a whole replicon.
ORDERED_LINEAGE_TARGETS = (*LINEAGE_TARGETS, "inversion", "transposition", "translocation",
                           "fission", "fusion", "chromosome_origination", "chromosome_loss")


def _preorder(tree) -> list[int]:
    """Species-tree node ids, parent before child — the order an inherited multiplier descends."""
    order: list[int] = []
    stack = [tree.root]
    while stack:
        i = stack.pop()
        order.append(i)
        stack.extend(tree.nodes[i].children)
    return order


def _objects(mods_by_target, targets) -> list:
    """Every distinct modifier **object** the rates carry, in written order, each appearing once.

    Identity, not equality: one `Random` written on duplication and on loss is one draw for the
    branch, so a fast lineage is fast at both; two separately built ones are two draws even with
    identical arguments (SPEC §5). That is the whole reason this is an identity walk rather than a
    set."""
    objects: list = []
    for target in targets:
        for m in mods_by_target.get(target, ()):
            if not any(m is o for o in objects):
                objects.append(m)
    return objects


def draw_lineage_multipliers(mods_by_target, tree, rng, targets=LINEAGE_TARGETS
                             ) -> "dict[str, dict[int, float]]":
    """One multiplier per species branch per rate: ``{target: {node id: multiplier}}``.

    ``mods_by_target`` is what each rate carries among lineages, from `Rate.carried_modifiers`. A
    `DRAWN` modifier gives every branch its own i.i.d. number; an `INHERITED` one (``Drift``) starts
    the root at `Inherited.initial` and perturbs it parent→child down the species tree, so close
    relatives keep similar rates. The two cannot be mixed on one rate (`check_one_memory`), but two
    rates may each take their own.

    The number belongs to the **branch**, not to a family: it multiplies the rate of every family on
    that branch, which is what lets one branch set the rate of a whole genome. On transfer the branch
    is the donor's — the rate says how often a lineage *donates* — and ``transfer_to`` still chooses
    the recipient.

    Every requested target gets a table, so a rate carrying nothing per lineage reads 1.0 on every
    branch and writes a column of ones rather than no column. ``{}`` when no rate carries anything,
    and then nothing is drawn and no randomness is consumed — so a run without a per-lineage draw is
    bit-identical to one from before this existed.

    Drawn over the **whole** species tree, extinct branches included, in one preorder walk: the
    table is the truth a method fitting a rate per branch is compared against, and a branch that
    left no descendants still had a rate. Factored out so the serial loop and the per-family engine
    draw them the same way — each hands in its own ``rng``, so the two give different though equally
    valid draws, as they do for everything else — and the per-family engine draws here in the
    parent, which is what keeps it worker-count invariant.
    """
    objects = tuple(_objects(mods_by_target, targets))
    if not objects:
        return {}
    values: dict[int, tuple[float, ...]] = {}
    for i in _preorder(tree):                                   # parent before child
        parent = tree.nodes[i].parent
        values[i] = (values_at_birth(objects, rng) if parent is None
                     else values_at_split(objects, values[parent], rng))
    where = {id(m): j for j, m in enumerate(objects)}
    return {target: {i: math.prod(v[where[id(m)]] for m in mods_by_target.get(target, ()))
                     for i, v in values.items()}
            for target in targets}


def lineage_multipliers_of(tables, labels, targets=LINEAGE_TARGETS
                           ) -> "dict[int, dict[str, float]]":
    """``{node: {target: multiplier}}``, one entry per species branch, from an engine's tables
    ``{target: {node: multiplier}}``.

    ``labels`` is the tree's node labels, and it decides **which** branches are listed: the table
    covers the tree the run was given, so a branch is in it whether or not a gene ever crossed it.
    Keyed by node id, as the rest of a result is; the file writes the label."""
    return {i: {target: float(tables[target][i]) for target in targets} for i in sorted(labels)}


def lineage_multipliers_header(targets=LINEAGE_TARGETS) -> str:
    """The header of ``lineage_multipliers.tsv``."""
    return "\t".join(("lineage", *targets))


def lineage_multipliers_row(label: str, row, targets=LINEAGE_TARGETS) -> str:
    """One branch's row of ``lineage_multipliers.tsv``."""
    return "\t".join((label, *(repr(float(row[target])) for target in targets)))


def lineage_multipliers_tsv(multipliers, labels, targets=LINEAGE_TARGETS) -> str:
    """``lineage_multipliers.tsv``: the header, then one row per species branch, in the order the
    tree numbers its nodes. The ``lineage`` column holds the branch's label, the name
    ``genomes.tsv`` and the species tree use, so a row can be matched to a branch without knowing
    the run's internal ids. A run whose rates do not vary among lineages writes the header alone."""
    rows = [lineage_multipliers_header(targets)]
    rows += [lineage_multipliers_row(labels[i], multipliers[i], targets)
             for i in sorted(multipliers)]
    return "\n".join(rows) + "\n"


def lineage_multipliers_from_tsv(text: str, labels=None) -> "dict[int | str, dict[str, float]]":
    """The multipliers `lineage_multipliers_tsv` wrote, read back.

    Keyed by node id when ``labels`` is given (the tree's ``labels()``, which the file's names are
    looked up in), and by the label itself otherwise — a table read on its own, without the tree it
    was written beside, still reads. With ``labels``, a name the tree does not hold raises: the two
    were written together, so a mismatch is a table read beside the wrong tree, and a row quietly
    keyed by its name would sit in the result looking like a node id."""
    lines = [line for line in text.splitlines() if line.strip()]
    header = lines[0].split("\t") if lines else []
    targets = tuple(header[1:])
    known = (LINEAGE_TARGETS, ORDERED_LINEAGE_TARGETS)
    if not header or header[0] != "lineage" or targets not in known:
        raise ValueError("lineage_multipliers.tsv must start with one of the headers "
                         + ", ".join(repr(lineage_multipliers_header(t)) for t in known))
    by_label = {label: i for i, label in (labels or {}).items()}
    out: dict[int | str, dict[str, float]] = {}
    for line in lines[1:]:
        label, *cells = line.split("\t")
        if labels is not None and label not in by_label:
            raise ValueError(f"lineage_multipliers.tsv has a row for branch {label!r}, which is "
                             f"not in the species tree read with it. The usual cause is a table and "
                             f"a tree taken from different runs.")
        out[by_label[label] if labels is not None else label] = {
            target: float(cell) for target, cell in zip(targets, cells)}
    return out
