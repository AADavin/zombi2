"""Each family's drawn rate multipliers, as a table: ``family_multipliers.tsv``.

A rate written with ``varying_among("families", law)`` draws a number for each family when the family
is born, and the family keeps that number for its whole life: the family's rate is the run's rate
multiplied by it. The run records the numbers, and the table is written with its other files, so a
method that fits a rate per family can be compared against the rates each family was simulated with.

Two levels draw such a number, so both write this table, in the same format and under the same name:
the genomes level, whose columns are its event rates, and the sequences level, whose one column is
the substitution rate (`zombi2.sequences.multipliers`). Each writes it into its own run directory,
because each drew its own numbers. Which columns a file has is what says which level wrote it.
"""

from __future__ import annotations

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
