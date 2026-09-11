"""The links a family genome run read from its own gene content, as a table: ``links.tsv``.

In a joint run a parameter can read the gene content the run is building. A run rate, a declared
family's own rate or a ``transfer_to`` written with ``scaled_by``, ``set_by`` or ``weighted_by`` on
``"genomes:<family>"``, ``"genomes:module:<group>"`` or ``"genomes:count"`` is one link. The run
records every link, and the table is written with its other files, so a method that looks for links
between gene families can be compared against it. A driver read from a trait or a file is not a link
here.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..params.connection import Driven
from ..params.mapping import Between, Table

#: the columns of ``links.tsv``, in order
COLUMNS = ("family", "target", "driver", "modifier", "mapping")

#: what the ``family`` column holds for a parameter of the run itself
RUN = "run"

#: the start of a live module driver, spelled as `zombi2.genomes.family` resolves it
_MODULE = "genomes:module:"

#: a declared family's own parameters that can carry a link, in the order its rows list them
_FAMILY_TARGETS = ("duplication", "transfer", "loss")


@dataclass(frozen=True)
class Link:
    """One row of ``links.tsv``: whose parameter it is (``family``, or ``"run"`` for the run's own),
    which parameter (``target``), what it reads (``driver``), the modifier that wrote the link
    (``modifier``), and ``mapping`` as text (see `mapping_text`)."""

    family: str
    target: str
    driver: str
    modifier: str
    mapping: str


def _number(x: object) -> str:
    """A table entry as text: the ``repr`` of the float, so ``1.0`` reads back as the same number."""
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return repr(float(x))
    return repr(x)


def mapping_text(mapping, driver: str, modules) -> str:
    """One link's mapping as text.

    A table is ``state=value`` pairs joined by ``;``, and a `Between` is ``donor>receiver=weight``
    pairs; both add ``default=`` when the default is not 1. A mapping read on a module is written at
    each of the module's completion levels, ``k/n=value``: a module of n families has exactly n + 1
    levels, so the text holds every number the run could have used. Any other mapping is written in
    its written form, which for a function is its name."""
    if isinstance(mapping, Table):
        parts = [f"{state}={_number(value)}" for state, value in mapping.per_state.items()]
        if mapping.default != 1.0:
            parts.append(f"default={_number(mapping.default)}")
        return ";".join(parts)
    if isinstance(mapping, Between):
        parts = [f"{a}>{b}={_number(w)}" for (a, b), w in mapping.per_pair.items()]
        if mapping.default != 1.0:
            parts.append(f"default={_number(mapping.default)}")
        return ";".join(parts)
    module = driver[len(_MODULE):] if driver.startswith(_MODULE) else None
    if module is not None and module in modules:
        n = len(modules[module])
        return ";".join(f"{k}/{n}={_number(mapping.multiplier(k / n))}" for k in range(n + 1))
    return repr(mapping)


def links_of(run_rates, run_transfer_to, declared, family_rates, family_transfer_to,
             modules) -> tuple[Link, ...]:
    """Every link a run reads from its own gene content, in a fixed order: the run's rates, the run's
    ``transfer_to``, then each declared family's rates and ``transfer_to``, in declaration order.

    ``run_rates`` is ``{target: Rate}``; ``family_rates`` is ``{target: {family index: Rate}}`` for
    the families whose own rate carries a modifier; ``family_transfer_to`` is
    ``{family index: rule}``; ``modules`` is ``{module name: family names}``."""
    rows: list[Link] = []

    def add(owner: str, target: str, modifiers) -> None:
        for m in modifiers:
            if isinstance(m, Driven) and isinstance(m.driver, str) and m.driver.startswith("genomes:"):
                rows.append(Link(owner, target, m.driver, m.verb,
                                 mapping_text(m.mapping, m.driver, modules)))

    for target, rate in run_rates.items():
        add(RUN, target, rate.modifiers)
    add(RUN, "transfer_to", (run_transfer_to,))
    for i, spec in enumerate(declared):
        for target in _FAMILY_TARGETS:
            rate = family_rates.get(target, {}).get(i)
            if rate is not None:
                add(spec.name, target, rate.modifiers)
        if i in family_transfer_to:
            add(spec.name, "transfer_to", (family_transfer_to[i],))
    return tuple(rows)


def links_tsv(links) -> str:
    """``links.tsv``: the header, then one row per link. A run with no link writes the header alone."""
    rows = ["\t".join(COLUMNS)]
    rows += ["\t".join((x.family, x.target, x.driver, x.modifier, x.mapping)) for x in links]
    return "\n".join(rows) + "\n"


def links_from_tsv(text: str) -> tuple[Link, ...]:
    """The links `links_tsv` wrote, read back."""
    header = "\t".join(COLUMNS)
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines or lines[0] != header:
        raise ValueError(f"links.tsv must start with the header {header!r}")
    return tuple(Link(*line.split("\t")) for line in lines[1:])
