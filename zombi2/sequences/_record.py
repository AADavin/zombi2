"""The sequence level's own history — every substitution, and every site gained or lost.

This is the one level that records nothing by default, and the reason is size. Three hundred sites on
a tree whose branches total thirty time units at rate 1.0 is nine thousand substitutions for **one**
family; a hundred families is close to a million rows, where a genome log for the same run holds a
few thousand. So ``record=True`` asks for it and nothing else turns it on.

**What it costs to record.** An ordinary run draws each branch's end from ``P(t) = exp(Q·bl)`` — one
matrix, one draw per site, and the path between the two ends is never simulated because nothing asks
what it was. A recorded run has to walk that path, so it runs a forward Gillespie over the sites: the
next event's waiting time, which site it lands on, and what it becomes. The two are the same process
and the same distribution at the branch's end; only one of them can say what happened in between.

**Naming a site.** A position in a lineage's own sequence shifts with every insertion above it, and
the alignment column index is stable only once the run has finished. So a site is named by an **id**,
minted where it arose and never reused, which is the frame the nucleotide level already uses for gene
events. An insertion row carries the id it follows, a deletion row the ids dropped; with those two a
reader rebuilds any lineage's column list at any moment, which is what turns an id back into a
position.

**Strand** is a column because the nucleotide level's blocks can sit reverse-complemented. Here it is
always ``+1``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SequenceEvent:
    """One row of the sequence level's log.

    ``kind`` is ``substitution``, ``insertion`` or ``deletion``. ``lineage`` is the species branch it
    happened on and ``gene`` the gene copy, so the pair names the same thing the alignment's
    ``n7_g12`` keys do. ``site`` is the site's id; ``after`` is the id it follows, on an insertion
    row, and ``-1`` for a run that landed before every column the family started with. ``from_state``
    and ``to_state`` are letters on a substitution and empty otherwise."""

    time: float
    kind: str
    lineage: str
    gene: int
    site: int
    strand: int = 1
    after: int = -1
    from_state: str = ""
    to_state: str = ""

    def row(self) -> tuple:
        return (f"{self.time:.6g}", self.kind, self.lineage, self.gene, self.site, self.strand,
                self.after, self.from_state, self.to_state)


HEADER = ("time", "event", "lineage", "gene", "site", "strand", "after", "from", "to")


def walk(states: np.ndarray, Q: np.ndarray, bl: float, rng, *, present=None):
    """A forward Gillespie down one branch: the end states, and every substitution on the way.

    Returns ``(end_states, [(position along the branch as a fraction, site, from, to)])``. ``bl`` is
    the branch's length in substitutions per site, so the waiting times are drawn in those units and
    the fraction is what maps a row back onto the branch.

    ``Q`` is normalised to one expected substitution per site per unit ``bl``, so the whole
    sequence's event rate is the sum of each site's own out-rate — which is what the draw below uses,
    updated at the one site each event touches rather than recomputed.

    ``present`` masks the columns this lineage actually carries. The engine evolves every column of
    the widened alphabet whether a lineage holds it or not, because that is what makes one alignment
    out of many; a *log* of a site the lineage does not have would be a claim about nothing."""
    out = states.copy()
    rates = -np.diagonal(Q)[out].astype(float)
    total = float(rates.sum())
    rows: list[tuple[float, int, int, int]] = []
    t = 0.0
    while total > 0.0:
        t += float(rng.exponential(1.0 / total))
        if t >= bl:
            break
        i = int(np.searchsorted(np.cumsum(rates), float(rng.random()) * total))
        i = min(i, out.size - 1)
        was = int(out[i])
        row = Q[was].copy()
        row[was] = 0.0
        new = int(rng.choice(row.size, p=row / row.sum()))
        out[i] = new
        if present is None or present[i]:
            rows.append((t / bl, i, was, new))
        total += float(-Q[new, new] + Q[was, was])
        rates[i] = float(-Q[new, new])
    return out, rows


@dataclass
class Recorder:
    """What a recorded run carries down the walk: where to put the rows, and how to name things.

    ``site_ids`` maps a column of the widened alignment to the id minted for it; ``None`` is a run
    with no indels, where a column *is* its id. ``present`` is the indel history's per-node mask, or
    ``None`` for the same reason."""

    events: list
    names: dict
    site_ids: "tuple[int, ...] | None" = None
    present: "dict | None" = None

    def mask(self, node):
        return None if self.present is None else self.present.get(id(node))

    def add(self, node, start: float, span: float, rows, alphabet: str) -> None:
        """Turn one **stretch** of a walk into rows, in the run's own frame.

        A position is in substitutions per site, and the time it maps to is that fraction of the
        stretch's own span. The stretch is a whole branch on an ordinary run, where the rate is
        constant down it. On a joint run it is one slice of one branch, or one segment of a slice —
        which is exactly the interval the rate *is* constant over, and why the mapping stays exact
        where the branch as a whole would not have been."""
        label = self.names.get(node.species, str(node.species))
        for frac, column, was, new in rows:
            site = column if self.site_ids is None else self.site_ids[column]
            self.events.append(SequenceEvent(
                start + frac * span, "substitution", label, node.copy, site,
                from_state=alphabet[was], to_state=alphabet[new]))

    def indel(self, kind: str, time: float, node, sites, after: int = -1) -> None:
        """One row per site gained or lost. A run of three is three rows, which keeps every row the
        same shape and makes the log greppable by site — the thing a reader actually wants."""
        label = self.names.get(node.species, str(node.species))
        for site in sites:
            self.events.append(SequenceEvent(time, kind, label, node.copy, site, after=after))


def recorder_for(record: bool, labels):
    """``(events, recorder)`` for a run, or ``(events, None)`` when it did not ask to record.

    Every engine that can record calls this rather than building a `Recorder` itself. Three of them
    built their own, and the third forgot: `record=True` on a joint run was accepted and logged
    nothing. One place to build it is one place to get it wrong."""
    events: list = []
    return events, (Recorder(events, labels) if record else None)
