"""The run summary — what came out, as JSON, beside the log that says what went in.

A run already records what it was *asked* for: the ``.log`` carries the command line, every resolved
parameter and each input's SHA-256. What it did not record is what *happened*, and that turned out to
be the thing everyone reconstructed by hand, each in their own way and one of them wrongly:

- the event log holds **one row per gene-tree edge**, so a duplication, a transfer and a speciation
  each write two rows. Counting by row inflates them exactly 2×. It is documented, with the right
  ``awk``, and it is still the mistake a reader makes first — one tester got the deduplicated numbers
  and another reported an asymmetry that was really in his own test;
- ``gene_trees/`` holds a file pair per family that *ever existed* while the summary line counts the
  families that **survived**, so "96 gene families" sat beside 213 files and neither number explained
  the other;
- the family-size cap is invisible. When it binds it discards duplications and arriving transfers, so
  realised rates fall below the declared ones, and nothing in the output said it had happened.

So this is a second file, not a longer log: the log is provenance and the summary is outcome, and a
script that wants the second should not have to skip the first. JSON because the people who asked for
it were writing collectors — the log's ``key<TAB>value``-with-comments shape cost two of them a
parser apiece.

Every count here is **deduplicated**: one number per event, not per edge.
"""

from __future__ import annotations

import json
import pathlib
import statistics


def _stats(values) -> dict:
    """``min``/``mean``/``max`` for a list of numbers, or nulls when there is nothing to describe."""
    vals = list(values)
    if not vals:
        return {"min": None, "mean": None, "max": None}
    return {"min": min(vals), "mean": round(statistics.fmean(vals), 6), "max": max(vals)}


def write_summary(path, payload: dict) -> None:
    """Write ``payload`` as indented JSON with a trailing newline.

    Key order is the dict's own, which is the order the level built it in — the readable order, and a
    stable one, so two runs of the same seed produce byte-identical files."""
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_summary(path) -> dict:
    """The counterpart, so a caller reading a run back does not hand-roll it."""
    return json.loads(pathlib.Path(path).read_text(encoding="utf-8"))


__all__ = ["write_summary", "read_summary"]
