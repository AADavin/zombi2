"""Where a result's files go inside the directory it was handed.

Most outputs are one file, and they land in the directory. A few are one file per gene family, per
block, or per node — a hundred families is two hundred Newick files, and a real genome times a real
tree is thousands of FASTAs — and those get a subdirectory of their own, or the handful of tables
they sit beside is buried under them.

Every ``Result.write`` groups them the same way, so a run written from Python has the layout the
manual describes and the one a ``zombi2`` command produces. ``flat=True`` is the escape hatch, for a
tool that wants one directory and nothing else; it is what ``--flat`` passes through.
"""

from __future__ import annotations

import pathlib
import shutil


def grouped_dir(base: pathlib.Path, name: str, flat: bool) -> pathlib.Path:
    """The directory a many-files-per-run output belongs in: ``base/name``, or ``base`` under
    ``flat``. Created on the way out, so the caller can write into it immediately.

    Safe to call repeatedly, and it is: a streamed run resolves its directory once per family, and
    two different outputs can share one (assembled genomes and the initial genome both land in
    ``genomes/``). Emptying the directory is `fresh_dirs`' job, once per write, precisely because
    this one cannot know whether it is being called for the first time."""
    d = base if flat else base / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def fresh_dirs(base: pathlib.Path, names, flat: bool) -> None:
    """Empty the ``names`` subdirectories of ``base`` — **once, before a write fills them**.

    These directories hold one file per *unit*: per gene family, per block, per node. The units are
    numbered by the run that made them, so writing a second run into the same place interleaves two
    sets — family 4 of the old run sits beside family 4 of the new one where the ids collide, and
    stays there alone where they do not. Nothing in the directory says which run a file came from.

    That was not hypothetical. A lecturer building a class practical re-ran the genome level at a
    higher loss rate, got a run with **zero** surviving families — announced clearly on the terminal —
    and ``zombi2 tools treedist`` then read a leftover ``gene_tree_fam4_extant.nwk`` and printed
    ``rf 0``, byte-identical to the previous run's answer. Separately, a student comparing transfer
    rates in one directory counted the leftovers and reported that gene trees were being written for
    families the run said had died out. Two people, two tasks, one cause, and neither saw an error:
    the symptom was a plausible number, not a crash.

    So a run's directory describes that run and nothing else. Called from the top of a ``write`` (and
    once when a streamed run opens its sink) rather than from `grouped_dir`, because the same
    directory is legitimately filled by more than one pass and clearing on each would leave only the
    last. Under ``flat`` this does nothing: there the caller shares one directory with every other
    output and every other level, so nothing in it can safely be called "ours" to remove."""
    if flat:
        return
    for name in names:
        d = base / name
        if not d.is_dir():
            continue
        for stale in d.iterdir():
            shutil.rmtree(stale) if stale.is_dir() else stale.unlink()


__all__ = ["grouped_dir", "fresh_dirs"]
