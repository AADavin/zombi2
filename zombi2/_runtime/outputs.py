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


def grouped_dir(base: pathlib.Path, name: str, flat: bool) -> pathlib.Path:
    """The directory a many-files-per-run output belongs in: ``base/name``, or ``base`` under
    ``flat``. Created on the way out, so the caller can write into it immediately."""
    d = base if flat else base / name
    d.mkdir(parents=True, exist_ok=True)
    return d


__all__ = ["grouped_dir"]
