"""Deprecated: the name `zombi2.tools.compact` had in 0.52.0, kept for one release.

`bundle_run` and `unbundle_run` run `compact_run` and `expand_run`, and warn that they will go.
"""
from __future__ import annotations

import warnings

from .compact import Converted, compact_run, expand_run


def _deprecated(old: str, new: str) -> None:
    warnings.warn(f"zombi2.tools.bundle.{old} is deprecated and will be removed in the next "
                  f"release; use zombi2.tools.compact.{new}", FutureWarning, stacklevel=3)


def bundle_run(directory, *, remove: bool = False) -> list[Converted]:
    """Deprecated: `zombi2.tools.compact.compact_run`."""
    _deprecated("bundle_run", "compact_run")
    return compact_run(directory, remove=remove)


def unbundle_run(directory, *, flat: bool = False, remove: bool = False) -> list[Converted]:
    """Deprecated: `zombi2.tools.compact.expand_run`."""
    _deprecated("unbundle_run", "expand_run")
    return expand_run(directory, flat=flat, remove=remove)


__all__ = ["Converted", "bundle_run", "unbundle_run"]
