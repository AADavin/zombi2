"""``reconparser`` — read the output of external reconciliation tools (ALE, AleRax).

Where :mod:`zombi2.tools.reconciliation` (ALElite) *computes* the reconciliation likelihood
of a ZOMBI2 gene family, this subpackage does the complementary **interop** job: it *reads*
what the established reconciliation programs already wrote. Point it at an ``ALEml`` /
``ALEml_undated`` result or an ``AleRax`` output directory and get back the reconciled trees,
the ML DTL rates, the log-likelihood, and the transfer / per-branch event tables as
``ete3`` trees and ``pandas`` DataFrames — the natural bridge for comparing a real
reconciliation against a ZOMBI2 simulation of the same system.

Public API::

    from zombi2.tools.reconparser import ALEParser, AleRaxRun, AleRaxFamily

- :class:`ALEParser` — classic ALE output (``.ucons_tree`` / ``.uTs`` / ``.uml_rec``; v0.4 and v1.0).
- :class:`AleRaxRun` / :class:`AleRaxFamily` — an AleRax run directory (v1.2+), fully lazy.

Optional dependency
-------------------
The parsers need ``ete3`` and ``pandas``, which are **not** part of the base install (``ete3``
in particular is heavy). Like the ``selection`` extra, they ship behind an opt-in extra::

    pip install 'zombi2[reconparser]'

Nothing here is re-exported into the top-level :mod:`zombi2` namespace, so plain
``import zombi2`` / ``import zombi2.tools`` never pulls in ``ete3`` — you only pay for it when
you import this subpackage explicitly.

Provenance
----------
Vendored from the standalone ``reconparser`` library (https://github.com/AADavin/reconparser,
commit ``a11e2b5``, MIT, same author) so ZOMBI2 users get it without a separate, not-yet-on-PyPI
install. The parser modules under :mod:`~zombi2.tools.reconparser.parsers` are kept close to
upstream to make future syncs easy.
"""

from __future__ import annotations

try:  # the parsers use ete3 + pandas pervasively at import time
    import ete3 as _ete3  # noqa: F401
    import pandas as _pandas  # noqa: F401
except ImportError as _exc:  # pragma: no cover - only hit without the optional extra
    raise ImportError(
        "zombi2.tools.reconparser needs the optional 'reconparser' extra (ete3, pandas). "
        "Install it with:  pip install 'zombi2[reconparser]'"
    ) from _exc

from .parsers.ale import ALEParser
from .parsers.alerax import AleRaxFamily, AleRaxRun

__all__ = ["ALEParser", "AleRaxRun", "AleRaxFamily"]
