"""``reconparser`` — read the output of external reconciliation tools (ALE, AleRax).

Where :mod:`zombi2.tools.reconciliation` (ALElite) *computes* the reconciliation likelihood
of a ZOMBI2 gene family, this subpackage does the complementary **interop** job: it *reads*
what the established reconciliation programs already wrote. Point it at an ``ALEml`` /
``ALEml_undated`` result or an ``AleRax`` output directory and get back the reconciled trees,
the ML DTL rates, the log-likelihood, and the transfer / per-branch event tables as native
ZOMBI2 :class:`~zombi2.tree.Tree` objects and ``pandas`` DataFrames — the natural bridge for
comparing a real reconciliation against a ZOMBI2 simulation of the same system.

Trees are parsed with ZOMBI2's own :func:`~zombi2.tree.read_newick`, so a parsed reconciliation
plugs straight into the rest of the package (the :mod:`~zombi2.tools.reconciliation` scorer, the
tree-distance tools, …) with no foreign tree type in between. Reconciliation annotations that ALE
and AleRax bake into node names (``.T@donor->recipient``, ``.D@…``, support values) are preserved
verbatim as :attr:`~zombi2.tree.TreeNode.name`.

Public API::

    from zombi2.tools.reconparser import ALEParser, AleRaxRun, AleRaxFamily

- :class:`ALEParser` — classic ALE output (``.ucons_tree`` / ``.uTs`` / ``.uml_rec``; v0.4 and v1.0).
- :class:`AleRaxRun` / :class:`AleRaxFamily` — an AleRax run directory (v1.2+), fully lazy.

Optional dependency
-------------------
The parsers need ``pandas`` for their table outputs, which is **not** part of the base install.
Like the ``selection`` extra, it ships behind an opt-in extra::

    pip install 'zombi2[reconparser]'

Nothing here is re-exported into the top-level :mod:`zombi2` namespace, so plain
``import zombi2`` / ``import zombi2.tools`` never pulls in ``pandas`` — you only pay for it when
you import this subpackage explicitly.

Provenance
----------
Vendored from the standalone ``reconparser`` library (https://github.com/AADavin/reconparser,
commit ``a11e2b5``, MIT, same author) so ZOMBI2 users get it without a separate, not-yet-on-PyPI
install. The parser modules under :mod:`~zombi2.tools.reconparser.parsers` are kept close to
upstream to make future syncs easy.
"""

from __future__ import annotations

try:  # the parsers use pandas pervasively at import time
    import pandas as _pandas  # noqa: F401
except ImportError as _exc:  # pragma: no cover - only hit without the optional extra
    raise ImportError(
        "zombi2.tools.reconparser needs the optional 'reconparser' extra (pandas). "
        "Install it with:  pip install 'zombi2[reconparser]'"
    ) from _exc

from .parsers.ale import ALEParser
from .parsers.alerax import AleRaxFamily, AleRaxRun

__all__ = ["ALEParser", "AleRaxRun", "AleRaxFamily"]
