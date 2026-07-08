"""ZOMBI2 experimental staging namespace.

Models exposed here are shipped so you can *use and iterate on them*, but they
have **not** yet cleared the core bar: they may lack full validation, their APIs
may change without notice, and each will eventually be either **promoted** into
the stable core (once validated and reviewed) or **removed**.

Import experimental models explicitly from ``zombi2.experimental`` -- nothing here
is re-exported into the top-level ``zombi2`` namespace, and nothing here is wired
into the ``zombi2`` command line or the model catalog. That separation is what
keeps the core surface small and stable.

See the model lifecycle in ``docs/contributing/model-lifecycle.md`` for how a
model moves ``idea -> experimental -> core``.
"""
from __future__ import annotations

import warnings

__all__ = ["ExperimentalWarning", "warn_experimental"]


class ExperimentalWarning(UserWarning):
    """A feature lives in :mod:`zombi2.experimental`: usable, but not yet
    validated or reviewed for the core, with an API that may change."""


# names already warned about, so a user sees the notice once per process rather
# than on every construction
_warned: set[str] = set()


def warn_experimental(name: str) -> None:
    """Emit an :class:`ExperimentalWarning` for ``name`` once per process.

    Experimental models call this from their constructor, so the first time a
    user builds one they are told it is not yet core-stable::

        class MyNewModel:
            def __init__(self, ...):
                warn_experimental("MyNewModel")
                ...
    """
    if name in _warned:
        return
    _warned.add(name)
    warnings.warn(
        f"{name} is experimental (zombi2.experimental): it is not yet validated "
        f"or reviewed for the core, and its API may change. It will be promoted "
        f"to the core or removed.",
        ExperimentalWarning,
        stacklevel=2,
    )


# ---------------------------------------------------------------------------
# Adding a model here: drop it in as a module under this package and re-export it,
#
#     from .my_new_model import MyNewModel
#     __all__.append("MyNewModel")
#
# and have its constructor call ``warn_experimental("MyNewModel")``.
#
# On promotion to the core a model LEAVES this package: it moves to a core module,
# is re-exported from ``zombi2``, drops its ``warn_experimental`` call, and gains a
# CLI surface + a catalog page. See docs/contributing/model-lifecycle.md.
#
# (Intra-genome gene conversion was promoted to the core in this way: it is now
# ``SharedRates(conversion=...)`` + ``zombi2.ConversionModel`` with a ``--conversion`` CLI flag.)
# ---------------------------------------------------------------------------

# Imported after warn_experimental is defined (the model imports it). ESM2Critic imports torch/esm
# lazily, so this line does NOT pull the optional zombi2[selection] dependencies.
from zombi2.experimental.selection import (  # noqa: E402
    Critic, ESM2Critic, FixedProfileCritic, PLMSelection,
)

__all__ += ["Critic", "ESM2Critic", "FixedProfileCritic", "PLMSelection"]

from zombi2.experimental.codon_selection import (  # noqa: E402
    CodonSelection, calibrate_beta, translate,
)

__all__ += ["CodonSelection", "calibrate_beta", "translate"]

from zombi2.experimental.genome_selection import CDS, GenomeSelection, read_cds_gff  # noqa: E402

__all__ += ["CDS", "GenomeSelection", "read_cds_gff"]

from zombi2.experimental.realism import frechet_esm_distance  # noqa: E402

__all__ += ["frechet_esm_distance"]
