"""Trait-evolution namespace (scikit-learn-style).

Re-exports the trait / biogeography public API so users can write::

    from zombi2.traits import OrnsteinUhlenbeck, DEC

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.traits.models` and
:mod:`zombi2.traits.biogeography`); it does not redefine anything.

Historical note: the trait *implementation* used to live in a flat module
(``zombi2/traits.py``). To make room for this scikit-learn-style namespace
without breaking ``from zombi2.traits import ...`` (which still works, because
every public trait name is re-exported below), the implementation was moved
verbatim to :mod:`zombi2.traits.models`. The private helper ``_expm`` is also
re-exported so existing callers keep working.
"""

from __future__ import annotations

from zombi2.traits.models import (
    BrownianMotion, OrnsteinUhlenbeck, MultivariateBrownian, MultivariateOU,
    MultiOptimumOU, ThresholdModel, EarlyBurst, Mk, CorrelatedBinary,
    CorrelatedBinaryK, HiddenStateMk, Cladogenesis, simulate_traits,
    replicate_traits, TraitResult,
    pagel_lambda, pagel_delta, pagel_kappa,
    _expm,  # noqa: F401 (re-exported for backward compatibility; see the module docstring)
)
from zombi2.traits.biogeography import DEC, simulate_biogeography

__all__ = [
    "BrownianMotion", "OrnsteinUhlenbeck", "MultivariateBrownian",
    "MultivariateOU", "MultiOptimumOU", "ThresholdModel", "EarlyBurst", "Mk",
    "CorrelatedBinary", "CorrelatedBinaryK", "HiddenStateMk", "simulate_traits",
    "replicate_traits",
    "TraitResult", "pagel_lambda", "pagel_delta", "pagel_kappa",
    "DEC", "simulate_biogeography", "Cladogenesis",
]
