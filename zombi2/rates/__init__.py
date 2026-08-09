"""The cross-level rate grammar (SPEC §5): ``effective rate = scope(base) × modifiers``.

Shared by every level, so it lives in one place. Reach the pieces as submodules::

    from zombi2.rates import scope, modifiers
    birth = scope.Global(1.0)
    birth = 1.0 * modifiers.OnTime({0: 1.0, 3: 0.3})

- ``scope`` — *per what?* ``PerCopy`` · ``PerLineage`` · ``PerSite`` · ``Global``
- ``modifiers`` — *depends on what?* ``OnTime`` · ``OnTotalDiversity`` · an inherited value
- ``values`` / ``verbs`` — the grammar's two halves, written as a grid rather than a list:
  ``Drawn(per="family", dist=LogNormal(0.0, 0.5))`` is what a value *is*, ``ScaledBy(habitat, {...})`` what a
  parameter *does* with it. They build the same objects ``modifiers`` does, so both spellings run
  identically; the grid is what makes the next cell free rather than a new class.
- ``rate`` — the internal ``Rate`` plumbing (users never build a ``Rate`` directly)
- ``parse`` — ``parse_rate("1.0 * OnTime({0: 1.0, 3: 0.3})")``: the same expression written as text,
  which is how the CLI and a ``--params`` file take a rate (SPEC §5, *one written form*)
- ``distributions`` — value/length distributions
"""

from .distributions import Exponential, Fixed, Gamma, Geometric, LogNormal, Uniform
from .modifiers import OnTime, OnTotalDiversity
from .values import UNITS, Clade, Drawn, Inherited, Measured, Time
from .verbs import ScaledBy, SetBy, Weights

#: The whole written vocabulary, from one place. `modifiers.WRITABLE` + `values.WRITABLE` +
#: `verbs.WRITABLE` + `distributions.WRITABLE` is the same list the text form whitelists, so what you
#: can import here and what you can write in a `--birth` flag are one surface rather than two. The
#: distributions are re-exported because `Drawn` and `Inherited` now require one: an argument you
#: must write should not need an import from a submodule nothing else asks you to know about.
__all__ = ["UNITS", "Clade", "Drawn", "Exponential", "Fixed", "Gamma", "Geometric", "Inherited",
           "LogNormal", "Measured", "OnTime", "OnTotalDiversity", "ScaledBy", "SetBy", "Time",
           "Uniform", "Weights"]
