"""The cross-level parameter grammar (SPEC §5): a parameter is what you set, a driver is what it
reads, a verb says what the reading does to it.

Shared by every level, so it lives in one place. A rate is written from its **scope**, and the verbs
chain onto it::

    from zombi2.params import PerCopy, PerLineage, LogNormal, Recipients

    birth       = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})
    loss        = PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))
    loss        = PerCopy().set_by("habitat.tsv", {'aquatic': 1.0, 'terrestrial': 0.25})
    transfer_to = Recipients().weighted_by("competence.tsv", {'competent': 3.0})

- **scopes** — *per what?* ``PerCopy`` · ``PerLineage`` · ``PerSite`` · ``PerChromosome`` · ``Global``
- **verbs** — *what does the reading do?* ``scaled_by`` multiplies, ``set_by`` replaces,
  ``weighted_by`` compares the candidates of a choice. Two drivers are written so often that each
  has a verb of its own, which is the only spelling for it: ``varying_among`` for `Random`,
  ``changing_at`` for the run's clock.
- **drivers** — *what is read?* `Random` · `TotalDiversity` · `Time` · `Clade` · a filename or a
  level name, which is the whole of conditioning and joining (SPEC §2)
- **laws** — *and then what happens to it?* a bare distribution is drawn and held; `Drift` carries
  it down the tree and perturbs it at each split
- **mappings** — *how does the driver's value become a number?* a dict is a `Table`, a function a
  `Curve`, a number a `Scalar`, and `Between` reads both ends of a transfer. Named only when the
  plain form will not do, since each is coerced from what you would have written anyway
- **choice rules** — *who receives?* ``"uniform"``, `Distance` (closer relatives likelier),
  `Clades` (weight the donor's clade against the recipient's), or a `Recipients()` rule
- ``Extent(...)`` and ``Recipients()`` — the entry points of the other two parameter kinds, written
  only when a verb is chained onto them

Reach the plumbing as submodules: ``parameter`` (the `Rate` and `Extent` the whole expression
evaluates to), ``parse`` (``parse_rate("PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})")``, which is
how the CLI and a ``--params`` file take a parameter — SPEC §5, *one written form*), ``connection``
(what a verb records), ``evaluate`` (what an engine calls), ``mapping``, ``distributions``.
"""

from .choice import Clades, Distance, Recipients
from .driver import Clade
from .distributions import Exponential, Fixed, Gamma, Geometric, LogNormal, Uniform
from .parameter import Extent
from .mapping import Between, Curve, Scalar, Schedule, Table
from .driver import Random, TotalDiversity
from .evaluate import UNITS
from .law import Drift
from .retired import RETIRED, name_message
from .scope import Global, PerChromosome, PerCopy, PerLineage, PerSite
from .driver import Time

#: The whole written vocabulary, from one place: **every name the text form may call is importable
#: from here**, so what you type in a ``--birth`` flag and what you import in Python are one surface
#: rather than two. The distributions and the mappings are re-exported for the same reason — an
#: argument you have to write should not need an import from a submodule nothing else asks you to
#: know about.
#:
#: Two entries are here and not in the text form's whitelist, and both for a stated reason. `Curve`
#: takes a function, which no flag can carry, so it is the one mapping that is Python-only (the
#: parser says so and points here). ``UNITS`` is data rather than a name you call — every unit the
#: grammar knows; an error quotes back the five a value may vary among (`evaluate.VARYING_UNITS`).
__all__ = ["UNITS", "Between", "Clade", "Clades", "Curve", "Distance", "Drift", "Exponential",
           "Extent", "Fixed", "Gamma", "Geometric", "Global", "LogNormal", "PerChromosome",
           "PerCopy", "PerLineage", "PerSite", "Random", "Recipients", "Scalar", "Schedule", "Table", "Time",
           "TotalDiversity", "Uniform"]


def __getattr__(name: str):
    """Answer a **retired** name with the sentence naming its replacement, rather than with
    "module has no attribute".

    The text form has carried that table since the grammar changed (`zombi2.params.retired`), and
    Python had nothing: ``from zombi2.params import ScaledBy`` said only that the name was absent,
    which is the one thing the reader already knew. This reads the same table, so a spelling
    retired once is answered the same way whether it was typed into a flag or into an import.

    It raises `ImportError` rather than `AttributeError` because ``from … import`` discards an
    ``AttributeError`` from here and substitutes its own generic "cannot import name", losing the
    sentence exactly where a port hits it first. The cost is that ``hasattr(zombi2.params,
    'ScaledBy')`` raises instead of answering ``False`` — worth it for a name nobody should be
    probing for.
    """
    if name in RETIRED:
        raise ImportError(name_message(name))
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
