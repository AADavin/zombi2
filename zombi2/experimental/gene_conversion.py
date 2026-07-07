"""Intra-genome (ectopic) gene conversion — **experimental**.

Within one genome, a gene copy of a family overwrites ("converts") another copy of the **same
family**. It is **non-reciprocal** (the donor copy is unchanged, the recipient is overwritten) and a
**replacement, not a duplication** — copy number does not change. It is the intra-genome analogue of
horizontal transfer: a transfer makes two *species'* copies coalesce more recently than the species
split; a conversion makes two *same-genome* copies coalesce more recently than their **duplication**.
Repeated, it homogenises a family's copies — **concerted evolution** — and correspondingly distorts
the reconstructed gene trees, pulling within-family coalescences toward the present.

Why it lives here (:mod:`zombi2.experimental`): the model runs and is validated (see
``tests/test_gene_conversion.py`` — a reduction, the structural invariants, and a homogenisation
oracle), but it has not yet been reviewed for the core, has no CLI surface or catalog page, and its
API may change. See ``docs/contributing/model-lifecycle.md``.

It is an unusual experimental model in one respect: intra-genome conversion is a new **event kind**,
not a self-contained rate model, so the *engine capability* to apply and reconstruct a conversion
(``EventType.CONVERSION``, :meth:`~zombi2.genomes.genome.UnorderedGenome.convert`, the simulator's
dispatch, the reconciliation) lives in the core engine, dormant — nothing fires a conversion unless a
rate model emits ``CONVERSION`` weights. This module is what activates it:

* :class:`GeneConversionRates` — a :class:`~zombi2.genomes.rates.SharedRates` that additionally emits
  the conversion events (so it is what turns the feature *on*);
* :class:`ConversionModel` — the conversion *mechanics* (donor directionality), passed to the
  simulator like a :class:`~zombi2.genomes.transfers.TransferModel`.

Usage::

    from zombi2.species import BirthDeath, simulate_species_tree
    from zombi2.genomes import simulate_genomes
    from zombi2.experimental import GeneConversionRates, ConversionModel

    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

    genomes = simulate_genomes(
        tree,
        GeneConversionRates(duplication=0.4, loss=0.1, conversion=1.0),
        conversions=ConversionModel(bias=0.0),   # bias=1.0 -> directional (toward the founder)
        initial_families=40, seed=42)

    trees = genomes.gene_trees()   # within-family copies coalesce toward their conversion times

On promotion to the core this module goes away: the rate becomes ``SharedRates(conversion=...)``, it
gains a ``--conversion`` CLI flag and a catalog page, drops the ``warn_experimental`` calls, and
``ConversionModel`` re-exports from ``zombi2``.
"""

from __future__ import annotations

from dataclasses import dataclass

from zombi2.experimental import warn_experimental
from zombi2.genomes.events import EventType
from zombi2.genomes.rates import EventWeight, SharedRates


@dataclass
class ConversionModel:
    """Intra-genome gene-conversion mechanics: donor directionality.

    The rate is set by :class:`GeneConversionRates`; this object says *what a conversion does* once
    it fires — specifically how the donor (template) copy is chosen — kept separate from the rate so
    the two compose freely, exactly like :class:`~zombi2.genomes.transfers.TransferModel`.

    Parameters
    ----------
    bias:
        Directionality of conversion, in ``[0, 1]``. The **recipient** (overwritten) copy is always
        chosen uniformly; ``bias`` controls the **donor**. ``0`` (default) draws the donor uniformly
        among the family's other copies — *unbiased*. ``1`` always picks the family's **founder /
        oldest lineage** (smallest ``Gene.origin_order``), homogenising the other copies *toward* it
        — directional. In between, with probability ``bias`` the oldest candidate donates, otherwise
        a uniform one does. ``bias`` is inert when a family holds exactly two copies (only one
        possible donor).
    """

    bias: float = 0.0

    def __post_init__(self):
        warn_experimental("ConversionModel")
        if not (0.0 <= self.bias <= 1.0):
            raise ValueError(f"bias must be in [0, 1], got {self.bias}")


class GeneConversionRates(SharedRates):
    """``SharedRates`` plus an intra-genome gene-**conversion** rate (**experimental**).

    Every gene family shares the same per-copy duplication / transfer / loss (and per-branch
    origination) rates, exactly as :class:`~zombi2.genomes.rates.SharedRates`, and additionally is
    subject to conversion at per-copy rate ``conversion``. A conversion needs both a donor and a
    recipient, so it fires only on a family holding **two or more copies**: a family with ``n``
    copies is converted at total rate ``conversion · n``.

    Being a ``SharedRates`` *subclass*, it is not the built-in model
    (``type(rates) is SharedRates`` is false), so it runs on the pure-Python engine automatically —
    the compiled Rust counts-only path is never asked to handle conversion (which is about gene-tree
    shape, not copy-number profiles).

    Pair with :class:`ConversionModel` (via ``simulate_genomes(..., conversions=...)``) to control
    donor directionality; without one, conversion is unbiased.
    """

    def __init__(self, duplication: float = 0.0, transfer: float = 0.0,
                 loss: float = 0.0, origination: float = 0.0, *, conversion: float = 0.0,
                 carrying_capacity: float | None = None):
        warn_experimental("GeneConversionRates")
        super().__init__(duplication, transfer, loss, origination,
                         carrying_capacity=carrying_capacity)
        if conversion < 0:
            raise ValueError(f"conversion rate must be >= 0, got {conversion}")
        self.conversion = float(conversion)

    def event_weights(self, genome, branch, time):
        out = super().event_weights(genome, branch, time)
        if self.conversion > 0 and genome.size() > 0:
            # per-copy rate, but a conversion needs a donor AND a recipient: only families with
            # >= 2 copies qualify, so emit one family-specific entry each (rate = conversion * cn).
            # Family-specific (never the family=None fast path), so the simulator receives the family.
            for family in genome.families():
                cn = genome.copy_number(family)
                if cn >= 2:
                    out.append(EventWeight(EventType.CONVERSION, family, self.conversion * cn))
        return out
