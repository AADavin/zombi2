"""How an intra-genome gene conversion resolves: donor directionality.

The rate says *how often* a conversion happens (``PerCopyRates(conversion=...)``); this object
says *what a conversion does* once it fires — specifically how the donor (template) copy is
chosen. Kept separate so the two concerns compose freely, exactly like
:class:`~zombi2.genomes.transfers.TransferModel`.

Gene conversion is the intra-genome analogue of horizontal transfer: within one genome a copy of
a family overwrites ("converts") another copy of the **same family**. It is **non-reciprocal**
(the donor is unchanged, the recipient is overwritten) and a **replacement, not a duplication** —
copy number does not change. Repeated, it homogenises a family's copies (**concerted evolution**)
and pulls their reconstructed within-family coalescences toward the present.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ConversionModel:
    """Intra-genome gene-conversion mechanics: donor directionality.

    The rate is set by ``PerCopyRates(conversion=...)``; this object says *what a conversion does*
    once it fires — how the donor (template) copy is chosen — kept separate from the rate so the
    two compose freely, exactly like :class:`~zombi2.genomes.transfers.TransferModel`.

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
        if not (0.0 <= self.bias <= 1.0):
            raise ValueError(f"bias must be in [0, 1], got {self.bias}")
