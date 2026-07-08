"""How a horizontal transfer resolves: recipient choice and additive vs replacement.

Rates say *how often* a transfer happens (the rate model); this object says *what a
transfer does* once it fires. Kept separate so the two concerns compose freely.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TransferModel:
    """Transfer mechanics.

    Parameters
    ----------
    replacement:
        Probability in [0, 1] that a transfer is a **replacement** rather than
        **additive**. A replacement adds the incoming copy and simultaneously removes one
        pre-existing copy of the same family in the recipient (a net-zero swap). It is only
        possible when the recipient already has a copy; otherwise the transfer is additive.
    distance_decay:
        Controls recipient choice by phylogenetic distance. ``None`` (default) picks a
        recipient uniformly among co-existing lineages. A positive value ``λ`` weights each
        candidate by ``exp(-λ · d)``, where ``d = 2·(t − t_MRCA)`` is the patristic distance
        between the donor and candidate at the transfer time ``t``. Larger ``λ`` = more
        local transfers; distant transfers are damped but never forbidden.
    allow_self:
        If ``True``, the donor lineage itself is an eligible recipient. A self-transfer
        creates a second copy in the same genome — mechanically a duplication — which lets
        one drop explicit duplications and run a transfer/loss-only model.
    receptivity:
        Per-branch **absorption** weights: a ``{branch_name: weight}`` map that biases which
        lineage *receives* a transfer (the counterpart to transfer *emission*, which is a rate and
        is scaled per branch via :class:`~zombi2.genomes.rates.BranchRates`). ``None`` (default)
        leaves recipient choice unweighted. Otherwise each candidate's selection weight is
        multiplied by its receptivity (branches not listed default to ``1.0``), composing with
        ``distance_decay`` if that is also set. A branch with weight ``0`` never receives; a branch
        with weight ``2`` is twice as likely as an unlisted one at the same distance. Reachable from
        the CLI via ``--branch-rates FILE`` (see
        :func:`~zombi2.genomes.read_rates.read_branch_rates`).
    """

    replacement: float = 0.0
    distance_decay: float | None = None
    allow_self: bool = False
    receptivity: dict | None = None

    def __post_init__(self):
        if not (0.0 <= self.replacement <= 1.0):
            raise ValueError(f"replacement must be in [0, 1], got {self.replacement}")
        if self.distance_decay is not None and self.distance_decay < 0:
            raise ValueError(f"distance_decay must be >= 0, got {self.distance_decay}")
        if self.receptivity is not None:
            recept = {str(k): float(v) for k, v in self.receptivity.items()}
            if any(v < 0 for v in recept.values()):
                raise ValueError("receptivity weights must be >= 0")
            self.receptivity = recept or None
