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
        is scaled per lineage via :class:`~zombi2.genomes.rates.LineageRates`). ``None`` (default)
        leaves recipient choice unweighted. Otherwise each candidate's selection weight is
        multiplied by its receptivity (branches not listed default to ``1.0``), composing with
        ``distance_decay`` if that is also set. A branch with weight ``0`` never receives; a branch
        with weight ``2`` is twice as likely as an unlisted one at the same distance. Reachable from
        the CLI via ``--branch-rates FILE`` (see
        :func:`~zombi2.genomes.read_rates.read_branch_rates`).
    pair:
        A :class:`PairModifier` biasing transfer by the **(donor, recipient) pair** — the mechanism
        behind transfer *highways*. ``None`` (default) leaves pair choice unbiased. Otherwise each
        candidate recipient's selection weight is additionally multiplied by the pair factor for
        ``(donor, recipient)``, composing with ``distance_decay`` and ``receptivity``. This is the
        *recipient-seam* counterpart of the emission-seam modifiers (per-lineage / per-family rates);
        see :class:`PairModifier` and ``docs/design/rate-modifiers.md``. Python API only for now, and
        it runs on the pure-Python engine.
    """

    replacement: float = 0.0
    distance_decay: float | None = None
    allow_self: bool = False
    receptivity: dict | None = None
    pair: "PairModifier | None" = None

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

    def bind(self, tree) -> None:
        """Resolve any donor→recipient pair modifier against the species tree (clade specs need it).
        Called once at the start of a simulation; a no-op when ``pair`` is ``None``."""
        if self.pair is not None:
            self.pair.bind(tree)


def _mrca(nodes):
    """The most recent common ancestor of ``nodes`` — their deepest (max-time) shared ancestor."""
    common = None
    for n in nodes:
        chain, x = set(), n
        while x is not None:
            chain.add(x)
            x = x.parent
        common = chain if common is None else (common & chain)
    if not common:
        raise ValueError("clade tips have no common ancestor (are they from this tree?)")
    return max(common, key=lambda x: x.time)


@dataclass
class PairModifier:
    """Donor→recipient multiplier on transfer — the recipient-seam counterpart of the emission-seam
    :class:`~zombi2.genomes.rates.Modifier`, and the mechanism behind transfer **highways**.

    Attach to a :class:`TransferModel` via ``pair=``. When a transfer fires, each candidate
    recipient's selection weight is multiplied by the factor for ``(donor, recipient)`` (default
    ``1.0``). Specify highways two combinable ways (overlapping specs multiply):

    * ``pairs`` — an explicit ``{(donor_name, recipient_name): factor}`` map of exact branch pairs.
    * ``blocks`` — a list of ``(from_clade, to_clade, factor)``: every donor in ``from_clade`` to
      every recipient in ``to_clade`` is scaled by ``factor``. A *clade* is a **node name** (str), or
      a **set/list of tip names** whose MRCA defines the clade (its whole subtree).

    Bound with the species tree at the start of a simulation (clade specs are resolved to subtrees
    then). Factors must be ``>= 0``; a ``0`` forbids that route. See ``docs/design/rate-modifiers.md``.
    """

    pairs: dict | None = None
    blocks: list | None = None

    def __post_init__(self):
        self._pairs = {(str(d), str(r)): float(f) for (d, r), f in (self.pairs or {}).items()}
        if any(f < 0 for f in self._pairs.values()):
            raise ValueError("pair factors must be >= 0")
        self._blocks = [(a, b, float(f)) for a, b, f in (self.blocks or [])]
        if any(f < 0 for _a, _b, f in self._blocks):
            raise ValueError("block factors must be >= 0")
        if self.pairs is None and self.blocks is None:
            raise ValueError("a PairModifier needs at least one of pairs or blocks")
        self._resolved: list = []  # (from_names:set, to_names:set, factor), filled by bind()

    def bind(self, tree) -> None:
        by_name = {n.name: n for n in tree.nodes_preorder()}
        self._resolved = [
            (self._subtree_names(self._clade(a, by_name)),
             self._subtree_names(self._clade(b, by_name)), f)
            for a, b, f in self._blocks
        ]

    @staticmethod
    def _subtree_names(node) -> set:
        names, stack = set(), [node]
        while stack:
            n = stack.pop()
            names.add(n.name)
            stack.extend(n.children)
        return names

    @staticmethod
    def _clade(spec, by_name):
        """Resolve a clade spec — a node name, or a set/list of tip names (their MRCA) — to a node."""
        if isinstance(spec, str):
            node = by_name.get(spec)
            if node is None:
                raise ValueError(f"unknown clade node name {spec!r}")
            return node
        tips = []
        for t in spec:
            node = by_name.get(str(t))
            if node is None:
                raise ValueError(f"unknown tip name {str(t)!r} in clade spec")
            tips.append(node)
        if not tips:
            raise ValueError("a clade tip-set must name at least one tip")
        return _mrca(tips)

    def factor(self, donor, recipient) -> float:
        """The multiplier for a transfer from ``donor`` to ``recipient`` (both tree nodes)."""
        f = self._pairs.get((donor.name, recipient.name), 1.0)
        if self._resolved:
            dn, rn = donor.name, recipient.name
            for from_names, to_names, bf in self._resolved:
                if dn in from_names and rn in to_names:
                    f *= bf
        return f
