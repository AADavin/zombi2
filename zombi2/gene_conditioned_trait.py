"""Gene-conditioned trait evolution — the ``genes:traits`` edge (gene content shapes a trait).

The reverse of [`traits:genes`](trait-linked-genomes) (where a trait drives gene content): here a
lineage's **gene content conditions how a phenotypic trait evolves**. A binary *modifier* gene comes
and goes along the tree (gain/loss), and its presence sets the trait's **optimum** — so a lineage
that acquires the gene is pulled toward a new adaptive peak, and one that loses it drifts back.
"Gene presence enables a trait shift": the ``genes:traits`` arrow, where *genes drive the trait*.

Gene content does not feed back into diversification here, so the tree is an *input* and this is an
**overlay**. The trait is an Ornstein–Uhlenbeck process whose optimum is ``theta_present`` while the
modifier is present and ``theta_absent`` while it is absent; the switch happens exactly when the gene
is gained or lost along a branch (its history is simulated first, as a two-state Markov chain).

    import zombi2 as z
    tree = z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=60, age=5, seed=1)
    m = z.GeneConditionedTrait(gene_gain=0.6, gene_loss=0.6, theta_absent=0.0, theta_present=5.0,
                               alpha=2.0, sigma2=0.5)
    res = z.simulate_gene_conditioned_trait(tree, m, seed=2)
    res.trait_values()        # trait at the extant tips (carriers sit near theta_present)
    res.gene_presence()       # which tips carry the modifier (0/1)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .traits import Mk, TraitResult, simulate_traits
from .tree import Tree, TreeNode


class GeneConditionedTrait:
    """A binary modifier gene that switches a continuous trait's optimum (the ``genes:traits`` edge).

    Parameters
    ----------
    gene_gain, gene_loss:
        Rates at which the modifier gene is gained (absent → present) and lost (present → absent).
    root_gene:
        Whether the modifier is present at the root (default ``False``).
    theta_absent, theta_present:
        The trait's OU optimum while the modifier is absent / present.
    alpha:
        OU mean-reversion strength (how fast the trait tracks its current optimum; ``>= 0``,
        ``0`` = Brownian motion with no pull).
    sigma2:
        Trait diffusion rate (variance per unit time; ``>= 0``).
    x0:
        Root trait value (default: ``theta_absent`` if the modifier starts absent, else
        ``theta_present``).
    """

    def __init__(self, *, gene_gain: float = 0.5, gene_loss: float = 0.5, root_gene: bool = False,
                 theta_absent: float = 0.0, theta_present: float = 5.0,
                 alpha: float = 1.0, sigma2: float = 1.0, x0: float | None = None):
        if gene_gain < 0 or gene_loss < 0:
            raise ValueError("gene_gain and gene_loss must be >= 0")
        if alpha < 0 or sigma2 < 0:
            raise ValueError("alpha and sigma2 must be >= 0")
        if gene_gain == 0 and gene_loss == 0 and not root_gene:
            raise ValueError("the modifier can never appear (gene_gain=gene_loss=0, root_gene=False)")
        self.gene_gain = float(gene_gain)
        self.gene_loss = float(gene_loss)
        self.root_gene = bool(root_gene)
        self.theta_absent = float(theta_absent)
        self.theta_present = float(theta_present)
        self.alpha = float(alpha)
        self.sigma2 = float(sigma2)
        self.x0 = (float(x0) if x0 is not None
                   else (self.theta_present if root_gene else self.theta_absent))

    def _ou_step(self, x: float, present: bool, dt: float, rng) -> float:
        """Exact OU transition over ``dt`` toward the optimum selected by the gene state."""
        theta = self.theta_present if present else self.theta_absent
        if self.alpha <= 0.0:                                   # Brownian-motion limit (no pull)
            std = (self.sigma2 * dt) ** 0.5
            return x + (rng.normal(0.0, std) if std > 0.0 else 0.0)
        e = np.exp(-self.alpha * dt)
        mean = theta + (x - theta) * e
        var = (self.sigma2 / (2.0 * self.alpha)) * (1.0 - e * e)
        return float(rng.normal(mean, var ** 0.5)) if var > 0.0 else float(mean)

    def __repr__(self):
        return (f"GeneConditionedTrait(theta_absent={self.theta_absent:g}, "
                f"theta_present={self.theta_present:g}, alpha={self.alpha:g})")


@dataclass
class GeneConditionedTraitResult:
    """The outcome of :func:`simulate_gene_conditioned_trait`.

    ``node_trait`` maps every node to its (continuous) trait value; ``gene`` is the modifier's own
    :class:`~zombi2.traits.TraitResult` (its presence 0/1 and per-branch history)."""

    tree: Tree
    model: GeneConditionedTrait
    node_trait: dict
    gene: TraitResult

    def trait_values(self) -> dict:
        """Trait value at each extant tip."""
        return {n: self.node_trait[n] for n in self.tree.extant_leaves()}

    def gene_presence(self) -> dict:
        """Modifier presence (0/1) at each extant tip."""
        return {n: int(self.gene.node_values[n]) for n in self.tree.extant_leaves()}

    def to_tsv(self, nodes: str = "all") -> str:
        """A ``node`` × (modifier, trait) table."""
        if nodes == "extant":
            selected = self.tree.extant_leaves()
        elif nodes == "leaves":
            selected = self.tree.leaves()
        elif nodes == "all":
            selected = self.tree.nodes()
        else:
            raise ValueError("nodes must be 'extant', 'leaves', or 'all'")
        lines = ["node\tmodifier\ttrait"]
        for n in selected:
            lines.append(f"{n.name}\t{int(self.gene.node_values[n])}\t{self.node_trait[n]:.6g}")
        return "\n".join(lines) + "\n"

    def to_newick(self) -> str:
        """Newick with each node's trait value in a ``[&trait=…]`` comment."""
        def rec(node: TreeNode) -> str:
            s = (f"({','.join(rec(c) for c in node.children)}){node.name}"
                 if node.children else node.name)
            s += f"[&trait={self.node_trait[node]:.6g}]"
            if node.parent is not None:
                s += f":{node.branch_length():.10g}"
            return s
        return rec(self.tree.root) + ";"


def simulate_gene_conditioned_trait(
    tree: Tree,
    model: GeneConditionedTrait,
    *,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
) -> GeneConditionedTraitResult:
    """Evolve a trait down ``tree`` whose optimum is set by a modifier gene's presence
    (the ``genes:traits`` edge — gene content drives the trait).

    The modifier gene's presence/absence history is simulated first (a two-state Markov chain with
    gain/loss rates), then the trait is evolved as an OU process whose optimum switches exactly when
    the gene is gained or lost along each branch. Returns a :class:`GeneConditionedTraitResult`.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    # 1) the modifier gene: a two-state (0 absent, 1 present) Markov chain along the tree
    gene_model = Mk([[0.0, model.gene_gain], [model.gene_loss, 0.0]], states=[0, 1],
                    root=(1 if model.root_gene else 0))
    gene = simulate_traits(tree, gene_model, rng=rng)          # .history[node] = [(state, dur), ...]

    # 2) the trait: OU whose optimum follows the gene state, switching within a branch exactly at
    #    each gain/loss event (walk the gene's per-branch stochastic map).
    node_trait = {tree.root: model.x0}
    for node in tree.nodes_preorder():
        if node.parent is None:
            continue
        x = node_trait[node.parent]
        for state, dur in gene.history[node]:
            x = model._ou_step(x, present=(state == 1), dt=dur, rng=rng)
        node_trait[node] = x

    return GeneConditionedTraitResult(tree=tree, model=model, node_trait=node_trait, gene=gene)
