"""Algorithm 1 — backward reconstructed birth-death, conditioned on N extant tips.

Under a constant-rate birth-death process, the reconstructed tree (extant lineages
only) is a *coalescent point process*: the internal node times are i.i.d. from a single
distribution, and the ranked topology is obtained by uniform coalescence
(Hartmann, Wong & Stadler 2010; Stadler 2011). We therefore

1. draw the internal node **ages** (time before present) by inverse-CDF sampling, and
2. assemble a random ranked tree from those times.

Node ages are drawn from the reconstructed-process CDF (complete sampling, ``ρ = 1``),
with ``r = λ − μ`` and tree age ``A``::

    g(a) = (1 − e^{−r a}) / (λ − μ e^{−r a})          # unnormalised
    F(a) = g(a) / g(A),        a ∈ (0, A)

Its density is highest near the present (``a → 0``) — the classic "pull of the present".
Inverting ``F(a) = u``::

    K = u · g(A)
    a = −(1/r) · ln[ (1 − λ K) / (1 − μ K) ]

with closed-form Yule (μ = 0) and critical (λ = μ) special cases.

A recommended additional validation (not runnable here) is a two-sample comparison
against TreeSim's ``sim.bd.taxa`` on identical (λ, μ, N, A).
"""

from __future__ import annotations

import math

import numpy as np

from .species_model import SpeciesTreeModel
from .tree import Tree, TreeNode


class SpeciesTreeSimulator:
    """Stateless sampler for reconstructed birth-death trees (survivors only)."""

    def simulate(self, model: SpeciesTreeModel, rng: np.random.Generator) -> Tree:
        model.validate()
        lam, mu, N = model.birth, model.death, model.n_tips
        A = float(model.age)
        total_age = A

        if model.age_type == "crown":
            # Root fixed at time 0; sample the other N-2 internal nodes.
            ages = [self._sample_age(rng.random(), lam, mu, A) for _ in range(N - 2)]
            internal_times = sorted([total_age - a for a in ages] + [0.0])
        else:  # "stem": all N-1 internal nodes are i.i.d.
            ages = [self._sample_age(rng.random(), lam, mu, A) for _ in range(N - 1)]
            internal_times = sorted(total_age - a for a in ages)

        return self._assemble(N, total_age, internal_times, rng)

    # --- inverse-CDF age sampler ------------------------------------------
    @staticmethod
    def _sample_age(u: float, lam: float, mu: float, A: float, tol: float = 1e-12) -> float:
        """Draw one internal-node age in (0, A) from the reconstructed-process CDF."""
        r = lam - mu
        if mu < tol:  # Yule (pure birth): F(a) = (1 - e^{-λa}) / (1 - e^{-λA})
            return -math.log1p(-u * (1.0 - math.exp(-lam * A))) / lam
        if abs(r) < tol * max(1.0, lam):  # critical, λ ≈ μ
            kp = u * (lam * A) / (1.0 + lam * A)
            return kp / (lam * (1.0 - kp))
        # general case
        e_a = math.exp(-r * A)
        g_a = (1.0 - e_a) / (lam - mu * e_a)
        k = u * g_a
        return -math.log((1.0 - lam * k) / (1.0 - mu * k)) / r

    # --- ranked-tree assembly ---------------------------------------------
    @staticmethod
    def _assemble(
        N: int, total_age: float, internal_times: list[float], rng: np.random.Generator
    ) -> Tree:
        """Coalesce N leaves at the given internal times (most recent first)."""
        active: list[TreeNode] = [
            TreeNode(name=f"n{i + 1}", time=total_age, is_extant=True) for i in range(N)
        ]
        for ft in sorted(internal_times, reverse=True):
            i, j = (int(x) for x in rng.choice(len(active), size=2, replace=False))
            parent = TreeNode(name="", time=ft)
            parent.add_child(active[i])
            parent.add_child(active[j])
            active = [x for k, x in enumerate(active) if k not in (i, j)]
            active.append(parent)

        root = active[0]
        tree = Tree(root, total_age)

        # deterministic internal-node names (preorder), root named "root"
        counter = 1
        for node in tree.nodes_preorder():
            if node.is_leaf():
                continue
            if node is root:
                node.name = "root"
            else:
                node.name = f"i{counter}"
                counter += 1
        return tree


def simulate_species_tree(model: SpeciesTreeModel, rng: np.random.Generator) -> Tree:
    """Convenience wrapper around :class:`SpeciesTreeSimulator`."""
    return SpeciesTreeSimulator().simulate(model, rng)
