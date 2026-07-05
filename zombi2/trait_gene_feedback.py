"""Trait<->gene feedback — the joint (bidirectional) traits<->genes model.

The joint of the two traits/genes edges, with the *same* coupled objects on both sides (as ClaSSE
couples one trait both ways):

* ``traits:genes`` — the trait sets a panel of families' **retention**: a present family is lost at
  ``base_loss·exp(-effect_loss·x)``, so a high trait value keeps the panel, a low one purges it;
* ``genes:traits`` — the panel sets the trait's **OU optimum**: the more of the panel a lineage
  carries, the higher its optimum ``theta(m) = theta_low + (theta_high - theta_low)·m/N``.

Together they close a loop: carrying the panel pulls the trait up, and a high trait keeps the panel —
so a lineage tends to fall into one of two self-reinforcing regimes (panel-rich & high-trait, or
panel-poor & low-trait), and the tips end up with the trait and the panel **correlated** even though
neither was imposed. Neither arrow points into S, so this is an **overlay** on a given tree.

Because the trait (a continuous OU diffusion) and the panel (a jump process) each depend on the
other's *current* value, they cannot be simulated one-then-the-other; they are integrated **together**
along each branch in ``steps`` small pieces (piecewise-constant coupling within a piece — the same
sub-segmenting :mod:`zombi2.trait_coupling` uses for a continuous trait covariate). Setting
``effect_loss = 0`` recovers pure ``genes:traits``; setting ``theta_high = theta_low`` recovers pure
``traits:genes`` — the joint model contains both single edges as limits.

    import zombi2 as z
    tree = z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=60, age=6, seed=1)
    m = z.TraitGeneFeedback(n_families=20, effect_loss=1.5, theta_low=-3, theta_high=3)
    res = z.simulate_trait_gene_feedback(tree, m, seed=2)
    res.trait_gene_correlation()     # tip trait vs panel occupancy — the emergent association
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .profiles import ProfileMatrix, _natkey
from .tree import Tree

_MAX_EXPONENT = 40.0                      # clamp on the loss exponent (matches trait_coupling)


class TraitGeneFeedback:
    """A continuous trait and a family panel coupled **both** ways (the traits<->genes joint model).

    Parameters
    ----------
    n_families:
        Panel size ``N`` (the coupled families; all respond to the trait with weight 1).
    effect_loss:
        Trait->gene strength: a present family's loss rate is ``base_loss·exp(-effect_loss·x)``
        (``0`` = the panel ignores the trait = pure ``genes:traits``).
    base_loss:
        Loss rate at trait value ``0``.
    gain:
        Trait-blind per-family gain rate (an absent family re-enters at this constant rate).
    theta_low, theta_high:
        Gene->trait coupling: the trait's OU optimum runs linearly from ``theta_low`` (empty panel)
        to ``theta_high`` (full panel). ``theta_high == theta_low`` = the trait ignores the genes =
        pure ``traits:genes``.
    alpha:
        OU mean-reversion strength (``0`` = Brownian, no pull toward the optimum).
    sigma2:
        Trait diffusion rate (variance per unit time).
    x0:
        Root trait value (default: midway between the two optima).
    root_fraction:
        Fraction of the panel present at the root (the first ``round(root_fraction·N)`` families).
    steps:
        Number of joint integration sub-steps per branch (higher = finer coupling; default 20).
    """

    def __init__(self, *, n_families: int = 20, effect_loss: float = 1.5, base_loss: float = 1.0,
                 gain: float = 1.0, theta_low: float = -3.0, theta_high: float = 3.0,
                 alpha: float = 1.0, sigma2: float = 0.5, x0: float | None = None,
                 root_fraction: float = 0.5, steps: int = 20):
        if n_families <= 0:
            raise ValueError("n_families must be positive")
        for name, val in (("effect_loss", effect_loss), ("base_loss", base_loss), ("gain", gain),
                          ("alpha", alpha), ("sigma2", sigma2)):
            if val < 0:
                raise ValueError(f"{name} must be >= 0")
        if not (0.0 <= root_fraction <= 1.0):
            raise ValueError(f"root_fraction must be in [0, 1], got {root_fraction}")
        if steps < 1:
            raise ValueError("steps must be >= 1")
        self.n_families = int(n_families)
        self.effect_loss = float(effect_loss)
        self.base_loss = float(base_loss)
        self.gain = float(gain)
        self.theta_low = float(theta_low)
        self.theta_high = float(theta_high)
        self.alpha = float(alpha)
        self.sigma2 = float(sigma2)
        self.x0 = float(x0) if x0 is not None else 0.5 * (self.theta_low + self.theta_high)
        self.root_fraction = float(root_fraction)
        self.steps = int(steps)
        self.panel_ids = [f"F{i}" for i in range(self.n_families)]

    def optimum(self, m: int) -> float:
        """The trait's OU optimum when ``m`` of the ``N`` panel families are present."""
        return self.theta_low + (self.theta_high - self.theta_low) * (m / self.n_families)

    def _ou_step(self, x: float, theta: float, dt: float, rng) -> float:
        """Exact OU transition over ``dt`` toward ``theta`` (Brownian if ``alpha == 0``)."""
        if self.alpha <= 0.0:
            std = (self.sigma2 * dt) ** 0.5
            return x + (rng.normal(0.0, std) if std > 0.0 else 0.0)
        e = math.exp(-self.alpha * dt)
        mean = theta + (x - theta) * e
        var = (self.sigma2 / (2.0 * self.alpha)) * (1.0 - e * e)
        return float(rng.normal(mean, var ** 0.5)) if var > 0.0 else float(mean)

    def __repr__(self):
        return (f"TraitGeneFeedback(n_families={self.n_families}, effect_loss={self.effect_loss:g}, "
                f"theta_low={self.theta_low:g}, theta_high={self.theta_high:g})")


@dataclass
class TraitGeneFeedbackResult:
    """The outcome of :func:`simulate_trait_gene_feedback`.

    ``node_trait`` maps every node to its trait value; ``node_presence`` to the boolean panel-presence
    array; ``profiles`` is the ``N × extant-species`` presence matrix.
    """

    tree: Tree
    model: TraitGeneFeedback
    node_trait: dict
    node_presence: dict
    profiles: ProfileMatrix

    def trait_values(self) -> dict:
        """Trait value at each extant tip."""
        return {n: self.node_trait[n] for n in self.tree.extant_leaves()}

    def panel_occupancy(self) -> dict:
        """Fraction of the panel present at each extant tip."""
        return {n: float(self.node_presence[n].mean()) for n in self.tree.extant_leaves()}

    def trait_gene_correlation(self) -> float:
        """Pearson correlation between tip trait value and tip panel occupancy — the emergent
        trait<->gene association the feedback writes into the data (``nan`` if either is constant)."""
        tips = self.tree.extant_leaves()
        if len(tips) < 2:
            return float("nan")
        x = np.array([self.node_trait[t] for t in tips])
        f = np.array([self.node_presence[t].mean() for t in tips])
        if x.std() == 0 or f.std() == 0:
            return float("nan")
        return float(np.corrcoef(x, f)[0, 1])


def _evolve_branch(x: float, present: np.ndarray, dt_total: float, model: TraitGeneFeedback,
                   rng) -> tuple[float, np.ndarray]:
    """Jointly integrate the trait (OU) and the panel (2-state per family) along one branch."""
    steps = model.steps
    dt = dt_total / steps if dt_total > 0 else 0.0
    if dt <= 0.0:
        return x, present
    n = model.n_families
    for _ in range(steps):
        m = int(present.sum())
        x = model._ou_step(x, model.optimum(m), dt, rng)          # gene -> trait (optimum from m)
        loss = model.base_loss * math.exp(max(-_MAX_EXPONENT, min(_MAX_EXPONENT,
                                                                   -model.effect_loss * x)))
        gain = model.gain                                          # trait -> gene (loss set by x)
        s = loss + gain
        if s > 0.0:
            decay = 1.0 - math.exp(-s * dt)                        # exact 2-state transition over dt
            p_lose = loss / s * decay
            p_gain = gain / s * decay
            u = rng.random(n)
            flip = (present & (u < p_lose)) | (~present & (u < p_gain))
            present = present ^ flip
    return x, present


def _panel_profile(node_presence, tips, model: TraitGeneFeedback) -> ProfileMatrix:
    """Full-panel families × extant-species presence matrix (every panel row kept)."""
    tips = sorted(tips, key=lambda n: _natkey(n.name))
    rows, cols, data = [], [], []
    for j, tip in enumerate(tips):
        pres = node_presence[tip]
        for i in np.nonzero(pres)[0]:
            rows.append(int(i)); cols.append(j); data.append(1)
    return ProfileMatrix(families=list(model.panel_ids), species=[t.name for t in tips],
                         coo=(rows, cols, data))


def simulate_trait_gene_feedback(
    tree: Tree,
    model: TraitGeneFeedback,
    *,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
) -> TraitGeneFeedbackResult:
    """Simulate the traits<->genes JOINT model down a **given** ``tree`` (the trait and the coupled
    panel evolve together, each modulating the other; see :class:`TraitGeneFeedback`).

    ``tree`` is any :class:`~zombi2.tree.Tree`. Returns a :class:`TraitGeneFeedbackResult` holding the
    trait and panel presence at every node, plus the extant profile matrix.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    n_root = int(round(model.root_fraction * model.n_families))
    root_present = np.zeros(model.n_families, dtype=bool)
    root_present[:n_root] = True

    node_trait = {tree.root: model.x0}
    node_presence = {tree.root: root_present}
    for node in tree.nodes_preorder():
        if node.parent is None:
            continue
        x, present = _evolve_branch(node_trait[node.parent],
                                    node_presence[node.parent].copy(),
                                    node.branch_length(), model, rng)
        node_trait[node] = x
        node_presence[node] = present

    profiles = _panel_profile(node_presence, tree.extant_leaves(), model)
    return TraitGeneFeedbackResult(tree=tree, model=model, node_trait=node_trait,
                                   node_presence=node_presence, profiles=profiles)
