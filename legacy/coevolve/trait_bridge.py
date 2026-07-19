"""Bridge: a state-target coupling on a trait's OU optimum.

The grammar's *rate*-target couplings compile onto the genome rate engine (a
:class:`~zombi2.coevolve.rate_bridge.CouplingModifier`); a **state**-target coupling on
``traits.optimum`` instead conditions the trait's **OU walk** — at each point the driver sets *where*
the trait is pulled (the optimum-shift JUMP of ``docs/design/coevolve-grammar.md`` §2.1). This module
walks that coupled OU trait, realizing the ``genes:traits`` edge (a modifier gene sets the optimum)
and the trait side of the ``traits↔genomes`` joint.

It runs on the OU / co-integrator path, **not** the genome rate engine, so it is independent of
``zombi2.genomes.rates``.
"""

from __future__ import annotations

import math

from zombi2.coevolve.grammar import DriverSignal, Response
from zombi2.tree import Tree


def _ou_step(x: float, theta: float, dt: float, alpha: float, sigma2: float, rng) -> float:
    """One exact Ornstein–Uhlenbeck transition over ``dt`` toward ``theta`` (Brownian if
    ``alpha <= 0``). Draws one normal per non-degenerate step — the same draw order the
    per-edge implementations used, so the walk is reproducible."""
    if alpha <= 0.0:                                            # Brownian-motion limit (no pull)
        std = (sigma2 * dt) ** 0.5
        return x + (rng.normal(0.0, std) if std > 0.0 else 0.0)
    e = math.exp(-alpha * dt)
    mean = theta + (x - theta) * e
    var = (sigma2 / (2.0 * alpha)) * (1.0 - e * e)
    return float(rng.normal(mean, var ** 0.5)) if var > 0.0 else float(mean)


def walk_optimum_coupled_trait(tree: Tree, driver: DriverSignal, response: Response, *,
                               alpha: float, sigma2: float, x0: float, rng) -> dict:
    """Walk an OU trait down ``tree`` whose optimum is set, at every point, by ``response`` applied
    to the ``driver``'s local value. Returns ``{node: trait_value}`` for every node.

    Along each branch the walk is sub-segmented at the driver's interior change points
    (:meth:`DriverSignal.refresh_times`), so the optimum switches **exactly** when the driver does —
    the exact stochastic-map coupling. This realizes a state-target coupling on ``traits.optimum``:
    ``response.state_offset(driver_value)`` is the OU optimum θ on each segment.
    """
    alpha = float(alpha)
    sigma2 = float(sigma2)
    node_trait = {tree.root: float(x0)}
    for node in tree.nodes_preorder():
        if node.parent is None:
            continue
        b0, b1 = node.parent.time, node.time
        x = node_trait[node.parent]
        # Cut the branch at the driver's interior change points on this lineage, then take one exact
        # OU step per constant-optimum segment (driver constant ⇒ θ constant within a segment).
        cuts = sorted(t for (t, br) in driver.refresh_times(b0, b1)
                      if br == node.name and b0 < t < b1)
        bounds = [b0, *cuts, b1]
        for s0, s1 in zip(bounds[:-1], bounds[1:]):
            dt = s1 - s0
            if dt <= 0.0:
                continue
            theta = response.state_offset(driver.value(node.name, s0))
            x = _ou_step(x, theta, dt, alpha, sigma2, rng)
        node_trait[node] = x
    return node_trait
