"""Species-tree models.

A model carries the process parameters (rates) and knows how to sample an internal-node
age under that process. Conditioning (number of tips, tree age) is supplied separately
to :func:`zombi2.simulate_species_tree`, keeping "what process" and "how many tips"
cleanly apart.
"""

from __future__ import annotations

import math


class BirthDeath:
    """Constant-rate birth-death process (speciation ``birth``, extinction ``death``)."""

    def __init__(self, birth: float, death: float = 0.0):
        self.birth = float(birth)
        self.death = float(death)

    def validate(self) -> None:
        if self.birth <= 0:
            raise ValueError(f"birth rate must be > 0, got {self.birth}")
        if self.death < 0:
            raise ValueError(f"death rate must be >= 0, got {self.death}")

    def sample_internal_age(self, u: float, A: float, tol: float = 1e-12) -> float:
        """Draw one internal-node age in (0, A) from the reconstructed-process CDF.

        With ``r = birth - death`` and ``F(a) = g(a)/g(A)``,
        ``g(a) = (1 - e^{-r a}) / (birth - death e^{-r a})``; inverting ``F(a) = u``
        gives the closed forms below (Yule / critical / general).
        """
        lam, mu = self.birth, self.death
        r = lam - mu
        if mu < tol:  # Yule (pure birth)
            return -math.log1p(-u * (1.0 - math.exp(-lam * A))) / lam
        if abs(r) < tol * max(1.0, lam):  # critical, birth ≈ death
            kp = u * (lam * A) / (1.0 + lam * A)
            return kp / (lam * (1.0 - kp))
        e_a = math.exp(-r * A)
        g_a = (1.0 - e_a) / (lam - mu * e_a)
        k = u * g_a
        return -math.log((1.0 - lam * k) / (1.0 - mu * k)) / r


class Yule(BirthDeath):
    """Pure-birth (Yule) process — a birth-death with no extinction."""

    def __init__(self, birth: float):
        super().__init__(birth, 0.0)
