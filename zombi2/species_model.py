"""Species-tree models.

A model carries the process parameters (rates) and knows how to sample an internal-node
age under that process. Conditioning (number of tips, tree age) is supplied separately
to :func:`zombi2.simulate_species_tree`, keeping "what process" and "how many tips"
cleanly apart.
"""

from __future__ import annotations

import math

import numpy as np


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

    def extinction_prob(self, tau: float, tol: float = 1e-12) -> float:
        """``E(τ)`` — probability a lineage present ``τ`` before the present leaves no sampled
        descendant. Constant-rate, complete sampling (ρ=1), so this is the probability of
        extinction by the present (used to place ghost lineages; see ``add_ghost_lineages``).
        """
        lam, mu = self.birth, self.death
        if mu <= 0.0:  # Yule: no extinction, complete sampling -> no ghosts
            return 0.0
        r = lam - mu
        if abs(r) < tol * max(1.0, lam):  # critical (birth ≈ death)
            x = lam * tau
            return x / (1.0 + x)
        e = math.exp(r * tau)
        return mu * (e - 1.0) / (lam * e - mu)


class Yule(BirthDeath):
    """Pure-birth (Yule) process — a birth-death with no extinction."""

    def __init__(self, birth: float):
        super().__init__(birth, 0.0)


class EpisodicBirthDeath:
    """Episodic (skyline) birth-death: piecewise-constant rates through time.

    Speciation and extinction rates are constant within each *epoch* and change at fixed
    ages (times before the present) — the model behind mass-extinction and shifting-regime
    scenarios. The reconstructed tree of a (time-varying) birth-death process is still a
    coalescent point process, so ``simulate_species_tree`` samples i.i.d. internal-node
    ages from this model's CDF and assembles a ranked tree exactly as for the constant-rate
    case. The CDF has no closed form here, so it is built numerically and inverted.

    Parameters
    ----------
    birth, death:
        Length-``K`` lists of speciation (λ) and extinction (μ) rates, one per epoch,
        ordered from the **present** backward (``birth[0]`` applies to the most recent
        epoch).
    shifts:
        The ``K-1`` epoch boundaries, as strictly increasing **ages** before the present.
        Epoch ``i`` covers ages ``(shifts[i-1], shifts[i])`` (with ``shifts[-1]=0`` and
        ``shifts[K]=∞`` implied).
    sampling_fraction:
        Probability ``ρ ∈ (0, 1]`` that an extant species is sampled (incomplete extant
        sampling). ``1.0`` = complete sampling.
    grid:
        Number of grid points for the numerical CDF (accuracy vs speed).

    Notes
    -----
    This models episodic *diversification* and incomplete *extant* sampling, both of which
    keep the tree ultrametric. Serial sampling *through time* (dated tips / fossils) is a
    separate, forward-simulation feature and is not modelled here.
    """

    def __init__(self, birth, death, shifts, *, sampling_fraction: float = 1.0, grid: int = 8000):
        self.birth = [float(x) for x in birth]
        self.death = [float(x) for x in death]
        self.shifts = [float(x) for x in shifts]
        self.rho = float(sampling_fraction)
        self.grid = int(grid)
        self._cache_A = None
        self._ages = None
        self._cdf = None

    def validate(self) -> None:
        k = len(self.birth)
        if len(self.death) != k:
            raise ValueError("birth and death must have the same length")
        if len(self.shifts) != k - 1:
            raise ValueError(f"need len(shifts) == len(birth) - 1, got {len(self.shifts)} and {k}")
        if k < 1:
            raise ValueError("need at least one epoch")
        if any(b < 0 for b in self.birth) or any(d < 0 for d in self.death):
            raise ValueError("rates must be >= 0")
        if not any(b > 0 for b in self.birth):
            raise ValueError("at least one epoch must have birth > 0")
        if any(s <= 0 for s in self.shifts) or list(self.shifts) != sorted(set(self.shifts)):
            raise ValueError("shifts must be strictly increasing positive ages")
        if not (0.0 < self.rho <= 1.0):
            raise ValueError(f"sampling_fraction must be in (0, 1], got {self.rho}")

    @staticmethod
    def _cumtrapz(y, x):
        out = np.zeros_like(y)
        out[1:] = np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(x))
        return out

    def _rates(self, ages):
        idx = np.searchsorted(self.shifts, ages, side="right")
        return np.array(self.birth)[idx], np.array(self.death)[idx]

    def _prepare(self, A: float) -> None:
        ages = np.linspace(0.0, A, self.grid + 1)
        h = ages[1] - ages[0]
        lam_g, mu_g = self._rates(ages)
        lam_m, mu_m = self._rates(ages[:-1] + 0.5 * h)  # midpoints for RK4

        def dE(l, m, e):
            return m - (l + m) * e + l * e * e

        # E(a) = probability a lineage at age a has no sampled extant descendant
        E = np.empty_like(ages)
        E[0] = 1.0 - self.rho
        for i in range(len(ages) - 1):
            e = E[i]
            k1 = dE(lam_g[i], mu_g[i], e)
            k2 = dE(lam_m[i], mu_m[i], e + 0.5 * h * k1)
            k3 = dE(lam_m[i], mu_m[i], e + 0.5 * h * k2)
            k4 = dE(lam_g[i + 1], mu_g[i + 1], e + h * k3)
            E[i + 1] = e + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)

        P = 1.0 - E                                   # survival (>=1 sampled descendant)
        cum_r = self._cumtrapz(lam_g - mu_g, ages)    # ∫ (λ - μ)
        density = lam_g * P**2 * np.exp(-cum_r)        # reconstructed node-age density
        cdf = self._cumtrapz(density, ages)
        cdf /= cdf[-1]
        self._cache_A, self._ages, self._cdf = A, ages, cdf

    def sample_internal_age(self, u: float, A: float) -> float:
        if self._cache_A != A:
            self._prepare(A)
        return float(np.interp(u, self._cdf, self._ages))
