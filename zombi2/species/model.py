"""Species-tree models.

A model carries the process parameters (rates) and knows how to sample an internal-node
age under that process. Conditioning (number of tips, tree age) is supplied separately
to :func:`zombi2.simulate_species_tree`, keeping "what process" and "how many tips"
cleanly apart.
"""

from __future__ import annotations

import math

import numpy as np

from zombi2.species._caps import FOSSILIZATION, REMOVAL, GrowthEngine, SpeciesCaps


def _finite(name: str, value: float) -> float:
    """Reject NaN / infinite parameters up front.

    NaN silently passes every ``> 0`` / ``>= 0`` comparison, so a plain range check lets it
    through and the process then produces a tree with NaN branch lengths. Infinities produce
    degenerate (inf/zero) trees. Catch both here with a clear message.
    """
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number, got {value}")
    return value


def _normalize_mass_extinctions(mass_extinctions):
    """Normalize ``mass_extinctions`` into a sorted list of ``(age, fraction)`` float pairs.

    A *mass extinction* is an instantaneous, tree-wide pulse: at ``age`` before the present
    every lineage alive at that instant independently **dies** with probability ``fraction``
    (equivalently, survives with probability ``1 - fraction``). Several pulses may be given.
    Times are ages before the present, matching :class:`EpisodicBirthDeath` ``shifts``.
    Returns ``[]`` when none are supplied. This is a *forward-simulation* feature — the pulses
    kill real lineages in time (see :func:`~zombi2.simulate_species_tree` with
    ``direction="forward"``).
    """
    if not mass_extinctions:
        return []
    return sorted((float(age), float(frac)) for age, frac in mass_extinctions)


def _validate_mass_extinctions(mes) -> None:
    for a, frac in mes:
        _finite("mass-extinction age", a)
        _finite("mass-extinction fraction", frac)
    ages = [a for a, _ in mes]
    if any(a <= 0 for a in ages):
        raise ValueError("mass-extinction times must be positive ages before the present")
    if len(set(ages)) != len(ages):
        raise ValueError("mass-extinction times must be distinct")
    for _, frac in mes:
        if not (0.0 < frac <= 1.0):
            raise ValueError(f"mass-extinction fraction must be in (0, 1], got {frac}")


class BirthDeath:
    """Constant-rate birth-death process (speciation ``birth`` λ, extinction ``death`` μ).

    Optional (forward-simulation) extras: serial ``fossilization`` (ψ, dated fossils), extant
    ``sampling_fraction`` (ρ), ``removal`` (r) on sampling, and ``mass_extinctions`` (tree-wide
    survival pulses) — see :func:`~zombi2.simulate_species_tree` with ``direction="forward"``.
    ``fossilization``, ``removal != 1`` and ``mass_extinctions`` are forward-only; the backward
    reconstructed sampler assumes ρ=1 and no mass extinctions.

    ``mass_extinctions`` is a list of ``(age, fraction)`` pairs: at each ``age`` before the
    present, every lineage independently dies with probability ``fraction`` (a cataclysm that
    wipes out that fraction of the standing diversity in one instant).
    """

    caps = SpeciesCaps(
        GrowthEngine.THINNING, supports_backward=True, supports_ghosts=True,
        supports_n_tips=True, forward_only_features=(FOSSILIZATION, REMOVAL),
    )

    def __init__(self, birth: float, death: float = 0.0, *,
                 fossilization: float = 0.0, sampling_fraction: float = 1.0,
                 removal: float = 1.0, mass_extinctions=None):
        self.birth = float(birth)
        self.death = float(death)
        self.fossilization = float(fossilization)
        self.sampling_fraction = float(sampling_fraction)
        self.removal = float(removal)
        self.mass_extinctions = _normalize_mass_extinctions(mass_extinctions)

    def validate(self) -> None:
        _finite("birth rate", self.birth)
        _finite("death rate", self.death)
        _finite("fossilization rate", self.fossilization)
        _finite("sampling_fraction", self.sampling_fraction)
        _finite("removal", self.removal)
        if self.birth <= 0:
            raise ValueError(f"birth rate must be > 0, got {self.birth}")
        if self.death < 0:
            raise ValueError(f"death rate must be >= 0, got {self.death}")
        if self.fossilization < 0:
            raise ValueError(f"fossilization rate must be >= 0, got {self.fossilization}")
        if not (0.0 < self.sampling_fraction <= 1.0):
            raise ValueError(f"sampling_fraction must be in (0, 1], got {self.sampling_fraction}")
        if not (0.0 <= self.removal <= 1.0):
            raise ValueError(f"removal must be in [0, 1], got {self.removal}")
        _validate_mass_extinctions(self.mass_extinctions)

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
        descendant, under sampling fraction ρ (``sampling_fraction``); ``E(0)=1-ρ``. Used to
        place ghost lineages (see :func:`~zombi2.add_ghost_lineages`).
        """
        lam, mu, rho = self.birth, self.death, self.sampling_fraction
        if mu <= 0.0 and rho >= 1.0:  # Yule + complete sampling -> no ghosts
            return 0.0
        r = lam - mu
        if abs(r) < tol * max(1.0, lam):  # critical (birth ≈ death)
            return 1.0 - rho / (1.0 + rho * lam * tau)
        e = math.exp(-r * tau)
        return 1.0 - r / (lam - (lam - r / rho) * e)


class Yule(BirthDeath):
    """Pure-birth (Yule) process — a birth-death with no (background) extinction.

    ``mass_extinctions`` still applies: a Yule radiation punctuated by cataclysms.
    """

    def __init__(self, birth: float, *, mass_extinctions=None):
        super().__init__(birth, 0.0, mass_extinctions=mass_extinctions)


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

    caps = SpeciesCaps(
        GrowthEngine.THINNING, supports_backward=True, supports_ghosts=True,
        supports_n_tips=False, incomplete_sampling_backward=True,
        forward_only_features=(FOSSILIZATION, REMOVAL),
    )

    def __init__(self, birth, death, shifts, *, fossilization=None,
                 sampling_fraction: float = 1.0, removal: float = 1.0,
                 mass_extinctions=None, grid: int = 8000):
        self.birth = [float(x) for x in birth]
        self.death = [float(x) for x in death]
        self.shifts = [float(x) for x in shifts]
        # per-epoch serial fossilization ψ (forward-only; None -> no fossils)
        self.fossilization = ([0.0] * len(self.birth) if fossilization is None
                              else [float(x) for x in fossilization])
        self.rho = float(sampling_fraction)
        self.removal = float(removal)
        # instantaneous tree-wide survival pulses layered on the episodic rates (forward-only)
        self.mass_extinctions = _normalize_mass_extinctions(mass_extinctions)
        self.grid = int(grid)
        self._cache_A = None
        self._ages = None
        self._cdf = None

    def validate(self) -> None:
        for x in (*self.birth, *self.death, *self.fossilization, *self.shifts):
            _finite("episodic rate/shift", x)
        _finite("sampling_fraction", self.rho)
        _finite("removal", self.removal)
        k = len(self.birth)
        if len(self.death) != k:
            raise ValueError("birth and death must have the same length")
        if len(self.fossilization) != k:
            raise ValueError("fossilization must have the same length as birth")
        if len(self.shifts) != k - 1:
            raise ValueError(f"need len(shifts) == len(birth) - 1, got {len(self.shifts)} and {k}")
        if k < 1:
            raise ValueError("need at least one epoch")
        if any(b < 0 for b in self.birth) or any(d < 0 for d in self.death) \
                or any(f < 0 for f in self.fossilization):
            raise ValueError("rates must be >= 0")
        if not any(b > 0 for b in self.birth):
            raise ValueError("at least one epoch must have birth > 0")
        if not (0.0 <= self.removal <= 1.0):
            raise ValueError(f"removal must be in [0, 1], got {self.removal}")
        if any(s <= 0 for s in self.shifts) or list(self.shifts) != sorted(set(self.shifts)):
            raise ValueError("shifts must be strictly increasing positive ages")
        if not (0.0 < self.rho <= 1.0):
            raise ValueError(f"sampling_fraction must be in (0, 1], got {self.rho}")
        _validate_mass_extinctions(self.mass_extinctions)

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
        self._E = E  # kept for ghost-lineage grafting (P[no sampled descendant] by age)

    def sample_internal_age(self, u: float, A: float) -> float:
        if self._cache_A != A:
            self._prepare(A)
        return float(np.interp(u, self._cdf, self._ages))


class ClaDS:
    """Cladogenetic diversification with rate shifts (ClaDS; Maliet, Hartig & Morlon 2019).

    Every lineage carries its **own** speciation rate. At each speciation the two daughters
    each inherit the parent's rate times an independent lognormal jump,
    ``λ_child = λ_parent · exp(N(log α, σ²))`` — a small multiplicative shift per branch, so
    rates drift lineage-by-lineage down the tree. ``α`` is the trend (``α<1`` = the typical
    slow-down of speciation toward the present); ``σ`` is the jump spread. Extinction is tied
    to speciation through a constant **turnover** ``ε = μ/λ`` (so ``μ_child = ε·λ_child``);
    ``ε=0`` is ClaDS0 (pure birth with shifts), ``ε>0`` is ClaDS2.

    This is a **forward-only** model — per-lineage rates have no closed-form reconstructed CDF,
    so the backward sampler cannot draw i.i.d. node ages. Use
    :func:`~zombi2.simulate_species_tree` with ``direction="forward"`` (``age`` or ``n_tips``
    mode; ``mass_extinctions`` require ``age`` mode). ``sampling_fraction`` (ρ) and
    ``mass_extinctions`` overlay exactly as for :class:`BirthDeath`.
    """

    caps = SpeciesCaps(GrowthEngine.GILLESPIE, supports_n_tips=True)

    def __init__(self, lambda_0: float, *, alpha: float = 0.9, sigma: float = 0.1,
                 turnover: float = 0.0, sampling_fraction: float = 1.0,
                 mass_extinctions=None):
        self.lambda_0 = float(lambda_0)
        self.alpha = float(alpha)
        self.sigma = float(sigma)
        self.turnover = float(turnover)
        self.sampling_fraction = float(sampling_fraction)
        self.mass_extinctions = _normalize_mass_extinctions(mass_extinctions)

    def validate(self) -> None:
        _finite("lambda_0 (root speciation rate)", self.lambda_0)
        _finite("alpha (rate-shift trend)", self.alpha)
        _finite("sigma (rate-shift spread)", self.sigma)
        _finite("turnover (ε = μ/λ)", self.turnover)
        _finite("sampling_fraction", self.sampling_fraction)
        if self.lambda_0 <= 0:
            raise ValueError(f"lambda_0 (root speciation rate) must be > 0, got {self.lambda_0}")
        if self.alpha <= 0:
            raise ValueError(f"alpha (rate-shift trend) must be > 0, got {self.alpha}")
        if self.sigma < 0:
            raise ValueError(f"sigma (rate-shift spread) must be >= 0, got {self.sigma}")
        if not (0.0 <= self.turnover < 1.0):
            raise ValueError(f"turnover (ε = μ/λ) must be in [0, 1), got {self.turnover}")
        if not (0.0 < self.sampling_fraction <= 1.0):
            raise ValueError(f"sampling_fraction must be in (0, 1], got {self.sampling_fraction}")
        _validate_mass_extinctions(self.mass_extinctions)


class DiversityDependent:
    """Diversity-dependent (density-dependent) birth–death: diversification slows as the tree
    fills its carrying capacity ``K`` (Rabosky & Lovette 2008; Etienne et al. 2012).

    The speciation rate declines linearly with the number of standing lineages ``n``,
    ``λ(n) = max(0, λ₀·(1 − n/K))``, while extinction ``μ`` is constant. The tree grows fast
    when small and saturates near ``K`` (with ``μ=0``) or near the equilibrium
    ``n* = K·(1 − μ/λ₀)``. All lineages share the current ``λ(n)``, so it is *homogeneous* but
    *time/diversity-varying* — the ecological counterpart of the per-family carrying capacity
    ZOMBI2 already offers for gene families.

    A **forward-only** model (the rate depends on the running lineage count). ``age`` or
    ``n_tips`` mode both work — but ``n_tips`` must be ``≤ K`` (the tree cannot grow past its
    capacity). ``sampling_fraction`` (ρ) and ``mass_extinctions`` overlay as for
    :class:`BirthDeath` (mass extinctions still require ``age`` mode).
    """

    caps = SpeciesCaps(GrowthEngine.GILLESPIE, supports_n_tips=True)

    def __init__(self, lambda_0: float, death: float = 0.0, *, carrying_capacity: float,
                 sampling_fraction: float = 1.0, mass_extinctions=None):
        self.lambda_0 = float(lambda_0)
        self.death = float(death)
        self.K = float(carrying_capacity)
        self.sampling_fraction = float(sampling_fraction)
        self.mass_extinctions = _normalize_mass_extinctions(mass_extinctions)

    def validate(self) -> None:
        _finite("lambda_0 (speciation rate)", self.lambda_0)
        _finite("death rate", self.death)
        _finite("carrying_capacity K", self.K)
        _finite("sampling_fraction", self.sampling_fraction)
        if self.lambda_0 <= 0:
            raise ValueError(f"lambda_0 (speciation rate) must be > 0, got {self.lambda_0}")
        if self.death < 0:
            raise ValueError(f"death rate must be >= 0, got {self.death}")
        if self.K <= 0:
            raise ValueError(f"carrying_capacity K must be > 0, got {self.K}")
        if not (0.0 < self.sampling_fraction <= 1.0):
            raise ValueError(f"sampling_fraction must be in (0, 1], got {self.sampling_fraction}")
        _validate_mass_extinctions(self.mass_extinctions)


class CladeShiftBirthDeath:
    """Constant-rate birth–death with a finite set of **clade-specific rate shifts**.

    Diversification runs at the background ``(birth, death)`` until a scheduled shift: at each
    ``age`` before the present, one lineage then alive — chosen uniformly at random, since
    contemporaneous lineages are exchangeable — and *all of its descendants* switch to a new
    ``(birth, death)`` regime. This is the discrete, hand-specified version of clade rate
    heterogeneity: "at this time, some clade starts diversifying under these rates" (a key
    innovation sparking a radiation, a lineage entering a slow-down). Give several shifts for
    several radiating/collapsing clades.

    A **forward-only** model in **age mode** (the shifts are scheduled as ages before a fixed
    present). ``sampling_fraction`` (ρ) and ``mass_extinctions`` overlay as for :class:`BirthDeath`.

    Parameters
    ----------
    birth, death:
        Background speciation/extinction rates (before any shift).
    clade_shifts:
        List of ``(age, birth, death)`` — at ``age`` before the present, a random extant lineage
        and its descendants adopt ``(birth, death)``.
    """

    caps = SpeciesCaps(GrowthEngine.GILLESPIE, supports_n_tips=False)

    def __init__(self, birth: float, death: float = 0.0, *, clade_shifts,
                 sampling_fraction: float = 1.0, mass_extinctions=None):
        self.birth = float(birth)
        self.death = float(death)
        self.clade_shifts = sorted((float(a), float(b), float(d)) for a, b, d in clade_shifts)
        self.sampling_fraction = float(sampling_fraction)
        self.mass_extinctions = _normalize_mass_extinctions(mass_extinctions)

    def validate(self) -> None:
        _finite("birth rate", self.birth)
        _finite("death rate", self.death)
        _finite("sampling_fraction", self.sampling_fraction)
        if self.birth <= 0:
            raise ValueError(f"birth rate must be > 0, got {self.birth}")
        if self.death < 0:
            raise ValueError(f"death rate must be >= 0, got {self.death}")
        if not self.clade_shifts:
            raise ValueError("CladeShiftBirthDeath needs at least one (age, birth, death) shift")
        for a, b, d in self.clade_shifts:
            _finite("clade-shift age", a)
            _finite("clade-shift birth rate", b)
            _finite("clade-shift death rate", d)
            if a <= 0:
                raise ValueError(f"clade-shift age must be a positive age before the present, got {a}")
            if b <= 0:
                raise ValueError(f"clade-shift birth rate must be > 0, got {b}")
            if d < 0:
                raise ValueError(f"clade-shift death rate must be >= 0, got {d}")
        if not (0.0 < self.sampling_fraction <= 1.0):
            raise ValueError(f"sampling_fraction must be in (0, 1], got {self.sampling_fraction}")
        _validate_mass_extinctions(self.mass_extinctions)
