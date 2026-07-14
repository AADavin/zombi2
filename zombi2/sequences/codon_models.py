"""Codon substitution models — GY94 and MG94, with ``dN/dS`` (``omega``) as a direct parameter.

Where :mod:`~zombi2.sequences.models` evolves nucleotides (``ACGT``) or amino acids independently,
a **codon** model evolves the coding sequence one *codon* at a time over the 61 sense codons, so the
distinction between **synonymous** (amino-acid preserving) and **non-synonymous** changes — the raw
material of ``dN/dS`` — is built into the state space. Selection enters through a single knob
``omega`` (``ω = dN/dS``): every non-synonymous instantaneous rate is multiplied by ``ω``, so ``ω<1``
is purifying selection, ``ω=1`` neutrality, ``ω>1`` positive selection. A transition/transversion
bias ``kappa`` acts on the underlying nucleotide change, exactly as in HKY/K80.

Two classical parameterisations, both time-reversible (so they reuse the scipy-free ``exp(Qt)``
eigendecomposition of :class:`~zombi2.sequences.models.SubstitutionModel` unchanged):

* **GY94** (Goldman & Yang 1994): for codons ``i≠j`` differing at exactly one nucleotide (a
  transition ``a↔g``/``c↔t`` or a transversion), ``Q_ij = π_j · κ^{ts} · ω^{nonsyn}``, where ``π_j``
  is the equilibrium frequency of the *target codon*. Changes touching >1 position or landing on a
  stop codon have rate 0. Codon frequencies come from a frequency model (``F1×4``/``F3×4``/``F61``).
* **MG94** (Muse & Gaut 1994): the mutation is written at the *nucleotide* level, so the rate uses the
  equilibrium frequency of the single *nucleotide* being introduced, ``Q_ij = π*_{b} · κ^{ts} ·
  ω^{nonsyn}``. Its stationary distribution is still the product-of-nucleotide codon frequency, so
  ``ω`` cleanly separates the neutral mutation process from selection.

Both are normalised so a branch of length ``t`` is ``t`` expected codon substitutions per codon site.
The models plug straight into :func:`~zombi2.sequences.models.evolve_on_tree` (their alphabet is the
61-codon tuple, so a codon-state array decodes to an in-frame DNA string) and, via
:func:`~zombi2.sequences.models.make_model`, into ``zombi2 sequence --subst-model gy94`` /
``mg94 --omega``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from zombi2.sequences.models import BASES, CODON_MODELS, SubstitutionModel, is_codon_model

# --------------------------------------------------------------------------- #
# The standard genetic code (NCBI translation table 1), in TCAG codon order.
# --------------------------------------------------------------------------- #
__all__ = [
    "gy94", "mg94", "make_codon_model", "is_codon_model", "CODON_MODELS",
    "GENETIC_CODE", "SENSE_CODONS", "STOP_CODONS", "translate",
    "f1x4", "f3x4", "f61", "expected_dnds",
    # codon site models (dN/dS varies among sites)
    "CodonSiteModel", "m1a", "m2a", "m3", "m7", "m8", "make_codon_site_model",
    "is_codon_site_model", "CODON_SITE_MODELS", "beta_category_omegas",
]

_TCAG = "TCAG"
CODONS64 = tuple(x + y + z for x in _TCAG for y in _TCAG for z in _TCAG)
# amino acid encoded by each codon in CODONS64 order ('*' = stop):
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
#: standard genetic code, ``{codon: amino-acid letter}`` (``'*'`` for the three stop codons)
GENETIC_CODE = dict(zip(CODONS64, _AAS))
#: the 61 sense codons, in TCAG order — the state alphabet of every codon model
SENSE_CODONS = tuple(c for c in CODONS64 if GENETIC_CODE[c] != "*")
#: the three stop codons (TAA, TAG, TGA in the standard code)
STOP_CODONS = frozenset(c for c in CODONS64 if GENETIC_CODE[c] == "*")

_NB = {b: i for i, b in enumerate(BASES)}          # A,C,G,T -> 0,1,2,3
_PURINES = frozenset("AG")


def translate(dna: str) -> str:
    """Translate an in-frame coding DNA string to amino acids (``'*'`` for a stop codon)."""
    if len(dna) % 3:
        raise ValueError(f"coding sequence length {len(dna)} is not a multiple of 3")
    dna = dna.upper()
    return "".join(GENETIC_CODE[dna[i:i + 3]] for i in range(0, len(dna), 3))


def _is_transition(a: str, b: str) -> bool:
    """True if the nucleotide change ``a->b`` is a transition (purine↔purine or pyrimidine↔pyrimidine)."""
    return (a in _PURINES) == (b in _PURINES)


def _single_diff(c1: str, c2: str):
    """If codons ``c1``/``c2`` differ at exactly one position, return ``(pos, a, b)`` else ``None``."""
    diffs = [k for k in range(3) if c1[k] != c2[k]]
    if len(diffs) != 1:
        return None
    k = diffs[0]
    return k, c1[k], c2[k]


# --------------------------------------------------------------------------- #
# Codon equilibrium frequencies
# --------------------------------------------------------------------------- #
def _validate_nt(nt) -> np.ndarray:
    nt = np.asarray(nt, dtype=float)
    if nt.shape != (4,) or nt.min() < 0 or not np.isclose(nt.sum(), 1.0):
        raise ValueError(f"nucleotide frequencies must be 4 non-negative values summing to 1, got {nt}")
    return nt / nt.sum()


def _validate_pos(nt3) -> np.ndarray:
    nt3 = np.asarray(nt3, dtype=float)
    if nt3.shape != (3, 4) or nt3.min() < 0 or not np.allclose(nt3.sum(axis=1), 1.0):
        raise ValueError("position-specific frequencies must be a 3×4 array whose rows sum to 1")
    return nt3 / nt3.sum(axis=1, keepdims=True)


def f1x4(base_freqs) -> np.ndarray:
    """``F1×4`` codon frequencies: ``π_codon ∝ Π f_nt`` from one set of base freqs ``(A,C,G,T)``."""
    f = _validate_nt(base_freqs)
    pi = np.array([f[_NB[c[0]]] * f[_NB[c[1]]] * f[_NB[c[2]]] for c in SENSE_CODONS])
    return pi / pi.sum()


def f3x4(position_freqs) -> np.ndarray:
    """``F3×4`` codon frequencies: ``π_codon ∝ Π f_nt`` from position-specific base freqs (3×4)."""
    f = _validate_pos(position_freqs)
    pi = np.array([f[0, _NB[c[0]]] * f[1, _NB[c[1]]] * f[2, _NB[c[2]]] for c in SENSE_CODONS])
    return pi / pi.sum()


def f61(codon_freqs=None) -> np.ndarray:
    """``F61`` codon frequencies: an explicit 61-vector over :data:`SENSE_CODONS` (default uniform)."""
    if codon_freqs is None:
        return np.full(len(SENSE_CODONS), 1.0 / len(SENSE_CODONS))
    pi = np.asarray(codon_freqs, dtype=float)
    if pi.shape != (len(SENSE_CODONS),) or pi.min() <= 0 or not np.isclose(pi.sum(), 1.0):
        raise ValueError(f"F61 needs {len(SENSE_CODONS)} positive frequencies summing to 1")
    return pi / pi.sum()


def _codon_pi_from(freqs) -> np.ndarray:
    """Resolve a GY94 ``freqs`` argument to a 61-codon stationary vector.

    ``None`` -> uniform F61; a length-4 ``(A,C,G,T)`` -> F1×4; a 3×4 array -> F3×4; a length-61
    array -> explicit F61.
    """
    if freqs is None:
        return f61()
    arr = np.asarray(freqs, dtype=float)
    if arr.shape == (4,):
        return f1x4(arr)
    if arr.shape == (3, 4):
        return f3x4(arr)
    if arr.shape == (len(SENSE_CODONS),):
        return f61(arr)
    raise ValueError("gy94 freqs must be None, a length-4 (A,C,G,T), a 3×4 array, or a length-61 array")


# --------------------------------------------------------------------------- #
# Matrix builders
# --------------------------------------------------------------------------- #
def _fill_diagonal(Q: np.ndarray) -> np.ndarray:
    """Set each diagonal to minus its off-diagonal row sum (a valid rate matrix); return ``Q``."""
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def _gy94_offdiag(kappa: float, omega: float, pi: np.ndarray) -> np.ndarray:
    """GY94 off-diagonal rates ``π_j · κ^{ts} · ω^{nonsyn}`` for single-nucleotide codon neighbours."""
    n = len(SENSE_CODONS)
    Q = np.zeros((n, n))
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            d = _single_diff(c1, c2) if i != j else None
            if d is None:
                continue
            _, a, b = d
            rate = pi[j]
            if _is_transition(a, b):
                rate *= kappa
            if GENETIC_CODE[c1] != GENETIC_CODE[c2]:
                rate *= omega
            Q[i, j] = rate
    return Q


def _mg94_nt3(freqs) -> np.ndarray:
    """Resolve an MG94 ``freqs`` argument to a 3×4 array of position-specific nucleotide frequencies."""
    if freqs is None:
        return np.full((3, 4), 0.25)
    arr = np.asarray(freqs, dtype=float)
    if arr.shape == (4,):
        return np.tile(_validate_nt(arr), (3, 1))
    if arr.shape == (3, 4):
        return _validate_pos(arr)
    raise ValueError("mg94 freqs must be None, a length-4 (A,C,G,T), or a 3×4 array")


def _mg94_offdiag(kappa: float, omega: float, nt3: np.ndarray) -> np.ndarray:
    """MG94 off-diagonal rates ``π*_b · κ^{ts} · ω^{nonsyn}`` (introduced-nucleotide frequency)."""
    n = len(SENSE_CODONS)
    Q = np.zeros((n, n))
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            d = _single_diff(c1, c2) if i != j else None
            if d is None:
                continue
            pos, a, b = d
            rate = nt3[pos, _NB[b]]
            if _is_transition(a, b):
                rate *= kappa
            if GENETIC_CODE[c1] != GENETIC_CODE[c2]:
                rate *= omega
            Q[i, j] = rate
    return Q


def _raw_codon_Q(base: str, kappa: float, omega: float, freqs):
    """Un-normalised codon rate matrix (diagonal filled) + stationary ``π`` for ``gy94``/``mg94``."""
    if kappa < 0 or omega < 0:
        raise ValueError("kappa and omega must be non-negative")
    base = base.lower()
    if base == "gy94":
        pi = _codon_pi_from(freqs)
        Q = _gy94_offdiag(kappa, omega, pi)
    elif base == "mg94":
        nt3 = _mg94_nt3(freqs)
        pi = f3x4(nt3)
        Q = _mg94_offdiag(kappa, omega, nt3)
    else:
        raise ValueError(f"unknown codon model {base!r} (choose from {CODON_MODELS})")
    return _fill_diagonal(Q), pi


def _mean_rate(Q: np.ndarray, pi: np.ndarray) -> float:
    """Expected substitution rate ``-Σ π_i Q_ii`` of a codon matrix."""
    return -(pi * np.diag(Q)).sum()


def gy94(kappa: float = 2.0, omega: float = 1.0, *, freqs=None) -> SubstitutionModel:
    """Goldman & Yang (1994) codon model — ``Q_ij = π_j · κ^{ts} · ω^{nonsyn}``.

    ``kappa`` is the transition/transversion rate ratio, ``omega`` (``= dN/dS``) scales every
    non-synonymous rate. ``freqs`` selects the codon frequency model: ``None`` = uniform ``F61``;
    a length-4 ``(A,C,G,T)`` = ``F1×4``; a 3×4 array = ``F3×4``; a length-61 vector = explicit ``F61``.
    """
    Q, pi = _raw_codon_Q("gy94", kappa, omega, freqs)
    scale = _mean_rate(Q, pi)
    if scale <= 0:
        raise ValueError("degenerate codon model (zero substitution rate)")
    return SubstitutionModel("GY94", Q / scale, pi, SENSE_CODONS)


def mg94(kappa: float = 2.0, omega: float = 1.0, *, freqs=None) -> SubstitutionModel:
    """Muse & Gaut (1994) codon model — ``Q_ij = π*_b · κ^{ts} · ω^{nonsyn}`` (mutation at the
    nucleotide level, so the rate uses the frequency of the introduced nucleotide ``b``).

    ``kappa``/``omega`` as in :func:`gy94`. ``freqs`` gives the nucleotide frequencies: ``None`` =
    uniform (0.25 each, i.e. ``F1×4``); a length-4 ``(A,C,G,T)`` = ``F1×4``; a 3×4 array = ``F3×4``
    (position-specific). The stationary distribution is the corresponding product codon frequency.
    """
    Q, pi = _raw_codon_Q("mg94", kappa, omega, freqs)
    scale = _mean_rate(Q, pi)
    if scale <= 0:
        raise ValueError("degenerate codon model (zero substitution rate)")
    return SubstitutionModel("MG94", Q / scale, pi, SENSE_CODONS)


def make_codon_model(name: str, *, kappa: float = 2.0, omega: float = 1.0,
                     freqs=None) -> SubstitutionModel:
    """Construct a codon model by name (``gy94`` or ``mg94``)."""
    name = name.lower()
    if name == "gy94":
        return gy94(kappa, omega, freqs=freqs)
    if name == "mg94":
        return mg94(kappa, omega, freqs=freqs)
    raise ValueError(f"unknown codon model {name!r} (choose from {CODON_MODELS})")


# --------------------------------------------------------------------------- #
# dN/dS bookkeeping (a validation oracle: recovers omega from the matrix itself)
# --------------------------------------------------------------------------- #
def _syn_masks():
    """Boolean 61×61 masks of single-nucleotide (synonymous, non-synonymous) codon neighbours."""
    n = len(SENSE_CODONS)
    syn = np.zeros((n, n), dtype=bool)
    non = np.zeros((n, n), dtype=bool)
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            if i == j or _single_diff(c1, c2) is None:
                continue
            if GENETIC_CODE[c1] == GENETIC_CODE[c2]:
                syn[i, j] = True
            else:
                non[i, j] = True
    return syn, non


def _fluxes(model: SubstitutionModel):
    """Stationary synonymous / non-synonymous substitution flux ``Σ_i π_i Σ_j Q_ij`` of a codon model."""
    if tuple(model.alphabet) != SENSE_CODONS:
        raise ValueError("expects a codon model over the 61 sense codons")
    syn, non = _syn_masks()
    w = model.stationary[:, None] * model.Q
    return float((w * syn).sum()), float((w * non).sum())


def expected_dnds(model: SubstitutionModel, neutral: SubstitutionModel) -> float:
    """The model's genome-wide ``dN/dS``, using ``neutral`` (its ``ω=1`` twin) to count opportunities.

    ``dN/dS = (N_flux/N_sites) / (S_flux/S_sites)`` where the "sites" (mutational opportunities) are
    the synonymous / non-synonymous fluxes of the neutral model. Per-model rate normalisation cancels,
    so for a single-``ω`` GY94/MG94 this returns ``omega`` exactly. It is an independent, matrix-level
    check that a codon model was assembled correctly (κ, the stop mask, and the ω weighting) — the
    simulation-based estimator in the tests confirms the same number from evolved sequences.
    """
    s_m, n_m = _fluxes(model)
    s_0, n_0 = _fluxes(neutral)
    if n_0 <= 0 or s_m <= 0 or s_0 <= 0:
        raise ValueError("degenerate fluxes; check the models share a mutation process")
    return (n_m / n_0) * (s_0 / s_m)


# --------------------------------------------------------------------------- #
# Codon site models — dN/dS varies among sites (Nielsen & Yang 1998; Yang et al. 2000)
#
# Selection is rarely uniform along a gene: most codons are conserved (ω<1), some evolve neutrally
# (ω=1), and a few may be under positive selection (ω>1). A *site model* draws each codon site's ω
# from a distribution and evolves it under the matching GY94/MG94 matrix. Because ω does not change
# the stationary distribution (it only reweights non-synonymous rates), every component matrix shares
# the same codon frequencies π and the same mutation process — they differ only in ω.
# --------------------------------------------------------------------------- #
def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function (Lentz's method, Numerical Recipes)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return h


def _betainc(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta ``I_x(a, b)`` — numpy/math only (no scipy)."""
    import math
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _beta_quantile(p: float, a: float, b: float) -> float:
    """Inverse of ``I_x(a, b)``: the ``x`` with ``I_x(a, b) = p`` (bisection). Numpy/math only."""
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _betainc(a, b, mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def beta_category_omegas(p: float, q: float, ncat: int) -> np.ndarray:
    """Discretise ``Beta(p, q)`` on ``[0, 1]`` into ``ncat`` equal-probability category **means**.

    The category values are the means of ``Beta(p, q)`` within each ``1/ncat`` probability bin (the
    codon analogue of the discrete-Gamma rates), so their average is exactly the beta mean
    ``p/(p+q)``. All values lie in ``[0, 1]`` — this is the ``ω`` distribution of the M7/M8 models.
    """
    if p <= 0 or q <= 0 or ncat < 1:
        raise ValueError("beta shapes p, q must be > 0 and ncat >= 1")
    cuts = [_beta_quantile(i / ncat, p, q) for i in range(1, ncat)]
    bounds = [0.0] + cuts + [1.0]
    mean = p / (p + q)
    omegas = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        # mean of Beta(p,q) truncated to [lo,hi] with mass 1/ncat = mean · ncat · [I_hi(p+1,q) - I_lo(...)]
        mass_hi = _betainc(p + 1.0, q, hi)
        mass_lo = _betainc(p + 1.0, q, lo)
        omegas.append(mean * ncat * (mass_hi - mass_lo))
    return np.asarray(omegas)


@dataclass(frozen=True)
class CodonSiteModel:
    """A mixture of codon models over sites: component matrices (one per ``ω`` class) + proportions.

    Every component is a GY94/MG94 matrix that shares the same mutation process and codon frequencies
    (hence the same stationary distribution) and differs only in ``ω``. A site draws a class from
    ``proportions`` and evolves under that component; :func:`~zombi2.sequences.models.evolve_on_tree`
    understands the mixture directly. The genome-wide ``dN/dS`` is the proportion-weighted mean ``ω``
    (:attr:`mean_omega`), because synonymous flux is shared and non-synonymous flux scales with ``ω``.
    """

    name: str
    components: "tuple[SubstitutionModel, ...]"
    proportions: np.ndarray
    omegas: np.ndarray

    def __post_init__(self):
        comps = tuple(self.components)
        object.__setattr__(self, "components", comps)
        object.__setattr__(self, "proportions", np.asarray(self.proportions, dtype=float))
        object.__setattr__(self, "omegas", np.asarray(self.omegas, dtype=float))
        n = len(comps)
        if n == 0 or self.proportions.shape != (n,) or self.omegas.shape != (n,):
            raise ValueError("components, proportions and omegas must have matching non-zero length")
        if self.proportions.min() < 0 or not np.isclose(self.proportions.sum(), 1.0):
            raise ValueError(f"proportions must be non-negative and sum to 1, got {self.proportions}")
        a0 = comps[0].alphabet
        if any(tuple(c.alphabet) != tuple(a0) for c in comps):
            raise ValueError("all components must share the same (codon) alphabet")
        if any(not np.allclose(c.stationary, comps[0].stationary) for c in comps):
            raise ValueError("all components must share the same stationary distribution")

    @property
    def stationary(self) -> np.ndarray:
        return self.components[0].stationary

    @property
    def alphabet(self):
        return self.components[0].alphabet

    @property
    def k(self) -> int:
        return self.components[0].k

    @property
    def unit(self) -> int:
        return self.components[0].unit

    @property
    def mean_omega(self) -> float:
        """The mixture's genome-wide ``dN/dS`` — the proportion-weighted mean of the class ``ω`` values."""
        return float((self.proportions * self.omegas).sum())


def _site_model(name, base, kappa, omegas, proportions, freqs) -> CodonSiteModel:
    """Assemble a site mixture: one raw matrix per ``ω`` class, all divided by a **single** shared
    scale (the proportion-weighted mean raw rate) so the mixture averages one substitution per codon
    per unit branch length. Sharing the scale — rather than normalising each class to rate 1 —
    preserves across-class rate heterogeneity (purifying classes evolve slower), which is what makes
    the mixture's genome-wide dN/dS equal the proportion-weighted mean ω."""
    omegas = np.asarray(omegas, dtype=float)
    props = np.asarray(proportions, dtype=float)
    raws = [_raw_codon_Q(base, kappa, float(w), freqs) for w in omegas]
    pi = raws[0][1]
    rates = np.array([_mean_rate(Q, pi) for Q, _ in raws])
    shared = float((props * rates).sum())
    if shared <= 0:
        raise ValueError("degenerate codon site model (zero mean substitution rate)")
    comps = tuple(SubstitutionModel(f"{name}#{i}", Q / shared, pi, SENSE_CODONS)
                  for i, (Q, _) in enumerate(raws))
    return CodonSiteModel(name, comps, props, omegas)


def m1a(kappa: float = 2.0, *, p0: float = 0.6, omega0: float = 0.1,
        freqs=None, base: str = "gy94") -> CodonSiteModel:
    """M1a "nearly neutral": a purifying class ``ω0 < 1`` (proportion ``p0``) and a neutral class
    ``ω = 1`` (proportion ``1 − p0``). No positive selection."""
    if not 0.0 <= p0 <= 1.0:
        raise ValueError("p0 must be in [0, 1]")
    if omega0 < 0 or omega0 > 1:
        raise ValueError("omega0 must be in [0, 1] (the purifying class)")
    return _site_model("M1a", base, kappa, [omega0, 1.0], [p0, 1.0 - p0], freqs)


def m2a(kappa: float = 2.0, *, p0: float = 0.6, omega0: float = 0.1, p1: float = 0.3,
        omega2: float = 2.0, freqs=None, base: str = "gy94") -> CodonSiteModel:
    """M2a "positive selection": M1a's purifying (``ω0``, ``p0``) and neutral (``ω=1``, ``p1``) classes
    plus a positive-selection class ``ω2 > 1`` with the remaining proportion ``1 − p0 − p1``."""
    p2 = 1.0 - p0 - p1
    if min(p0, p1, p2) < 0:
        raise ValueError("p0 + p1 must be <= 1 (p2 = 1 - p0 - p1 is the positive class)")
    if omega0 > 1 or omega2 < 1:
        raise ValueError("expected omega0 <= 1 (purifying) and omega2 >= 1 (positive)")
    return _site_model("M2a", base, kappa, [omega0, 1.0, omega2], [p0, p1, p2], freqs)


def m3(kappa: float = 2.0, *, omegas, proportions, freqs=None,
       base: str = "gy94") -> CodonSiteModel:
    """M3 "discrete": ``K`` site classes with free ``ω`` values and proportions (renormalised)."""
    omegas = np.asarray(omegas, dtype=float)
    props = np.asarray(proportions, dtype=float)
    if omegas.shape != props.shape or omegas.ndim != 1 or omegas.size == 0:
        raise ValueError("omegas and proportions must be equal-length 1-D arrays")
    if props.min() < 0 or props.sum() <= 0 or omegas.min() < 0:
        raise ValueError("proportions must be non-negative (and not all zero); omegas non-negative")
    return _site_model("M3", base, kappa, omegas, props / props.sum(), freqs)


def m7(kappa: float = 2.0, *, beta_p: float, beta_q: float, ncat: int = 4,
       freqs=None, base: str = "gy94") -> CodonSiteModel:
    """M7 "beta": ``ω ~ Beta(beta_p, beta_q)`` on ``[0, 1]``, discretised into ``ncat`` equal-weight
    classes. A flexible purifying/neutral null with no positive selection."""
    omegas = beta_category_omegas(beta_p, beta_q, ncat)
    return _site_model("M7", base, kappa, omegas, np.full(ncat, 1.0 / ncat), freqs)


def m8(kappa: float = 2.0, *, beta_p: float, beta_q: float, p0: float = 0.9,
       omega_s: float = 2.0, ncat: int = 4, freqs=None, base: str = "gy94") -> CodonSiteModel:
    """M8 "beta & ω": M7's ``Beta(beta_p, beta_q)`` classes with total proportion ``p0``, plus one
    positive-selection class ``ω_s ≥ 1`` with proportion ``1 − p0``. The standard positive-selection
    test against M7."""
    if not 0.0 <= p0 <= 1.0:
        raise ValueError("p0 must be in [0, 1]")
    if omega_s < 1:
        raise ValueError("omega_s must be >= 1 (the positive-selection class)")
    beta_omegas = beta_category_omegas(beta_p, beta_q, ncat)
    omegas = np.concatenate([beta_omegas, [omega_s]])
    props = np.concatenate([np.full(ncat, p0 / ncat), [1.0 - p0]])
    return _site_model("M8", base, kappa, omegas, props, freqs)


CODON_SITE_MODELS = ("m1a", "m2a", "m3", "m7", "m8")


def is_codon_site_model(name: str) -> bool:
    """True if ``name`` is a codon site model (``m1a``/``m2a``/``m3``/``m7``/``m8``)."""
    return name.lower() in CODON_SITE_MODELS


def make_codon_site_model(name: str, *, kappa: float = 2.0, base: str = "gy94", freqs=None,
                          p0=None, omega0=None, p1=None, omega2=None,
                          beta_p=None, beta_q=None, omega_s=None, ncat: int = 4,
                          omegas=None, proportions=None) -> CodonSiteModel:
    """Construct a codon site model by name, passing only the parameters that model uses.

    Parameters left ``None`` fall back to each constructor's default; ``base`` is the underlying
    codon matrix (``gy94``/``mg94``). Raises if a required parameter for the chosen model is missing.
    """
    name = name.lower()

    def need(**kw):
        missing = [k for k, v in kw.items() if v is None]
        if missing:
            raise ValueError(f"codon site model {name!r} requires {missing}")

    if name == "m1a":
        kw = dict(kappa=kappa, freqs=freqs, base=base)
        if p0 is not None:
            kw["p0"] = p0
        if omega0 is not None:
            kw["omega0"] = omega0
        return m1a(**kw)
    if name == "m2a":
        kw = dict(kappa=kappa, freqs=freqs, base=base)
        for k, v in (("p0", p0), ("omega0", omega0), ("p1", p1), ("omega2", omega2)):
            if v is not None:
                kw[k] = v
        return m2a(**kw)
    if name == "m3":
        need(omegas=omegas, proportions=proportions)
        return m3(kappa=kappa, omegas=omegas, proportions=proportions, freqs=freqs, base=base)
    if name == "m7":
        need(beta_p=beta_p, beta_q=beta_q)
        return m7(kappa=kappa, beta_p=beta_p, beta_q=beta_q, ncat=ncat, freqs=freqs, base=base)
    if name == "m8":
        need(beta_p=beta_p, beta_q=beta_q)
        kw = dict(kappa=kappa, beta_p=beta_p, beta_q=beta_q, ncat=ncat, freqs=freqs, base=base)
        if p0 is not None:
            kw["p0"] = p0
        if omega_s is not None:
            kw["omega_s"] = omega_s
        return m8(**kw)
    raise ValueError(f"unknown codon site model {name!r} (choose from {CODON_SITE_MODELS})")
