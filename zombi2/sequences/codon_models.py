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

import numpy as np

from zombi2.sequences.models import BASES, CODON_MODELS, SubstitutionModel, is_codon_model

# --------------------------------------------------------------------------- #
# The standard genetic code (NCBI translation table 1), in TCAG codon order.
# --------------------------------------------------------------------------- #
__all__ = [
    "gy94", "mg94", "make_codon_model", "is_codon_model", "CODON_MODELS",
    "GENETIC_CODE", "SENSE_CODONS", "STOP_CODONS", "translate",
    "f1x4", "f3x4", "f61", "expected_dnds",
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
def _normalise(name: str, Q: np.ndarray, pi: np.ndarray) -> SubstitutionModel:
    """Fill the diagonal, scale to one expected substitution per codon site, and wrap as a model."""
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    scale = -(pi * np.diag(Q)).sum()
    if scale <= 0:
        raise ValueError("degenerate codon model (zero substitution rate)")
    return SubstitutionModel(name, Q / scale, pi, SENSE_CODONS)


def gy94(kappa: float = 2.0, omega: float = 1.0, *, freqs=None) -> SubstitutionModel:
    """Goldman & Yang (1994) codon model — ``Q_ij = π_j · κ^{ts} · ω^{nonsyn}``.

    ``kappa`` is the transition/transversion rate ratio, ``omega`` (``= dN/dS``) scales every
    non-synonymous rate. ``freqs`` selects the codon frequency model: ``None`` = uniform ``F61``;
    a length-4 ``(A,C,G,T)`` = ``F1×4``; a 3×4 array = ``F3×4``; a length-61 vector = explicit ``F61``.
    """
    if kappa < 0 or omega < 0:
        raise ValueError("kappa and omega must be non-negative")
    pi = _codon_pi_from(freqs)
    n = len(SENSE_CODONS)
    Q = np.zeros((n, n))
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            if i == j:
                continue
            d = _single_diff(c1, c2)
            if d is None:
                continue
            _, a, b = d
            rate = pi[j]
            if _is_transition(a, b):
                rate *= kappa
            if GENETIC_CODE[c1] != GENETIC_CODE[c2]:
                rate *= omega
            Q[i, j] = rate
    return _normalise("GY94", Q, pi)


def mg94(kappa: float = 2.0, omega: float = 1.0, *, freqs=None) -> SubstitutionModel:
    """Muse & Gaut (1994) codon model — ``Q_ij = π*_b · κ^{ts} · ω^{nonsyn}`` (mutation at the
    nucleotide level, so the rate uses the frequency of the introduced nucleotide ``b``).

    ``kappa``/``omega`` as in :func:`gy94`. ``freqs`` gives the nucleotide frequencies: ``None`` =
    uniform (0.25 each, i.e. ``F1×4``); a length-4 ``(A,C,G,T)`` = ``F1×4``; a 3×4 array = ``F3×4``
    (position-specific). The stationary distribution is the corresponding product codon frequency.
    """
    if kappa < 0 or omega < 0:
        raise ValueError("kappa and omega must be non-negative")
    if freqs is None:
        nt3 = np.full((3, 4), 0.25)
    else:
        arr = np.asarray(freqs, dtype=float)
        if arr.shape == (4,):
            nt3 = np.tile(_validate_nt(arr), (3, 1))
        elif arr.shape == (3, 4):
            nt3 = _validate_pos(arr)
        else:
            raise ValueError("mg94 freqs must be None, a length-4 (A,C,G,T), or a 3×4 array")
    pi = f3x4(nt3)
    n = len(SENSE_CODONS)
    Q = np.zeros((n, n))
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            if i == j:
                continue
            d = _single_diff(c1, c2)
            if d is None:
                continue
            pos, a, b = d
            rate = nt3[pos, _NB[b]]
            if _is_transition(a, b):
                rate *= kappa
            if GENETIC_CODE[c1] != GENETIC_CODE[c2]:
                rate *= omega
            Q[i, j] = rate
    return _normalise("MG94", Q, pi)


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
