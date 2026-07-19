"""Substitution models — the **menu** (``sequence-api.md``).

A substitution model is the *chemistry* of a sequence: a ``K×K`` rate matrix ``Q`` (normalised to
one expected substitution per site per unit branch length) and its stationary frequencies ``π``.
Different models are genuinely different matrices — Jukes–Cantor, K80, HKY85, GTR differ in their
transition/transversion structure and base composition — so, unlike the clock, they do **not**
collapse to one grammar: they stay a menu of constructors, each taking its own physical parameters
(``SPEC §4`` — "faking a grammar over the matrices would be worse than a menu").

This slice ships the four **nucleotide** models (4 states, ``ACGT``). Amino-acid and codon models,
and across-site ``+Γ`` heterogeneity, are named later slices; adding them is a pure extension of this
menu, no refactor.

Every model here is time-reversible, so the transition matrix over a branch of length ``t`` (in
substitutions/site), ``P(t) = exp(Q·t)``, is computed by eigendecomposition of the *symmetric*
matrix ``B = diag(√π)·Q·diag(1/√π)`` (numpy only, no scipy):
``P(t) = diag(1/√π)·V·exp(Λt)·Vᵀ·diag(√π)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: the nucleotide alphabet, in the order ``Q`` and ``stationary`` follow.
BASES = "ACGT"


@dataclass(frozen=True)
class SubstitutionModel:
    """A ``K``-state reversible model: a normalised ``K×K`` rate matrix ``Q``, its stationary
    frequencies, and the ordered ``alphabet`` whose order ``Q`` / ``stationary`` follow.

    Built through the menu constructors (:func:`jc69`, :func:`k80`, :func:`hky85`, :func:`gtr`),
    never directly. The reversible eigendecomposition behind :meth:`p_matrix` is precomputed once
    in ``__post_init__``.
    """

    name: str
    Q: np.ndarray
    stationary: np.ndarray
    alphabet: str = BASES

    def __post_init__(self) -> None:
        # Precompute the reversible eigendecomposition once (numpy only) for fast, scipy-free exp(Qt).
        pi = self.stationary
        sq = np.sqrt(pi)
        B = (sq[:, None] * self.Q) / sq[None, :]   # symmetric similarity transform of Q
        B = (B + B.T) / 2.0                         # kill round-off asymmetry before eigh
        w, V = np.linalg.eigh(B)
        # the pieces of P(t) = diag(1/√π) · V · exp(Λt) · Vᵀ · diag(√π)
        object.__setattr__(self, "_eigvals", w)
        object.__setattr__(self, "_left", V / sq[:, None])       # diag(1/√π) · V
        object.__setattr__(self, "_right", (V * sq[:, None]).T)   # Vᵀ · diag(√π)

    @property
    def k(self) -> int:
        """Number of states in the alphabet."""
        return self.Q.shape[0]

    def p_matrix(self, t: float) -> np.ndarray:
        """Transition probabilities over branch length ``t`` (substitutions/site).

        ``P(t) = exp(Qt)`` via the reversible eigendecomposition; clipped to ``[0, ∞)`` to scrub tiny
        negative round-off so every row is a valid probability distribution.
        """
        # The BLAS matmul kernel can raise spurious FP flags on larger matrices even when every input
        # is finite and the result is a valid stochastic matrix; silence them — the clip is the guard.
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            P = (self._left * np.exp(self._eigvals * t)) @ self._right
        return np.clip(P, 0.0, None)


def _reversible_model(name: str, S: np.ndarray, pi, alphabet: str = BASES) -> SubstitutionModel:
    """Build a normalised reversible model from a symmetric exchangeability matrix ``S`` and ``pi``.

    ``Q_ij = S_ij · pi_j`` (i≠j), scaled so the expected rate ``-Σ pi_i Q_ii = 1`` (branch lengths
    are then in substitutions/site). ``S`` must be symmetric with a zero diagonal.
    """
    S = np.asarray(S, dtype=float)
    pi = np.asarray(pi, dtype=float)
    k = pi.shape[0]
    if pi.shape != (k,) or pi.min() <= 0 or not np.isclose(pi.sum(), 1.0):
        raise ValueError(f"stationary frequencies must be strictly positive and sum to 1, got {pi} "
                         "(a zero-frequency state makes the rate matrix degenerate)")
    if S.shape != (k, k) or (S < 0).any() or not np.allclose(S, S.T):
        raise ValueError("exchangeabilities must be a symmetric non-negative K×K matrix")
    pi = pi / pi.sum()          # renormalise (published freqs round to 1 only to ~1e-6)
    Q = S * pi[None, :]
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    scale = -(pi * np.diag(Q)).sum()
    if scale <= 0:
        raise ValueError("degenerate substitution model (zero rate)")
    return SubstitutionModel(name, Q / scale, pi, alphabet)


def _gtr_model(name: str, exch, pi) -> SubstitutionModel:
    """Build a GTR-family nucleotide model from 6 exchangeabilities ``[AC,AG,AT,CG,CT,GT]`` and freqs.

    ``Q_ij = exch_ij · pi_j`` (i≠j), scaled so ``-Σ pi_i Q_ii = 1`` (branch lengths in subs/site).
    """
    pi = np.asarray(pi, dtype=float)
    if pi.shape != (4,):
        raise ValueError(f"stationary frequencies must be 4 values, got {pi}")
    ac, ag, at, cg, ct, gt = (float(x) for x in exch)
    if min(ac, ag, at, cg, ct, gt) < 0:
        raise ValueError("exchangeability rates must be non-negative")
    S = np.array([[0, ac, ag, at], [ac, 0, cg, ct], [ag, cg, 0, gt], [at, ct, gt, 0]], dtype=float)
    return _reversible_model(name, S, pi, BASES)


def jc69() -> SubstitutionModel:
    """Jukes–Cantor (1969): equal rates, equal base frequencies. No free parameters."""
    return _gtr_model("JC69", [1, 1, 1, 1, 1, 1], [0.25] * 4)


def k80(kappa: float = 2.0) -> SubstitutionModel:
    """Kimura 2-parameter (1980): transition/transversion ratio ``kappa``, equal frequencies."""
    return _gtr_model("K80", [1, kappa, 1, 1, kappa, 1], [0.25] * 4)


def hky85(kappa: float = 2.0, freqs=(0.25, 0.25, 0.25, 0.25)) -> SubstitutionModel:
    """HKY85 (Hasegawa–Kishino–Yano 1985): transition bias ``kappa`` with unequal base ``freqs`` (A,C,G,T)."""
    return _gtr_model("HKY85", [1, kappa, 1, 1, kappa, 1], freqs)


def gtr(rates=(1, 1, 1, 1, 1, 1), freqs=(0.25, 0.25, 0.25, 0.25)) -> SubstitutionModel:
    """General time-reversible: 6 exchangeabilities ``[AC,AG,AT,CG,CT,GT]`` and freqs (A,C,G,T)."""
    return _gtr_model("GTR", rates, freqs)


def decode(states: np.ndarray, alphabet: str = BASES) -> str:
    """Map an array of integer states back to a string over ``alphabet`` (default ``ACGT``)."""
    return "".join(alphabet[i] for i in np.asarray(states))


__all__ = ["SubstitutionModel", "jc69", "k80", "hky85", "gtr", "decode", "BASES"]
