"""Simulate DNA or protein sequences along a gene tree under a substitution model.

ZOMBI2 gene trees carry *branch lengths in substitutions per site* (the phylograms produced by
:class:`~zombi2.SequenceEvolution`). This module evolves an actual sequence down such a tree —
the classic "simulate sequences across a gene tree" of ZOMBI 1 — recording the sequence at **every**
node (internal and tip). Tip sequences are the observable gene alignment; the internal-node
sequences let the nucleotide model reconstruct the genome sequence at every ancestor.

The engine is generic over the alphabet size ``K``: nucleotides (``K=4``, ``ACGT``) and amino
acids (``K=20``, the 20-letter protein alphabet). A :class:`SubstitutionModel` carries a normalised
``K×K`` rate matrix ``Q`` (one expected substitution per site per unit branch length) plus the
stationary frequencies.

* **Nucleotide** (continuous-time Markov, 4 states ``ACGT``): :func:`jc69`, :func:`k80`,
  :func:`hky85`, :func:`gtr`.
* **Amino acid** (20 states): :func:`poisson` (equal rates, uniform frequencies) and the empirical
  models :func:`lg`, :func:`wag`, :func:`jtt`, :func:`dayhoff`, each built from published
  exchangeabilities ``S`` (symmetric) and stationary frequencies ``π`` as ``Q_ij = S_ij·π_j``.

Optional across-site rate heterogeneity via a discrete-Gamma (:class:`GammaRates`), for any ``K``.

The transition matrix over a branch of length ``t`` is ``P(t) = exp(Q·t)``. Because every model here
is time-reversible, ``exp(Qt)`` is computed by eigendecomposition of the *symmetric* matrix
``B = diag(√π) Q diag(1/√π)`` (numpy only — no scipy): ``P(t) = diag(1/√π) V exp(Λt) Vᵀ diag(√π)``.
A child state is drawn per site from ``P(t)[parent_state]``; everything is vectorised over sites.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._aa_models import (
    _DAYHOFF_EXCH, _DAYHOFF_PI, _JTT_EXCH, _JTT_PI, _LG_EXCH, _LG_PI, _WAG_EXCH, _WAG_PI,
)

BASES = "ACGT"
_CODE = {b: i for i, b in enumerate(BASES)}

# The 20-letter amino-acid alphabet, in the PAML column order used by every empirical matrix below.
AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"


# --------------------------------------------------------------------------- #
# Substitution models
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SubstitutionModel:
    """A K-state model: a normalised ``K×K`` rate matrix ``Q``, stationary freqs, and its alphabet.

    ``alphabet`` is the ordered string of state symbols (``"ACGT"`` for nucleotides, the 20-letter
    protein alphabet for amino acids); ``Q`` and ``stationary`` follow that order.
    """

    name: str
    Q: np.ndarray
    stationary: np.ndarray
    alphabet: str = BASES

    def __post_init__(self):
        # Precompute the reversible eigendecomposition once (numpy only) for fast, scipy-free exp(Qt).
        pi = self.stationary
        sq = np.sqrt(pi)
        B = (sq[:, None] * self.Q) / sq[None, :]      # symmetric similarity transform of Q
        B = (B + B.T) / 2.0                            # kill round-off asymmetry before eigh
        w, V = np.linalg.eigh(B)
        # store the pieces of P(t) = diag(1/√π) · V · exp(Λt) · Vᵀ · diag(√π)
        object.__setattr__(self, "_eigvals", w)
        object.__setattr__(self, "_left", V / sq[:, None])     # diag(1/√π) · V
        object.__setattr__(self, "_right", (V * sq[:, None]).T)  # Vᵀ · diag(√π)

    @property
    def k(self) -> int:
        """Number of states in the alphabet."""
        return self.Q.shape[0]

    def p_matrix(self, t: float) -> np.ndarray:
        """Transition probabilities over branch length ``t`` (substitutions/site).

        ``P(t) = exp(Qt)`` via the reversible eigendecomposition; clipped to ``[0, ∞)`` to scrub
        tiny negative round-off so every row is a valid probability distribution.
        """
        P = (self._left * np.exp(self._eigvals * t)) @ self._right
        return np.clip(P, 0.0, None)


def _reversible_model(name, S, pi, alphabet) -> SubstitutionModel:
    """Build a normalised reversible model from a symmetric exchangeability matrix ``S`` and ``pi``.

    ``Q_ij = S_ij · pi_j`` (i≠j), scaled so the expected rate ``-Σ pi_i Q_ii = 1`` (branch lengths
    are then in substitutions/site). ``S`` must be symmetric with a zero diagonal.
    """
    S = np.asarray(S, dtype=float)
    pi = np.asarray(pi, dtype=float)
    k = pi.shape[0]
    if pi.shape != (k,) or pi.min() < 0 or not np.isclose(pi.sum(), 1.0):
        raise ValueError(f"stationary frequencies must be non-negative and sum to 1, got {pi}")
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


def _tri_to_symmetric(tri, k: int) -> np.ndarray:
    """Expand a flat lower-triangular list (row 1..k-1, PAML order) into a symmetric ``k×k`` matrix."""
    S = np.zeros((k, k))
    it = iter(tri)
    for i in range(1, k):
        for j in range(i):
            S[i, j] = S[j, i] = next(it)
    return S


def _gtr_model(name, exch, pi) -> SubstitutionModel:
    """Build a GTR-family model from 6 exchangeabilities ``[AC,AG,AT,CG,CT,GT]`` and freqs ``pi``.

    ``Q_ij = exch_ij · pi_j`` (i≠j), scaled so the expected rate ``-Σ pi_i Q_ii = 1`` (branch
    lengths are then in substitutions/site).
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
    """Jukes--Cantor: equal rates, equal base frequencies."""
    return _gtr_model("JC69", [1, 1, 1, 1, 1, 1], [0.25] * 4)


def k80(kappa: float = 2.0) -> SubstitutionModel:
    """Kimura 2-parameter: transition/transversion ratio ``kappa``, equal frequencies."""
    return _gtr_model("K80", [1, kappa, 1, 1, kappa, 1], [0.25] * 4)


def hky85(kappa: float = 2.0, freqs=(0.25, 0.25, 0.25, 0.25)) -> SubstitutionModel:
    """HKY85: transition bias ``kappa`` with unequal base frequencies ``freqs`` (A,C,G,T)."""
    return _gtr_model("HKY85", [1, kappa, 1, 1, kappa, 1], freqs)


def gtr(rates=(1, 1, 1, 1, 1, 1), freqs=(0.25, 0.25, 0.25, 0.25)) -> SubstitutionModel:
    """General time-reversible: 6 exchangeabilities ``[AC,AG,AT,CG,CT,GT]`` and freqs (A,C,G,T)."""
    return _gtr_model("GTR", rates, freqs)


# --------------------------------------------------------------------------- #
# Amino-acid (protein) models — 20 states, empirical exchangeabilities + frequencies
# --------------------------------------------------------------------------- #
def _empirical_aa(name, tri, pi) -> SubstitutionModel:
    """Build a 20-state protein model from PAML lower-triangular exchangeabilities and freqs."""
    S = _tri_to_symmetric(tri, 20)
    return _reversible_model(name, S, pi, AMINO_ACIDS)


def poisson() -> SubstitutionModel:
    """Poisson (Felsenstein-81 for proteins): equal exchangeabilities, uniform frequencies."""
    S = np.ones((20, 20)) - np.eye(20)
    return _reversible_model("Poisson", S, np.full(20, 1.0 / 20.0), AMINO_ACIDS)


def lg() -> SubstitutionModel:
    """LG (Le & Gascuel 2008) — the modern default protein model."""
    return _empirical_aa("LG", _LG_EXCH, _LG_PI)


def wag() -> SubstitutionModel:
    """WAG (Whelan & Goldman 2001)."""
    return _empirical_aa("WAG", _WAG_EXCH, _WAG_PI)


def jtt() -> SubstitutionModel:
    """JTT (Jones, Taylor & Thornton 1992)."""
    return _empirical_aa("JTT", _JTT_EXCH, _JTT_PI)


def dayhoff() -> SubstitutionModel:
    """Dayhoff (Dayhoff et al. 1978, PAML values)."""
    return _empirical_aa("Dayhoff", _DAYHOFF_EXCH, _DAYHOFF_PI)


_NT_MODELS = {"jc69": jc69, "k80": k80, "hky85": hky85, "gtr": gtr}
_AA_MODELS = {"poisson": poisson, "lg": lg, "wag": wag, "jtt": jtt, "dayhoff": dayhoff}
_MODELS = {**_NT_MODELS, **_AA_MODELS}

#: names of the DNA (nucleotide) substitution models
DNA_MODELS = tuple(_NT_MODELS)
#: names of the protein (amino-acid) substitution models
PROTEIN_MODELS = tuple(_AA_MODELS)


def is_protein_model(name: str) -> bool:
    """True if ``name`` is one of the amino-acid models (used to auto-detect DNA vs protein)."""
    return name.lower() in _AA_MODELS


def make_model(name: str, *, kappa: float = 2.0, freqs=None, rates=None) -> SubstitutionModel:
    """Construct a model by name — DNA (``jc69``/``k80``/``hky85``/``gtr``) or protein
    (``poisson``/``lg``/``wag``/``jtt``/``dayhoff``) — with the relevant parameters.

    Protein models take no parameters (empirical); the DNA parameters are ignored for them.
    """
    name = name.lower()
    if name == "jc69":
        return jc69()
    if name == "k80":
        return k80(kappa)
    if name == "hky85":
        return hky85(kappa, freqs or (0.25, 0.25, 0.25, 0.25))
    if name == "gtr":
        return gtr(rates or (1, 1, 1, 1, 1, 1), freqs or (0.25, 0.25, 0.25, 0.25))
    if name in _AA_MODELS:
        return _AA_MODELS[name]()
    raise ValueError(f"unknown substitution model {name!r} (choose from {sorted(_MODELS)})")


# --------------------------------------------------------------------------- #
# Across-site rate heterogeneity (discrete Gamma)
# --------------------------------------------------------------------------- #
def _gammainc_lower(a: float, x: float) -> float:
    """Regularised lower incomplete gamma ``P(a, x)`` — numpy/math only (no scipy).

    Series expansion for ``x < a+1``, continued fraction otherwise (Numerical Recipes).
    """
    import math
    if x <= 0.0:
        return 0.0
    gln = math.lgamma(a)
    if x < a + 1.0:                       # series
        ap, s, term = a, 1.0 / a, 1.0 / a
        for _ in range(1000):
            ap += 1.0
            term *= x / ap
            s += term
            if abs(term) < abs(s) * 1e-14:
                break
        return s * math.exp(-x + a * math.log(x) - gln)
    # continued fraction (Lentz)
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return 1.0 - math.exp(-x + a * math.log(x) - gln) * h


def _gamma_quantile(p: float, shape: float) -> float:
    """Inverse of the regularised lower incomplete gamma: ``x`` with ``P(shape, x) = p``.

    Bisection on ``_gammainc_lower`` (a mean-1 Gamma uses ``scale = 1/shape``, applied by the
    caller). Numpy/math only.
    """
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return float("inf")
    lo, hi = 0.0, 1.0
    while _gammainc_lower(shape, hi) < p:
        hi *= 2.0
        if hi > 1e12:
            break
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _gammainc_lower(shape, mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass
class GammaRates:
    """Discrete-Gamma across-site rates: ``k`` equal-probability categories, mean rate 1.

    Category rates are the **means** within each equal-probability bin of a mean-1 Gamma(shape)
    (Yang 1994), computed numpy-only (no scipy). Works for any alphabet size.
    """

    shape: float
    k: int = 4

    def __post_init__(self):
        if self.shape <= 0 or self.k < 1:
            raise ValueError("Gamma shape must be > 0 and k >= 1")
        a = float(self.shape)
        # bin boundaries on a mean-1 Gamma (scale = 1/a): quantiles at i/k
        cuts = [_gamma_quantile(i / self.k, a) / a for i in range(1, self.k)]
        # mean of Gamma(a, scale=1/a) truncated to each bin = [F_{a+1}(hi)-F_{a+1}(lo)] * k
        # where F_{a+1} is the regularised lower incomplete gamma at shape a+1 (in units of scale=1/a)
        bounds = [0.0] + cuts + [float("inf")]
        rates = []
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            plo = _gammainc_lower(a + 1.0, lo * a) if lo > 0 else 0.0
            phi = _gammainc_lower(a + 1.0, hi * a) if hi != float("inf") else 1.0
            rates.append((phi - plo) * self.k)
        self.rates = np.asarray(rates) / np.mean(rates)


# --------------------------------------------------------------------------- #
# Sequence <-> integer state coding
# --------------------------------------------------------------------------- #
def encode(seq: str, rng: np.random.Generator | None = None,
           pi: np.ndarray | None = None, alphabet: str = BASES) -> np.ndarray:
    """Map a sequence string to integer states over ``alphabet``.

    Ambiguous/other letters draw from ``pi`` (or uniform). Defaults to the ``ACGT`` alphabet, so
    existing nucleotide callers are unchanged.
    """
    k = len(alphabet)
    code = {b: i for i, b in enumerate(alphabet)}
    states = np.empty(len(seq), dtype=np.int8)
    unknown = []
    for i, ch in enumerate(seq.upper()):
        c = code.get(ch, -1)
        if c < 0:
            unknown.append(i)
            c = 0
        states[i] = c
    if unknown:
        rng = rng or np.random.default_rng()
        p = pi if pi is not None else np.full(k, 1.0 / k)
        states[unknown] = rng.choice(k, size=len(unknown), p=p)
    return states


def decode(states: np.ndarray, alphabet: str = BASES) -> str:
    """Map integer states back to a string over ``alphabet`` (default ``ACGT``)."""
    arr = np.asarray(states)
    return "".join(alphabet[i] for i in arr)


_COMP = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string (for segments on the − strand)."""
    return seq.translate(_COMP)[::-1]


# --------------------------------------------------------------------------- #
# Evolve a sequence down a tree
# --------------------------------------------------------------------------- #
def evolve_on_tree(root, subst: dict, model: SubstitutionModel,
                   rng: np.random.Generator, *, root_seq=None, length: int | None = None,
                   gamma: GammaRates | None = None) -> dict:
    """Evolve a sequence over a gene tree; return ``{node.gid: end_sequence}`` for every node.

    ``root`` is a tree node (``reconciliation._Node``: has ``.gid`` and ``.children``); ``subst`` maps
    each node to the substitution length of the branch **ending at it** (as
    :func:`~zombi2.sequence_evolution._annotate` computes). The root's incoming sequence is
    ``root_seq`` (a DNA string) or, if ``None``, a fresh draw of ``length`` sites from the model's
    stationary frequencies. Each node's recorded sequence is its state at the end of its branch (a
    speciation/duplication/transfer/loss time, or the present for a tip); with ``subst=0`` a node
    simply copies its parent.
    """
    pi = model.stationary
    alphabet = model.alphabet
    k = model.k
    if root_seq is not None:
        root_states = encode(root_seq, rng, pi, alphabet)
        length = root_states.shape[0]
    else:
        if length is None:
            raise ValueError("give either root_seq or length")
        root_states = rng.choice(k, size=length, p=pi).astype(np.int8)

    if gamma is not None:
        site_cat = rng.integers(gamma.k, size=length)
    out: dict = {}
    pcache: dict = {}

    def p_for(t: float) -> np.ndarray:
        key = round(float(t), 12)
        P = pcache.get(key)
        if P is None:
            P = model.p_matrix(key)
            pcache[key] = P
        return P

    def sample(parent_states: np.ndarray, cum: np.ndarray) -> np.ndarray:
        r = rng.random(parent_states.shape[0])
        return (r[:, None] < cum[parent_states]).argmax(1).astype(np.int8)

    def visit(node, parent_states: np.ndarray) -> None:
        t = float(subst.get(node, 0.0))
        if t <= 0.0:
            states = parent_states
        elif gamma is None:
            states = sample(parent_states, p_for(t).cumsum(1))
        else:
            states = np.empty(length, dtype=np.int8)
            for c in range(gamma.k):
                mask = site_cat == c
                if mask.any():
                    cum = p_for(t * gamma.rates[c]).cumsum(1)
                    states[mask] = sample(parent_states[mask], cum)
        out[node.gid] = states
        for child in node.children:
            visit(child, states)

    visit(root, root_states)
    return {gid: decode(s, alphabet) for gid, s in out.items()}


# --------------------------------------------------------------------------- #
# FASTA I/O
# --------------------------------------------------------------------------- #
def read_fasta(path) -> dict:
    """Read a (optionally ``.gz``) FASTA into ``{seqid: sequence}`` (seqid = first token of header)."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    out: dict = {}
    name = None
    chunks: list[str] = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        out[name] = "".join(chunks)
    return out


def write_fasta(path, records: dict, *, gzip_out: bool = False, width: int = 70) -> None:
    """Write ``{name: sequence}`` to a FASTA file (gzipped when ``gzip_out`` or path ends ``.gz``)."""
    import gzip
    gz = gzip_out or str(path).endswith(".gz")
    opener = gzip.open if gz else open
    with opener(path, "wt") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
