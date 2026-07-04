"""Simulate DNA sequences along a gene tree under a nucleotide substitution model.

ZOMBI2 gene trees carry *branch lengths in substitutions per site* (the phylograms produced by
:class:`~zombi2.SequenceEvolution`). This module evolves an actual DNA sequence down such a tree —
the classic "simulate sequences across a gene tree" of ZOMBI 1 — recording the sequence at **every**
node (internal and tip). Tip sequences are the observable gene alignment; the internal-node
sequences let the nucleotide model reconstruct the genome sequence at every ancestor.

Models (continuous-time Markov, 4 states ``ACGT``): :func:`jc69`, :func:`k80`, :func:`hky85`,
:func:`gtr`, each a :class:`SubstitutionModel` carrying a rate matrix ``Q`` normalised to one
expected substitution per site per unit branch length, plus the stationary base frequencies.
Optional across-site rate heterogeneity via a discrete-Gamma (:class:`GammaRates`).

The transition matrix over a branch of length ``t`` is ``P(t) = expm(Q·t)``; a child state is drawn
per site from ``P(t)[parent_state]``. Everything is vectorised over sites with numpy.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm

BASES = "ACGT"
_CODE = {b: i for i, b in enumerate(BASES)}


# --------------------------------------------------------------------------- #
# Substitution models
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SubstitutionModel:
    """A nucleotide model: a normalised 4×4 rate matrix ``Q`` and stationary frequencies."""

    name: str
    Q: np.ndarray
    stationary: np.ndarray

    def p_matrix(self, t: float) -> np.ndarray:
        """Transition probabilities over branch length ``t`` (substitutions/site)."""
        return expm(self.Q * t)


def _gtr_model(name, exch, pi) -> SubstitutionModel:
    """Build a GTR-family model from 6 exchangeabilities ``[AC,AG,AT,CG,CT,GT]`` and freqs ``pi``.

    ``Q_ij = exch_ij · pi_j`` (i≠j), scaled so the expected rate ``-Σ pi_i Q_ii = 1`` (branch
    lengths are then in substitutions/site).
    """
    pi = np.asarray(pi, dtype=float)
    if pi.shape != (4,) or pi.min() < 0 or not np.isclose(pi.sum(), 1.0):
        raise ValueError(f"stationary frequencies must be 4 non-negative values summing to 1, got {pi}")
    ac, ag, at, cg, ct, gt = (float(x) for x in exch)
    if min(ac, ag, at, cg, ct, gt) < 0:
        raise ValueError("exchangeability rates must be non-negative")
    S = np.array([[0, ac, ag, at], [ac, 0, cg, ct], [ag, cg, 0, gt], [at, ct, gt, 0]])
    Q = S * pi[None, :]
    np.fill_diagonal(Q, -Q.sum(axis=1))
    scale = -(pi * np.diag(Q)).sum()
    if scale <= 0:
        raise ValueError("degenerate substitution model (zero rate)")
    return SubstitutionModel(name, Q / scale, pi)


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


_MODELS = {"jc69": jc69, "k80": k80, "hky85": hky85, "gtr": gtr}


def make_model(name: str, *, kappa: float = 2.0, freqs=None, rates=None) -> SubstitutionModel:
    """Construct a model by name (``jc69``/``k80``/``hky85``/``gtr``) with the relevant parameters."""
    name = name.lower()
    if name == "jc69":
        return jc69()
    if name == "k80":
        return k80(kappa)
    if name == "hky85":
        return hky85(kappa, freqs or (0.25, 0.25, 0.25, 0.25))
    if name == "gtr":
        return gtr(rates or (1, 1, 1, 1, 1, 1), freqs or (0.25, 0.25, 0.25, 0.25))
    raise ValueError(f"unknown substitution model {name!r} (choose from {sorted(_MODELS)})")


# --------------------------------------------------------------------------- #
# Across-site rate heterogeneity (discrete Gamma)
# --------------------------------------------------------------------------- #
@dataclass
class GammaRates:
    """Discrete-Gamma across-site rates: ``k`` equal-probability categories, mean rate 1."""

    shape: float
    k: int = 4

    def __post_init__(self):
        if self.shape <= 0 or self.k < 1:
            raise ValueError("Gamma shape must be > 0 and k >= 1")
        from scipy.stats import gamma as _g
        mids = _g.ppf((np.arange(self.k) + 0.5) / self.k, a=self.shape, scale=1.0 / self.shape)
        self.rates = mids / mids.mean()


# --------------------------------------------------------------------------- #
# Sequence <-> integer state coding
# --------------------------------------------------------------------------- #
def encode(seq: str, rng: np.random.Generator | None = None,
           pi: np.ndarray | None = None) -> np.ndarray:
    """Map a DNA string to integer states; ambiguous/other letters draw from ``pi`` (or uniform)."""
    states = np.empty(len(seq), dtype=np.int8)
    unknown = []
    for i, ch in enumerate(seq.upper()):
        c = _CODE.get(ch, -1)
        if c < 0:
            unknown.append(i)
            c = 0
        states[i] = c
    if unknown:
        rng = rng or np.random.default_rng()
        p = pi if pi is not None else np.full(4, 0.25)
        states[unknown] = rng.choice(4, size=len(unknown), p=p)
    return states


def decode(states: np.ndarray) -> str:
    """Map integer states back to a DNA string."""
    arr = np.asarray(states)
    return "".join(BASES[i] for i in arr)


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
    Q, pi = model.Q, model.stationary
    if root_seq is not None:
        root_states = encode(root_seq, rng, pi)
        length = root_states.shape[0]
    else:
        if length is None:
            raise ValueError("give either root_seq or length")
        root_states = rng.choice(4, size=length, p=pi).astype(np.int8)

    if gamma is not None:
        site_cat = rng.integers(gamma.k, size=length)
    out: dict = {}
    pcache: dict = {}

    def p_for(t: float) -> np.ndarray:
        key = round(float(t), 12)
        P = pcache.get(key)
        if P is None:
            P = expm(Q * key)
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
    return {gid: decode(s) for gid, s in out.items()}


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
