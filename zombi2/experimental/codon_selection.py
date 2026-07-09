"""Codon-level Halpern--Bruno mutation-selection (**P2: frozen, codon**).

P1 evolves amino-acid sequences. This module evolves the underlying **coding DNA** one codon at a
time: mutation happens at the *nucleotide* level (a nucleotide substitution model), while selection
acts on the *amino acid* the codon encodes (a protein language model's per-site preference). That
split is biologically correct and buys the headline result for free:

* **synonymous** changes leave the amino acid unchanged, so their fixation factor is ``h(0) = 1`` --
  they evolve neutrally;
* **non-synonymous** changes are scrutinised by the language model;

so **dN/dS emerges as an output** of the model rather than being imposed -- exactly the signal the
ADH1/dN-dS work cares about. :meth:`CodonSelection.dnds` returns the model's expected genome-wide
``omega`` (which drops below 1 as ``beta`` rises).

Each codon site is a 61-state (sense-codon) Halpern--Bruno mutation-selection process, built exactly
like P1 but over codons: ``Q(c1->c2) = mu(c1->c2) · h(F_c2 - F_c1)`` for single-nucleotide codon
neighbours, ``F(codon) = beta · ln(preference[aa(codon)])``, stationary ``pi_target ∝ pi_mut · ...``.
Nonsense mutations (to a stop codon) are treated as lethal and excluded, per standard codon models;
the coding sequence is taken 5'->3', in frame, stop-free (strand/GFF wiring is the next slice).

Like P1, nothing here imports torch/esm at module load.
"""
from __future__ import annotations

import math

import numpy as np

from zombi2.experimental import warn_experimental
from zombi2.experimental.selection import (
    Critic, FixedProfileCritic, _ExpmSite, _evolve_frozen, _hb_fixation,
)
from zombi2.sequences.models import AMINO_ACIDS, BASES, SubstitutionModel, hky85

__all__ = ["CodonSelection", "GENETIC_CODE", "SENSE_CODONS", "calibrate_beta", "translate"]

# --------------------------------------------------------------------------- #
# The standard genetic code
# --------------------------------------------------------------------------- #
_TCAG = "TCAG"
_CODON_ORDER = [x + y + z for x in _TCAG for y in _TCAG for z in _TCAG]
# NCBI transl_table 1, in TCAG codon order:
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = dict(zip(_CODON_ORDER, _AAS))
SENSE_CODONS = [c for c in _CODON_ORDER if GENETIC_CODE[c] != "*"]      # 61
STOP_CODONS = frozenset(c for c in _CODON_ORDER if GENETIC_CODE[c] == "*")
_CODON_INDEX = {c: i for i, c in enumerate(SENSE_CODONS)}
_AA_INDEX = {a: i for i, a in enumerate(AMINO_ACIDS)}
# amino-acid index (into AMINO_ACIDS) encoded by each sense codon
_CODON_AA = np.array([_AA_INDEX[GENETIC_CODE[c]] for c in SENSE_CODONS])


def translate(dna: str) -> str:
    """Translate an in-frame coding DNA string to its amino-acid sequence (``*`` for stop)."""
    if len(dna) % 3:
        raise ValueError(f"coding sequence length {len(dna)} is not a multiple of 3")
    return "".join(GENETIC_CODE[dna[i:i + 3].upper()] for i in range(0, len(dna), 3))


def _encode_codons(dna: str) -> np.ndarray:
    if len(dna) % 3:
        raise ValueError(f"coding sequence length {len(dna)} is not a multiple of 3")
    out = np.empty(len(dna) // 3, dtype=np.int8)
    for i in range(0, len(dna), 3):
        c = dna[i:i + 3].upper()
        idx = _CODON_INDEX.get(c)
        if idx is None:
            what = "stop codon" if c in STOP_CODONS else "invalid codon"
            raise ValueError(f"coding sequence contains a {what} {c!r} at position {i}")
        out[i // 3] = idx
    return out


def _decode_codons(states: np.ndarray) -> str:
    return "".join(SENSE_CODONS[i] for i in states)


# --------------------------------------------------------------------------- #
# Neutral codon mutation from a nucleotide model
# --------------------------------------------------------------------------- #
def _codon_mutation(nuc: SubstitutionModel):
    """From a nucleotide model, the neutral codon exchange matrix ``mu`` (61x61, single-nt neighbours
    only) and neutral codon stationary ``pi_mut`` (product of nucleotide frequencies)."""
    nb = {b: i for i, b in enumerate(BASES)}
    n = len(SENSE_CODONS)
    mu = np.zeros((n, n))
    for i, c1 in enumerate(SENSE_CODONS):
        for j, c2 in enumerate(SENSE_CODONS):
            if i == j:
                continue
            diffs = [k for k in range(3) if c1[k] != c2[k]]
            if len(diffs) == 1:
                k = diffs[0]
                mu[i, j] = nuc.Q[nb[c1[k]], nb[c2[k]]]
    pi = np.array([nuc.stationary[nb[c[0]]] * nuc.stationary[nb[c[1]]] * nuc.stationary[nb[c[2]]]
                   for c in SENSE_CODONS])
    return mu, pi / pi.sum()


def _codon_targets(pi_mut: np.ndarray, F: np.ndarray) -> np.ndarray:
    """Codon target stationary ``pi_target ∝ pi_mut · e^F``, in log space (max-shifted).

    Entries for very unfit codons may underflow to exactly 0 -- that is the *correct* near-delta
    stationary and it is safe here: the ``expm`` kernel never needs ``pi > 0`` (unlike the P1 eigh),
    and the ``scale`` and dN/dS weightings only ever *sum over* ``pi_target``, so zeros drop out.
    (An earlier 1e-12 floor -- a P1 holdover -- made :meth:`CodonSelection.dnds` bottom out and then
    grow with ``beta`` instead of decaying to 0; removing it restores the monotone decay.)
    """
    logt = np.log(pi_mut) + F
    logt -= logt.max()
    t = np.exp(logt)
    return t / t.sum()


def _codon_site_model(mu: np.ndarray, pi_mut: np.ndarray, aa_pref: np.ndarray,
                      beta: float) -> _ExpmSite:
    """One codon Halpern--Bruno model for a site with amino-acid preference ``aa_pref`` (20,).

    All sites are scaled by the **same** global neutral rate, so a branch of length ``t`` is ``t``
    expected substitutions per **nucleotide site** under neutrality. The scale is the stationary-mean
    codon exit rate ``Σ pi_mut · (neutral exit rate)`` divided by 3 (three nucleotide sites per codon):
    this makes the codon clock a *per-nucleotide-site* clock, so coding and non-coding regions that
    share one branch-length scale (e.g. the nucleotide genome model) diverge at the same neutral rate
    -- without the /3 a codon block under-diverges ~3x versus a nucleotide block. A conserved site
    therefore accrues *fewer* real substitutions (its own mean rate is below the global rate) --
    emergent site-rate heterogeneity, and the visible face of dN/dS < 1. (A per-site 'mean rate 1'
    scaling would erase that heterogeneity and, for near-delta stationaries, underflow to 0.)
    """
    F = beta * np.log(np.clip(aa_pref, 1e-12, None))[_CODON_AA]      # fitness per codon
    pi_target = _codon_targets(pi_mut, F)
    dF = F[None, :] - F[:, None]                                     # dF[a, b] = F_b - F_a
    Q = mu * _hb_fixation(dF)                                        # off-diagonal (mu diagonal is 0)
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(1))
    # /3: convert the per-codon neutral exit rate to a per-nucleotide-site rate (3 nt per codon), so the
    # codon clock matches the nucleotide clock when both share a branch-length scale (genome models).
    scale = float((pi_mut * mu.sum(1)).sum()) / 3.0                 # neutral mean rate per nt site (global)
    if scale <= 0:
        raise ValueError("degenerate mutation model (zero substitution rate)")
    return _ExpmSite(Q / scale, pi_target)


# --------------------------------------------------------------------------- #
# CodonSelection: the user-facing overlay
# --------------------------------------------------------------------------- #
class CodonSelection:
    """Overlay language-model-guided selection on **coding-DNA** evolution (codon-level, frozen).

    Mutation comes from ``nuc_model`` (a nucleotide model; default :func:`~zombi2.sequences.models.hky85`);
    selection comes from ``critic`` (a :class:`~zombi2.experimental.selection.Critic`, which scores the
    *translated* protein). ``beta`` is the selection strength; ``0`` is a neutral codon model.
    """

    def __init__(self, critic: Critic, *, beta: float = 1.0,
                 nuc_model: SubstitutionModel | None = None):
        warn_experimental("CodonSelection")
        try:
            beta = float(beta)
        except (TypeError, ValueError):
            raise ValueError(f"beta must be a number, got {beta!r}") from None
        if not math.isfinite(beta) or beta < 0:
            raise ValueError(f"beta must be a finite value >= 0, got {beta}")
        self.critic = critic
        self.beta = beta
        self.nuc = nuc_model if nuc_model is not None else hky85()
        if self.nuc.alphabet != BASES:
            raise ValueError("nuc_model must be a nucleotide (ACGT) model, e.g. hky85()")
        self.mu, self.pi_mut = _codon_mutation(self.nuc)
        # single-nt neighbours split into synonymous / non-synonymous (for dN/dS)
        neighbour = self.mu > 0
        same_aa = _CODON_AA[:, None] == _CODON_AA[None, :]
        self._syn = neighbour & same_aa
        self._nonsyn = neighbour & ~same_aa

    def _profile(self, protein: str) -> np.ndarray:
        prof = self.critic.profile(protein)
        if prof.shape != (len(protein), 20):
            raise ValueError(f"critic profile shape {prof.shape} != ({len(protein)}, 20)")
        return prof

    def _site_models(self, protein: str) -> list[_ExpmSite]:
        prof = self._profile(protein)
        return [_codon_site_model(self.mu, self.pi_mut, prof[i], self.beta)
                for i in range(len(protein))]

    def evolve_coding_family(self, root, subst: dict, root_dna: str, *,
                             rng: np.random.Generator) -> dict:
        """Evolve one gene family's **coding DNA** down its node tree. ``root_dna`` is 5'->3', in
        frame, stop-free. Returns ``{node.gid: coding_dna}`` for every node."""
        site_models = self._site_models(translate(root_dna))
        root_states = _encode_codons(root_dna)
        evolved = _evolve_frozen(root, subst, site_models, root_states, rng)
        return {gid: _decode_codons(st) for gid, st in evolved.items()}

    def evolve_coding_families(self, node_trees: dict, root_dnas: dict, *,
                               seed: int | None = None) -> dict:
        """Evolve many independent families (see
        :meth:`~zombi2.experimental.selection.PLMSelection.evolve_families`); ``root_dnas`` maps
        ``family -> coding DNA``. Returns ``{family: {gid: coding_dna}}``."""
        rng = np.random.default_rng(seed)
        out: dict = {}
        for fam, root_dna in root_dnas.items():
            trees = node_trees.get(fam)
            if not trees:
                continue
            entry = trees.get("extant") or trees.get("complete")
            if entry is None:
                continue
            root_node, subst = entry
            out[fam] = self.evolve_coding_family(root_node, subst, root_dna, rng=rng)
        return out

    def dnds(self, protein: str) -> float:
        """The model's expected genome-wide **dN/dS** (``omega``) for a protein of this length and
        preference. Synonymous flux is neutral (``dS == 1``), so ``omega == dN`` = the stationary-
        averaged non-synonymous fixation factor. ``beta == 0`` gives ``omega == 1``; larger ``beta``
        drives it below 1 (purifying selection)."""
        prof = self._profile(protein)
        num = den = 0.0
        for i in range(len(protein)):
            F = self.beta * np.log(np.clip(prof[i], 1e-12, None))[_CODON_AA]
            pi_target = _codon_targets(self.pi_mut, F)
            h = _hb_fixation(F[None, :] - F[:, None])
            w = pi_target[:, None] * self.mu                 # pi_i · mu_ij (neutral nonsyn flux)
            num += float((w * h * self._nonsyn).sum())
            den += float((w * self._nonsyn).sum())
        return num / den

    def dnds_syn_check(self, protein: str) -> float:
        """dS (should be 1.0 for any beta -- synonymous changes are neutral). Exposed for validation."""
        prof = self._profile(protein)
        num = den = 0.0
        for i in range(len(protein)):
            F = self.beta * np.log(np.clip(prof[i], 1e-12, None))[_CODON_AA]
            pi_target = _codon_targets(self.pi_mut, F)
            h = _hb_fixation(F[None, :] - F[:, None])
            w = pi_target[:, None] * self.mu
            num += float((w * h * self._syn).sum())
            den += float((w * self._syn).sum())
        return num / den


def calibrate_beta(critic: Critic, protein: str, target_dnds: float, *,
                   nuc_model: SubstitutionModel | None = None, hi: float = 64.0,
                   tol: float = 1e-4, max_iter: int = 64) -> float:
    """Find the selection strength ``beta`` whose codon model gives expected dN/dS ≈ ``target_dnds``
    (a value in ``(0, 1)``) for ``protein`` -- the usable inverse of :meth:`CodonSelection.dnds`, so a
    user can ask for a target ω instead of guessing ``beta``. The critic is queried **once** (its
    profile is reused across the search); ``dnds`` is monotone decreasing in ``beta``, so a bisection
    on ``[0, hi]`` converges to the *smallest* beta reaching the target. Raises if ``target_dnds`` is
    unreachable below ``hi``, or if the search does not reach ``tol`` within ``max_iter`` iterations
    (e.g. a target so small it sits on the dN/dS underflow plateau -- raise ``tol`` or lower it)."""
    if not 0.0 < target_dnds < 1.0:
        raise ValueError(f"target_dnds must be in (0, 1), got {target_dnds}")
    fixed = FixedProfileCritic(critic.profile(protein))          # one critic call, reused for every beta

    def omega(b: float) -> float:
        return CodonSelection(fixed, beta=b, nuc_model=nuc_model).dnds(protein)

    if omega(hi) > target_dnds:
        raise ValueError(f"target dN/dS {target_dnds} needs beta > {hi}; pass a larger hi=")
    lo, high = 0.0, float(hi)
    for _ in range(max_iter):
        mid = 0.5 * (lo + high)
        w = omega(mid)
        # accept only a genuine (positive) match: w underflowed to 0 means beta is too high, so keep
        # lowering it -- otherwise a target below the underflow floor returns an arbitrary plateau beta.
        if w > 0.0 and abs(w - target_dnds) <= tol:
            return mid
        lo, high = (mid, high) if w > target_dnds else (lo, mid)   # omega decreasing: too high -> raise beta
    mid = 0.5 * (lo + high)
    w = omega(mid)
    if abs(w - target_dnds) > tol:
        raise ValueError(f"calibrate_beta did not reach tol={tol} within max_iter={max_iter} "
                         f"(best dN/dS {w:.4g} at beta {mid:.4g}); raise max_iter or tol")
    return mid
