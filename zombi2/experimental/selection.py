"""Protein-language-model-guided selection over gene-family sequence evolution (**P1: frozen**).

ZOMBI2 evolves protein-coding sequences under a neutral substitution model, so a long branch
drifts a protein into non-protein noise. This module overlays **selection** on that neutral
process using a protein language model (e.g. ESM2) as an empirical fitness landscape, framed as a
**mutation--selection** model:

    substitution  =  mutation (the neutral model)  x  selection (the language model)

`beta` is the one knob: the strength of purifying selection (0 = neutral drift; large = strong
constraint). It plays the role of the population-scaled selection coefficient and is what drives
emergent dN/dS -- which we measure as an *output*, not calibrate to.

Two modes share one engine (see the design doc):

* **frozen** (this P1): call the critic **once** on the root protein to read a per-site amino-acid
  preference, bake it into each site's own substitution process, then evolve in one ordinary pass
  down the gene tree. Sites are independent -- no epistasis, no particle filter, cheap and
  closed-form. Each site is a **Halpern--Bruno / Sella--Hirsh mutation--selection** process: its
  rate is the base (neutral) model's rate times the fixation-probability factor
  ``h(F_b - F_a)``, with ``h(x) = x / (1 - e^{-x})`` and scaled fitness ``F = beta·ln(preference)``.
  The stationary distribution is ``pi_target ∝ pi_mut · preference**beta``, and at ``beta == 0`` the
  kernel reduces **exactly** to the base model (``h(0) = 1``). Here the amino-acid base model plays
  the role of the neutral process; the literal codon-level Halpern--Bruno with a nucleotide mutation
  backbone (and the dN/dS it yields) is P2.
* **live** (P3): a particle filter that re-scores the *current* sequence at time-slices to capture
  epistasis. Not implemented here.

Nothing in this module imports ``torch``/``esm`` at module load -- only :class:`ESM2Critic`
imports them, lazily, in its constructor. So ``import zombi2.experimental.selection`` and the whole
frozen path work without the optional ``zombi2[selection]`` dependencies.
"""
from __future__ import annotations

import math
from abc import ABC, abstractmethod

import numpy as np

from zombi2.experimental import warn_experimental
from zombi2.sequences.models import AMINO_ACIDS, SubstitutionModel, decode, encode, lg

__all__ = ["Critic", "FixedProfileCritic", "ESM2Critic", "PLMSelection"]


# --------------------------------------------------------------------------- #
# Critic: the pluggable language-model interface
# --------------------------------------------------------------------------- #
class Critic(ABC):
    """A pluggable protein-language-model interface. ESM2 is the first implementation;
    ProtT5 / ESM-C / a hand-supplied profile all slot in behind the same two methods.

    Both methods speak the 20-letter amino-acid alphabet in ZOMBI2's PAML order
    (:data:`~zombi2.sequences.models.AMINO_ACIDS`).
    """

    @abstractmethod
    def profile(self, seq: str) -> np.ndarray:
        """Per-site amino-acid preference for ``seq``: an ``(len(seq), 20)`` array whose rows are
        probability distributions over the 20 amino acids. Used by the **frozen** mode."""

    def score(self, seqs: list[str]) -> np.ndarray:
        """Per-sequence fitness (e.g. mean pseudo-log-likelihood per residue), one value per input.
        A whole-sequence readout for analysis/validation; neither the frozen nor the live overlay
        needs it (both drive selection through :meth:`profile`). Optional for profile-only critics."""
        raise NotImplementedError("this Critic implements profile() only, not score()")


class FixedProfileCritic(Critic):
    """A critic that returns a *fixed* per-site preference, ignoring the input sequence.

    Useful for injecting a known landscape (tests, ablations, or user-supplied preferences without a
    language model). ``profile`` is an ``(L, 20)`` array over :data:`AMINO_ACIDS`; rows are
    normalised on construction.
    """

    def __init__(self, profile: np.ndarray):
        warn_experimental("FixedProfileCritic")
        p = np.asarray(profile, dtype=float)
        if p.ndim != 2 or p.shape[1] != 20:
            raise ValueError(f"profile must be (L, 20) over the 20 amino acids, got {p.shape}")
        if (p < 0).any() or not np.isfinite(p).all():
            raise ValueError("profile entries must be finite and non-negative")
        rowsum = p.sum(1, keepdims=True)
        if (rowsum <= 0).any():
            raise ValueError("every profile row must have positive total mass")
        self._profile = p / rowsum

    def profile(self, seq: str) -> np.ndarray:
        if len(seq) != self._profile.shape[0]:
            raise ValueError(f"sequence length {len(seq)} != fixed profile length "
                             f"{self._profile.shape[0]}")
        return self._profile


class ESM2Critic(Critic):
    """ESM2 as the critic. Imports ``torch`` and ``esm`` **lazily** in ``__init__`` so the rest of
    this module has no heavy dependency; install them with ``pip install 'zombi2[selection]'``.

    Weights are downloaded and cached by ``esm`` on first use.
    """

    def __init__(self, model_name: str = "esm2_t6_8M_UR50D", device: str = "cpu"):
        warn_experimental("ESM2Critic")
        try:
            import esm
            import torch
        except ImportError as exc:  # pragma: no cover - exercised only without the extra installed
            raise ImportError(
                "ESM2Critic needs the optional protein-language-model dependencies; install them "
                "with:  pip install 'zombi2[selection]'"
            ) from exc
        try:
            loader = getattr(esm.pretrained, model_name)
        except AttributeError as exc:
            raise ValueError(f"unknown ESM2 model {model_name!r}") from exc
        self._torch = torch
        self._model, self._alphabet = loader()
        self._model.eval().to(device)
        self._bc = self._alphabet.get_batch_converter()
        self._device = device
        # ESM2 vocab columns for the 20 amino acids, in ZOMBI2's AMINO_ACIDS order
        self._cols = [self._alphabet.get_idx(a) for a in AMINO_ACIDS]

    def profile(self, seq: str) -> np.ndarray:
        torch = self._torch
        with torch.no_grad():
            _, _, toks = self._bc([("x", seq)])
            logits = self._model(toks.to(self._device))["logits"][0]        # (L+2, vocab)
            p = torch.softmax(logits, dim=-1)[1:len(seq) + 1][:, self._cols]  # (L, 20)
            arr = np.array(p.tolist())                                        # avoid torch<->numpy bridge
        return arr / arr.sum(1, keepdims=True)

    def score(self, seqs: list[str]) -> np.ndarray:
        torch = self._torch
        if not seqs:
            return np.zeros(0)
        with torch.no_grad():
            _, _, toks = self._bc([(str(i), s) for i, s in enumerate(seqs)])
            logp = torch.log_softmax(self._model(toks.to(self._device))["logits"], dim=-1)
            out = [logp[i, torch.arange(1, len(s) + 1), toks[i, 1:len(s) + 1]].mean().item()
                   for i, s in enumerate(seqs)]
        return np.array(out)


# --------------------------------------------------------------------------- #
# Frozen Halpern-Bruno mutation-selection kernel
# --------------------------------------------------------------------------- #
def _site_targets(profile: np.ndarray, model: SubstitutionModel, beta: float) -> np.ndarray:
    """Per-site target stationary distributions ``pi_target_i ∝ pi_mut · preference_i**beta``.

    Computed in **log space** with a per-row max shift so a whole row can never underflow to 0 even
    at large ``beta`` (which would otherwise give NaN). ``pi_mut`` is the base model's stationary. At
    ``beta == 0`` every row is ``pi_mut``. A small floor keeps every entry strictly positive (the
    reversible eigendecomposition needs ``pi > 0``). Returns ``(L, 20)``.
    """
    pref = np.clip(np.asarray(profile, dtype=float), 1e-12, None)
    logt = np.log(model.stationary)[None, :] + float(beta) * np.log(pref)
    logt -= logt.max(1, keepdims=True)
    t = np.exp(logt)
    t /= t.sum(1, keepdims=True)
    t = np.clip(t, 1e-12, None)
    return t / t.sum(1, keepdims=True)


def _hb_fixation(dF: np.ndarray) -> np.ndarray:
    """Halpern--Bruno / Sella--Hirsh relative fixation factor ``h(x) = x / (1 - e^{-x})``, with
    ``h(0) = 1``. Evaluated stably for both signs (rewritten as ``x·e^x/(e^x - 1)`` where ``x < 0``)."""
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        pos = dF / (1.0 - np.exp(-dF))                 # stable for dF > 0
        ex = np.exp(dF)
        neg = dF * ex / (ex - 1.0)                      # stable for dF < 0 (no overflow)
    return np.where(np.abs(dF) < 1e-9, 1.0, np.where(dF > 0.0, pos, neg))


def _site_model(base_Q: np.ndarray, pi_target: np.ndarray, pi_mut: np.ndarray,
                name: str, alphabet: str) -> SubstitutionModel:
    """One per-site Halpern--Bruno model: the neutral (base) rate ``mu_ab`` times the fixation factor
    ``h(F_b - F_a)`` for scaled fitness ``F = ln(pi_target / pi_mut)``. Reversible with stationary
    ``pi_target`` (detailed balance holds); renormalised to mean rate 1. At ``F == 0`` (beta=0) it is
    exactly the base model, since ``h(0) = 1``."""
    logF = np.log(pi_target) - np.log(pi_mut)          # scaled fitness, up to an additive constant
    dF = logF[None, :] - logF[:, None]                 # dF[a, b] = F_b - F_a
    Q = np.array(base_Q, dtype=float)
    np.fill_diagonal(Q, 0.0)                           # keep only the off-diagonal neutral rates mu_ab
    Q = Q * _hb_fixation(dF)
    np.fill_diagonal(Q, -Q.sum(1))
    scale = -(pi_target * np.diag(Q)).sum()
    if scale <= 0:
        raise ValueError("degenerate per-site model (zero substitution rate)")
    return SubstitutionModel(name, Q / scale, pi_target, alphabet)


def _site_models(targets: np.ndarray, base_Q: np.ndarray, pi_mut: np.ndarray,
                 alphabet: str) -> list[SubstitutionModel]:
    """One Halpern--Bruno model per site (stationary ``targets[i]``, neutral backbone ``base_Q``)."""
    return [_site_model(base_Q, targets[i], pi_mut, f"HB[{i}]", alphabet)
            for i in range(targets.shape[0])]


def _jump(site_models: list[SubstitutionModel], parent: np.ndarray, t: float,
          rng: np.random.Generator) -> np.ndarray:
    """Sample end states after branch length ``t`` -- each site under its own per-site model."""
    r = rng.random(parent.shape[0])
    out = np.empty(parent.shape[0], dtype=np.int8)
    for i in range(parent.shape[0]):
        c = np.cumsum(site_models[i].p_matrix(t)[parent[i]])
        c /= c[-1]              # renormalise: P(t) rows can be ~1e-5 off for a near-delta stationary
        out[i] = np.searchsorted(c, r[i], side="right")
    return out


def _evolve_frozen(root, subst: dict, site_models: list[SubstitutionModel],
                   root_states: np.ndarray, rng: np.random.Generator) -> dict:
    """Evolve down the gene tree with per-site models; return ``{node.gid: states}`` for every node.

    Mirrors :func:`~zombi2.sequences.models.evolve_on_tree` (same split rule: both children of a node
    inherit its committed state), but each branch applies each site's own substitution process.
    """
    out: dict = {}

    def visit(node, parent_states):
        t = float(subst.get(node, 0.0))
        states = parent_states if t <= 1e-12 else _jump(site_models, parent_states, t, rng)
        out[node.gid] = states
        for child in node.children:
            visit(child, states)

    visit(root, root_states)
    return out


# --------------------------------------------------------------------------- #
# PLMSelection: the user-facing overlay
# --------------------------------------------------------------------------- #
class PLMSelection:
    """Overlay language-model-guided selection on gene-family sequence evolution.

    Parameters
    ----------
    critic:
        A :class:`Critic` (e.g. :class:`ESM2Critic`, or :class:`FixedProfileCritic`).
    mode:
        ``"frozen"`` -- read the per-site preference **once** on the root, so sites are independent
        (no epistasis). ``"live"`` -- the **epistatic** mode: re-read the preference from the *current*
        sequence every ``refresh`` substitutions/site as it drifts, so each site feels the others'
        changes. Both use the same Halpern--Bruno kernel; live just refreshes it.
    beta:
        Selection strength (a finite value ``>= 0``). ``0`` reduces to the neutral base model.
    model:
        The neutral (mutation) amino-acid model; defaults to :func:`~zombi2.sequences.models.lg`.
        Must be a 20-state protein model.
    refresh:
        Live mode only: how often (in substitutions/site, *along each lineage*) to re-evaluate the
        critic on the current sequence. Smaller = finer epistasis + more critic calls; ``inf`` never
        refreshes (identical to ``frozen``). Default ``0.25``.
    """

    def __init__(self, critic: Critic, *, mode: str = "frozen", beta: float = 1.0,
                 model: SubstitutionModel | None = None, refresh: float = 0.25):
        warn_experimental("PLMSelection")
        if mode not in ("frozen", "live"):
            raise ValueError(f"mode must be 'frozen' or 'live', got {mode!r}")
        try:
            beta = float(beta)
        except (TypeError, ValueError):
            raise ValueError(f"beta must be a number, got {beta!r}") from None
        if not math.isfinite(beta) or beta < 0:
            raise ValueError(f"beta must be a finite value >= 0, got {beta}")
        try:                               # validated for both modes so a bad value can never linger
            refresh = float(refresh)
        except (TypeError, ValueError):
            raise ValueError(f"refresh must be a number, got {refresh!r}") from None
        if not refresh > 1e-9:             # allows inf (= never refresh); rejects <=0, nan, sub-epsilon
            raise ValueError(f"refresh must be > 1e-9 substitutions/site, got {refresh}")
        self.critic = critic
        self.mode = mode
        self.beta = beta
        self.refresh = refresh
        self.model = model if model is not None else lg()
        if self.model.k != 20 or self.model.alphabet != AMINO_ACIDS:
            raise ValueError("PLMSelection needs a 20-state amino-acid model (e.g. lg())")

    def site_targets(self, root_protein: str) -> np.ndarray:
        """The per-site target stationary distributions the frozen model would evolve toward for
        ``root_protein`` -- ``(len, 20)`` over :data:`AMINO_ACIDS`. Handy for inspection/plots."""
        prof = self.critic.profile(root_protein)
        self._check_profile(prof, root_protein)
        return _site_targets(prof, self.model, self.beta)

    def _build_models(self, protein: str) -> list[SubstitutionModel]:
        """The per-site Halpern--Bruno models for ``protein`` (validates the critic's profile shape)."""
        prof = self.critic.profile(protein)
        self._check_profile(prof, protein)
        targets = _site_targets(prof, self.model, self.beta)
        return _site_models(targets, self.model.Q, self.model.stationary, self.model.alphabet)

    def evolve_family(self, root, subst: dict, root_protein: str, *,
                      rng: np.random.Generator) -> dict:
        """Evolve one gene family's protein down its node tree.

        ``root`` is a tree node (any object with ``.gid`` and ``.children`` -- e.g. a reconciliation
        ``_Node``); ``subst`` maps each node to the substitution length of the branch ending at it
        (as :func:`~zombi2.sequences.evolution._annotate` produces). In ``frozen`` mode the per-site
        kernel is read once on the root; in ``live`` mode it is re-read from the current sequence
        every ``refresh`` substitutions/site (epistasis). Returns ``{node.gid: protein}``.
        """
        alphabet = self.model.alphabet
        root_states = encode(root_protein, rng, self.model.stationary, alphabet)
        # build the kernel from the *concrete* root sequence actually placed at the root (resolving any
        # ambiguity codes), so frozen and live see the same sequence -- keeping refresh=inf exactly frozen
        root_seq = decode(root_states, alphabet)
        if self.mode == "frozen":
            evolved = _evolve_frozen(root, subst, self._build_models(root_seq), root_states, rng)
        else:
            evolved = self._evolve_live(root, subst, root_states, rng)
        return {gid: decode(st, alphabet) for gid, st in evolved.items()}

    def _evolve_live(self, root, subst: dict, root_states: np.ndarray,
                     rng: np.random.Generator) -> dict:
        """Live/epistatic descent: evolve under the Halpern--Bruno kernel, re-reading the critic on the
        current sequence every ``self.refresh`` substitutions/site along each lineage (the profile
        therefore tracks the drifting context). Children inherit the parent's committed sequence and
        its current kernel; each then refreshes on its own schedule."""
        alphabet = self.model.alphabet
        out: dict = {}

        def build(states: np.ndarray) -> list:
            return self._build_models(decode(states, alphabet))

        def visit(node, parent_states, models, accrued):
            states = parent_states
            remaining = float(subst.get(node, 0.0))
            while remaining > 1e-12:
                step = min(self.refresh - accrued, remaining)
                states = _jump(models, states, step, rng)
                accrued += step
                remaining -= step
                if accrued >= self.refresh - 1e-12:        # drifted a full refresh interval -> re-read
                    models = build(states)
                    accrued = 0.0
            out[node.gid] = states
            for child in node.children:
                visit(child, states, models, accrued)

        visit(root, root_states, build(root_states), 0.0)
        return out

    def evolve_families(self, node_trees: dict, root_seqs: dict, *,
                        seed: int | None = None) -> dict:
        """Evolve many independent families. ``node_trees`` is a mapping ``family -> {"complete":
        (root, subst) | None, "extant": (root, subst) | None}`` (as
        :meth:`~zombi2.sequences.evolution.SequenceEvolution.scale_families_trees` returns), and
        ``root_seqs`` maps ``family -> root protein``. Returns ``{family: {gid: protein}}``.

        A single RNG is threaded across families (in ``root_seqs`` iteration order), so the result is
        reproducible for a given ``seed`` but families are not order-independent draws.
        """
        rng = np.random.default_rng(seed)
        out: dict = {}
        for fam, root_protein in root_seqs.items():
            trees = node_trees.get(fam)
            if not trees:
                continue
            entry = trees.get("extant") or trees.get("complete")
            if entry is None:
                continue
            root_node, subst = entry
            out[fam] = self.evolve_family(root_node, subst, root_protein, rng=rng)
        return out

    def _check_profile(self, prof: np.ndarray, seq: str) -> None:
        if prof.shape != (len(seq), self.model.k):
            raise ValueError(f"critic profile shape {prof.shape} != ({len(seq)}, {self.model.k})")
