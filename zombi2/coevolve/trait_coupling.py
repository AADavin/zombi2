"""Trait-conditioned gene families — couple a phenotypic trait to gene-family dynamics.

This is the *forward* generator behind studies that read gene-family patterns as a record of
a trait's history (Davin 2025: timing the bacterial tree from the Great Oxidation Event). A
trait — aerobic/anaerobic, a habitat, an oxygen-tolerance level — evolves down the species
tree (:mod:`zombi2.traits`); a designated set of gene families then gains or loses membership
at rates that **depend on the local trait value**, so the resulting phylogenetic profile
carries a known, trait-linked signal.

**Direction.** This module implements *genes conditioned on a trait* (the trait is simulated
first, then fed in as a per-branch, per-time covariate). It is the general building block —
the reverse coupling (a trait read off gene content) and the fully joint co-evolutionary
model are limiting cases to layer on top later.

**Mechanism (retention).** For a *responsive* family ``i`` with coupling weight ``w_i`` and
local trait value ``s(branch, t)`` the loss rate is

    loss_i = base_loss · exp(-effect_loss · w_i · s),

so where the trait favours it (``w_i·s`` large) the family is retained, and where it does not
the family is purged. Gain is horizontal transfer (a trait-blind influx): a family flows in
and the trait-modulated loss then *selectively retains* it, which is what writes the
trait↔gene association into the
profiles. So the **net** gene content of a lineage tracks its trait even though the influx
itself is trait-blind. ``effect_gain`` optionally scales a lineage's transfer (HGT) activity
by its trait too; it is off by default.

**Trait as a covariate in time.** A *discrete* trait contributes its exact stochastic
character map — the per-branch ``(state, duration)`` segments — so the coupling is exact: the
segment boundaries become refresh points (:meth:`RateModel.refresh_times`) at which the
simulator recomputes the branch's rates. A *continuous* trait (Brownian/OU) is sub-segmented
into ``trait_steps`` pieces per branch with the value interpolated between the node endpoints
(a piecewise-constant covariate; a diffusion-bridge refinement can replace the interpolation
later).

Only the family-side rate model is new; the simulator, genome and output are untouched, so a
coupled run reuses the whole pipeline (profiles, gene trees, reconciliations) via
:meth:`TraitGeneResult.genomes`.
"""

from __future__ import annotations

import copy
import math
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field

import numpy as np

from zombi2.genomes.events import EventType
from zombi2.genomes.genome import Gene, UnorderedGenome
from zombi2.genomes.genome_sim import GenomeSimulator
from zombi2.genomes.profiles import ProfileMatrix, _natkey
from zombi2.genomes.rates import ModifiedRates, Rates
from zombi2.traits.models import TraitResult, simulate_traits
from zombi2.genomes.transfers import TransferModel
from zombi2.tree import Tree
from zombi2.coevolve.grammar import Scalar
from zombi2.coevolve.rate_bridge import CouplingModifier

#: Clamp on the loss/gain exponent so an extreme ``effect·w·s`` never overflows ``exp`` (or
#: drives the Gillespie loop into an instant loss/gain hot cycle).
_MAX_EXPONENT = 40.0


def _clamp(x: float) -> float:
    return max(-_MAX_EXPONENT, min(_MAX_EXPONENT, x))


def _is_scalar(value) -> bool:
    arr = np.asarray(value)
    return arr.ndim == 0 or arr.size == 1


# ═══════════════════════════════════════════════════════════════════════════════
# The trait trajectory: a (branch, time) -> value lookup, plus its change points
# ═══════════════════════════════════════════════════════════════════════════════
class TraitTrajectory:
    """The trait's value along every species branch, as a piecewise-constant function of time.

    Built from a :class:`~zombi2.traits.TraitResult` by :meth:`from_result`. Each branch (keyed
    by node name) holds a list of segment start times and the trait value on each segment; the
    value is looked up right-continuously (:meth:`value`), and the interior segment boundaries
    are the refresh points the simulator honours (:meth:`refresh_times`).
    """

    __slots__ = ("_starts", "_vals", "_boundaries", "_btimes", "_default")

    def __init__(self, starts, vals, boundaries, default):
        self._starts = starts          # branch name -> ascending list of segment start times
        self._vals = vals              # branch name -> aligned list of trait values
        self._boundaries = boundaries  # ascending [(time, branch_name)] interior jump points
        self._btimes = [b[0] for b in boundaries]
        self._default = float(default)

    def value(self, branch: str, t: float) -> float:
        """Trait value on ``branch`` at time ``t`` (right-continuous at segment boundaries)."""
        starts = self._starts.get(branch)
        if not starts:
            return self._default
        i = bisect_right(starts, t) - 1
        return self._vals[branch][i if i >= 0 else 0]

    def refresh_times(self, t0: float, t1: float) -> list[tuple[float, str]]:
        """The interior trait-change points strictly inside ``(t0, t1)`` — ``(time, branch)``."""
        lo = bisect_right(self._btimes, t0)
        hi = bisect_left(self._btimes, t1)
        return self._boundaries[lo:hi]

    @classmethod
    def from_result(cls, result: TraitResult, *, steps: int = 16, state_values=None) -> "TraitTrajectory":
        """Discretise a trait history into a per-branch trajectory.

        Discrete traits use their exact stochastic map (``result.history``); continuous traits
        are sub-segmented into ``steps`` pieces per branch, the value linearly interpolated
        between the parent and child node values and held constant across each piece.
        ``state_values`` optionally maps a discrete state index to a numeric value (default: the
        index itself, so a binary trait is 0/1).
        """
        if steps < 1:
            raise ValueError("trait_steps must be >= 1")
        tree = result.tree
        if not _is_scalar(result.node_values[tree.root]):
            raise ValueError("trait–gene coupling needs a univariate (scalar) trait; got a "
                             "multivariate one")

        def state_num(idx) -> float:
            return float(idx) if state_values is None else float(state_values[int(idx)])

        discrete = result.kind == "discrete"
        starts: dict[str, list[float]] = {}
        vals: dict[str, list[float]] = {}
        boundaries: list[tuple[float, str]] = []

        root_val = state_num(result.node_values[tree.root]) if discrete \
            else float(result.node_values[tree.root])

        for node in tree.nodes_preorder():
            if node.parent is None:
                continue
            b0, b1 = node.parent.time, node.time
            name = node.name
            if discrete and result.history is not None:
                segs = result.history.get(node) or []
                seg_starts, seg_vals, t = [], [], b0
                for state, dur in segs:
                    seg_starts.append(t)
                    seg_vals.append(state_num(state))
                    t += dur
                if not seg_starts:  # a zero-length branch, or no recorded history
                    seg_starts, seg_vals = [b0], [state_num(result.node_values[node])]
            else:  # continuous: sub-segment with interpolated, piecewise-constant values
                start_v = float(result.node_values[node.parent])
                end_v = float(result.node_values[node])
                dt = (b1 - b0) / steps
                seg_starts = [b0 + k * dt for k in range(steps)]
                seg_vals = [start_v + (end_v - start_v) * ((k + 0.5) / steps) for k in range(steps)]
            starts[name] = seg_starts
            vals[name] = seg_vals
            boundaries.extend((st, name) for st in seg_starts[1:])  # interior jumps only

        boundaries.sort(key=lambda x: x[0])
        return cls(starts, vals, boundaries, root_val)


# ═══════════════════════════════════════════════════════════════════════════════
# The coupling specification
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class TraitGeneCoupling:
    """A panel of ``n_families`` families, each with a trait-coupling weight, plus base rates.

    ``weights[i]`` is family ``i``'s coupling weight: ``0`` = inert (evolves at the base rates),
    positive = favoured by a high trait value, negative = favoured by a low one. Build the panel
    with :meth:`build` (which populates the weights from an explicit id list, a count, or a
    fraction); the families are named ``F0 .. F{n-1}`` (override with ``prefix``).

    Parameters
    ----------
    n_families  : panel size ``N``.
    weights     : length-``N`` array of per-family coupling weights.
    effect_loss : coupling strength on retention — loss scales by ``exp(-effect_loss·w_i·s)``.
    effect_gain : coupling strength on a lineage's transfer (HGT) activity — the transfer rate
                  scales by ``exp(effect_gain·s)`` (0 = field-blind gain, the default).
    base_loss   : baseline per-copy loss rate (the loss at ``w_i·s = 0``).
    transfer    : per-copy horizontal-transfer rate — the (field-blind) gain channel.
    duplication : per-copy duplication rate (trait-independent).
    origination : background rate of brand-new, uncoupled families (0 → closed panel).
    prefix      : family-id prefix; family ``i`` is ``f"{prefix}{i}"``.
    state_values: optional map from a discrete state index to a numeric trait value.
    """

    n_families: int
    weights: np.ndarray
    effect_loss: float = 1.0
    effect_gain: float = 0.0
    base_loss: float = 1.0
    transfer: float = 0.5
    duplication: float = 0.0
    origination: float = 0.0
    prefix: str = "F"
    state_values: np.ndarray | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if self.n_families <= 0:
            raise ValueError("n_families must be positive")
        self.weights = np.asarray(self.weights, dtype=float)
        if self.weights.shape != (self.n_families,):
            raise ValueError(f"weights must have shape ({self.n_families},), got {self.weights.shape}")
        for name in ("base_loss", "transfer", "duplication", "origination"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be >= 0")
        self.panel_ids: list[str] = [f"{self.prefix}{i}" for i in range(self.n_families)]
        self.index: dict[str, int] = {fam: i for i, fam in enumerate(self.panel_ids)}
        # id -> weight for the responsive families only (the hot loop reads this)
        self.weights_by_id: dict[str, float] = {
            self.panel_ids[i]: float(w) for i, w in enumerate(self.weights) if w != 0.0
        }
        if self.state_values is not None:
            self.state_values = np.asarray(self.state_values, dtype=float)

    def null(self, kind="neutral", **kwargs):
        """Decoupled **null** for the ``traits:genes`` arrow (trait → gene retention).

        ``"neutral"`` cuts both coupling channels (``effect_loss = effect_gain = 0``): every
        family, responsive or not, evolves at its base rates, so the panel no longer tracks the
        trait. The character-independent (``"cid"``) null is a **workflow** — drive the panel with
        a *hidden* trait while observing a *second, independent neutral trait* — not a model
        transform; the CLI ``--null cid`` builds it. See
        :doc:`the null-models guide </guide/coevolution_nulls>`.
        """
        kind = kind.lower()
        if kind == "neutral":
            m = copy.copy(self)
            m.effect_loss = 0.0
            m.effect_gain = 0.0
            return m
        if kind == "cid":
            raise TypeError(
                "the traits:genes CID null is a workflow (a hidden trait drives the panel, an "
                "independent neutral trait is observed), not a model transform; use the CLI "
                "`--null cid` — see docs/guide/coevolution_nulls.md")
        if kind == "timing":
            raise ValueError("traits:genes has no 'timing' null; use kind='neutral'")
        raise ValueError(f"unknown null kind {kind!r}; expected 'neutral'")

    @property
    def responsive_ids(self) -> list[str]:
        """The ids of the families with a non-zero coupling weight."""
        return list(self.weights_by_id)

    @property
    def n_responsive(self) -> int:
        return len(self.weights_by_id)

    # --- constructor -------------------------------------------------------
    @classmethod
    def build(cls, n_families: int, responsive, *, weight: float = 1.0, signed: bool = False,
              seed: int | None = None, rng: np.random.Generator | None = None,
              **rate_kw) -> "TraitGeneCoupling":
        """Build a coupling, choosing which families respond.

        ``responsive`` selects the responsive set:

        * an **int** ``k`` — ``k`` families chosen uniformly at random;
        * a **float** ``0 < f <= 1`` — a fraction ``f`` of the panel, chosen at random;
        * an **iterable of ids/indices** — those exact families (``"F3"`` or ``3``).

        Each responsive family gets weight ``+weight``; with ``signed=True`` the sign is
        randomised (so some families are favoured by a high trait value and others by a low
        one). Remaining rate parameters are passed through to the dataclass.
        """
        if rng is None:
            rng = np.random.default_rng(seed)
        w = np.zeros(n_families, dtype=float)
        idxs = _resolve_responsive(n_families, responsive, rng)
        for i in idxs:
            sign = 1.0 if not signed else (1.0 if rng.random() < 0.5 else -1.0)
            w[i] = sign * weight
        return cls(n_families=n_families, weights=w, **rate_kw)


def _resolve_responsive(n: int, responsive, rng) -> list[int]:
    """Turn the ``responsive`` selector (count / fraction / id list) into panel indices."""
    if isinstance(responsive, bool):
        raise TypeError("responsive must be an int, float, or iterable of ids — not bool")
    if isinstance(responsive, (int, np.integer)):
        k = min(int(responsive), n)
        return sorted(int(i) for i in rng.choice(n, size=k, replace=False)) if k > 0 else []
    if isinstance(responsive, float):
        if not 0.0 <= responsive <= 1.0:
            raise ValueError(f"a fractional responsive set must be in [0, 1], got {responsive}")
        k = int(round(responsive * n))
        return sorted(int(i) for i in rng.choice(n, size=k, replace=False)) if k > 0 else []
    # an explicit iterable of ids ("F3") or indices (3)
    out: list[int] = []
    for tok in responsive:
        if isinstance(tok, (int, np.integer)):
            idx = int(tok)
        else:
            s = str(tok).strip()
            digits = "".join(ch for ch in s if ch.isdigit())
            if not digits:
                raise ValueError(f"cannot read a family index from {tok!r} (want e.g. 'F3' or 3)")
            idx = int(digits)
        if not 0 <= idx < n:
            raise ValueError(f"responsive family index {idx} out of range for n_families={n}")
        out.append(idx)
    return sorted(set(out))


# ═══════════════════════════════════════════════════════════════════════════════
# The rate model
# ═══════════════════════════════════════════════════════════════════════════════
class TraitGeneRates(ModifiedRates):
    """The ``traits:genomes`` edge, compiled onto the grammar.

    A :class:`~zombi2.genomes.rates.ModifiedRates` over a :class:`~zombi2.genomes.rates.Rates`
    (``per="copy"``) base, with one :class:`~zombi2.coevolve.rate_bridge.CouplingModifier` per trait-coupled channel:

    * **retention** — a responsive family's per-copy loss is scaled by ``exp(-effect_loss·w_i·s)``
      (``s`` the local trait value), so a favoured family is kept and a disfavoured one purged; an
      inert family (weight 0) keeps ``base_loss·copies``;
    * **gain** — an optional field-blind ``TRANSFER`` scaled by ``exp(effect_gain·s)`` (off by
      default).

    Duplication and background origination come straight from the base. The trajectory's breakpoints
    are each modifier's :meth:`~zombi2.genomes.rates.Modifier.refresh_times`, so the simulator
    refreshes a branch exactly at a gene event or a trait change — all inherited from
    ``ModifiedRates`` (nothing re-implemented here). Emission is per-family on every channel via the
    ``ModifiedRates`` expansion — distributionally identical to the previous mixed
    per-family-loss / aggregate-gain emission (see ``docs/design/coevolve-grammar.md`` §6).
    """

    def __init__(self, coupling: TraitGeneCoupling, trajectory: TraitTrajectory):
        self.c = coupling
        self.traj = trajectory
        base = Rates(loss=coupling.base_loss, transfer=coupling.transfer,
                     duplication=coupling.duplication, origination=coupling.origination)  # per="copy"
        # retention: a responsive family sees w_i·s, so loss ×= exp(-effect_loss·w_i·s); inert
        # families (absent from weights_by_id) are unmodified.
        modifiers = [CouplingModifier(trajectory, Scalar(-coupling.effect_loss), EventType.LOSS,
                                      weights=coupling.weights_by_id)]
        if coupling.effect_gain != 0.0:      # optional field-blind gain: transfer ×= exp(effect_gain·s)
            modifiers.append(
                CouplingModifier(trajectory, Scalar(coupling.effect_gain), EventType.TRANSFER))
        super().__init__(base, modifiers)


# ═══════════════════════════════════════════════════════════════════════════════
# Driving a trait-linked simulation
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class TraitGeneResult:
    """Output of :func:`simulate_trait_conditioned_genomes`.

    ``profiles`` is the ``N × extant-species`` panel matrix (all panel rows kept, including
    all-absent ones); ``trait`` is the :class:`~zombi2.traits.TraitResult` the genes were
    conditioned on. :meth:`genomes` promotes the run to a standard :class:`~zombi2.Genomes` so
    the whole output pipeline (gene trees, reconciliations, ``write()``) applies unchanged.
    """

    species_tree: Tree
    profiles: ProfileMatrix
    trait: TraitResult
    leaf_genomes: dict
    event_log: object
    coupling: TraitGeneCoupling = field(repr=False, default=None)

    def genomes(self):
        """Promote to a :class:`~zombi2.Genomes` (shares the standard output writers)."""
        from zombi2.genomes.simulation import Genomes
        return Genomes(species_tree=self.species_tree, leaf_genomes=self.leaf_genomes,
                       event_log=self.event_log, profiles=self.profiles)


def _seed_panel_factory(seed_families):
    """A ``genome_factory`` that seeds a genome with exactly ``seed_families`` present."""
    def factory(ids):
        g = UnorderedGenome(ids)
        for fam in seed_families:
            g._add(Gene(ids.new_gene(), fam))
        return g
    return factory


def _panel_profile(leaf_genomes, coupling: TraitGeneCoupling) -> ProfileMatrix:
    """Profile over the *full* panel, keeping every panel row (even all-absent ones)."""
    species_nodes = sorted(leaf_genomes, key=lambda n: _natkey(n.name))
    species = [n.name for n in species_nodes]
    index = coupling.index
    rows, cols, data = [], [], []
    for j, node in enumerate(species_nodes):
        genome = leaf_genomes[node]
        for fam in genome.families():
            k = index.get(fam)
            if k is None:
                continue  # a background (originated) family — not part of the panel
            cn = genome.copy_number(fam)
            if cn:
                rows.append(k); cols.append(j); data.append(cn)
    return ProfileMatrix(families=list(coupling.panel_ids), species=species, coo=(rows, cols, data))


def simulate_trait_conditioned_genomes(
    tree: Tree,
    trait,
    coupling: TraitGeneCoupling,
    *,
    trait_steps: int = 16,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    transfers: TransferModel | None = None,
    initial_presence=None,
) -> TraitGeneResult:
    """Simulate gene families along ``tree`` conditioned on a trait.

    ``trait`` is either a trait **model** (e.g. :class:`~zombi2.Mk`, :class:`~zombi2.BrownianMotion`)
    — evolved here with :func:`~zombi2.simulate_traits` — or an already-simulated
    :class:`~zombi2.traits.TraitResult` to reuse. ``coupling`` names the responsive panel and the
    effect sizes; ``trait_steps`` is the within-branch resolution for a continuous trait (ignored
    for a discrete one, whose exact stochastic map is used).

    The panel is seeded present at the root (pass ``initial_presence`` as a length-``N`` 0/1 mask
    for a different start). Transfers default to full replacement so a re-acquired family does not
    stack copies — keeping the panel cleanly presence/absence.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    result = trait if isinstance(trait, TraitResult) else simulate_traits(tree, trait, rng=rng)
    trajectory = TraitTrajectory.from_result(result, steps=trait_steps,
                                             state_values=coupling.state_values)

    if initial_presence is None:
        seed_families = list(coupling.panel_ids)
    else:
        mask = np.asarray(initial_presence)
        if mask.shape != (coupling.n_families,):
            raise ValueError(f"initial_presence must have shape ({coupling.n_families},)")
        seed_families = [fam for fam, on in zip(coupling.panel_ids, mask) if on]

    rates = TraitGeneRates(coupling, trajectory)
    tm = transfers if transfers is not None else TransferModel(replacement=1.0)
    gres = GenomeSimulator().simulate(
        tree, rates, rng, initial_size=0, transfers=tm,
        genome_factory=_seed_panel_factory(seed_families),
    )
    return TraitGeneResult(
        species_tree=tree,
        profiles=_panel_profile(gres.leaf_genomes, coupling),
        trait=result,
        leaf_genomes=gres.leaf_genomes,
        event_log=gres.event_log,
        coupling=coupling,
    )


# Backwards-compatible aliases — the ``traits:genomes`` edge standardised on the ``TraitGene*`` stem
# (matching its config :class:`TraitGeneCoupling` and the joint :class:`TraitGeneFeedback`); the old
# ``TraitLinked*`` names named the same edge. They remain the **same objects**, so existing deep
# imports (``from zombi2.coevolve.trait_coupling import TraitLinkedRates``) keep working silently;
# the package surfaces (``zombi2`` / ``zombi2.coevolve``) resolve them with a DeprecationWarning.
TraitLinkedRates = TraitGeneRates
TraitLinkedResult = TraitGeneResult
simulate_trait_linked_genomes = simulate_trait_conditioned_genomes
