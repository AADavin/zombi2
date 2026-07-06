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
the family is purged. Gain is horizontal transfer (a field-blind influx), exactly as in the
:class:`~zombi2.PottsRates` coupling model: a family flows in and the trait-modulated loss
then *selectively retains* it, which is what writes the trait↔gene association into the
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
:meth:`TraitLinkedResult.genomes`.
"""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field

import numpy as np

from zombi2.genomes.events import EventType
from zombi2.genomes.genome import Gene, UnorderedGenome
from zombi2.genomes.genome_sim import GenomeSimulator
from zombi2.genomes.profiles import ProfileMatrix, _natkey
from zombi2.genomes.rates import EventWeight, RateModel
from zombi2.traits.models import TraitResult, simulate_traits
from zombi2.genomes.transfers import TransferModel
from zombi2.tree import Tree

#: Clamp on the loss/gain exponent so an extreme ``effect·w·s`` never overflows ``exp`` (or
#: drives the Gillespie loop into an instant loss/gain hot cycle). Matches :mod:`zombi2.coupling`.
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
class TraitLinkedRates(RateModel):
    """Trait-conditioned loss over a fixed family panel, plus field-blind transfer gain.

    For every present family it emits a per-family loss weight — ``base_loss·copies`` for an
    inert family, and ``base_loss·copies·exp(-effect_loss·w_i·s)`` for a responsive one, with
    ``s`` the local trait value from the :class:`TraitTrajectory`. Gain is a single field-blind
    ``TRANSFER`` channel (optionally scaled by ``exp(effect_gain·s)``), plus optional
    duplication and background origination.

    The weights change only at gene events and at the trajectory's own breakpoints, so
    ``time_dependent`` stays ``False`` and the simulator refreshes a branch exactly when it
    fires an event or crosses a trait change (:meth:`refresh_times`) — no blunt full refresh.
    """

    def __init__(self, coupling: TraitGeneCoupling, trajectory: TraitTrajectory):
        self.c = coupling
        self.traj = trajectory
        self._w = coupling.weights_by_id

    def event_weights(self, genome, branch, time):
        c = self.c
        s = self.traj.value(branch, time)
        out: list[EventWeight] = []

        for fam in genome.families():
            cn = genome.copy_number(fam)
            w = self._w.get(fam, 0.0)
            if w != 0.0 and c.effect_loss != 0.0 and s != 0.0:
                rate = c.base_loss * cn * math.exp(_clamp(-c.effect_loss * w * s))
            else:
                rate = c.base_loss * cn
            if rate > 0.0:
                out.append(EventWeight(EventType.LOSS, fam, rate))

        n = genome.size()
        if c.duplication > 0.0 and n > 0:
            out.append(EventWeight(EventType.DUPLICATION, None, c.duplication * n))
        if c.transfer > 0.0 and n > 0:
            rate = c.transfer * n
            if c.effect_gain != 0.0 and s != 0.0:
                rate *= math.exp(_clamp(c.effect_gain * s))
            out.append(EventWeight(EventType.TRANSFER, None, rate))
        if c.origination > 0.0:
            out.append(EventWeight(EventType.ORIGINATION, None, c.origination))
        return out

    def refresh_times(self, t0, t1):
        return self.traj.refresh_times(t0, t1)


# ═══════════════════════════════════════════════════════════════════════════════
# Driving a trait-linked simulation
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class TraitLinkedResult:
    """Output of :func:`simulate_trait_linked_genomes`.

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


def simulate_trait_linked_genomes(
    tree: Tree,
    trait,
    coupling: TraitGeneCoupling,
    *,
    trait_steps: int = 16,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    transfers: TransferModel | None = None,
    initial_presence=None,
) -> TraitLinkedResult:
    """Simulate gene families along ``tree`` conditioned on a trait.

    ``trait`` is either a trait **model** (e.g. :class:`~zombi2.Mk`, :class:`~zombi2.BrownianMotion`)
    — evolved here with :func:`~zombi2.simulate_traits` — or an already-simulated
    :class:`~zombi2.traits.TraitResult` to reuse. ``coupling`` names the responsive panel and the
    effect sizes; ``trait_steps`` is the within-branch resolution for a continuous trait (ignored
    for a discrete one, whose exact stochastic map is used).

    The panel is seeded present at the root (pass ``initial_presence`` as a length-``N`` 0/1 mask
    for a different start). Transfers default to full replacement so a re-acquired family does not
    stack copies — keeping the panel cleanly presence/absence, as in :func:`~zombi2.simulate_coupled`.
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

    rates = TraitLinkedRates(coupling, trajectory)
    tm = transfers if transfers is not None else TransferModel(replacement=1.0)
    gres = GenomeSimulator().simulate(
        tree, rates, rng, initial_size=0, transfers=tm,
        genome_factory=_seed_panel_factory(seed_families),
    )
    return TraitLinkedResult(
        species_tree=tree,
        profiles=_panel_profile(gres.leaf_genomes, coupling),
        trait=result,
        leaf_genomes=gres.leaf_genomes,
        event_log=gres.event_log,
        coupling=coupling,
    )
