"""Forward-in-time species-tree simulation.

The default :func:`~zombi2.simulate_species_tree` runs *backward* and yields the reconstructed
tree (surviving lineages only). This module runs the birth–death process *forward* and returns
the **complete** tree — extinct lineages included natively (``is_extant=False`` leaves at their
death times), no ghost-grafting needed. It complements :func:`~zombi2.add_ghost_lineages`: two
routes to the same complete-tree object, which the forward gene simulator then uses with
transfers from the dead for free.

Conventions (matching the backward crown tree): the tree is rooted at the **crown** (two
lineages at ``time == 0``) and the present is at ``total_age``; ``age`` is therefore the crown
age. It is conditioned on ≥2 sampled survivors.

Supported models: constant-rate :class:`~zombi2.BirthDeath`/:class:`~zombi2.Yule` (either stop
mode), and :class:`~zombi2.EpisodicBirthDeath` (time-varying λ/μ and incomplete sampling
``ρ<1``) in **age mode only** — the present must be fixed for "age before present" to be defined
forward.

Any of these may carry ``mass_extinctions`` — instantaneous, tree-wide survival pulses at
specified ages before the present, where every standing lineage independently dies with a given
fraction. Because they place their times as ages before the present, they too require age mode.
The killed lineages become ordinary extinct leaves, so the forward gene simulator treats them as
ghost transfer partners with no extra work.

Three further models have rates that are *constant between events* (not time-varying), so they are
grown by an **exact Gillespie** loop (``_grow_gillespie``) instead of the thinning loop above:
:class:`~zombi2.ClaDS` (per-lineage rates that shift at each speciation),
:class:`~zombi2.DiversityDependent` (rates that decline as the tree fills a carrying capacity), and
:class:`~zombi2.CladeShiftBirthDeath` (scheduled clade-specific rate shifts). All are forward-only;
the first two support ``age`` or ``n_tips`` mode, the clade-shift model ``age`` mode only.
"""

from __future__ import annotations

import bisect
import math

import numpy as np

from zombi2.species.model import (
    BirthDeath, CladeShiftBirthDeath, ClaDS, DiversityDependent, EpisodicBirthDeath,
    SharedBirthDeath,
)
from zombi2.species._caps import GrowthEngine, species_caps
from zombi2.tree import Tree, TreeNode


class _ForwardRates:
    """Rates as functions of tree-time ``t`` (0 = crown, present at ``present``). Episodic
    rates map tree-time to age-before-present ``present - t`` (age mode, so ``present`` is
    fixed). Provides ``rates(t) -> (λ, μ, ψ)`` (ψ = serial/fossil sampling rate), a thinning
    bound, the extant sampling fraction ρ, the removal probability ``r`` on sampling, and
    ``mass_extinctions`` — a list of ``(tree_time, survival)`` pulses (each surviving lineage
    is kept with probability ``survival`` at that instant), sorted by increasing tree-time."""

    __slots__ = ("rates", "rate_bound", "rho", "removal", "mass_extinctions")

    def __init__(self, model, present):
        if isinstance(model, EpisodicBirthDeath):
            model.validate()
            shifts, births, deaths, foss = (model.shifts, model.birth, model.death,
                                            model.fossilization)

            def rates(t):
                i = bisect.bisect_right(shifts, present - t)
                return births[i], deaths[i], foss[i]

            self.rates = rates
            self.rate_bound = max(b + d + f for b, d, f in zip(births, deaths, foss))
            self.rho = model.rho
            self.removal = model.removal
        elif isinstance(model, BirthDeath):  # Yule is a subclass
            model.validate()
            b, d, psi = model.birth, model.death, model.fossilization
            self.rates = lambda t: (b, d, psi)
            self.rate_bound = b + d + psi
            self.rho = model.sampling_fraction
            self.removal = model.removal
        else:
            raise NotImplementedError(
                f"forward simulation supports BirthDeath/Yule and EpisodicBirthDeath, "
                f"not {type(model).__name__}"
            )
        # ages before present -> tree-time; death fraction -> survival probability
        self.mass_extinctions = sorted(
            (present - age, 1.0 - frac) for age, frac in getattr(model, "mass_extinctions", ())
        )


def _new_crown():
    """A fresh crown for a forward trial: a root at time 0 with two live children at time 0.
    Returns ``(root, live)``. Shared by both forward growth loops."""
    root = TreeNode(name="", time=0.0)
    live = []
    for _ in range(2):
        child = TreeNode(name="", time=0.0)
        root.add_child(child)
        live.append(child)
    return root, live


def _grow(view, age, n_tips, rng, max_lineages):
    """One forward trial from a crown of two lineages (thinning handles time-varying rates).
    Returns ``(crown_node, end_time)`` or ``None`` to reject (extinct / <2 sampled survivors)."""
    root, live = _new_crown()
    bound = view.rate_bound
    mass_ext = view.mass_extinctions
    me_idx = 0
    t = 0.0
    end = None
    while True:
        n = len(live)
        if n == 0:
            return None
        if n_tips is not None and n == n_tips:
            # place the present strictly after the N-th lineage appeared: the tree age is the last
            # speciation time plus a memoryless waiting time to the next event (the standard
            # birth-death-conditioned-on-N convention). Using end = t would give the two newest tips
            # and their parent zero-length pendant edges — degenerate to an age-0 tree at n_tips == 2.
            lam, mu, psi = view.rates(t)
            R = n * (lam + mu + psi)
            end = t + rng.exponential(1.0 / (R if R > 0.0 else n * bound))
            break
        if n > max_lineages:
            raise RuntimeError(
                f"forward tree exceeded max_lineages={max_lineages}; explosive parameters "
                "(birth >> death over this age) — lower the age/rates or raise max_lineages"
            )
        dt = rng.exponential(1.0 / (n * bound))
        # a scheduled mass extinction before the next candidate event fires first (the
        # thinning process is memoryless, so we advance to it and redraw the waiting time)
        if me_idx < len(mass_ext) and mass_ext[me_idx][0] <= t + dt:
            me_t, survival = mass_ext[me_idx]
            me_idx += 1
            t = me_t
            kept = []
            for node in live:
                if survival >= 1.0 or rng.random() < survival:
                    kept.append(node)
                else:  # struck down by the mass extinction -> an extinct leaf at this instant
                    node.time = me_t
                    node.is_extant = False
            live = kept
            continue
        if age is not None and t + dt >= age:
            end = age
            break
        t += dt
        lam, mu, psi = view.rates(t)
        total = lam + mu + psi
        if total <= 0.0 or rng.random() >= total / bound:
            continue  # thinned out (or an epoch with no events)
        i = int(rng.integers(n))
        node = live[i]
        node.time = t
        live[i] = live[-1]
        live.pop()
        r = rng.random() * total
        if r < lam:  # speciation
            a = TreeNode(name="", time=t)
            b = TreeNode(name="", time=t)
            node.add_child(a)
            node.add_child(b)
            live.append(a)
            live.append(b)
        elif r < lam + mu:  # extinction
            node.is_extant = False
        else:  # serial (through-time) sampling
            node.is_extant = False
            node.sampled = True
            if rng.random() >= view.removal:  # not removed -> sampled ancestor (lineage continues)
                cont = TreeNode(name="", time=t)
                node.add_child(cont)
                live.append(cont)
            # else: removed -> a dated fossil tip (node stays a leaf)

    if _finalize_present(live, end, view.rho, rng) < 2:
        return None
    return root, end


def _finalize_present(live, end, rho, rng) -> int:
    """Survivors reach the present at ``end``; each is sampled (kept extant) with probability
    ``rho``. Unsampled survivors are marked ``is_extant=False`` (ghost tips). Returns the number
    of sampled tips (the run is conditioned on ≥2)."""
    n_sampled = 0
    for node in live:
        node.time = end
        if rho >= 1.0 or rng.random() < rho:
            node.is_extant = True
            node.sampled = True
            n_sampled += 1
        else:
            node.is_extant = False
    return n_sampled


# --- exact-Gillespie forward growth (rates constant between events) ------------
# For models whose per-lineage (λ, μ) do not vary with time between events, we sample the exact
# next-event time from the summed rate — no thinning bound needed. Each live lineage carries an
# opaque *state* (aligned with ``live``) that the view maps to (λ, μ); the loop also interleaves a
# merged, time-ordered timeline of scheduled tree-wide events (mass extinctions, clade shifts).
# Covers ClaDS (state = the lineage's own λ), diversity-dependent BD (rates depend only on the
# live count n), and clade-shift BD (state = the lineage's current (λ, μ) regime).

def _build_schedule(present, mass_extinctions=(), clade_shifts=()):
    """Merge scheduled tree-wide events into one list of ``(tree_time, kind, payload)`` sorted by
    tree-time. Ages before present map to tree-time ``present - age``. ``kind`` is
    ``"mass_extinction"`` (payload = survival probability) or ``"clade_shift"`` (payload = the new
    ``(λ, μ)`` regime)."""
    events = [(present - a, "mass_extinction", 1.0 - f) for a, f in mass_extinctions]
    events += [(present - a, "clade_shift", (b, d)) for a, b, d in clade_shifts]
    return sorted(events, key=lambda e: e[0])


class _ClaDSView:
    """Per-lineage ClaDS rates. Each lineage's state is its own λ; μ = ε·λ (constant turnover).
    Daughters jump multiplicatively at each speciation."""

    __slots__ = ("initial_state", "eps", "log_alpha", "sigma", "rho", "scheduled")

    def __init__(self, model: ClaDS, present):
        self.initial_state = model.lambda_0
        self.eps = model.turnover
        self.log_alpha = math.log(model.alpha)
        self.sigma = model.sigma
        self.rho = model.sampling_fraction
        self.scheduled = _build_schedule(present, model.mass_extinctions)

    def lineage_rates(self, lam_i, n):
        return lam_i, self.eps * lam_i

    def split(self, lam_i, rng):
        m1 = math.exp(self.log_alpha + self.sigma * rng.normal())
        m2 = math.exp(self.log_alpha + self.sigma * rng.normal())
        return lam_i * m1, lam_i * m2


class _DDView:
    """Diversity-dependent rates: all lineages share λ(n) = max(0, λ₀·(1 − n/K)), constant μ.
    The per-lineage state is an unused dummy (rates come from the count n, not the lineage)."""

    __slots__ = ("initial_state", "lambda_0", "mu", "K", "rho", "scheduled")

    def __init__(self, model: DiversityDependent, present):
        self.initial_state = model.lambda_0
        self.lambda_0 = model.lambda_0
        self.mu = model.death
        self.K = model.K
        self.rho = model.sampling_fraction
        self.scheduled = _build_schedule(present, model.mass_extinctions)

    def lineage_rates(self, state, n):
        return self.lambda_0 * max(0.0, 1.0 - n / self.K), self.mu

    def split(self, state, rng):
        return state, state


class _SharedView:
    """Shared-clock birth–death. The *total* birth and death rates are fixed, so each of the ``n``
    live lineages carries ``(birth/n, death/n)`` — the summed rate is a constant ``(birth, death)``
    regardless of ``n`` (equivalently, per-lineage ``λ(n) = birth/n``). All lineages share it
    equally, so the actor is chosen uniformly. The per-lineage state is an unused dummy."""

    __slots__ = ("initial_state", "birth", "death", "rho", "scheduled")

    def __init__(self, model: SharedBirthDeath, present):
        self.initial_state = None
        self.birth = model.birth
        self.death = model.death
        self.rho = model.sampling_fraction
        self.scheduled = _build_schedule(present, model.mass_extinctions)

    def lineage_rates(self, state, n):
        return self.birth / n, self.death / n

    def split(self, state, rng):
        return state, state


class _ShiftView:
    """Clade-specific rate shifts. Each lineage's state is its current ``(λ, μ)`` regime; daughters
    inherit the parent's regime. Shifts arrive as scheduled ``clade_shift`` events that reassign one
    random live lineage's regime (its descendants inherit it)."""

    __slots__ = ("initial_state", "rho", "scheduled")

    def __init__(self, model: CladeShiftBirthDeath, present):
        self.initial_state = (model.birth, model.death)
        self.rho = model.sampling_fraction
        self.scheduled = _build_schedule(present, model.mass_extinctions, model.clade_shifts)

    def lineage_rates(self, state, n):
        return state  # (λ, μ)

    def split(self, state, rng):
        return state, state


#: Gillespie models -> their per-lineage rate view. Data-driven replacement for the isinstance
#: view ladder in ``simulate_forward`` (keyed by exact type; these classes have no subclasses).
_GILLESPIE_VIEWS = {ClaDS: _ClaDSView, DiversityDependent: _DDView, CladeShiftBirthDeath: _ShiftView,
                    SharedBirthDeath: _SharedView}


def _weighted_index(weights, total, rng) -> int:
    """Index sampled proportional to ``weights`` (which sum to ``total``)."""
    x = rng.random() * total
    cum = 0.0
    for i, w in enumerate(weights):
        cum += w
        if x < cum:
            return i
    return len(weights) - 1


def _cull(live, state, me_t, survival, rng):
    """Apply a mass extinction: keep each lineage with probability ``survival``; the rest become
    extinct leaves at ``me_t``. Returns the pruned ``(live, state)`` lists (kept in lockstep)."""
    kept_live, kept_state = [], []
    for node, st in zip(live, state):
        if survival >= 1.0 or rng.random() < survival:
            kept_live.append(node)
            kept_state.append(st)
        else:
            node.time = me_t
            node.is_extant = False
    return kept_live, kept_state


def _apply_scheduled(event, live, state, rng):
    """Apply one scheduled tree-wide event, returning the (possibly new) ``(live, state, t)``."""
    t, kind, payload = event
    if kind == "mass_extinction":
        live, state = _cull(live, state, t, payload, rng)
    else:  # "clade_shift": a uniformly chosen live lineage (and its descendants) adopts (λ, μ)
        if live:
            state[int(rng.integers(len(live)))] = payload
    return live, state, t


def _grow_gillespie(view, age, n_tips, rng, max_lineages):
    """One forward trial for a rates-constant-between-events model (ClaDS / diversity-dependent /
    clade-shift). ``state`` runs in lockstep with ``live``, holding each lineage's opaque rate
    state. Returns ``(crown_node, end_time)`` or ``None`` to reject (extinct / <2 survivors)."""
    root, live = _new_crown()
    state = [view.initial_state for _ in live]
    scheduled = view.scheduled
    s_idx = 0
    t = 0.0
    end = None
    while True:
        n = len(live)
        if n == 0:
            return None
        if n > max_lineages:
            raise RuntimeError(
                f"forward tree exceeded max_lineages={max_lineages}; explosive parameters — "
                "lower the age/rates or raise max_lineages"
            )
        rates = [view.lineage_rates(state[i], n) for i in range(n)]  # (λ, μ) per lineage
        totals = [b + d for b, d in rates]
        total_rate = math.fsum(totals)
        if n_tips is not None and n == n_tips:
            # present strictly after the N-th birth: last event + Exp(total rate). See _grow.
            end = t + rng.exponential(1.0 / total_rate) if total_rate > 0.0 else t
            break
        if total_rate <= 0.0:
            # nothing stochastic can happen (e.g. diversity-dependent at capacity with μ=0): jump
            # to the next scheduled event, or coast to the present
            if s_idx < len(scheduled):
                live, state, t = _apply_scheduled(scheduled[s_idx], live, state, rng)
                s_idx += 1
                continue
            if age is None:
                raise RuntimeError(
                    "the process stalled below the requested --tips (diversity-dependent tree at "
                    "carrying capacity with no extinction); use --age, a larger K, or death > 0"
                )
            end = age
            break
        dt = rng.exponential(1.0 / total_rate)
        if s_idx < len(scheduled) and scheduled[s_idx][0] <= t + dt:
            live, state, t = _apply_scheduled(scheduled[s_idx], live, state, rng)
            s_idx += 1
            continue
        if age is not None and t + dt >= age:
            end = age
            break
        t += dt
        i = _weighted_index(totals, total_rate, rng)
        node, st_i = live[i], state[i]
        b, d = rates[i]
        node.time = t
        live[i] = live[-1]
        live.pop()
        state[i] = state[-1]
        state.pop()
        if rng.random() * (b + d) < b:  # speciation
            sa, sb = view.split(st_i, rng)
            child_a = TreeNode(name="", time=t)
            child_b = TreeNode(name="", time=t)
            node.add_child(child_a)
            node.add_child(child_b)
            live.append(child_a)
            state.append(sa)
            live.append(child_b)
            state.append(sb)
        else:  # extinction
            node.is_extant = False

    if _finalize_present(live, end, view.rho, rng) < 2:
        return None
    return root, end


def _at_present(t: float, present: float) -> bool:
    """True if a leaf reached the present (an unsampled-extant *ghost*), rather than dying before it."""
    return abs(t - present) <= 1e-9 * max(1.0, abs(present))


def _name(tree: Tree) -> None:
    # Naming convention: sampled-extant leaves n*, unsampled-extant "ghost" leaves u* (alive at the
    # present but not sampled under ρ<1), extinct/fossil leaves e* (gone before the present),
    # internal nodes i*. Extinct and unsampled are different fates, so they get different letters.
    present = tree.total_age
    extant = unsampled = extinct = internal = 0
    for node in tree.nodes_preorder():
        if node is tree.root:
            node.name = "root"
        elif node.is_leaf():
            if node.is_extant:
                extant += 1
                node.name = f"n{extant}"
            elif not node.sampled and _at_present(node.time, present):
                unsampled += 1
                node.name = f"u{unsampled}"
            else:
                extinct += 1
                node.name = f"e{extinct}"
        else:
            internal += 1
            node.name = f"i{internal}"


def simulate_forward(
    model,
    *,
    age: float | None = None,
    n_tips: int | None = None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_attempts: int = 10_000,
    max_lineages: int = 1_000_000,
) -> Tree:
    """Forward implementation behind ``simulate_species_tree(..., direction="forward")``:
    a **complete** species tree (extinct/fossil lineages included).

    Provide exactly one stopping condition:

    * ``age`` — grow for this crown age; the number of extant tips is random.
    * ``n_tips`` — grow until this many extant lineages first coexist; the age is random.
      (Not for ``EpisodicBirthDeath``, whose present must be fixed.)

    The tree is rooted at the crown (``time == 0``), the present is at ``total_age``, extinct
    leaves carry ``is_extant=False`` at their death times, and extant leaves ``is_extant=True``
    at the present. Under incomplete sampling (``ρ<1``), extant but unsampled lineages are marked
    ``is_extant=False`` too. So :func:`~zombi2.simulate_genomes` treats extinct/unsampled lineages
    as ghost transfer partners automatically. The run is conditioned on ≥2 sampled survivors.

    ``ClaDS`` and ``DiversityDependent`` (rates constant between events) are grown by an exact
    Gillespie loop; the other models by thinning. Both routes share the mass-extinction, present-
    sampling, naming and conditioning logic.
    """
    if (age is None) == (n_tips is None):
        raise ValueError("provide exactly one of `age` or `n_tips`")
    caps = species_caps(model)           # loud TypeError for an unregistered model type
    heterogeneous = caps.growth is GrowthEngine.GILLESPIE
    if n_tips is not None and not caps.supports_n_tips:
        raise NotImplementedError(
            f"{type(model).__name__} is defined against a fixed present (time-varying rates or "
            "scheduled shifts), so it requires `age` mode, not `n_tips`"
        )
    model.validate()
    from zombi2.species.sim import _check_age, _check_n_tips
    if age is not None:
        _check_age(age)
    if n_tips is not None:
        _check_n_tips(n_tips, max_lineages)
        n_tips = int(n_tips)
        # n_tips mode stops on the count of *standing* lineages; present-day sampling (ρ<1) then
        # subsamples them, so the run would return ~ρ·n_tips sampled tips, not n_tips. Reject the
        # combination (backward simulation rejects ρ<1 too) — use age mode for incomplete sampling.
        rho = getattr(model, "sampling_fraction", 1.0)
        if rho < 1.0:
            raise ValueError(
                f"n_tips mode requires complete sampling, but sampling_fraction={rho:g}; it counts "
                "standing lineages, so ρ<1 would yield ~ρ·n_tips sampled tips rather than n_tips. "
                "Use age mode with sampling_fraction<1, or set sampling_fraction=1 for n_tips mode."
            )
    if isinstance(model, DiversityDependent) and n_tips is not None and n_tips > model.K:
        raise ValueError(
            f"n_tips ({n_tips}) exceeds the carrying capacity K ({model.K:g}); a diversity-"
            "dependent tree cannot grow past K, so use age mode or a smaller n_tips"
        )
    if isinstance(model, CladeShiftBirthDeath) and any(a >= age for a, _, _ in model.clade_shifts):
        raise ValueError(
            f"every clade-shift age must be < the crown age ({age}); "
            f"got ages {[a for a, _, _ in model.clade_shifts]}"
        )
    mes = getattr(model, "mass_extinctions", None)
    if mes:
        if n_tips is not None:
            raise NotImplementedError(
                "mass extinctions are placed at an age before the present, so they require "
                "`age` mode (a fixed present), not `n_tips`"
            )
        if any(a >= age for a, _ in mes):
            raise ValueError(
                f"every mass-extinction age must be < the crown age ({age}); "
                f"got ages {[a for a, _ in mes]}"
            )
    if rng is None:
        rng = np.random.default_rng(seed)

    present = age if age is not None else 0.0
    if heterogeneous:
        view = _GILLESPIE_VIEWS[type(model)](model, present)
        grow = _grow_gillespie
    else:
        view = _ForwardRates(model, present=present)
        grow = _grow
    for _ in range(max_attempts):
        result = grow(view, age, n_tips, rng, max_lineages)
        if result is not None:
            root, end_time = result
            tree = Tree(root, end_time)
            _name(tree)
            return tree

    raise RuntimeError(
        f"forward simulation produced no surviving tree in {max_attempts} attempts "
        "(the process kept going extinct); raise max_attempts, lower death, or use the "
        "backward simulate_species_tree"
    )
