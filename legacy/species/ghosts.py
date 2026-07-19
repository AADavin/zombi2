"""Un-pruning: graft dead (ghost) lineages onto a reconstructed species tree.

The backward :func:`~zombi2.simulate_species_tree` yields only lineages with sampled extant
descendants. This module adds the dead ones back using the exact conditional law of the
complete tree given the reconstructed one: along each reconstructed edge, ghost lineages
attach as an inhomogeneous Poisson process with intensity ``λ(t)·E(t)``, where ``E`` is the
probability of leaving no sampled descendant (Nee, May & Harvey 1994; Stadler 2009;
Lambert & Stadler 2013). Each ghost roots a birth–death subtree conditioned on no sampled
descendant. See ``docs/ghost_lineages.md``.

Supports the constant-rate :class:`~zombi2.BirthDeath`/:class:`~zombi2.Yule` and the
:class:`~zombi2.EpisodicBirthDeath` (time-varying rates, incomplete sampling ``ρ<1``). With
``ρ=1`` ghosts are lineages that go extinct before the present; with ``ρ<1`` they also include
lineages alive today but unsampled. Once grafted, the forward gene-family simulator uses them
automatically as transfer donors/recipients ("transfer from the dead").
"""

from __future__ import annotations

import bisect
import math

import numpy as np

from zombi2.species.model import BirthDeath, EpisodicBirthDeath
from zombi2.tree import Tree, TreeNode


class _ProcessView:
    """Model-agnostic access to the rates the grafting needs, as functions of age ``τ`` (time
    before the present): ``E(τ)`` (no-sampled-descendant prob), ``rates(τ) -> (λ, μ)``, the
    sampling fraction ``ρ``, and thinning bounds ``birth_bound`` (max λ) / ``rate_bound``
    (max λ+μ)."""

    __slots__ = ("E", "rates", "rho", "birth_bound", "rate_bound")

    def __init__(self, model, total_age):
        if isinstance(model, EpisodicBirthDeath):
            model.validate()
            if getattr(model, "_cache_A", None) != total_age:
                model._prepare(total_age)
            ages, egrid = model._ages, model._E
            shifts, births, deaths = model.shifts, model.birth, model.death

            def E(tau):
                return float(np.interp(tau, ages, egrid))

            def rates(tau):
                i = bisect.bisect_right(shifts, tau)
                return births[i], deaths[i]

            self.E = E
            self.rates = rates
            self.rho = model.rho
            self.birth_bound = max(births)
            self.rate_bound = max(b + d for b, d in zip(births, deaths))
        elif isinstance(model, BirthDeath):  # Yule is a subclass -> supported (μ=0 -> no ghosts)
            model.validate()
            b, d = model.birth, model.death
            self.E = model.extinction_prob  # ρ-aware (uses model.sampling_fraction)
            self.rates = lambda tau: (b, d)
            self.rho = model.sampling_fraction
            self.birth_bound = b
            self.rate_bound = b + d
        else:
            raise NotImplementedError(
                f"add_ghost_lineages supports BirthDeath/Yule and EpisodicBirthDeath, "
                f"not {type(model).__name__}"
            )


def _next_event(bt, total_age, view, rng):
    """First event of a lineage born at absolute time ``bt``, via thinning (exact for the
    piecewise-constant episodic rates). Returns ``("present"|"birth"|"death", time)``."""
    t = bt
    bound = view.rate_bound
    while True:
        t += rng.exponential(1.0 / bound)
        if t >= total_age:
            return "present", total_age
        lam, mu = view.rates(total_age - t)
        r = lam + mu
        if r > 0.0 and rng.random() < r / bound:
            return ("birth", t) if rng.random() < lam / r else ("death", t)


def _grow_ghost_rejection(t_start, total_age, view, rng, max_size):
    """Grow one birth–death subtree forward from ``t_start``, conditioned (by rejection) on
    leaving no sampled descendant. Returns the subtree root, or ``None`` to reject."""
    root = TreeNode(name="", time=t_start, is_extant=False)
    root._birth = t_start
    stack = [root]
    size = 1
    while stack:
        node = stack.pop()
        bt = node._birth
        del node._birth
        kind, et = _next_event(bt, total_age, view, rng)
        node.time = et
        if kind == "present":
            if rng.random() < view.rho:  # sampled tip -> violates the conditioning
                return None
            # else: an unsampled extant ghost (ρ<1); is_extant already False
        elif kind == "birth":
            for _ in range(2):
                child = TreeNode(name="", time=et, is_extant=False)
                child._birth = et
                node.add_child(child)
                stack.append(child)
            size += 2
            if size > max_size:  # runaway -> almost surely leaves a sampled tip
                return None
        # "death": a dead-before-present leaf; nothing more to grow
    return root


class _HazardSampler:
    """Direct (rejection-free) sampler for a ghost subtree via Doob's h-transform.

    Conditioning a lineage on leaving no sampled descendant turns the birth–death process into
    one with per-lineage birth ``λ(τ)·E(τ)`` and death ``μ(τ)/E(τ)`` (τ = age before present).
    The death rate diverges as ``τ→0`` when ``ρ=1`` (``E→0``), which forces extinction before
    the present — exactly as it should. We sample the next-event age by inverting the cumulative
    hazard ``T(τ)=∫_τ^A g``, tabulated on a grid, with the log-singular first cell handled
    analytically (there ``g ≈ 1/τ``, so ``∫_{τ_e}^{τ1} g = ln(τ1/τ_e)``).
    """

    def __init__(self, view, total_age, grid=4000):
        self.view = view
        taus = np.linspace(0.0, total_age, grid + 1)
        E = np.array([view.E(t) for t in taus])
        lam = np.empty_like(taus)
        mu = np.empty_like(taus)
        for i, t in enumerate(taus):
            lam[i], mu[i] = view.rates(t)
        with np.errstate(divide="ignore", invalid="ignore"):
            g = lam * E + np.where(E > 0.0, mu / np.where(E > 0.0, E, 1.0), np.inf)
        # cumulative hazard from the top: T[i] = ∫_{taus[i]}^{A} g dτ (trapezoid, skip cell 0
        # when it is singular). T is decreasing in τ.
        T = np.zeros_like(taus)
        for i in range(grid - 1, 0, -1):
            T[i] = T[i + 1] + 0.5 * (g[i] + g[i + 1]) * (taus[i + 1] - taus[i])
        self.singular = E[0] <= 0.0  # ρ=1 → E(0)=0 → log singularity in cell [0, taus[1]]
        if self.singular:
            T[0] = math.inf
        else:
            T[0] = T[1] + 0.5 * (g[0] + g[1]) * (taus[1] - taus[0])
        self.taus, self.T = taus, T
        self.tau1, self.T1, self.T0 = taus[1], T[1], T[0]
        self.rho = view.rho
        # reversed regular arrays for inverse lookup (T increasing)
        self._Tr = T[1:][::-1]
        self._taur = taus[1:][::-1]

    def _T_at(self, tau):
        if self.singular and tau < self.tau1:
            return self.T1 + math.log(self.tau1 / tau)
        return float(np.interp(tau, self.taus, self.T))

    def sample_event_age(self, tau_b, rng):
        """First event of a lineage born at age ``tau_b``; ``("present"|"birth"|"death", τ_e)``."""
        v = self._T_at(tau_b) + rng.exponential(1.0)
        if v > self.T0:  # only reachable when ρ<1 (T0 finite): the lineage reaches the present
            return "present", 0.0
        if v > self.T1:  # first cell (0, tau1)
            if self.singular:
                tau_e = self.tau1 * math.exp(-(v - self.T1))
            else:
                tau_e = self.tau1 * (v - self.T0) / (self.T1 - self.T0)
        else:
            tau_e = float(np.interp(v, self._Tr, self._taur))
        ee = self.view.E(tau_e)
        le, me = self.view.rates(tau_e)
        ge = le * ee + me / ee
        kind = "birth" if rng.random() * ge < le * ee else "death"
        return kind, tau_e


def _grow_ghost_htransform(t_start, total_age, hazard, rng):
    """Grow one ghost subtree via the h-transform — every realisation already leaves no sampled
    descendant, so there is no rejection."""
    root = TreeNode(name="", time=t_start, is_extant=False)
    root._tau = total_age - t_start
    stack = [root]
    while stack:
        node = stack.pop()
        tau_b = node._tau
        del node._tau
        kind, tau_e = hazard.sample_event_age(tau_b, rng)
        if kind == "present":
            node.time = total_age  # unsampled extant ghost
            continue
        node.time = total_age - tau_e
        if kind == "birth":
            for _ in range(2):
                child = TreeNode(name="", time=node.time, is_extant=False)
                child._tau = tau_e
                node.add_child(child)
                stack.append(child)
        # "death": a dead-before-present leaf
    return root


def add_ghost_lineages(
    tree: Tree,
    model,
    *,
    method: str = "rejection",
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_subtree_size: int = 4096,
    max_attempts: int = 100_000,
) -> Tree:
    """Graft dead (ghost) lineages onto a reconstructed ``tree`` **in place**, and return it.

    ``model`` should be the same species-tree model used to build the tree
    (:class:`~zombi2.BirthDeath`/:class:`~zombi2.Yule`, or
    :class:`~zombi2.EpisodicBirthDeath`). Ghost lineages attach along each edge as a Poisson
    process with intensity ``λ(t)·E(t)`` (τ = time before present) and each roots a
    birth–death subtree conditioned on leaving no sampled descendant; every ghost node gets
    ``is_extant=False``, and the new leaves are named ``e*`` (extinct), the new internal nodes
    ``i*`` — the same convention as the rest of ZOMBI2. The sampled (extant) leaves are left
    untouched, so pruning back to them recovers the original tree.

    Parameters
    ----------
    method:
        How each ghost subtree is grown, conditioned on leaving no sampled descendant:
        ``"rejection"`` (default; simple and exact) or ``"htransform"`` (rejection-free via
        Doob's h-transform, faster when extinction/sampling makes rejection retry a lot).
        Both are statistically equivalent.
    max_subtree_size:
        Reject (and retry) a ghost subtree that exceeds this many lineages — a cheap guard
        against runaway supercritical growth (bias is negligible). ``"rejection"`` only.
    max_attempts:
        Safety cap on rejection retries per ghost. ``"rejection"`` only.
    """
    if method not in ("rejection", "htransform"):
        raise ValueError(f"method must be 'rejection' or 'htransform', got {method!r}")
    if rng is None:
        rng = np.random.default_rng(seed)

    view = _ProcessView(model, tree.total_age)
    total_age = tree.total_age
    birth_bound = view.birth_bound
    if birth_bound <= 0.0:  # no speciation -> no ghosts possible
        return tree
    hazard = _HazardSampler(view, total_age) if method == "htransform" else None

    # snapshot the reconstructed edges before mutating the tree
    edges = [n for n in tree.nodes_preorder() if n.parent is not None]

    for node in edges:
        parent = node.parent  # current upper endpoint of the (shrinking) segment
        t = parent.time
        t1 = node.time
        # inhomogeneous Poisson via thinning: propose at rate birth_bound, accept w.p.
        # λ(τ)·E(τ) / birth_bound (≤ 1 since λ ≤ birth_bound and E ≤ 1)
        while True:
            t += rng.exponential(1.0 / birth_bound)
            if t >= t1:
                break
            tau = total_age - t
            lam, _ = view.rates(tau)
            if rng.random() >= lam * view.E(tau) / birth_bound:
                continue  # thinned out
            if hazard is not None:
                ghost = _grow_ghost_htransform(t, total_age, hazard, rng)
            else:
                ghost = None
                for _ in range(max_attempts):
                    ghost = _grow_ghost_rejection(t, total_age, view, rng, max_subtree_size)
                    if ghost is not None:
                        break
                if ghost is None:
                    continue  # gave up on this attachment (extremely rare)
            # splice a binary junction between `parent` and `node` at time t
            junction = TreeNode(name="", time=t, is_extant=False)
            parent.children[parent.children.index(node)] = junction
            junction.parent = parent
            junction.add_child(node)   # the reconstructed lineage continues
            junction.add_child(ghost)  # the dead sibling
            parent = junction          # next ghost on this edge sits above `node`, below M

    # Name the new nodes with the shared convention (sampled-extant n*, unsampled-extant ghosts u*,
    # extinct e*, internal i*), continuing past the indices already used on the reconstructed tree.
    def _max_index(prefix: str) -> int:
        used = [int(n.name[len(prefix):]) for n in tree.nodes_preorder()
                if n.name.startswith(prefix) and n.name[len(prefix):].isdigit()]
        return max(used, default=0)

    n_ctr, u_ctr, e_ctr, i_ctr = _max_index("n"), _max_index("u"), _max_index("e"), _max_index("i")
    present = tree.total_age
    for node in tree.nodes_preorder():
        if node.name:
            continue
        if node.is_leaf():
            if node.is_extant:
                n_ctr += 1
                node.name = f"n{n_ctr}"
            elif not node.sampled and abs(node.time - present) <= 1e-9 * max(1.0, abs(present)):
                u_ctr += 1
                node.name = f"u{u_ctr}"
            else:
                e_ctr += 1
                node.name = f"e{e_ctr}"
        else:
            i_ctr += 1
            node.name = f"i{i_ctr}"
    return tree
