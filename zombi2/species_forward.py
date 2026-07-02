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
"""

from __future__ import annotations

import bisect

import numpy as np

from .species_model import (
    BirthDeath,
    EpisodicBirthDeath,
    EpisodicFossilizedBirthDeath,
    FossilizedBirthDeath,
)
from .tree import Tree, TreeNode


class _ForwardRates:
    """Rates as functions of tree-time ``t`` (0 = crown, present at ``present``). Episodic
    rates map tree-time to age-before-present ``present - t`` (age mode, so ``present`` is
    fixed). Provides ``rates(t) -> (λ, μ, ψ)`` (ψ = serial/fossil sampling rate), a thinning
    bound, the extant sampling fraction ρ, and the removal probability ``r`` on sampling."""

    __slots__ = ("rates", "rate_bound", "rho", "removal")

    def __init__(self, model, present):
        self.removal = 1.0
        if isinstance(model, EpisodicFossilizedBirthDeath):
            model.validate()
            shifts, births, deaths, foss = (model.shifts, model.birth, model.death,
                                            model.fossilization)

            def rates(t):
                i = bisect.bisect_right(shifts, present - t)
                return births[i], deaths[i], foss[i]

            self.rates = rates
            self.rate_bound = max(b + d + f for b, d, f in zip(births, deaths, foss))
            self.rho = model.sampling
            self.removal = model.removal
        elif isinstance(model, FossilizedBirthDeath):
            model.validate()
            b, d, psi = model.birth, model.death, model.fossilization
            self.rates = lambda t: (b, d, psi)
            self.rate_bound = b + d + psi
            self.rho = model.sampling
            self.removal = model.removal
        elif isinstance(model, EpisodicBirthDeath):
            model.validate()
            shifts, births, deaths = model.shifts, model.birth, model.death

            def rates(t):
                i = bisect.bisect_right(shifts, present - t)
                return births[i], deaths[i], 0.0

            self.rates = rates
            self.rate_bound = max(b + d for b, d in zip(births, deaths))
            self.rho = model.rho
        elif isinstance(model, BirthDeath):  # Yule is a subclass
            b, d = model.birth, model.death
            self.rates = lambda t: (b, d, 0.0)
            self.rate_bound = b + d
            self.rho = 1.0
        else:
            raise NotImplementedError(
                f"forward simulation supports BirthDeath/Yule, EpisodicBirthDeath, "
                f"FossilizedBirthDeath and EpisodicFossilizedBirthDeath, not "
                f"{type(model).__name__}"
            )


def _grow(view, age, n_tips, rng, max_lineages):
    """One forward trial from a crown of two lineages (thinning handles time-varying rates).
    Returns ``(crown_node, end_time)`` or ``None`` to reject (extinct / <2 sampled survivors)."""
    root = TreeNode(name="", time=0.0)
    live = []
    for _ in range(2):
        child = TreeNode(name="", time=0.0)
        root.add_child(child)
        live.append(child)
    bound = view.rate_bound
    t = 0.0
    end = None
    while True:
        n = len(live)
        if n == 0:
            return None
        if n_tips is not None and n == n_tips:
            end = t
            break
        if n > max_lineages:
            raise RuntimeError(
                f"forward tree exceeded max_lineages={max_lineages}; explosive parameters "
                "(birth >> death over this age) — lower the age/rates or raise max_lineages"
            )
        dt = rng.exponential(1.0 / (n * bound))
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

    n_sampled = 0
    for node in live:  # survivors reach the present; sample each with probability ρ
        node.time = end
        if view.rho >= 1.0 or rng.random() < view.rho:
            node.is_extant = True
            node.sampled = True
            n_sampled += 1
        else:
            node.is_extant = False  # alive but unsampled -> a ghost tip
    if n_sampled < 2:
        return None
    return root, end


def _name(tree: Tree) -> None:
    leaf = internal = 1
    for node in tree.nodes_preorder():
        if node is tree.root:
            node.name = "root"
        elif node.is_leaf():
            node.name = f"n{leaf}"
            leaf += 1
        else:
            node.name = f"i{internal}"
            internal += 1


def simulate_species_tree_forward(
    model,
    *,
    age: float | None = None,
    n_tips: int | None = None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_attempts: int = 10_000,
    max_lineages: int = 1_000_000,
) -> Tree:
    """Simulate a **complete** species tree forward in time (extinct lineages included).

    Provide exactly one stopping condition:

    * ``age`` — grow for this crown age; the number of extant tips is random.
    * ``n_tips`` — grow until this many extant lineages first coexist; the age is random.
      (Constant-rate models only — the present must be fixed for episodic rates.)

    The tree is rooted at the crown (``time == 0``), the present is at ``total_age``, extinct
    leaves carry ``is_extant=False`` at their death times, and extant leaves ``is_extant=True``
    at the present. Under incomplete sampling (``EpisodicBirthDeath`` with ``ρ<1``), extant but
    unsampled lineages are marked ``is_extant=False`` too. So :func:`~zombi2.simulate_genomes`
    treats extinct/unsampled lineages as ghost transfer partners automatically. The run is
    conditioned on ≥2 sampled survivors.
    """
    if (age is None) == (n_tips is None):
        raise ValueError("provide exactly one of `age` or `n_tips`")
    if isinstance(model, (EpisodicBirthDeath, EpisodicFossilizedBirthDeath)):
        if n_tips is not None:
            raise NotImplementedError(
                "episodic forward simulation requires `age` (the present must be fixed to map "
                "age-before-present); n_tips mode is constant-rate only"
            )
    elif not isinstance(model, (BirthDeath, FossilizedBirthDeath)):
        raise NotImplementedError(
            f"forward simulation supports BirthDeath/Yule, EpisodicBirthDeath, "
            f"FossilizedBirthDeath and EpisodicFossilizedBirthDeath, not {type(model).__name__}"
        )
    model.validate()
    if age is not None and age <= 0:
        raise ValueError(f"age must be > 0, got {age}")
    if n_tips is not None and n_tips < 2:
        raise ValueError(f"n_tips must be >= 2, got {n_tips}")
    if rng is None:
        rng = np.random.default_rng(seed)

    view = _ForwardRates(model, present=(age if age is not None else 0.0))
    for _ in range(max_attempts):
        result = _grow(view, age, n_tips, rng, max_lineages)
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
