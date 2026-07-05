"""Algorithm 1 — backward reconstructed birth-death, conditioned on N extant tips.

Under a constant-rate birth-death process the reconstructed tree (extant lineages only)
is a *coalescent point process*: internal node times are i.i.d. from a single
distribution and the ranked topology is a uniform coalescence (Hartmann, Wong & Stadler
2010; Stadler 2011). We therefore (1) draw internal node ages via the model's
inverse-CDF sampler and (2) assemble a random ranked tree from those times.

A recommended extra validation (not runnable here) is a two-sample comparison against
TreeSim's ``sim.bd.taxa`` on identical (birth, death, N, age).
"""

from __future__ import annotations

import math

import numpy as np

from .species_model import BirthDeath, CladeShiftBirthDeath, ClaDS, DiversityDependent
from .tree import Tree, TreeNode


def _check_age(age: float) -> None:
    """Common ``age`` validation (finite + positive). NaN silently passes ``age <= 0``."""
    if not math.isfinite(age):
        raise ValueError(f"age must be a finite number, got {age}")
    if age <= 0:
        raise ValueError(f"age must be > 0, got {age}")


def _check_n_tips(n_tips, max_lineages: int) -> None:
    """Common ``n_tips`` validation (whole number, >= 2, within the lineage cap)."""
    if isinstance(n_tips, float):
        if not n_tips.is_integer():
            raise ValueError(f"n_tips must be a whole number, got {n_tips}")
        n_tips = int(n_tips)
    if not isinstance(n_tips, (int, np.integer)):
        raise ValueError(f"n_tips must be an integer, got {n_tips!r}")
    if n_tips < 2:
        raise ValueError(f"n_tips must be >= 2, got {n_tips}")
    if n_tips > max_lineages:
        raise ValueError(
            f"n_tips ({n_tips}) exceeds max_lineages ({max_lineages}); assembling a tree this "
            "large would be extremely slow — lower n_tips or raise max_lineages"
        )


def simulate_species_tree(
    model,
    *,
    n_tips: int | None = None,
    age: float | None = None,
    direction: str = "backward",
    age_type: str = "crown",
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_attempts: int = 10_000,
    max_lineages: int = 1_000_000,
) -> Tree:
    """Simulate a species tree under a birth–death ``model``.

    Parameters
    ----------
    direction:
        ``"backward"`` (default) samples the **reconstructed** tree (survivors only) as a
        coalescent point process, conditioned on **both** ``n_tips`` and ``age``.
        ``"forward"`` simulates the process forward and returns the **complete** tree (extinct
        and fossil lineages included), conditioned on **exactly one** of ``age`` or ``n_tips``.
    n_tips, age:
        Number of extant species ``N`` (>= 2) and the crown age. Backward needs both; forward
        needs one (the other is then random). ``age_type`` (``"crown"``/``"stem"``) applies to
        backward only.
    seed / rng:
        Provide a seed (a fresh generator is made) or an explicit numpy Generator.
    max_attempts, max_lineages:
        Forward-only safety bounds (rejection retries / lineage-count cap).

    Notes
    -----
    Forward-only model features (``fossilization`` > 0, ``removal`` != 1) and constant-rate
    incomplete sampling (``sampling_fraction`` < 1) are rejected under ``direction="backward"``.
    """
    if direction not in ("backward", "forward"):
        raise ValueError(f"direction must be 'backward' or 'forward', got {direction!r}")
    if rng is None:
        rng = np.random.default_rng(seed)

    if direction == "forward":
        from .species_forward import simulate_forward
        return simulate_forward(model, age=age, n_tips=n_tips, rng=rng,
                                max_attempts=max_attempts, max_lineages=max_lineages)

    # --- backward: reconstructed tree conditioned on (n_tips, age) ---------
    if isinstance(model, (ClaDS, DiversityDependent, CladeShiftBirthDeath)):
        raise ValueError(
            f"{type(model).__name__} has per-lineage/diversity-dependent rates with no closed-form "
            "reconstructed CDF; it is forward-only — use direction='forward'"
        )
    model.validate()
    if n_tips is None or age is None:
        raise ValueError("backward simulation needs both `n_tips` and `age`")
    foss = model.fossilization
    if (sum(foss) if isinstance(foss, list) else foss) > 0 or getattr(model, "removal", 1.0) != 1.0:
        raise ValueError("fossilization / removal are forward-only; use direction='forward'")
    if getattr(model, "mass_extinctions", None):
        raise ValueError(
            "mass_extinctions kill real lineages forward in time and are not represented in the "
            "backward reconstructed tree; use direction='forward'"
        )
    if isinstance(model, BirthDeath) and model.sampling_fraction < 1.0:
        raise ValueError(
            "constant-rate backward sampling assumes complete sampling (ρ=1); use "
            "EpisodicBirthDeath for incomplete sampling, or direction='forward'"
        )
    _check_n_tips(n_tips, max_lineages)
    _check_age(age)
    n_tips = int(n_tips)
    if age_type not in ("crown", "stem"):
        raise ValueError(f"age_type must be 'crown' or 'stem', got {age_type!r}")

    A = float(age)
    if age_type == "crown":  # root fixed at time 0; sample the other N-2 internal nodes
        ages = [model.sample_internal_age(rng.random(), A) for _ in range(n_tips - 2)]
        internal_times = sorted([A - a for a in ages] + [0.0])
    else:  # stem: all N-1 internal nodes i.i.d.
        ages = [model.sample_internal_age(rng.random(), A) for _ in range(n_tips - 1)]
        internal_times = sorted(A - a for a in ages)

    return _assemble(n_tips, A, internal_times, rng)


def _assemble(n_tips: int, total_age: float, internal_times: list[float], rng) -> Tree:
    """Coalesce N leaves at the given internal times (most recent first)."""
    active: list[TreeNode] = [
        TreeNode(name=f"n{i + 1}", time=total_age, is_extant=True) for i in range(n_tips)
    ]
    for ft in sorted(internal_times, reverse=True):
        i, j = sorted(int(x) for x in rng.choice(len(active), size=2, replace=False))  # i < j
        parent = TreeNode(name="", time=ft)
        parent.add_child(active[i])
        parent.add_child(active[j])
        # remove the two coalesced lineages in O(1) each (swap-with-last, larger index first)
        active[j] = active[-1]
        active.pop()
        active[i] = active[-1]
        active.pop()
        active.append(parent)

    root = active[0]
    tree = Tree(root, total_age)

    counter = 1
    for node in tree.nodes_preorder():
        if node.is_leaf():
            continue
        if node is root:
            node.name = "root"
        else:
            node.name = f"i{counter}"
            counter += 1
    return tree
