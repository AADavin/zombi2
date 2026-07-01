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

import numpy as np

from .species_model import BirthDeath
from .tree import Tree, TreeNode


def simulate_species_tree(
    model: BirthDeath,
    *,
    n_tips: int,
    age: float,
    age_type: str = "crown",
    seed: int | None = None,
    rng: np.random.Generator | None = None,
) -> Tree:
    """Simulate a reconstructed species tree (survivors only).

    Parameters
    ----------
    model:
        A species-tree model, e.g. :class:`~zombi2.BirthDeath` or :class:`~zombi2.Yule`.
    n_tips:
        Number of extant species ``N`` (>= 2) to condition on.
    age:
        Tree age. ``age_type="crown"`` (default): age of the root/MRCA (root at time 0,
        leaves at ``age``). ``age_type="stem"``: time of origin (a stem precedes the crown).
    seed / rng:
        Provide a seed (a fresh generator is made) or an explicit numpy Generator.
    """
    if rng is None:
        rng = np.random.default_rng(seed)
    model.validate()
    if n_tips < 2:
        raise ValueError(f"n_tips must be >= 2, got {n_tips}")
    if age <= 0:
        raise ValueError(f"age must be > 0, got {age}")
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
        i, j = (int(x) for x in rng.choice(len(active), size=2, replace=False))
        parent = TreeNode(name="", time=ft)
        parent.add_child(active[i])
        parent.add_child(active[j])
        active = [x for k, x in enumerate(active) if k not in (i, j)]
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
