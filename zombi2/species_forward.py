"""Forward-in-time species-tree simulation.

The default :func:`~zombi2.simulate_species_tree` runs *backward* and yields the reconstructed
tree (surviving lineages only). This module runs the birth–death process *forward* and returns
the **complete** tree — extinct lineages included natively (``is_extant=False`` leaves at their
death times), no ghost-grafting needed. It complements :func:`~zombi2.add_ghost_lineages`: two
routes to the same complete-tree object, which the forward gene simulator then uses with
transfers from the dead for free.

Conventions (matching the backward crown tree): the returned tree is rooted at the **crown**
(first speciation) at ``time == 0``, with the present at ``total_age``. ``age`` is therefore the
crown age. The process is conditioned on the origin lineage speciating and on ≥2 extant
survivors. v1 supports the constant-rate :class:`~zombi2.BirthDeath`/:class:`~zombi2.Yule`.
"""

from __future__ import annotations

import numpy as np

from .species_model import BirthDeath
from .tree import Tree, TreeNode


def _grow(lam, mu, age, n_tips, rng, max_lineages):
    """One forward birth–death trial from a single origin lineage. Returns
    ``(crown_node, crown_time, end_time)`` or ``None`` to reject (extinct / <2 survivors)."""
    rate = lam + mu
    p_birth = lam / rate
    origin = TreeNode(name="", time=0.0)
    live = [origin]
    t = 0.0
    crown = None
    end = None
    while True:
        n = len(live)
        if n == 0:
            return None  # died out
        if crown is None and n >= 2:
            crown = t  # the origin has speciated → crown established
        if n_tips is not None and crown is not None and n == n_tips:
            end = t
            break
        if n > max_lineages:
            raise RuntimeError(
                f"forward tree exceeded max_lineages={max_lineages}; explosive parameters "
                "(birth >> death over this age) — lower the age/rates or raise max_lineages"
            )
        dt = rng.exponential(1.0 / (n * rate))
        if age is not None and crown is not None and t + dt >= crown + age:
            end = crown + age
            break
        t += dt
        i = int(rng.integers(n))
        node = live[i]
        node.time = t
        live[i] = live[-1]
        live.pop()
        if rng.random() < p_birth:  # speciation
            a = TreeNode(name="", time=t)
            b = TreeNode(name="", time=t)
            node.add_child(a)
            node.add_child(b)
            live.append(a)
            live.append(b)
        else:  # extinction
            node.is_extant = False

    if len(live) < 2:  # need a non-trivial surviving tree (also guarantees a crown)
        return None
    for node in live:  # survivors become extant tips at the present
        node.time = end
        node.is_extant = True
    return origin, crown, end


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
    model: BirthDeath,
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

    The tree is rooted at the crown (``time == 0``), the present is at ``total_age``, extinct
    leaves carry ``is_extant=False`` at their death times, and extant leaves ``is_extant=True``
    at the present — so :func:`~zombi2.simulate_genomes` treats extinct lineages as ghost
    transfer donors/recipients automatically. The run is conditioned on ≥2 extant survivors;
    high extinction may need more ``max_attempts``.

    v1 supports the constant-rate :class:`~zombi2.BirthDeath`/:class:`~zombi2.Yule`.
    """
    if (age is None) == (n_tips is None):
        raise ValueError("provide exactly one of `age` or `n_tips`")
    if not isinstance(model, BirthDeath):  # Yule is a subclass -> supported
        raise NotImplementedError(
            "forward simulation currently supports the constant-rate BirthDeath/Yule model; "
            "episodic forward rates are a planned follow-up"
        )
    model.validate()
    if age is not None and age <= 0:
        raise ValueError(f"age must be > 0, got {age}")
    if n_tips is not None and n_tips < 2:
        raise ValueError(f"n_tips must be >= 2, got {n_tips}")
    if rng is None:
        rng = np.random.default_rng(seed)

    lam, mu = model.birth, model.death
    for _ in range(max_attempts):
        result = _grow(lam, mu, age, n_tips, rng, max_lineages)
        if result is None:
            continue
        crown_node, crown_time, end_time = result
        tree = Tree(crown_node, end_time - crown_time)
        for node in tree.nodes_preorder():  # shift so the crown sits at time 0
            node.time -= crown_time
        _name(tree)
        return tree

    raise RuntimeError(
        f"forward simulation produced no surviving tree in {max_attempts} attempts "
        "(the process kept going extinct); raise max_attempts, lower death, or use the "
        "backward simulate_species_tree"
    )
