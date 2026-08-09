"""Traits — helpers shared by the continuous and discrete engines (tree preorder, the correlation matrix and its symmetric square root, and the driver resolution a driven rate needs)."""

from __future__ import annotations


import numpy as np

from .._runtime.progress import track
from ..params.modifiers import Driven
from ..tree import Tree


def _preorder(tree: Tree, progress: bool = False):
    """Node ids in an order that visits every node **after its parent** (a valid preorder). The
    forward engine always gives a child a higher id than its parent, so ascending id order suffices
    — the same monotonic-id fact ``genomes.prune`` relies on in reverse. No recursion needed.

    Every trait engine walks the tree exactly this way, so ``progress`` is handled here once rather
    than around each of their loops."""
    return track(sorted(tree.nodes), "traits", unit="node", enabled=progress)



def _symmetric_sqrt(matrix: np.ndarray) -> np.ndarray:
    """A symmetric square root ``L`` of a symmetric PSD matrix (``L @ L == matrix``), via its
    eigendecomposition; tiny negative eigenvalues from round-off are clipped to zero."""
    w, V = np.linalg.eigh((matrix + matrix.T) / 2.0)
    return (V * np.sqrt(np.clip(w, 0.0, None))) @ V.T



def _correlation_matrix(traits: list, correlation) -> np.ndarray:
    """The ``k×k`` correlation matrix from pairwise ``{(a, b): ρ}`` — 1 on the diagonal, ρ off it
    (symmetric), 0 for unspecified pairs. Validates each ρ ∈ [−1, 1] and that the matrix is
    positive-semidefinite (the ρ values must be jointly consistent)."""
    idx = {t: i for i, t in enumerate(traits)}
    R = np.eye(len(traits))
    for pair, rho in (correlation or {}).items():
        if not (isinstance(pair, tuple) and len(pair) == 2):
            raise ValueError(f"correlation keys must be (trait_a, trait_b) pairs, got {pair!r}")
        a, b = pair
        if a not in idx or b not in idx:
            raise ValueError(f"correlation key {pair!r} names a trait not in {traits}")
        if a == b:
            raise ValueError(f"correlation key {pair!r} is a self-correlation")
        if isinstance(rho, bool) or not isinstance(rho, (int, float)) or not -1.0 <= rho <= 1.0:
            raise ValueError(f"correlation for {pair!r} must be a number in [−1, 1], got {rho!r}")
        R[idx[a], idx[b]] = R[idx[b], idx[a]] = float(rho)
    if float(np.linalg.eigvalsh(R).min()) < -1e-9:
        raise ValueError(
            "the correlation matrix is not positive-semidefinite — the given ρ values are jointly "
            "inconsistent (e.g. three traits cannot all be strongly negatively correlated)."
        )
    return R



def _driven_mods(rate) -> list:
    """The **driven** modifiers a rate carries — the ones written with ``scaled_by(driver,
    mapping)`` — or ``[]`` when it carries none. A non-empty list
    means the rate reads another level on each lineage, so the engine must thread a ``drivers``
    value and step where the driver switches."""
    return [m for m in rate.modifiers if isinstance(m, Driven)]



def _resolve_drivers(mods: list, tree: Tree, level: str) -> dict:
    """Resolve a rate's driven modifiers into one `~zombi2.params.conditioned.DriverTrajectory` per
    driver, keyed by the modifier's ``key`` — the per-lineage lookup the engine reads as it walks the
    tree. The genome level's shape (``genomes/family.py``): dedupe by ``key`` so a driver shared
    across rates resolves once, resolve each driver (a written trait log, a grown ``TraitsResult``
    handed over in memory, or a genome's ``presence`` / ``completion``), then check that the mapping
    can actually fire.

    A mapping whose states never occur in the driver would leave every lineage at the default factor
    — the run would be the undriven model wearing a driven rate — so it is refused here, naming the
    driver, rather than passed over in silence.

    ``level`` is which trait engine is asking (``"traits.continuous"`` / ``"traits.discrete"``), which
    is what lets a driver refuse a level that may not read it at all — a trait may read anything, but
    the resolver is one and the name is what makes that a decision rather than an omission.

    No driver ⇒ an empty dict, and the engine's loop stays exactly what it was."""
    if not mods:
        return {}
    from ..params.conditioned import check_mapping_fires, resolve_driver

    by_key: dict = {}
    for m in mods:
        by_key.setdefault(m.key, m)
    trajs = {key: resolve_driver(m.driver, tree, step=m.step, level=level)
             for key, m in by_key.items()}
    for m in mods:
        label = m.driver if isinstance(m.driver, str) else f"<{type(m.driver).__name__}>"
        check_mapping_fires(m.mapping, trajs[m.key].states(), driver_label=label)
    return trajs


