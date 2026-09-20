"""Traits — helpers shared by the continuous and discrete engines (tree preorder, the correlation matrix and its symmetric square root, and the driver resolution a driven rate needs)."""

from __future__ import annotations

import math

import numpy as np

from .._runtime.progress import track
from ..params.connection import Driven
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



def _expm(M: np.ndarray) -> np.ndarray:
    """The matrix exponential ``e^M``, in numpy alone (ZOMBI2 does not depend on scipy).

    Scaling and squaring with a [6/6] Padé approximant (Moler & Van Loan 2003, "method 3"): halve
    ``M`` until its norm is at most 0.5, where the approximant is accurate to double precision, then
    square the result back up as many times. It needs no eigendecomposition, so a drift matrix with
    complex or repeated eigenvalues is no special case."""
    M = np.asarray(M, dtype=float)
    norm = float(np.linalg.norm(M, np.inf))
    s = max(0, int(math.ceil(math.log2(norm / 0.5)))) if norm > 0.5 else 0
    A = M / (2.0 ** s)
    q = 6
    c = 1.0
    X = np.eye(len(A))
    N = np.eye(len(A))
    D = np.eye(len(A))
    for k in range(1, q + 1):
        c = c * (q - k + 1) / (k * (2 * q - k + 1))
        X = A @ X
        N = N + c * X
        D = D + (-1) ** k * c * X
    E = np.linalg.solve(D, N)
    for _ in range(s):
        E = E @ E
    return E



def _ou_transition(P: np.ndarray, Sigma: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """The exact branch transition of a multivariate OU ``dx = −P(x − θ)dt + Σ^{1/2}dW``: the decay
    ``e^{−P·dt}``, which takes the start's distance from θ to the mean's, and the covariance
    ``∫₀^dt e^{−Ps} Σ e^{−Pᵀs} ds`` accrued over the branch.

    Both come from one exponential of a block matrix (Van Loan 1978): the exponential of
    ``[[P, Σ], [0, −Pᵀ]]·dt`` holds ``e^{−Pᵀdt}`` in its lower-right block, whose transpose is the
    decay, and the decay times that block's upper-right one is the covariance. That needs no Lyapunov
    solve, so it holds for any ``P`` — one whose traits never settle, or settle while circling."""
    k = len(P)
    C = np.zeros((2 * k, 2 * k))
    C[:k, :k] = P
    C[:k, k:] = Sigma
    C[k:, k:] = -P.T
    E = _expm(C * dt)
    decay = E[k:, k:].T                     # e^{−P·dt}
    cov = decay @ E[:k, k:]
    return decay, (cov + cov.T) / 2.0      # symmetric in exact arithmetic; drop the round-off



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


