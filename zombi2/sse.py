"""State-dependent speciation and extinction — trait-dependent diversification.

Unlike :func:`~zombi2.simulate_traits`, which drops a trait onto a *fixed* tree, here the
trait and the tree **co-evolve**: a lineage's discrete state sets its speciation and
extinction rates, so the character shapes the shape of the tree (Maddison, Midford & Otto
2007). This is the BiSSE / MuSSE family (``diversitree``).

The process runs **forward** from a crown of two lineages (the same convention as
:func:`~zombi2.simulate_species_tree(..., direction="forward")`): each live lineage carries a
state ``i`` and competes to

* **speciate** at rate ``λ_i`` — both daughters inherit its current state;
* **go extinct** at rate ``μ_i`` — the lineage dies (a leaf with ``is_extant=False``);
* **change state** ``i → j`` at rate ``Q[i, j]`` — anagenetic change along the branch.

The result is a **complete** tree (extinct lineages included) plus the realized character
history, returned as a :class:`~zombi2.traits.TraitResult` whose ``.tree`` is the simulated
tree. Its observable ``.values`` are the extant tips' states; ``z.prune(result.tree)`` gives
the reconstructed (survivors-only) tree for downstream analysis.

    import zombi2 as z
    res = z.simulate_sse(z.BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.2,
                                 q01=0.1, q10=0.1), age=4.0, seed=1)
    res.values             # {extant leaf: state}
    res.tree               # complete species tree (states drove its shape)
"""

from __future__ import annotations

import numpy as np

from .tree import Tree, TreeNode
from .species_forward import _name
from .traits import TraitResult


# --------------------------------------------------------------------------- models
class MuSSE:
    """Multi-state speciation and extinction: state-dependent birth/death over ``k`` states.

    Parameters
    ----------
    birth:
        Length-``k`` speciation rates ``λ_i`` (``>= 0``).
    death:
        Length-``k`` extinction rates ``μ_i`` (``>= 0``).
    Q:
        ``k x k`` anagenetic transition-rate matrix (off-diagonals ``>= 0``; the diagonal is
        recomputed so rows sum to zero), exactly as in :class:`~zombi2.Mk`.
    states:
        Optional labels for the ``k`` states (default ``0 .. k-1``).
    """

    kind = "discrete"  # so TraitResult treats states as discrete labels

    def __init__(self, birth, death, Q, states=None):
        self.lambdas = np.asarray(birth, dtype=float)
        self.mus = np.asarray(death, dtype=float)
        k = self.lambdas.shape[0]
        if self.lambdas.shape != (k,) or self.mus.shape != (k,):
            raise ValueError("birth and death must be equal-length 1-D rate vectors")
        if np.any(self.lambdas < 0) or np.any(self.mus < 0):
            raise ValueError("birth and death rates must be >= 0")
        Q = np.asarray(Q, dtype=float)
        if Q.shape != (k, k):
            raise ValueError(f"Q must be {k} x {k} to match the {k} rate states")
        Q = Q.copy()
        offdiag = Q.copy()
        np.fill_diagonal(offdiag, 0.0)
        if np.any(offdiag < 0):
            raise ValueError("off-diagonal transition rates Q[i, j] must be >= 0")
        np.fill_diagonal(Q, -offdiag.sum(axis=1))
        self.Q = Q
        self.k = k
        self.states = list(range(k)) if states is None else list(states)
        if len(self.states) != k:
            raise ValueError(f"states has length {len(self.states)}, expected {k}")

    def stationary_distribution(self) -> np.ndarray:
        """Stationary distribution of the character's transition matrix ``Q`` (``π Q = 0``)."""
        A = np.vstack([self.Q.T[:-1], np.ones(self.k)])
        b = np.zeros(self.k)
        b[-1] = 1.0
        pi = np.clip(np.linalg.solve(A, b), 0.0, None)
        return pi / pi.sum()

    def __repr__(self) -> str:
        return f"MuSSE(k={self.k})"


class BiSSE(MuSSE):
    """Binary-state speciation and extinction (Maddison, Midford & Otto 2007).

    Two states (``0`` and ``1``), each with its own speciation and extinction rate, and
    asymmetric transitions between them.

    Parameters
    ----------
    lambda0, lambda1:
        Speciation rates in states 0 and 1.
    mu0, mu1:
        Extinction rates in states 0 and 1.
    q01, q10:
        Transition rates ``0 -> 1`` and ``1 -> 0``.
    states:
        Optional labels for the two states (default ``(0, 1)``).
    """

    def __init__(self, lambda0, lambda1, mu0, mu1, q01, q10, states=(0, 1)):
        if q01 < 0 or q10 < 0:
            raise ValueError("transition rates q01, q10 must be >= 0")
        super().__init__(birth=[lambda0, lambda1], death=[mu0, mu1],
                         Q=[[0.0, q01], [q10, 0.0]], states=list(states))

    def __repr__(self) -> str:
        return (f"BiSSE(lambda0={self.lambdas[0]:g}, lambda1={self.lambdas[1]:g}, "
                f"mu0={self.mus[0]:g}, mu1={self.mus[1]:g}, "
                f"q01={self.Q[0, 1]:g}, q10={self.Q[1, 0]:g})")


# --------------------------------------------------------------------------- engine
class _Lineage:
    """A live lineage: its growing node, current state, and the segments of its branch so far."""

    __slots__ = ("node", "state", "seg_start", "segs")

    def __init__(self, node, state, seg_start):
        self.node = node
        self.state = state
        self.seg_start = seg_start
        self.segs = []

    def close(self, t, node_values, history):
        """End the branch at time ``t``: record the final segment and the node's state."""
        self.segs.append((self.state, t - self.seg_start))
        node_values[self.node] = self.state
        history[self.node] = self.segs


def _simulate_sse(model, age, n_tips, root_state, rng, max_lineages):
    """One forward trial from a crown of two lineages sharing the root state.

    Returns ``(root, end_time, node_values, history)`` or ``None`` to reject (whole clade
    extinct / fewer than two extant survivors).
    """
    lambdas, mus, Q, k = model.lambdas, model.mus, model.Q, model.k
    out_rate = -np.diag(Q)                       # total transition rate leaving each state
    bound = float(np.max(lambdas + mus + out_rate))
    if bound <= 0.0:
        raise ValueError("all rates are zero; nothing can happen")

    root = TreeNode(name="", time=0.0)
    node_values = {root: root_state}
    history = {root: []}
    live = []
    for _ in range(2):                           # crown: two lineages at time 0, root state
        child = TreeNode(name="", time=0.0)
        root.add_child(child)
        live.append(_Lineage(child, root_state, 0.0))

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
                f"SSE tree exceeded max_lineages={max_lineages}; explosive parameters — "
                "lower the age/rates or raise max_lineages"
            )
        dt = rng.exponential(1.0 / (n * bound))
        if age is not None and t + dt >= age:
            end = age
            break
        t += dt
        idx = int(rng.integers(n))
        lin = live[idx]
        i = lin.state
        lam, mu, ri = lambdas[i], mus[i], out_rate[i]
        total = lam + mu + ri
        if rng.random() >= total / bound:        # thinned out
            continue
        r = rng.random() * total
        if r < lam:                              # speciation: both daughters inherit state i
            lin.node.time = t
            lin.close(t, node_values, history)
            live[idx] = live[-1]
            live.pop()
            for _ in range(2):
                d = TreeNode(name="", time=t)
                lin.node.add_child(d)
                live.append(_Lineage(d, i, t))
        elif r < lam + mu:                       # extinction
            lin.node.time = t
            lin.node.is_extant = False
            lin.close(t, node_values, history)
            live[idx] = live[-1]
            live.pop()
        else:                                    # anagenetic transition i -> j
            probs = Q[i].copy()
            probs[i] = 0.0
            probs /= ri
            j = int(rng.choice(k, p=probs))
            lin.segs.append((i, t - lin.seg_start))
            lin.seg_start = t
            lin.state = j

    for lin in live:                             # survivors reach the present
        lin.node.time = end
        lin.node.is_extant = True
        lin.node.sampled = True
        lin.close(end, node_values, history)
    if len(live) < 2:
        return None
    return root, end, node_values, history


def simulate_sse(
    model,
    *,
    age: float | None = None,
    n_tips: int | None = None,
    root_state=None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_attempts: int = 10_000,
    max_lineages: int = 1_000_000,
) -> TraitResult:
    """Simulate a tree jointly with a trait that drives its diversification (BiSSE / MuSSE).

    Provide exactly one stopping condition:

    * ``age`` — grow for this crown age; the number of extant tips is random.
    * ``n_tips`` — grow until this many extant lineages first coexist; the age is random.

    The run starts from a crown of two lineages sharing the root state (``root_state``, an
    integer state index; if ``None``, drawn from the character's stationary distribution) and
    is conditioned on at least two extant survivors. Returns a
    :class:`~zombi2.traits.TraitResult` whose ``.tree`` is the **complete** simulated tree
    (extinct leaves carry ``is_extant=False``), ``.values`` are the extant tips' states, and
    ``.history`` is the realized character map along every branch.
    """
    if (age is None) == (n_tips is None):
        raise ValueError("provide exactly one of `age` or `n_tips`")
    if not isinstance(model, MuSSE):
        raise TypeError("model must be a BiSSE or MuSSE instance")
    if age is not None and age <= 0:
        raise ValueError(f"age must be > 0, got {age}")
    if n_tips is not None and n_tips < 2:
        raise ValueError(f"n_tips must be >= 2, got {n_tips}")
    if root_state is not None and not (0 <= int(root_state) < model.k):
        raise ValueError(f"root_state must be a state index in [0, {model.k - 1}]")
    if root_state is None and np.allclose(model.Q - np.diag(np.diag(model.Q)), 0.0):
        raise ValueError("with no state transitions (Q = 0) the stationary root state is "
                         "undefined; pass an explicit root_state")
    if rng is None:
        rng = np.random.default_rng(seed)

    for _ in range(max_attempts):
        r0 = (int(root_state) if root_state is not None
              else int(rng.choice(model.k, p=model.stationary_distribution())))
        result = _simulate_sse(model, age, n_tips, r0, rng, max_lineages)
        if result is not None:
            root, end_time, node_values, history = result
            tree = Tree(root, end_time)
            _name(tree)
            return TraitResult(tree=tree, model=model, node_values=node_values,
                               history=history, kind="discrete")

    raise RuntimeError(
        f"SSE simulation produced no surviving tree in {max_attempts} attempts "
        "(the clade kept going extinct); raise max_attempts or lower the extinction rates"
    )
