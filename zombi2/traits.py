"""Simulation of phenotypic **traits** on a timetree.

A trait is a value that evolves *along the branches* of a fixed tree — a body size, a
gene-expression level, a discrete character state (habitat, chromosome number, presence of
a structure). This is the classic phylogenetic-comparative-methods setting (Felsenstein
1985): given a tree, drop a value at the root and let it wander down every branch to the
tips.

The engine is the same preorder walk that :mod:`~zombi2.rate_variation` already uses — each
node inherits its parent's end-state and then evolves across its own branch — so a trait
model needs only to answer *"starting from state ``x``, where do I end up after ``dt`` units
of time?"*. Two families share that one driver:

* **continuous** traits (e.g. :class:`BrownianMotion`) draw the branch **endpoint** directly
  from the exact transition distribution — no path simulation is needed, and the node-by-node
  walk reproduces the exact multivariate-normal law over the tips;
* **discrete** traits (e.g. :class:`Mk`) simulate the continuous-time Markov jumps **exactly**
  along each branch (Gillespie), so the per-branch list of ``(state, duration)`` segments *is*
  the realized character history — a stochastic character map (Huelsenbeck 2003; Nielsen 2002)
  for free.

The result (:class:`TraitResult`) exposes the observable tip values, the (exact, not inferred)
ancestral node states, and — for discrete traits — the full branch histories.

    import zombi2 as z
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=30, age=5.0, seed=1)

    bm = z.simulate_traits(tree, z.BrownianMotion(sigma2=0.5), seed=1)
    bm.values                      # {extant leaf: float}

    mk = z.simulate_traits(tree, z.Mk.equal_rates(3, 0.4), seed=1)
    mk.values                      # {extant leaf: state index}
    mk.history                     # {node: [(state, duration), ...]} — the stochastic map
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .tree import Tree, TreeNode


# --------------------------------------------------------------------------- linear algebra
def _matrix_sqrt(S: np.ndarray) -> np.ndarray:
    """A symmetric square root of a symmetric PSD matrix (``M`` with ``M @ M == S``)."""
    w, V = np.linalg.eigh(S)
    return V @ np.diag(np.sqrt(np.clip(w, 0.0, None))) @ V.T


def _lyapunov(A: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Solve the continuous Lyapunov equation ``A V + V Aᵀ = Q`` for ``V``."""
    k = A.shape[0]
    K = np.kron(np.eye(k), A) + np.kron(A, np.eye(k))
    v = np.linalg.solve(K, Q.reshape(-1, order="F"))
    return v.reshape((k, k), order="F")


def _mvn(mean: np.ndarray, cov: np.ndarray, rng) -> np.ndarray:
    """Draw one multivariate normal, robust to tiny negative eigenvalues from round-off."""
    cov = (cov + cov.T) / 2.0
    w, V = np.linalg.eigh(cov)
    return mean + V @ (np.sqrt(np.clip(w, 0.0, None)) * rng.standard_normal(len(mean)))


# --------------------------------------------------------------------------- models
class BrownianMotion:
    """Brownian-motion evolution of a continuous trait (Felsenstein 1985).

    Along a branch of duration ``t`` the trait performs a random walk

        ``dX = trend·dt + σ·dW``,

    so the endpoint given the start ``x`` is normal with mean ``x + trend·t`` and variance
    ``sigma2·t``. Simulated node-by-node in preorder this reproduces the exact tip law: the
    tips are multivariate-normal with mean ``x0 + trend·(root-to-tip depth)`` and covariance
    ``sigma2 · C``, where ``C`` is the shared-path-length matrix of the tree.

    Parameters
    ----------
    sigma2:
        Diffusion rate (variance accrued per unit time). Must be ``>= 0``.
    x0:
        Trait value at the root (default ``0.0``). Override per run with
        ``simulate_traits(..., root_state=)``.
    trend:
        Directional drift per unit time (default ``0.0``); a non-zero value makes the walk
        biased (a "trend" model, e.g. ``phytools``).
    """

    kind = "continuous"

    def __init__(self, sigma2: float, x0: float = 0.0, trend: float = 0.0):
        if sigma2 < 0:
            raise ValueError(f"sigma2 must be >= 0, got {sigma2}")
        self.sigma2 = float(sigma2)
        self.x0 = float(x0)
        self.trend = float(trend)

    def root_value(self, rng):
        return self.x0

    def evolve(self, x, dt, t0, rng):
        """Endpoint of the branch: ``Normal(x + trend·dt, sigma2·dt)``. No segments."""
        std = (self.sigma2 * dt) ** 0.5
        end = x + self.trend * dt + (rng.normal(0.0, std) if std > 0.0 else 0.0)
        return end, None

    def __repr__(self) -> str:
        extra = f", trend={self.trend:g}" if self.trend else ""
        return f"BrownianMotion(sigma2={self.sigma2:g}, x0={self.x0:g}{extra})"


class Mk:
    """The Mk model: a continuous-time Markov chain over ``k`` discrete states (Lewis 2001).

    Transitions between states follow an instantaneous rate matrix ``Q`` (off-diagonals
    ``Q[i, j] >= 0`` are the rate ``i -> j``; each diagonal is set to ``-Σ_{j≠i} Q[i, j]`` so
    rows sum to zero). Along every branch the jumps are simulated **exactly** by the Gillespie
    algorithm, so the recorded ``(state, duration)`` segments are the realized history.

    This generalizes :class:`~zombi2.RateVariation` (whose ``Q`` is the banded nearest-neighbour
    walk over ordered rate bins) to an arbitrary transition structure. Convenience constructors
    cover the standard sub-models:

    * :meth:`equal_rates` — one shared rate (``ER``);
    * :meth:`symmetric` — ``Q[i, j] == Q[j, i]`` (``SYM``);
    * the raw constructor — any matrix, including all-rates-different (``ARD``) and ordered
      (meristic) characters.

    Parameters
    ----------
    Q:
        ``k x k`` rate matrix. Off-diagonals must be ``>= 0``; the diagonal is (re)computed.
    states:
        Optional labels for the ``k`` states (default ``0 .. k-1``); used for output only.
    root:
        Root-state policy: ``"uniform"`` (default, pick a state uniformly at random),
        ``"stationary"`` (draw from the chain's stationary distribution), an integer state
        index, or a length-``k`` probability vector. Override per run with
        ``simulate_traits(..., root_state=)``.
    """

    kind = "discrete"

    def __init__(self, Q, states=None, root="uniform"):
        Q = np.asarray(Q, dtype=float)
        if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
            raise ValueError(f"Q must be a square matrix, got shape {Q.shape}")
        k = Q.shape[0]
        if k < 2:
            raise ValueError("Mk needs at least 2 states")
        Q = Q.copy()
        offdiag = Q.copy()
        np.fill_diagonal(offdiag, 0.0)
        if np.any(offdiag < 0):
            raise ValueError("off-diagonal rates Q[i, j] must be >= 0")
        np.fill_diagonal(Q, -offdiag.sum(axis=1))  # rows sum to zero
        self.Q = Q
        self.k = k
        self.states = list(range(k)) if states is None else list(states)
        if len(self.states) != k:
            raise ValueError(f"states has length {len(self.states)}, expected {k}")
        self._root = self._parse_root(root)

    @classmethod
    def equal_rates(cls, k: int, rate: float = 1.0, states=None, root="uniform") -> "Mk":
        """``ER``: every allowed transition happens at the same ``rate``."""
        if k < 2:
            raise ValueError("Mk needs at least 2 states")
        if rate < 0:
            raise ValueError(f"rate must be >= 0, got {rate}")
        Q = np.full((k, k), float(rate))
        np.fill_diagonal(Q, 0.0)
        return cls(Q, states=states, root=root)

    @classmethod
    def symmetric(cls, rates, states=None, root="uniform") -> "Mk":
        """``SYM``: a symmetric rate matrix (``i -> j`` and ``j -> i`` share a rate).

        ``rates`` is the ``k x k`` matrix of off-diagonal rates; it is symmetrized
        (``(R + Rᵀ)/2``) and the diagonal is ignored.
        """
        R = np.asarray(rates, dtype=float)
        if R.ndim != 2 or R.shape[0] != R.shape[1]:
            raise ValueError("rates must be a square matrix")
        R = (R + R.T) / 2.0
        np.fill_diagonal(R, 0.0)
        return cls(R, states=states, root=root)

    @classmethod
    def ordered(cls, k: int, rate: float = 1.0, states=None, root="uniform") -> "Mk":
        """An **ordered** (meristic) character: transitions only between **adjacent** states
        (``i <-> i±1``), each at ``rate`` — a tridiagonal ``Q``. This is the character-state
        analogue of :class:`~zombi2.RateVariation`'s nearest-neighbour walk over ordered bins.
        """
        if k < 2:
            raise ValueError("Mk needs at least 2 states")
        if rate < 0:
            raise ValueError(f"rate must be >= 0, got {rate}")
        Q = np.zeros((k, k))
        for i in range(k - 1):
            Q[i, i + 1] = rate
            Q[i + 1, i] = rate
        return cls(Q, states=states, root=root)

    def _parse_root(self, root):
        """Return either an int state index, the string ``"stationary"``, or a prob vector."""
        if isinstance(root, str):
            if root == "uniform":
                return np.full(self.k, 1.0 / self.k)
            if root == "stationary":
                return "stationary"
            raise ValueError("root string must be 'uniform' or 'stationary'")
        if np.isscalar(root):
            r = int(root)
            if not (0 <= r < self.k):
                raise ValueError(f"root index must be in [0, {self.k - 1}], got {root}")
            return r
        p = np.asarray(root, dtype=float)
        if p.shape != (self.k,) or np.any(p < 0) or not np.isclose(p.sum(), 1.0):
            raise ValueError("root vector must be a length-k probability distribution")
        return p / p.sum()

    def stationary_distribution(self) -> np.ndarray:
        """The chain's stationary distribution ``π`` (``π Q = 0``, ``Σ π = 1``)."""
        A = np.vstack([self.Q.T[:-1], np.ones(self.k)])
        b = np.zeros(self.k)
        b[-1] = 1.0
        pi = np.linalg.solve(A, b)
        return np.clip(pi, 0.0, None) / np.clip(pi, 0.0, None).sum()

    def transition_matrix(self, t: float) -> np.ndarray:
        """``P(t) = exp(Q·t)`` — the probability of state ``j`` after time ``t`` from state ``i``.

        Computed by eigendecomposition (valid for the diagonalizable ``Q`` of these models);
        provided for users and for checking simulations.
        """
        vals, vecs = np.linalg.eig(self.Q * t)
        P = (vecs @ np.diag(np.exp(vals)) @ np.linalg.inv(vecs)).real
        return np.clip(P, 0.0, 1.0)

    def root_value(self, rng):
        r = self._root
        if isinstance(r, (int, np.integer)):
            return int(r)
        p = self.stationary_distribution() if isinstance(r, str) else r
        return int(rng.choice(self.k, p=p))

    def evolve(self, state, dt, t0, rng):
        """Exact Gillespie simulation of the CTMC along a branch of duration ``dt``.

        Returns ``(end_state, segments)`` where ``segments`` is a list of ``(state, duration)``
        pieces whose durations sum to ``dt`` — the realized character history on this branch.
        """
        segments = []
        elapsed = 0.0
        current = int(state)
        while True:
            rate_out = -self.Q[current, current]
            if rate_out <= 0.0:  # absorbing state
                segments.append((current, dt - elapsed))
                return current, segments
            wait = rng.exponential(1.0 / rate_out)
            if elapsed + wait >= dt:
                segments.append((current, dt - elapsed))
                return current, segments
            segments.append((current, wait))
            elapsed += wait
            probs = self.Q[current].copy()
            probs[current] = 0.0
            probs /= rate_out
            current = int(rng.choice(self.k, p=probs))

    def __repr__(self) -> str:
        return f"Mk(k={self.k})"


class CorrelatedBinary(Mk):
    """Correlated evolution of two binary characters (Pagel 1994).

    Two binary traits **X** and **Y** evolve jointly over the four states ``(X, Y)`` in
    ``{0,1}²``, changing **one trait at a time** (simultaneous double changes have rate 0).
    Each trait's gain/loss rate may depend on the *other* trait's current state — that
    dependence is exactly *correlated evolution*. When X's rates do not depend on Y and vice
    versa (see :meth:`independent`) the two traits evolve independently: Pagel's null model,
    against which the dependent fit is tested.

    It is a 4-state :class:`Mk` with the two double-transitions zeroed, so it simulates with
    :func:`simulate_traits` like any discrete trait; each node's value is the ``(X, Y)`` pair,
    read from ``result.values``. Decompose with ``x = state[0]``, ``y = state[1]``.

    The eight rates name the changing trait, its direction, and the *other* trait's state:

    Parameters
    ----------
    x_gain_y0, x_gain_y1:
        Rate of ``X: 0 → 1`` while ``Y`` is 0 / 1.
    x_loss_y0, x_loss_y1:
        Rate of ``X: 1 → 0`` while ``Y`` is 0 / 1.
    y_gain_x0, y_gain_x1:
        Rate of ``Y: 0 → 1`` while ``X`` is 0 / 1.
    y_loss_x0, y_loss_x1:
        Rate of ``Y: 1 → 0`` while ``X`` is 0 / 1.
    root:
        Root policy as in :class:`Mk`, or an ``(X, Y)`` tuple pinning the root pair.
    """

    # state index = 2*X + Y  ->  0:(0,0)  1:(0,1)  2:(1,0)  3:(1,1)
    STATES = [(0, 0), (0, 1), (1, 0), (1, 1)]

    def __init__(self, x_gain_y0, x_gain_y1, x_loss_y0, x_loss_y1,
                 y_gain_x0, y_gain_x1, y_loss_x0, y_loss_x1, root="uniform"):
        rates = [x_gain_y0, x_gain_y1, x_loss_y0, x_loss_y1,
                 y_gain_x0, y_gain_x1, y_loss_x0, y_loss_x1]
        if any(r < 0 for r in rates):
            raise ValueError("all transition rates must be >= 0")
        Q = np.zeros((4, 4))
        Q[0, 1] = y_gain_x0   # (0,0) -> (0,1): Y gains, X=0
        Q[0, 2] = x_gain_y0   # (0,0) -> (1,0): X gains, Y=0
        Q[1, 0] = y_loss_x0   # (0,1) -> (0,0): Y loses, X=0
        Q[1, 3] = x_gain_y1   # (0,1) -> (1,1): X gains, Y=1
        Q[2, 0] = x_loss_y0   # (1,0) -> (0,0): X loses, Y=0
        Q[2, 3] = y_gain_x1   # (1,0) -> (1,1): Y gains, X=1
        Q[3, 1] = x_loss_y1   # (1,1) -> (0,1): X loses, Y=1
        Q[3, 2] = y_loss_x1   # (1,1) -> (1,0): Y loses, X=1
        if isinstance(root, tuple):
            root = 2 * int(root[0]) + int(root[1])
        super().__init__(Q, states=self.STATES, root=root)

    @classmethod
    def independent(cls, x_gain, x_loss, y_gain, y_loss, root="uniform"):
        """Pagel's null model: X and Y evolve independently (each rate ignores the other trait)."""
        return cls(x_gain, x_gain, x_loss, x_loss,
                   y_gain, y_gain, y_loss, y_loss, root=root)

    def __repr__(self) -> str:
        return "CorrelatedBinary(Pagel 1994)"


class HiddenStateMk(Mk):
    """A discrete character with hidden rate classes — a hidden Markov model (corHMM).

    The observed character evolves over ``O`` states, but its transition rates depend on an
    unobserved **hidden class** (a rate category) drawn from ``H`` classes that themselves switch
    along the tree. The full state is the ``(observed, hidden)`` pair; only the observed part is
    reported by :meth:`TraitResult.labeled_values`, while :meth:`TraitResult.changes` shows the
    full history including hidden switches. This captures rate heterogeneity a plain :class:`Mk`
    cannot (Beaulieu et al. 2013).

    Parameters
    ----------
    observed_rates:
        One ``O x O`` rate matrix **per hidden class** (a length-``H`` sequence): the observed
        character's transition rates while in that class (e.g. a slow class and a fast class).
    hidden_rate:
        ``H x H`` matrix of switch rates between hidden classes (applied within any observed
        state), or a scalar for a symmetric all-to-all rate.
    observed_states, hidden_states:
        Optional labels (defaults ``0 .. O-1`` and ``0 .. H-1``).
    root:
        Root policy as in :class:`Mk`, over the ``O*H`` product states.
    """

    def __init__(self, observed_rates, hidden_rate, observed_states=None,
                 hidden_states=None, root="uniform"):
        mats = [np.asarray(R, dtype=float) for R in observed_rates]
        H = len(mats)
        if H < 1:
            raise ValueError("need at least one hidden class")
        O = mats[0].shape[0]
        if any(R.shape != (O, O) for R in mats):
            raise ValueError("each observed_rates matrix must be O x O with the same O")
        obs = list(range(O)) if observed_states is None else list(observed_states)
        hid = list(range(H)) if hidden_states is None else list(hidden_states)
        if len(obs) != O or len(hid) != H:
            raise ValueError("observed_states / hidden_states lengths must match the matrices")
        if np.isscalar(hidden_rate):
            HR = np.full((H, H), float(hidden_rate))
            np.fill_diagonal(HR, 0.0)
        else:
            HR = np.asarray(hidden_rate, dtype=float)
            if HR.shape != (H, H):
                raise ValueError("hidden_rate must be an H x H matrix or a scalar")

        n = O * H
        Q = np.zeros((n, n))
        for h in range(H):
            R = mats[h]
            for o in range(O):
                for o2 in range(O):
                    if o2 != o:
                        Q[o * H + h, o2 * H + h] = R[o, o2]
        for o in range(O):
            for h in range(H):
                for h2 in range(H):
                    if h2 != h:
                        Q[o * H + h, o * H + h2] = HR[h, h2]

        self._O, self._H, self._obs = O, H, obs
        states = [(obs[o], hid[h]) for o in range(O) for h in range(H)]
        super().__init__(Q, states=states, root=root)

    def discretize(self, index):
        """The observed state of a product-state ``index`` (hidden class collapsed)."""
        return self._obs[int(index) // self._H]

    def __repr__(self) -> str:
        return f"HiddenStateMk(observed={self._O}, hidden={self._H})"


class OrnsteinUhlenbeck:
    """Ornstein–Uhlenbeck evolution of a continuous trait (Hansen 1997; Butler & King 2004).

    The trait is pulled toward an optimum ``theta`` with strength ``alpha`` while it diffuses,

        ``dX = alpha·(theta − X)·dt + σ·dW``,

    a model of stabilizing selection / adaptation. The exact transition over a branch of
    duration ``t`` from ``x`` is normal with mean ``theta + (x − theta)·e^{−alpha·t}`` and
    variance ``sigma2/(2·alpha)·(1 − e^{−2·alpha·t})``; the stationary law is
    ``N(theta, sigma2/(2·alpha))``. As ``alpha → 0`` this approaches :class:`BrownianMotion`.

    Parameters
    ----------
    sigma2:
        Diffusion rate (``>= 0``).
    alpha:
        Mean-reversion / selection strength (``> 0``; use :class:`BrownianMotion` for ``0``).
    theta:
        Optimum value the trait is pulled toward.
    x0:
        Root value (default: ``theta`` — start at the optimum). Override per run with
        ``simulate_traits(..., root_state=)``.
    """

    kind = "continuous"

    def __init__(self, sigma2: float, alpha: float, theta: float, x0: float | None = None):
        if sigma2 < 0:
            raise ValueError(f"sigma2 must be >= 0, got {sigma2}")
        if alpha <= 0:
            raise ValueError(f"alpha must be > 0 (use BrownianMotion for alpha=0), got {alpha}")
        self.sigma2 = float(sigma2)
        self.alpha = float(alpha)
        self.theta = float(theta)
        self.x0 = self.theta if x0 is None else float(x0)

    def root_value(self, rng):
        return self.x0

    def evolve(self, x, dt, t0, rng):
        e = np.exp(-self.alpha * dt)
        mean = self.theta + (x - self.theta) * e
        var = self.sigma2 / (2.0 * self.alpha) * (1.0 - e * e)
        end = mean + (rng.normal(0.0, var ** 0.5) if var > 0.0 else 0.0)
        return end, None

    def __repr__(self) -> str:
        return (f"OrnsteinUhlenbeck(sigma2={self.sigma2:g}, alpha={self.alpha:g}, "
                f"theta={self.theta:g})")


class MultivariateBrownian:
    """Brownian motion of a **vector-valued** trait with a rate (covariance) matrix ``R``.

    A length-``k`` trait diffuses so the increment over a branch of duration ``t`` is
    ``MVN(trend·t, R·t)``. The off-diagonal ``R[a, b]`` couples dimensions ``a`` and ``b``, so
    this is the model of **correlated** continuous-trait evolution (``mvMORPH``, ``Rphylopars``):
    the tips are jointly multivariate-normal with covariance ``R ⊗ C`` (``C`` = the tree's
    shared-path-length matrix). Each node's value is a length-``k`` :class:`numpy.ndarray`.

    Parameters
    ----------
    R:
        ``k x k`` symmetric positive-semidefinite rate matrix.
    x0:
        Root vector (default: zeros).
    trend:
        Per-dimension directional drift (default: zeros).
    """

    kind = "continuous"

    def __init__(self, R, x0=None, trend=None):
        R = np.asarray(R, dtype=float)
        if R.ndim != 2 or R.shape[0] != R.shape[1]:
            raise ValueError(f"R must be a square matrix, got shape {R.shape}")
        k = R.shape[0]
        R = (R + R.T) / 2.0
        if np.linalg.eigvalsh(R).min() < -1e-9:
            raise ValueError("R must be positive semidefinite")
        self.k = k
        self.R = R
        self._sqrtR = _matrix_sqrt(R)
        self.x0 = np.zeros(k) if x0 is None else np.asarray(x0, dtype=float)
        self.trend = np.zeros(k) if trend is None else np.asarray(trend, dtype=float)
        if self.x0.shape != (k,) or self.trend.shape != (k,):
            raise ValueError(f"x0 and trend must be length-{k} vectors")

    def root_value(self, rng):
        return self.x0.copy()

    def evolve(self, x, dt, t0, rng):
        inc = self.trend * dt + (dt ** 0.5) * (self._sqrtR @ rng.standard_normal(self.k))
        return x + inc, None

    def __repr__(self) -> str:
        return f"MultivariateBrownian(k={self.k})"


class MultivariateOU:
    """Ornstein–Uhlenbeck evolution of a **vector-valued** trait (multivariate OU; ``mvMORPH``).

        ``dX = A·(theta − X)·dt + Σ^{1/2}·dW``,

    with mean-reversion matrix ``A`` (``alpha``), optimum vector ``theta``, and diffusion
    covariance ``R`` (``Σ``). The exact branch transition has mean
    ``theta + e^{−A·t}·(x − theta)`` and covariance ``V − e^{−A·t}·V·e^{−Aᵀ·t}``, where the
    stationary covariance ``V`` solves the Lyapunov equation ``A·V + V·Aᵀ = R``. For a scalar
    ``alpha`` (``A = alpha·I``) this gives ``V = R/(2·alpha)``, matching per-dimension OU with
    correlated diffusion.

    Parameters
    ----------
    R:
        ``k x k`` symmetric positive-semidefinite diffusion covariance.
    alpha:
        Mean reversion, as a scalar (isotropic ``alpha·I``), a length-``k`` vector (a diagonal
        ``A``), or a ``k x k`` matrix. Its eigenvalues must have positive real part (stable).
    theta:
        Optimum vector (length ``k``).
    x0:
        Root vector (default: ``theta``).
    """

    kind = "continuous"

    def __init__(self, R, alpha, theta, x0=None):
        R = np.asarray(R, dtype=float)
        if R.ndim != 2 or R.shape[0] != R.shape[1]:
            raise ValueError(f"R must be a square matrix, got shape {R.shape}")
        k = R.shape[0]
        R = (R + R.T) / 2.0
        if np.linalg.eigvalsh(R).min() < -1e-9:
            raise ValueError("R must be positive semidefinite")

        a = np.asarray(alpha, dtype=float)
        if a.ndim == 0:
            A = float(a) * np.eye(k)
        elif a.ndim == 1:
            if a.shape != (k,):
                raise ValueError(f"alpha vector must have length {k}")
            A = np.diag(a)
        elif a.shape == (k, k):
            A = a
        else:
            raise ValueError("alpha must be a scalar, a length-k vector, or a k x k matrix")
        if np.linalg.eigvals(A).real.min() <= 0:
            raise ValueError("alpha (A) must have eigenvalues with positive real part "
                             "(a stable, mean-reverting process)")

        theta = np.asarray(theta, dtype=float)
        if theta.shape != (k,):
            raise ValueError(f"theta must be a length-{k} vector")
        self.k = k
        self.R = R
        self.A = A
        self.theta = theta
        self.x0 = theta.copy() if x0 is None else np.asarray(x0, dtype=float)
        if self.x0.shape != (k,):
            raise ValueError(f"x0 must be a length-{k} vector")
        self.V = _lyapunov(A, R)  # stationary covariance
        w, Vec = np.linalg.eig(A)  # for a fast exp(-A·dt)
        self._wA, self._VA, self._VAinv = w, Vec, np.linalg.inv(Vec)

    def root_value(self, rng):
        return self.x0.copy()

    def _E(self, dt):
        """``exp(-A·dt)`` via the cached eigendecomposition of ``A``."""
        return (self._VA @ np.diag(np.exp(-self._wA * dt)) @ self._VAinv).real

    def evolve(self, x, dt, t0, rng):
        E = self._E(dt)
        mean = self.theta + E @ (x - self.theta)
        cov = self.V - E @ self.V @ E.T
        return _mvn(mean, cov, rng), None

    def __repr__(self) -> str:
        return f"MultivariateOU(k={self.k})"


class MultiOptimumOU:
    """OU with a different optimum on each painted regime of the tree (``OUwie``, ``ouch``).

    Different parts of the tree adapt toward different optima: each branch belongs to a
    **regime**, and the trait follows an Ornstein–Uhlenbeck process pulled toward that regime's
    optimum ``theta``. The regimes come from a **discrete stochastic map** — typically an
    :class:`Mk` trait simulated on the *same* tree — so a regime may even switch partway along a
    branch, and the OU is integrated exactly, piece by piece.

    Optionally ``alpha`` and ``sigma2`` also vary by regime (the ``OUMV`` / ``OUMA`` / ``OUMVA``
    variants); by default they are shared while only ``theta`` differs (``OUM``).

    Parameters
    ----------
    regimes:
        A **discrete** :class:`TraitResult` (e.g. from ``simulate_traits(tree, Mk...)``) whose
        per-branch history paints the regimes. Simulate this model on that *same* tree.
    theta:
        Optima, one per regime state (indexed by the regime's state index).
    alpha:
        Mean-reversion strength (``> 0``): a scalar shared by all regimes, or one per regime.
    sigma2:
        Diffusion rate (``>= 0``): a scalar, or one per regime.
    x0:
        Root value (default: the optimum of the regime at the root).
    """

    kind = "continuous"

    def __init__(self, regimes, theta, alpha, sigma2, x0=None):
        if getattr(regimes, "kind", None) != "discrete" or regimes.history is None:
            raise ValueError("regimes must be a discrete TraitResult carrying a stochastic map")
        self.regimes = regimes
        r = len(regimes.model.states)
        self.theta = np.asarray(theta, dtype=float)
        if self.theta.shape != (r,):
            raise ValueError(f"theta must give one optimum per regime ({r} regimes)")
        self.alpha = np.full(r, float(alpha)) if np.isscalar(alpha) else np.asarray(alpha, float)
        self.sigma2 = np.full(r, float(sigma2)) if np.isscalar(sigma2) else np.asarray(sigma2, float)
        if self.alpha.shape != (r,) or self.sigma2.shape != (r,):
            raise ValueError(f"alpha and sigma2 must be scalars or length-{r} per-regime vectors")
        if np.any(self.alpha <= 0):
            raise ValueError("alpha must be > 0")
        if np.any(self.sigma2 < 0):
            raise ValueError("sigma2 must be >= 0")
        root_regime = int(regimes.node_values[regimes.tree.root])
        self.x0 = float(self.theta[root_regime]) if x0 is None else float(x0)

    def root_value(self, rng):
        return self.x0

    def evolve_branch(self, node, x, rng):
        """Integrate the OU across this branch's regime segments (each piece exact)."""
        segs = self.regimes.history.get(node)
        if segs is None:
            raise ValueError("MultiOptimumOU must be simulated on the same tree its regimes "
                             "were painted on")
        for regime, dt in segs:
            th, al, s2 = self.theta[regime], self.alpha[regime], self.sigma2[regime]
            e = np.exp(-al * dt)
            mean = th + (x - th) * e
            var = s2 / (2.0 * al) * (1.0 - e * e)
            x = mean + (rng.normal(0.0, var ** 0.5) if var > 0.0 else 0.0)
        return x, None

    def __repr__(self) -> str:
        return f"MultiOptimumOU(n_regimes={len(self.theta)})"


class ThresholdModel:
    """Felsenstein's (2012) threshold model: a discrete state from an underlying continuous liability.

    An unobserved **liability** evolves by Brownian motion; the observed discrete state is the
    interval the liability currently falls in, cut by an ordered set of ``thresholds`` (``k-1``
    thresholds give ``k`` states). This links continuous and discrete evolution and naturally
    produces correlated / polymorphic discrete characters. Only the ratio of the thresholds to
    the diffusion scale is identifiable, so ``sigma2`` is fixed to ``1`` by default.

    The evolving value at each node is the liability (``result.values`` / ``ancestral_states``);
    the observed discrete state comes from :meth:`TraitResult.labeled_values` (or ``to_tsv`` /
    ``to_newick``, which report the discrete state).

    Parameters
    ----------
    thresholds:
        Strictly increasing cut points on the liability axis. For a binary trait, ``[0.0]``.
    sigma2:
        Liability diffusion rate (default ``1.0``).
    x0:
        Root liability (default ``0.0``).
    states:
        Optional labels for the ``len(thresholds)+1`` states (default ``0 .. k-1``).
    """

    kind = "continuous"  # the liability is continuous; the observed state is derived

    def __init__(self, thresholds, sigma2: float = 1.0, x0: float = 0.0, states=None):
        self.thresholds = np.asarray(thresholds, dtype=float)
        if self.thresholds.ndim != 1 or self.thresholds.size < 1:
            raise ValueError("thresholds must be a non-empty 1-D array of cut points")
        if np.any(np.diff(self.thresholds) <= 0):
            raise ValueError("thresholds must be strictly increasing")
        if sigma2 < 0:
            raise ValueError(f"sigma2 must be >= 0, got {sigma2}")
        self.sigma2 = float(sigma2)
        self.x0 = float(x0)
        k = self.thresholds.size + 1
        self.states = list(range(k)) if states is None else list(states)
        if len(self.states) != k:
            raise ValueError(f"states has length {len(self.states)}, expected {k}")

    def root_value(self, rng):
        return self.x0

    def evolve(self, x, dt, t0, rng):
        std = (self.sigma2 * dt) ** 0.5
        return x + (rng.normal(0.0, std) if std > 0.0 else 0.0), None

    def discretize(self, liability):
        """The observed discrete state: how many thresholds the liability exceeds."""
        return self.states[int(np.searchsorted(self.thresholds, liability))]

    def __repr__(self) -> str:
        return f"ThresholdModel(k={self.thresholds.size + 1})"


class EarlyBurst:
    """Early-burst / ACDC: Brownian motion whose rate changes exponentially through time
    (Blomberg et al. 2003; Harmon et al. 2010).

    The diffusion rate at absolute time ``t`` (root at 0) is ``σ²(t) = sigma2 · e^{rate·t}``.
    With ``rate < 0`` the rate **decays** — most divergence happens early, the signature of an
    adaptive radiation (an *early burst*); with ``rate > 0`` it **accelerates** (the AC of
    ACDC); ``rate = 0`` is plain Brownian motion. The variance accrued over a branch spanning
    ``[t1, t2]`` is the exact integral ``sigma2·(e^{rate·t2} − e^{rate·t1})/rate``, so the tips
    stay multivariate-normal.

    Parameters
    ----------
    sigma2:
        Diffusion rate **at the root** (``>= 0``).
    rate:
        Exponential rate of change of ``σ²`` through time; ``< 0`` = early burst, ``> 0`` =
        accelerating, ``0`` = Brownian motion.
    x0:
        Root value (default ``0.0``).
    trend:
        Optional directional drift per unit time (default ``0.0``).
    """

    kind = "continuous"

    def __init__(self, sigma2: float, rate: float, x0: float = 0.0, trend: float = 0.0):
        if sigma2 < 0:
            raise ValueError(f"sigma2 must be >= 0, got {sigma2}")
        self.sigma2 = float(sigma2)
        self.rate = float(rate)
        self.x0 = float(x0)
        self.trend = float(trend)

    def root_value(self, rng):
        return self.x0

    def evolve(self, x, dt, t0, rng):
        t1, t2 = t0, t0 + dt
        if self.rate == 0.0:
            var = self.sigma2 * dt
        else:
            var = self.sigma2 * (np.exp(self.rate * t2) - np.exp(self.rate * t1)) / self.rate
        end = x + self.trend * dt + (rng.normal(0.0, var ** 0.5) if var > 0.0 else 0.0)
        return end, None

    def __repr__(self) -> str:
        return f"EarlyBurst(sigma2={self.sigma2:g}, rate={self.rate:g})"


# --------------------------------------------------------------------------- result
@dataclass
class TraitResult:
    """The outcome of :func:`simulate_traits`.

    ``node_values`` maps **every** node to its (exact) trait value — internal nodes are the
    true ancestral states, not an inference. ``history`` is the per-branch stochastic character
    map for discrete models (``{node: [(state, duration), ...]}``) and ``None`` for continuous
    ones. Discrete states are stored as integer indices; :meth:`label` maps an index to its
    user-supplied state label.
    """

    tree: Tree
    model: object
    node_values: dict
    history: dict | None
    kind: str

    # --- observable / derived views ---------------------------------------
    @property
    def values(self) -> dict:
        """Observable tip values — the **extant** leaves only (the comparative-data vector).

        For discrete models these are raw state **indices**; call :meth:`labeled_values` (or
        :meth:`label`) to decode them to the model's state labels.
        """
        return {n: self.node_values[n] for n in self.tree.extant_leaves()}

    def labeled_values(self) -> dict:
        """Observable extant-tip values with discrete state indices decoded to their labels
        (identical to :attr:`values` for continuous traits)."""
        return {n: self.label(v) for n, v in self.values.items()}

    def leaf_values(self) -> dict:
        """Values at every leaf (including extinct/fossil tips)."""
        return {n: self.node_values[n] for n in self.tree.leaves()}

    def ancestral_states(self) -> dict:
        """Values at every internal node (exact ancestral states)."""
        return {n: self.node_values[n] for n in self.tree.internal_nodes()}

    def label(self, value):
        """Map a stored node value to its **observable** label: a discrete state index to its
        state label, a hidden-state or threshold model's value to the *observed* state (hidden
        classes / liabilities collapsed), else the value unchanged."""
        discretize = getattr(self.model, "discretize", None)
        if discretize is not None:
            return discretize(value)
        if self.kind == "discrete":
            return self.model.states[int(value)]
        return value

    def full_label(self, value):
        """The **complete** discrete state label (e.g. a hidden-state model's ``(observed,
        hidden)`` pair), without collapsing hidden classes; the value itself for continuous."""
        if self.kind == "discrete":
            return self.model.states[int(value)]
        return value

    def changes(self) -> list:
        """Discrete only: realized transitions as ``(node, time, from_label, to_label)``, in time
        order along the tree (the events of the stochastic map). Uses the **full** state label, so
        hidden-state transitions are visible."""
        if self.kind != "discrete":
            raise ValueError("changes() is only defined for discrete-trait models")
        out = []
        for node in self.tree.nodes_preorder():
            segs = self.history.get(node)
            if not segs or node.parent is None:
                continue
            t = node.parent.time
            for (state, dur), (nxt, _) in zip(segs, segs[1:]):
                t += dur
                out.append((node, t, self.full_label(state), self.full_label(nxt)))
        return out

    # --- I/O --------------------------------------------------------------
    def to_tsv(self, nodes: str = "extant") -> str:
        """Two-column ``node<TAB>trait`` table.

        ``nodes`` selects the rows: ``"extant"`` (default, observable tips), ``"leaves"``
        (all tips), or ``"all"`` (every node, i.e. tips + ancestral states).
        """
        if nodes == "extant":
            selected = self.tree.extant_leaves()
        elif nodes == "leaves":
            selected = self.tree.leaves()
        elif nodes == "all":
            selected = self.tree.nodes()
        else:
            raise ValueError("nodes must be 'extant', 'leaves', or 'all'")
        lines = ["node\ttrait"]
        for n in selected:
            v = self.node_values[n]
            cell = repr(self.label(v)) if isinstance(self.label(v), str) else _fmt(self.label(v))
            lines.append(f"{n.name}\t{cell}")
        return "\n".join(lines) + "\n"

    def to_newick(self, annotate: bool = True) -> str:
        """Newick with each node's trait value in a BEAST/FigTree comment ``[&trait=…]``.

        With ``annotate=False`` this is just the tree's own Newick.
        """
        if not annotate:
            return self.tree.to_newick()

        def rec(node: TreeNode) -> str:
            if node.children:
                inner = ",".join(rec(c) for c in node.children)
                s = f"({inner}){node.name}"
            else:
                s = node.name
            s += f"[&trait={_fmt(self.label(self.node_values[node]))}]"
            if node.parent is not None:
                s += f":{node.branch_length():.10g}"
            return s

        return rec(self.tree.root) + ";"


def _fmt(v) -> str:
    if isinstance(v, str):
        return v
    arr = np.asarray(v)
    if arr.ndim >= 1:  # a multivariate trait (floats) or a range tuple (labels) -> {a,b,c}
        return "{" + ",".join(_fmt(x) for x in arr.ravel()) + "}"
    if arr.dtype.kind == "f":
        return f"{float(v):.6g}"
    return str(v)


# --------------------------------------------------------------------------- driver
def simulate_traits(
    tree: Tree,
    model,
    *,
    root_state=None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
) -> TraitResult:
    """Evolve a single trait down ``tree`` under ``model`` and return a :class:`TraitResult`.

    The trait starts at the root (from ``model``'s own root policy, or ``root_state`` if given)
    and is evolved branch by branch in preorder — each node inherits its parent's end-state.
    Works on any :class:`~zombi2.tree.Tree`: a simulated species tree, or a gene tree loaded via
    :func:`~zombi2.read_newick`.

    Parameters
    ----------
    tree:
        The timetree to evolve the trait on.
    model:
        A trait model, e.g. :class:`BrownianMotion` (continuous) or :class:`Mk` (discrete).
    root_state:
        Optional explicit value/state to pin at the root, overriding the model's default.
    seed / rng:
        Reproducibility, as elsewhere in ZOMBI2 (pass a seed, or your own numpy ``Generator``).
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    node_values: dict = {}
    history: dict | None = {} if model.kind == "discrete" else None

    root = tree.root
    node_values[root] = model.root_value(rng) if root_state is None else root_state
    if history is not None:
        history[root] = []

    # branch-aware models (e.g. MultiOptimumOU, whose optimum changes along a branch) implement
    # evolve_branch(node, start, rng); the rest use the plain evolve(state, dt, t0, rng).
    branch_hook = getattr(model, "evolve_branch", None)
    for node in tree.nodes_preorder():
        if node.parent is None:
            continue
        start = node_values[node.parent]
        if branch_hook is not None:
            end, segs = branch_hook(node, start, rng)
        else:
            end, segs = model.evolve(start, node.branch_length(), node.parent.time, rng)
        node_values[node] = end
        if history is not None:
            history[node] = segs

    return TraitResult(tree=tree, model=model, node_values=node_values,
                       history=history, kind=model.kind)


def replicate_traits(tree: Tree, model, n: int, *, seed: int | None = None,
                     rng: np.random.Generator | None = None) -> list:
    """Simulate the **same** trait ``n`` times on ``tree`` under ``model`` — independent draws
    with identical parameters — returning a list of ``n`` :class:`TraitResult`.

    Reproducible given ``seed`` / ``rng`` (a single generator is advanced across the replicates,
    so the draws are independent). Useful for building a comparative dataset or an empirical null
    distribution: each replicate is a fresh realization of the trait on the fixed tree.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if rng is None:
        rng = np.random.default_rng(seed)
    return [simulate_traits(tree, model, rng=rng) for _ in range(n)]


# --------------------------------------------------------------------------- Pagel tree transforms
def _rebuild(tree: Tree, time_fn) -> Tree:
    """Copy ``tree`` with each non-root node's time set by ``time_fn(node, parent_new_time)``."""

    def rec(node: TreeNode, parent_new_time) -> TreeNode:
        nt = 0.0 if node.parent is None else time_fn(node, parent_new_time)
        copy = TreeNode(name=node.name, time=nt,
                        is_extant=node.is_extant, sampled=node.sampled)
        for child in node.children:
            copy.add_child(rec(child, nt))
        return copy

    new = Tree(rec(tree.root, None), tree.total_age)
    new.total_age = max((leaf.time for leaf in new.leaves()), default=0.0)
    return new


def pagel_lambda(tree: Tree, lam: float) -> Tree:
    """Pagel's **λ** (1999): scale phylogenetic signal, returning a new tree to evolve on.

    Internal (shared) branch depths are multiplied by ``lam`` while tip depths are held fixed,
    so a Brownian trait's between-species covariance is scaled by ``lam`` without changing its
    per-species variance. ``lam = 1`` is the original tree (full signal); ``lam = 0`` collapses
    every internal node to the root (a star tree — independent tips, no signal).
    """
    if not (0.0 <= lam <= 1.0):
        raise ValueError(f"lambda must be in [0, 1], got {lam}")
    return _rebuild(tree, lambda node, _pnt: (lam * node.time if node.children else node.time))


def pagel_delta(tree: Tree, delta: float) -> Tree:
    """Pagel's **δ** (1999): raise node depths to the power ``delta`` (root and tips fixed).

    ``delta > 1`` pushes divergence toward the present (late, species-specific evolution);
    ``delta < 1`` toward the root (early evolution); ``delta = 1`` is the original tree.
    """
    if delta <= 0:
        raise ValueError(f"delta must be > 0, got {delta}")
    T = tree.total_age
    if T <= 0:
        raise ValueError("tree has non-positive age")
    return _rebuild(tree, lambda node, _pnt: T * (node.time / T) ** delta)


def pagel_kappa(tree: Tree, kappa: float) -> Tree:
    """Pagel's **κ** (1999): raise each branch length to the power ``kappa``.

    ``kappa = 1`` is the original tree; ``kappa = 0`` sets every branch to length 1 (a
    *speciational* / punctuational model — change accrues per speciation event, not per unit
    time); intermediate values interpolate. Tip depths are not preserved.
    """
    if kappa < 0:
        raise ValueError(f"kappa must be >= 0, got {kappa}")
    return _rebuild(tree, lambda node, pnt: pnt + node.branch_length() ** kappa)
