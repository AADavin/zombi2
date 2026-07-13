"""State-dependent speciation and extinction — trait-dependent diversification.

Unlike :func:`~zombi2.simulate_traits`, which drops a trait onto a *fixed* tree, here the
trait and the tree **co-evolve**: a lineage's discrete state sets its speciation and
extinction rates, so the character shapes the shape of the tree (Maddison, Midford & Otto
2007). This is the BiSSE / MuSSE family (``diversitree``), the ``traits:species`` edge of the
:doc:`coevolution model <models/coevolution>`.

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

from zombi2.tree import Tree, TreeNode
from zombi2.species.forward import _name
from zombi2.traits.models import Cladogenesis, TraitResult


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
        try:
            pi = np.clip(np.linalg.solve(A, b), 0.0, None)
        except np.linalg.LinAlgError:
            # a reducible character (a frozen/absorbing state, or several absorbing states) has no
            # unique stationary distribution; fall back to uniform. This vector only seeds the root
            # state, so a degenerate character need not crash the simulation with a LinAlgError.
            pi = np.ones(self.k)
        s = pi.sum()
        return pi / s if s > 0 else np.full(self.k, 1.0 / self.k)

    def null(self, kind="neutral", *, n_hidden=2, hidden_switch=None, **kwargs):
        """Return a decoupled **null** model — the same process with the ``traits:species``
        arrow (character → diversification) cut. See :doc:`the null-models guide
        </guide/coevolution_nulls>`.

        ``kind``:

        * ``"neutral"`` — the character no longer sets the rates: every state gets the same
          speciation and extinction rate (the mean of the per-state rates), while the
          anagenetic ``Q`` is left untouched so the character still evolves. This is the
          naive constant-rate null.
        * ``"cid"`` — a **k-state Character-Independent Diversification** null: the observed
          character's rate spread is reassigned to ``n_hidden`` hidden classes (spanning the
          observed speciation range, with a shared extinction), uncorrelated with the observed
          state, which keeps its own transitions but no longer sets the rates. For a **binary**
          character this returns the dedicated :class:`CID`; use :meth:`BiSSE.null`.

        There is no ``"timing"`` null for this edge (no at-speciation event to redistribute).
        """
        kind = kind.lower()
        if kind == "neutral":
            lam = float(self.lambdas.mean())
            mu = float(self.mus.mean())
            return MuSSE(birth=[lam] * self.k, death=[mu] * self.k, Q=self.Q,
                         states=list(self.states))
        if kind == "cid":
            H = int(n_hidden)
            if H < 2:
                raise ValueError(f"a CID null needs n_hidden >= 2, got {H}")
            offdiag = self.Q - np.diag(np.diag(self.Q))
            switch = (float(hidden_switch) if hidden_switch is not None
                      else max(float(offdiag.max()), 1e-9))
            lam_levels = np.linspace(float(self.lambdas.min()), float(self.lambdas.max()), H)
            mu_const = float(self.mus.mean())
            k, n = self.k, self.k * H
            birth = np.zeros(n)
            death = np.zeros(n)
            Q = np.zeros((n, n))
            for o in range(k):
                for h in range(H):
                    idx = o * H + h
                    birth[idx] = lam_levels[h]                  # rate set by the HIDDEN class only
                    death[idx] = mu_const
                    for o2 in range(k):                          # observed transitions within class
                        if o2 != o and self.Q[o, o2] != 0.0:
                            Q[idx, o2 * H + h] = self.Q[o, o2]
                    for h2 in range(H):                          # hidden switches within obs. state
                        if h2 != h:
                            Q[idx, o * H + h2] = switch
            return _CharacterIndependentMuSSE(birth, death, Q, list(self.states), H)
        if kind == "timing":
            raise ValueError(
                "traits:species has no 'timing' null (there is no at-speciation event to "
                "decouple); use kind='neutral' or 'cid'")
        raise ValueError(f"unknown null kind {kind!r}; expected 'neutral' or 'cid'")

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

    def null(self, kind="cid", *, n_hidden=2, hidden_switch=None, **kwargs):
        """Return a decoupled **null** for the ``traits:species`` arrow.

        ``kind``:

        * ``"cid"`` *(default)* — a :class:`CID` null: the observed rate spread is preserved
          but reassigned to a **hidden** class uncorrelated with the observed trait, which
          keeps its own transitions (``q01``/``q10``) but no longer sets the rates. With
          ``n_hidden=2`` (the canonical **CID-2**) the two hidden regimes carry this model's
          two observed rate pairs exactly; ``n_hidden>2`` spreads that many regimes across the
          same rate range (e.g. ``n_hidden=4`` = **CID-4**). ``hidden_switch`` sets how fast
          lineages move between regimes (default: the larger observed transition rate).
        * ``"neutral"`` — the constant-rate null (see :meth:`MuSSE.null`).
        """
        kind = kind.lower()
        if kind == "cid":
            q01, q10 = float(self.Q[0, 1]), float(self.Q[1, 0])
            switch = float(hidden_switch) if hidden_switch is not None else max(q01, q10, 1e-9)
            l0, l1 = float(self.lambdas[0]), float(self.lambdas[1])
            m0, m1 = float(self.mus[0]), float(self.mus[1])
            if n_hidden < 2:
                raise ValueError(f"a CID null needs n_hidden >= 2, got {n_hidden}")
            if n_hidden == 2:                       # CID-2: the two observed rate pairs exactly
                rates = [(l0, m0), (l1, m1)]
            else:                                   # CID-H: H regimes bracketing the rate range
                lam = np.linspace(min(l0, l1), max(l0, l1), n_hidden)
                mu = np.linspace(min(m0, m1), max(m0, m1), n_hidden)
                rates = list(zip(lam.tolist(), mu.tolist()))
            return CID(rates, switch, q01=q01, q10=q10)
        return super().null(kind, **kwargs)

    def __repr__(self) -> str:
        return (f"BiSSE(lambda0={self.lambdas[0]:g}, lambda1={self.lambdas[1]:g}, "
                f"mu0={self.mus[0]:g}, mu1={self.mus[1]:g}, "
                f"q01={self.Q[0, 1]:g}, q10={self.Q[1, 0]:g})")


class HiSSE(MuSSE):
    """Hidden-State Speciation and Extinction (Beaulieu & O'Meara 2016).

    Extends :class:`BiSSE` with unobserved **hidden classes**: each of the two observed states
    comes in ``H`` hidden variants with their own diversification rates, so rate heterogeneity is
    not falsely pinned on the observed character (the classic failure mode of a raw BiSSE fit).
    Give one :class:`BiSSE` per hidden class — a diversification "regime" — plus the rates at
    which lineages switch between classes. Simulate with :func:`simulate_sse`; the tips report the
    **observed** state (the hidden class is collapsed by :meth:`~zombi2.TraitResult.labeled_values`,
    but visible via the full ``(observed, hidden)`` node values).

    Parameters
    ----------
    classes:
        A sequence of :class:`BiSSE` models, one per hidden class.
    hidden_transition:
        ``H x H`` matrix of switch rates between hidden classes (applied within an observed
        state), or a scalar for a symmetric all-to-all rate.
    hidden_states:
        Optional labels for the ``H`` hidden classes (default ``0 .. H-1``).
    """

    def __init__(self, classes, hidden_transition, hidden_states=None):
        classes = list(classes)
        H = len(classes)
        if H < 1:
            raise ValueError("need at least one hidden class")
        if any(not isinstance(c, BiSSE) for c in classes):
            raise TypeError("classes must be BiSSE instances (one diversification regime each)")
        hid = list(range(H)) if hidden_states is None else list(hidden_states)
        if len(hid) != H:
            raise ValueError("hidden_states length must match the number of classes")
        if np.isscalar(hidden_transition):
            HR = np.full((H, H), float(hidden_transition))
            np.fill_diagonal(HR, 0.0)
        else:
            HR = np.asarray(hidden_transition, dtype=float)
            if HR.shape != (H, H):
                raise ValueError("hidden_transition must be an H x H matrix or a scalar")

        n = 2 * H
        birth = np.zeros(n)
        death = np.zeros(n)
        Q = np.zeros((n, n))
        for h, c in enumerate(classes):
            for o in (0, 1):
                birth[o * H + h] = c.lambdas[o]
                death[o * H + h] = c.mus[o]
            Q[0 * H + h, 1 * H + h] = c.Q[0, 1]     # observed 0->1 within this class
            Q[1 * H + h, 0 * H + h] = c.Q[1, 0]     # observed 1->0 within this class
        for o in (0, 1):
            for h in range(H):
                for h2 in range(H):
                    if h2 != h:
                        Q[o * H + h, o * H + h2] = HR[h, h2]

        self._H = H
        states = [(o, hid[h]) for o in (0, 1) for h in range(H)]
        super().__init__(birth=birth, death=death, Q=Q, states=states)

    def discretize(self, index):
        """The observed state (0 or 1) of a product-state ``index`` (hidden class collapsed)."""
        return int(index) // self._H

    def null(self, kind="neutral", **kwargs):
        """The ``neutral`` null of a HiSSE **collapses the hidden classes** — they existed only to
        carry rate heterogeneity, which the neutral null removes — to a plain constant-rate
        :class:`BiSSE` on the observed binary character (so the output carries no hidden state).
        ``cid``/``timing`` are not defined (a HiSSE is already a hidden-state model)."""
        if kind.lower() == "neutral":
            lam = float(self.lambdas.mean())
            mu = float(self.mus.mean())
            q01 = float(self.Q[0, self._H])          # observed 0->1 within hidden class 0
            q10 = float(self.Q[self._H, 0])          # observed 1->0 within hidden class 0
            return BiSSE(lam, lam, mu, mu, q01, q10)
        raise ValueError("a HiSSE is already a hidden-state model; its only null is 'neutral' "
                         "(a constant-rate BiSSE)")

    def __repr__(self) -> str:
        return f"HiSSE(hidden={self._H})"


class CID(HiSSE):
    """Character-Independent Diversification — the honest null for the ``traits:species`` edge.

    A :class:`HiSSE` constrained so that **within each hidden class the two observed states
    share the same rates**: diversification varies across the unobserved *hidden* classes only,
    never with the observed character. The tree therefore carries genuine rate heterogeneity
    (fast and slow clades) that a raw :class:`BiSSE` fit cannot legitimately pin on the observed
    trait — the null of Beaulieu & O'Meara (2016), and the model behind
    :meth:`BiSSE.null(kind="cid") <BiSSE.null>`.

    ``CID-2`` (two hidden classes) matches a :class:`BiSSE`'s rate spread; ``CID-4`` (four)
    matches a full :class:`HiSSE`. See :doc:`the null-models guide </guide/coevolution_nulls>`.

    Parameters
    ----------
    class_rates:
        One ``(lambda, mu)`` pair per hidden class — the speciation and extinction rate of
        *both* observed states in that class.
    hidden_transition:
        ``H x H`` switch-rate matrix between hidden classes, or a scalar for a symmetric
        rate — how fast lineages move between diversification regimes.
    q01, q10:
        Anagenetic transition rates of the **observed** binary character (shared across
        classes). At least one must be ``> 0`` so the observed character still evolves; a
        frozen character would be trivially, uninterestingly independent of the rates.
    hidden_states:
        Optional labels for the ``H`` hidden classes (default ``0 .. H-1``).

    Examples
    --------
    >>> CID.two(lambda_slow=0.5, lambda_fast=2.0, mu=0.2, switch=0.15)   # CID-2
    CID(classes=2)
    """

    def __init__(self, class_rates, hidden_transition, *, q01=0.1, q10=0.1, hidden_states=None):
        rates = [tuple(float(v) for v in r) for r in class_rates]
        if len(rates) < 1:
            raise ValueError("need at least one hidden class")
        if any(len(r) != 2 for r in rates):
            raise ValueError("each class_rates entry must be a (lambda, mu) pair")
        if q01 < 0 or q10 < 0:
            raise ValueError("observed transition rates q01, q10 must be >= 0")
        if q01 == 0 and q10 == 0:
            raise ValueError(
                "a CID null needs the observed character to still evolve (q01 or q10 > 0); "
                "with both 0 the observed state is frozen and trivially independent")
        classes = [BiSSE(lam, lam, mu, mu, q01, q10) for (lam, mu) in rates]
        super().__init__(classes, hidden_transition, hidden_states=hidden_states)
        self._class_rates = rates
        self._q = (float(q01), float(q10))

    @classmethod
    def two(cls, lambda_slow, lambda_fast, *, mu=0.0, mu_slow=None, mu_fast=None,
            switch=0.1, q01=0.1, q10=0.1, hidden_states=None):
        """**CID-2**: two hidden diversification regimes (slow, fast). ``mu`` is a shared
        extinction rate; override per class with ``mu_slow`` / ``mu_fast``."""
        ms = mu if mu_slow is None else mu_slow
        mf = mu if mu_fast is None else mu_fast
        return cls([(lambda_slow, ms), (lambda_fast, mf)], switch,
                   q01=q01, q10=q10, hidden_states=hidden_states)

    @classmethod
    def four(cls, rates, *, switch=0.1, q01=0.1, q10=0.1, hidden_states=None):
        """**CID-4**: four hidden regimes. ``rates`` is four ``(lambda, mu)`` pairs (or four
        bare lambdas, taking a shared extinction of 0)."""
        norm = [tuple(r) if isinstance(r, (tuple, list)) else (r, 0.0) for r in rates]
        if len(norm) != 4:
            raise ValueError("CID.four needs exactly four (lambda, mu) classes")
        return cls(norm, switch, q01=q01, q10=q10, hidden_states=hidden_states)

    def __repr__(self) -> str:
        return f"CID(classes={self._H})"


class _CharacterIndependentMuSSE(MuSSE):
    """A k-state MuSSE whose diversification depends only on a **hidden** class — the k-state
    generalisation of :class:`CID`, returned by :meth:`MuSSE.null(kind="cid") <MuSSE.null>`.

    Product states are ``(observed, hidden)``; :meth:`discretize` collapses them to the observed
    label so outputs never leak the hidden class (as :class:`HiSSE` does for the binary case).
    """

    def __init__(self, birth, death, Q, observed_states, n_hidden):
        H = int(n_hidden)
        obs = list(observed_states)
        states = [(obs[o], h) for o in range(len(obs)) for h in range(H)]
        super().__init__(birth=birth, death=death, Q=Q, states=states)
        self._H = H
        self._observed = obs

    def discretize(self, index):
        """The observed state of a product-state ``index`` (hidden class collapsed)."""
        return self._observed[int(index) // self._H]

    def __repr__(self) -> str:
        return f"MuSSE-CID(observed={len(self._observed)}, hidden={self._H})"


class QuaSSE:
    """Quantitative-trait Speciation and Extinction (FitzJohn 2010).

    A **continuous** trait diffuses (Brownian motion) along every lineage, and the speciation
    and extinction rates are functions of its current value — so a quantitative character shapes
    the tree. Simulate with :func:`simulate_sse`; the result is a *continuous*
    :class:`~zombi2.traits.TraitResult` whose ``.tree`` is the complete tree and whose ``.values``
    are the extant tips' trait values.

    The rate functions must be **bounded**: an unbounded ``λ(x)`` under a diffusing, unbounded
    ``x`` has no valid thinning bound. Pass ``rate_bound`` = an upper bound on ``λ(x) + μ(x)`` over
    all ``x``; :meth:`sigmoid` builds a convenient bounded rate function.

    Parameters
    ----------
    speciation, extinction:
        Callables ``x -> rate`` returning non-negative, bounded rates.
    sigma2:
        Diffusion rate of the trait (Brownian motion).
    rate_bound:
        An upper bound on ``speciation(x) + extinction(x)`` over all ``x`` (for exact thinning).
    x0:
        Root trait value (default ``0.0``); override per run with ``simulate_sse(..., root_state=)``.
    drift:
        Optional Brownian drift of the trait (default ``0.0``).
    """

    kind = "continuous"

    def __init__(self, speciation, extinction, sigma2, *, rate_bound, x0=0.0, drift=0.0):
        if not callable(speciation) or not callable(extinction):
            raise TypeError("speciation and extinction must be callables x -> rate")
        if sigma2 < 0:
            raise ValueError("sigma2 must be >= 0")
        if rate_bound <= 0:
            raise ValueError("rate_bound must be > 0")
        self.speciation = speciation
        self.extinction = extinction
        self.sigma2 = float(sigma2)
        self.rate_bound = float(rate_bound)
        self.x0 = float(x0)
        self.drift = float(drift)

    @staticmethod
    def sigmoid(low, high, center=0.0, slope=1.0):
        """A bounded rate function ``low + (high-low)/(1 + e^{-slope·(x-center)})`` in ``[low, high]``."""
        low, high, center, slope = float(low), float(high), float(center), float(slope)
        return lambda x: low + (high - low) / (1.0 + np.exp(-slope * (x - center)))

    def null(self, kind="neutral", *, rate=None, **kwargs):
        """Return a decoupled **null** for the ``traits:species`` arrow.

        ``"neutral"`` makes speciation **trait-independent** — a constant λ (so the diffusing trait
        no longer sets the rate), keeping the diffusion and extinction unchanged. The constant is
        the coupled speciation at the root value ``x0`` by default; override with ``rate=``. There
        is no ``"cid"`` (that is a discrete-character null) or ``"timing"`` null for QuaSSE.
        """
        kind = kind.lower()
        if kind == "neutral":
            lam0 = float(self.speciation(self.x0)) if rate is None else float(rate)
            return QuaSSE(lambda x: lam0, self.extinction, sigma2=self.sigma2,
                          rate_bound=self.rate_bound, x0=self.x0, drift=self.drift)
        if kind == "cid":
            raise TypeError("QuaSSE has no 'cid' null (CID is a discrete-character null); use "
                            "kind='neutral', or a BiSSE/MuSSE model for a cid null")
        if kind == "timing":
            raise ValueError("traits:species has no 'timing' null; use kind='neutral'")
        raise ValueError(f"unknown null kind {kind!r}; expected 'neutral'")

    def __repr__(self) -> str:
        return f"QuaSSE(sigma2={self.sigma2:g}, rate_bound={self.rate_bound:g})"


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


def _simulate_sse(model, age, n_tips, root_state, rng, max_lineages, clado=None):
    """One forward trial from a crown of two lineages sharing the root state.

    With ``clado`` (a :class:`~zombi2.traits.Cladogenesis` kernel) each daughter's state jumps at
    its birth (the ``species:traits`` arrow); combined with the state-dependent rates this is the
    full ClaSSE feedback.

    Returns ``(root, end_time, node_values, history)`` or ``None`` to reject (whole clade
    extinct / fewer than two extant survivors).
    """
    lambdas, mus, Q, k = model.lambdas, model.mus, model.Q, model.k
    out_rate = -np.diag(Q)                       # total transition rate leaving each state
    bound = float(np.max(lambdas + mus + out_rate))
    if bound <= 0.0:
        raise ValueError("all rates are zero; nothing can happen")

    def born(state):                             # a daughter's starting state (cladogenetic jump)
        return clado.apply(state, model, rng) if clado is not None else state

    root = TreeNode(name="", time=0.0)
    node_values = {root: root_state}
    history = {root: []}
    live = []
    for _ in range(2):                           # crown: two lineages at time 0, root state
        child = TreeNode(name="", time=0.0)
        root.add_child(child)
        live.append(_Lineage(child, born(root_state), 0.0))

    t = 0.0
    end = None
    while True:
        n = len(live)
        if n == 0:
            return None
        if n_tips is not None and n == n_tips:
            # present strictly after the N-th birth: last event + Exp(total event rate over the N
            # live lineages), the standard SSE-conditioned-on-N convention (avoids zero-length tips)
            R = sum(lambdas[s.state] + mus[s.state] + out_rate[s.state] for s in live)
            end = t + rng.exponential(1.0 / (R if R > 0.0 else n * bound))
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
        if r < lam:                              # speciation: daughters inherit i (jumped by clado)
            lin.node.time = t
            lin.close(t, node_values, history)
            live[idx] = live[-1]
            live.pop()
            for _ in range(2):
                d = TreeNode(name="", time=t)
                lin.node.add_child(d)
                live.append(_Lineage(d, born(i), t))
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


def _simulate_quasse(model, age, n_tips, x0, rng, max_lineages, clado=None):
    """One forward QuaSSE trial: a crown of two lineages carrying a diffusing continuous trait
    whose value sets each lineage's speciation/extinction rate (exact thinning against
    ``rate_bound``). With ``clado`` each daughter's value also jumps at its birth (the
    ``species:traits`` arrow — cladogenetic bursts). Returns ``(root, end, node_values)`` or
    ``None`` to reject."""
    spec, ext = model.speciation, model.extinction
    sigma2, drift, bound = model.sigma2, model.drift, model.rate_bound

    def born(x):                                 # a daughter's starting value (cladogenetic jump)
        return clado.apply(x, model, rng) if clado is not None else x

    root = TreeNode(name="", time=0.0)
    node_values = {root: x0}
    live = []  # [node, trait value]
    for _ in range(2):
        child = TreeNode(name="", time=0.0)
        root.add_child(child)
        live.append([child, born(x0)])

    def diffuse(interval):
        std = (sigma2 * interval) ** 0.5
        for L in live:
            L[1] += drift * interval + (rng.normal(0.0, std) if std > 0.0 else 0.0)

    t = 0.0
    end = None
    while True:
        n = len(live)
        if n == 0:
            return None
        if n_tips is not None and n == n_tips:
            # present strictly after the N-th birth: last event + Exp(total rate); see _simulate_sse
            R = sum(spec(L[1]) + ext(L[1]) for L in live)
            end = t + rng.exponential(1.0 / (R if R > 0.0 else n * bound))
            break
        if n > max_lineages:
            raise RuntimeError(
                f"QuaSSE tree exceeded max_lineages={max_lineages}; explosive parameters — "
                "lower the age/rates or raise max_lineages"
            )
        dt = rng.exponential(1.0 / (n * bound))
        if age is not None and t + dt >= age:
            diffuse(age - t)
            end = age
            break
        t += dt
        diffuse(dt)
        idx = int(rng.integers(n))
        node, x = live[idx]
        lam, mu = spec(x), ext(x)
        total = lam + mu
        if total > bound * (1.0 + 1e-9):
            # A rate_bound below the true speciation+extinction breaks rejection sampling:
            # total/bound > 1 accepts (almost) every candidate, silently biasing the tree.
            # Fail loudly rather than return a wrong tree.
            raise ValueError(
                f"rate_bound={bound:g} is below speciation+extinction={total:g} at trait "
                f"x={x:g}; it must be a true upper bound over all reachable trait values"
            )
        if total <= 0.0 or rng.random() >= total / bound:       # thinned out
            continue
        node.time = t
        node_values[node] = x
        live[idx] = live[-1]
        live.pop()
        if rng.random() * total < lam:                          # speciation: daughters inherit x (jumped)
            for _ in range(2):
                d = TreeNode(name="", time=t)
                node.add_child(d)
                live.append([d, born(x)])
        else:                                                   # extinction
            node.is_extant = False

    for node, x in live:                                        # survivors reach the present
        node.time = end
        node.is_extant = True
        node.sampled = True
        node_values[node] = x
    if len(live) < 2:
        return None
    return root, end, node_values


def simulate_sse(
    model,
    *,
    age: float | None = None,
    n_tips: int | None = None,
    root_state=None,
    cladogenesis: Cladogenesis | None = None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_attempts: int = 10_000,
    max_lineages: int = 1_000_000,
) -> TraitResult:
    """Simulate a tree jointly with a trait that drives its diversification (BiSSE / MuSSE).

    Provide exactly one stopping condition:

    * ``age`` — grow for this crown age; the number of extant tips is random.
    * ``n_tips`` — grow until this many extant lineages first coexist; the age is random.

    The run starts from a crown of two lineages sharing the root state and is conditioned on at
    least two extant survivors. Returns a :class:`~zombi2.traits.TraitResult` whose ``.tree`` is
    the **complete** simulated tree (extinct leaves carry ``is_extant=False``) and ``.values``
    are the extant tips' states.

    For a **discrete** model (:class:`BiSSE`/:class:`MuSSE`/:class:`HiSSE`) ``root_state`` is an
    integer state index (default: the character's stationary distribution) and ``.history`` is
    the realized character map. For :class:`QuaSSE` the trait is **continuous**: ``root_state``
    is the initial trait value (default the model's ``x0``) and ``.history`` is ``None``.

    ``cladogenesis`` (a :class:`~zombi2.traits.Cladogenesis` kernel) additionally jumps each
    daughter's state **at speciation** — the ``species:traits`` arrow. Combining it with the
    state-dependent rates gives the full **ClaSSE** feedback (trait ↔ tree): the trait both shapes
    the tree *and* is kicked by its branching.
    """
    if (age is None) == (n_tips is None):
        raise ValueError("provide exactly one of `age` or `n_tips`")
    if age is not None and age <= 0:
        raise ValueError(f"age must be > 0, got {age}")
    if n_tips is not None and n_tips < 2:
        raise ValueError(f"n_tips must be >= 2, got {n_tips}")
    if rng is None:
        rng = np.random.default_rng(seed)
    clado = cladogenesis if (cladogenesis is not None and cladogenesis.is_active()) else None

    if isinstance(model, QuaSSE):
        x0 = model.x0 if root_state is None else float(root_state)
        for _ in range(max_attempts):
            result = _simulate_quasse(model, age, n_tips, x0, rng, max_lineages, clado)
            if result is not None:
                root, end_time, node_values = result
                tree = Tree(root, end_time)
                _name(tree)
                return TraitResult(tree=tree, model=model, node_values=node_values,
                                   history=None, kind="continuous")
        raise RuntimeError(
            f"QuaSSE simulation produced no surviving tree in {max_attempts} attempts; "
            "raise max_attempts or lower the extinction rate"
        )

    if not isinstance(model, MuSSE):
        raise TypeError("model must be a BiSSE, MuSSE, HiSSE, or QuaSSE instance")
    if root_state is not None and not (0 <= int(root_state) < model.k):
        raise ValueError(f"root_state must be a state index in [0, {model.k - 1}]")
    # With no anagenetic transitions (Q = 0) the stationary root state is undefined — but if
    # cladogenesis is active it supplies the state dynamics, so a uniform root is fine there.
    q_is_zero = np.allclose(model.Q - np.diag(np.diag(model.Q)), 0.0)
    if root_state is None and q_is_zero and clado is None:
        raise ValueError("with no state transitions (Q = 0) the stationary root state is "
                         "undefined; pass an explicit root_state (or add cladogenesis)")

    for _ in range(max_attempts):
        if root_state is not None:
            r0 = int(root_state)
        elif q_is_zero:                          # clado active (else we raised): uniform root
            r0 = int(rng.integers(model.k))
        else:
            r0 = int(rng.choice(model.k, p=model.stationary_distribution()))
        result = _simulate_sse(model, age, n_tips, r0, rng, max_lineages, clado)
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
