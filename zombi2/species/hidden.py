"""Hidden lineages for a fixed sampled tree (issue #480).

A real dated tree shows the lineages that were sampled. The lineages that went extinct, and the
survivors nobody sampled, are **hidden**: they belong to the complete tree and leave no sampled
descendant. `add_hidden_lineages` simulates them under a birth–death process with sampling,
conditioned on the given tree, so that pruning the result to the sampled tips gives the input back
exactly. A benchmark of transfers from hidden lineages on a real tree is what it is for.

The model, with age ``t`` measured before the present, constant rates λ (birth) and μ (death), and
each present-day lineage sampled with probability ρ. ``E(t)`` is the probability that a lineage alive
at age ``t`` leaves no sampled descendant::

    A(t) = e^(−rt) + ρλ(1 − e^(−rt)) / r,   r = λ − μ        (A(t) = 1 + ρλt when λ = μ)
    E(t) = 1 − ρ / A(t)

Along every sampled branch, hidden lineages attach as a Poisson process of rate ``2λE(t)``: a
speciation whose other daughter leaves no sampled descendant, given that the branch itself does.
Each hidden lineage then evolves conditioned on leaving none, which is the Doob transform with
``h = E``: birth at ``λE(t)``, death at ``μ/E(t)``. Its next event is drawn exactly by inverting the
cumulative hazard ``H(a, u) = (λ + μ)(a − u) + log E(a) − log E(u)`` between the ages ``a`` and
``u < a``. A hidden lineage that reaches the present is unsampled, with no further draw.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..rng import stream
from ..tree import Node, Tree, _depths, _preorder


@dataclass(frozen=True)
class _Model:
    """λ, μ, ρ, and the functions of age the sampler reads."""

    birth: float
    death: float
    sampling: float

    def A(self, t: float) -> float:
        r = self.birth - self.death
        if abs(r) < 1e-12:
            return 1.0 + self.sampling * self.birth * t
        e = math.exp(-r * t)
        return e + self.sampling * self.birth * (1.0 - e) / r

    def E(self, t: float) -> float:
        """The probability that a lineage alive at age ``t`` leaves no sampled descendant."""
        return max(0.0, 1.0 - self.sampling / self.A(t))

    def hazard(self, a: float, u: float) -> float:
        """The cumulative hazard of a hidden lineage's next event between ages ``a`` and ``u < a``;
        infinite when ``E(u)`` is 0, which is the present under full sampling."""
        eu = self.E(u)
        if eu <= 0.0:
            return math.inf
        return (self.birth + self.death) * (a - u) + math.log(self.E(a)) - math.log(eu)

    def next_event_age(self, a: float, rng) -> float | None:
        """The age of a hidden lineage's next event, below ``a``, or ``None`` when it reaches the
        present: one exponential draw, inverted through `hazard` by bisection."""
        x = float(rng.exponential(1.0))
        if self.hazard(a, 0.0) <= x:
            return None
        lo, hi = 0.0, a                         # hazard(a, hi) = 0 < x <= hazard(a, lo)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self.hazard(a, mid) < x:
                hi = mid
            else:
                lo = mid
            if hi - lo <= 1e-15 * max(a, 1.0):
                break
        return 0.5 * (lo + hi)

    def birth_share(self, u: float) -> float:
        """The probability that a hidden lineage's event at age ``u`` is a birth rather than a death."""
        e = self.E(u)
        b, d = self.birth * e, (self.death / e if e > 0.0 else math.inf)
        return 0.0 if math.isinf(d) else b / (b + d)


def _check(birth, death, sampling) -> _Model:
    for name, value in (("birth", birth), ("death", death), ("sampling", sampling)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{name} must be a finite number, got {value!r}")
    if birth <= 0:
        raise ValueError(f"birth must be positive, got {birth!r}")
    if death < 0:
        raise ValueError(f"death must be non-negative, got {death!r}")
    if not 0.0 < sampling <= 1.0:
        raise ValueError(f"sampling is the probability that a surviving lineage was sampled, in "
                         f"(0, 1], got {sampling!r}")
    return _Model(float(birth), float(death), float(sampling))


def _check_sampled(sampled: Tree) -> float:
    """The input is a sampled tree: every tip extant, every tip at the present. Returns the crown
    age, the depth of the tips below the root's end."""
    tips = sampled.leaves()
    if not tips:
        raise ValueError("the sampled tree has no tips")
    bad = sorted(sampled.nodes[i].fate for i in tips if sampled.nodes[i].fate != "extant")
    if bad:
        raise ValueError(
            f"every tip of the sampled tree must be extant, since the hidden lineages are what the "
            f"tree does not show; found tips with fate {sorted(set(bad))}. A tree read from a Newick "
            f"with every tip at the present has that already.")
    depth = _depths(sampled)                                 # from the crown, the root's end
    tip_depths = [depth[i] for i in tips]
    lo, hi = min(tip_depths), max(tip_depths)
    if hi > 0 and (hi - lo) > 1e-6 * hi:
        raise ValueError(
            f"the sampled tree is not ultrametric (tip depths differ by {hi - lo:.3g}); every sampled "
            f"tip is at the present. zombi2.tree.make_ultrametric snaps rounding noise.")
    return hi


def add_hidden_lineages(sampled, *, birth, death=0.0, sampling=1.0, stem=None, seed=None,
                        names=None, progress=False):
    """Add the lineages a sampled tree does not show, under a birth–death process with sampling,
    conditioned on the tree (#480). Returns a `SpeciesResult` whose complete tree holds the input's
    lineages and the hidden ones, and whose ``extant_tree`` is the input again.

    ``sampled`` is a dated tree whose tips are all at the present, as `zombi2.tree.read_newick`
    gives for a published tree; ``names`` is the name map it returned, written to ``names.tsv``.
    ``birth``, ``death`` and ``sampling`` are the rates and the sampling probability the hidden
    lineages are simulated under; `fit_birth_death_rates` estimates the rates from the tree.
    ``stem`` sets the length of the root branch above the crown, where hidden lineages attach too;
    without it the input's own root branch is kept, which for a tree read from a Newick without a
    root length is zero. Time in the result runs forward from the origin, the top of the stem.

    Sampled tips keep the fate ``extant``. A hidden lineage alive at the present is ``unsampled``;
    one that died is ``extinct``. Every attachment is a speciation node on the sampled branch, and
    which side the hidden daughter takes is drawn, so a drawing shows no systematic side. Pruning the
    result to its extant tips recovers the input exactly, node order included.

    The function simulates hidden lineages for a benchmark; it does not estimate the true ones. Under
    ``death=0`` and ``sampling=1`` there are none to add, and the result is the input with the fates
    on. Transfer rates on the completed tree need recalibration, since recipients now include the
    hidden lineages.
    """
    from . import Event, SpeciesResult

    model = _check(birth, death, sampling)
    crown_age = _check_sampled(sampled)
    rng, seed = stream("species", seed)
    root_in = sampled.nodes[sampled.root]
    stem_len = root_in.end_time - root_in.birth_time if stem is None else float(stem)
    if not math.isfinite(stem_len) or stem_len < 0:
        raise ValueError(f"stem must be a non-negative length, got {stem!r}")
    present = crown_age + stem_len                       # the output's time axis: origin at 0
    depth = _depths(sampled)                             # from the crown, the root's end

    def age_of_end(i: int) -> float:
        """The age at which sampled node ``i``'s branch ends."""
        return present - (depth[i] + stem_len)

    nodes: dict[int, Node] = {}
    events: list[Event] = []
    counter = 0

    def new_node(parent: int | None, birth_age: float, end_age: float, fate: str) -> int:
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, present - birth_age, present - end_age, (), fate)
        return i

    def grow_hidden(parent: int, age: float) -> int:
        """A hidden lineage born at ``age`` under ``parent``: its whole subtree, conditioned on
        leaving no sampled descendant. Returns the id of its first node."""
        first = new_node(parent, age, age, "alive")
        stack = [(first, age)]
        while stack:
            i, a = stack.pop()
            u = model.next_event_age(a, rng)
            nd = nodes[i]
            if u is None:                                    # reaches the present, unsampled
                nodes[i] = Node(i, nd.parent, nd.birth_time, present, (), "unsampled")
                continue
            if float(rng.random()) < model.birth_share(u):   # a split into two hidden lineages
                c1 = new_node(i, u, u, "alive")
                c2 = new_node(i, u, u, "alive")
                nodes[i] = Node(i, nd.parent, nd.birth_time, present - u, (c1, c2), "speciation")
                events.append(Event(present - u, "speciation", i, (c1, c2)))
                stack.append((c2, u))
                stack.append((c1, u))
            else:                                            # dies
                nodes[i] = Node(i, nd.parent, nd.birth_time, present - u, (), "extinct")
                events.append(Event(present - u, "extinction", i))
        return first

    def attachments(a_old: float, a_young: float) -> list[float]:
        """The ages at which hidden lineages attach along a branch from age ``a_old`` down to
        ``a_young``: a Poisson process of rate ``2λE(t)``, thinned from rate ``2λ``."""
        out: list[float] = []
        bound = 2.0 * model.birth
        t = a_old
        while True:
            t -= float(rng.exponential(1.0 / bound))
            if t <= a_young:
                return sorted(out, reverse=True)
            if float(rng.random()) < model.E(t):
                out.append(t)

    # walk the sampled tree, parent before child; each sampled branch becomes a chain of
    # attachment nodes ending in the sampled node itself, which keeps the input's children order
    def link(node: int, hidden_first: int, hidden_left: bool, continuation: int) -> None:
        """An attachment node's two children, in the drawn order, and its speciation event."""
        nd = nodes[node]
        kids = (hidden_first, continuation) if hidden_left else (continuation, hidden_first)
        nodes[node] = Node(node, nd.parent, nd.birth_time, nd.end_time, kids, "speciation")
        events.append(Event(nd.end_time, "speciation", node, kids))

    new_id: dict[int, int] = {}            # sampled node -> its own node, the last of its branch's chain
    first_id: dict[int, int] = {}          # sampled node -> the first node of that chain
    for s in _preorder(sampled):
        sn = sampled.nodes[s]
        parent_new = None if sn.parent is None else new_id[sn.parent]
        a_start = present if sn.parent is None else age_of_end(sn.parent)
        a_end = age_of_end(s)
        pending = None                     # (attachment node, its hidden daughter, hidden on the left)
        for t in attachments(a_start, a_end):
            node = new_node(parent_new, a_start, t, "speciation")
            first_id.setdefault(s, node)
            if pending is not None:
                link(*pending, continuation=node)
            hidden_first = grow_hidden(node, t)
            pending = (node, hidden_first, float(rng.random()) < 0.5)
            parent_new, a_start = node, t
        me = new_node(parent_new, a_start, a_end, "extant" if not sn.children else "speciation")
        first_id.setdefault(s, me)
        if pending is not None:
            link(*pending, continuation=me)
        new_id[s] = me
    for s in _preorder(sampled):
        sn = sampled.nodes[s]
        if sn.children:
            me = new_id[s]
            kids = tuple(first_id[c] for c in sn.children)   # a child branch starts at its chain's first node
            nd = nodes[me]
            nodes[me] = Node(me, nd.parent, nd.birth_time, nd.end_time, kids, "speciation")
            events.append(Event(nd.end_time, "speciation", me, kids))
    events.sort(key=lambda e: (e.time, e.kind != "speciation", e.node))
    # the root is the topmost node: an attachment on the stem sits above the sampled root
    tree = Tree(nodes, next(i for i, n in nodes.items() if n.parent is None))
    labels = {new_id[i]: name for i, name in (names or {}).items() if i in new_id}
    return SpeciesResult(tree, events, seed, names=labels)


# --- the rates, from the tree ---------------------------------------------------------------------

def _log_likelihood(ages: list[float], crown: float, model: _Model) -> float:
    """The log-likelihood of the sampled tree's internal node ages given its crown age, under the
    reconstructed birth–death process with sampling (Nee et al. 1994; Stadler 2009): each age below
    the crown is an independent draw with density ``λ p1(x) / Z``, where
    ``p1(t) = ρ r² e^(−rt) / (ρλ + c e^(−rt))²`` with ``c = λ(1−ρ) − μ``, and ``Z`` is that density's
    integral from the present to the crown age, so the ages are a truncated distribution on
    ``[0, crown]``. In the pure-birth case with full sampling this is the truncated exponential."""
    lam, mu, rho = model.birth, model.death, model.sampling
    r = lam - mu
    c = lam * (1.0 - rho) - mu
    if abs(r) < 1e-12:
        return -math.inf

    def log_p1(t: float) -> float:
        denom = rho * lam + c * math.exp(-r * t)
        if denom <= 0:
            return -math.inf
        return math.log(rho) + 2.0 * math.log(abs(r)) - r * t - 2.0 * math.log(denom)

    # Z = ∫_0^crown λ p1(x) dx. With g(x) = ρ r / (ρλ + c e^(−rx)), dg/dx = c p1(x), so
    # Z = (λ / c)(g(crown) − g(0)) when c ≠ 0, and λ(1 − e^(−λ crown)) / λ... in the pure-birth
    # limit c = 0, where p1 = e^(−rt): Z = (λ / r)(1 − e^(−r crown)).
    if abs(c) < 1e-12:
        z = (lam / r) * (1.0 - math.exp(-r * crown))
    else:
        g_crown = rho * r / (rho * lam + c * math.exp(-r * crown))
        z = (lam / c) * (g_crown - rho)
    if not z > 0 or not math.isfinite(z):
        return -math.inf
    return sum(math.log(lam) + log_p1(x) - math.log(z) for x in ages)


def _branching_ages(sampled: Tree) -> tuple[list[float], float]:
    depth = _depths(sampled)
    tips = sampled.leaves()
    crown = max(depth[i] for i in tips)                   # the crown age: depth is from the crown
    internal = [i for i, n in sampled.nodes.items() if n.children and i != sampled.root]
    ages = [crown - depth[i] for i in internal]
    return ages, crown


def fit_birth_death_rates(sampled, *, sampling=1.0, death=None) -> tuple[float, float]:
    """Birth and death rates fitted to a sampled tree by maximum likelihood, given the sampling
    probability (#480): the reconstructed-tree likelihood conditional on the crown age. With
    ``death`` given, the birth rate alone is fitted; with ``death=None`` both are, by a
    Nelder–Mead search on their logarithms. With all three free only ``ρλ`` and ``λ − μ`` are
    identifiable, which is why ``sampling`` is an input. Returns ``(birth, death)``."""
    _check(1.0, 0.0 if death is None else death, sampling)
    _check_sampled(sampled)
    ages, crown = _branching_ages(sampled)
    if len(ages) < 2:
        raise ValueError("fitting the rates needs a tree with at least four tips")

    def objective(params):
        lam = math.exp(params[0])
        mu = death if death is not None else math.exp(params[1])
        if lam <= 0 or lam > 1e4 or (mu is not None and mu > lam * 50):
            return math.inf
        return -_log_likelihood(ages, crown, _Model(lam, float(mu), float(sampling)))

    n = 1 if death is not None else 2
    start = [math.log(max(len(ages), 2) / crown)] + ([math.log(0.1 * len(ages) / crown)] if n == 2 else [])
    best = _nelder_mead(objective, start)
    lam = math.exp(best[0])
    mu = float(death) if death is not None else math.exp(best[1])
    return lam, mu


def _nelder_mead(f, x0, step: float = 0.5, iterations: int = 400, tol: float = 1e-10) -> list[float]:
    """A small Nelder–Mead for one or two parameters: numpy only, no scipy."""
    n = len(x0)
    simplex = [list(x0)] + [[x0[j] + (step if j == i else 0.0) for j in range(n)] for i in range(n)]
    values = [f(p) for p in simplex]
    for _ in range(iterations):
        order = sorted(range(n + 1), key=lambda k: values[k])
        simplex = [simplex[k] for k in order]
        values = [values[k] for k in order]
        if abs(values[-1] - values[0]) < tol and max(
                abs(simplex[-1][j] - simplex[0][j]) for j in range(n)) < 1e-8:
            break
        centroid = [sum(p[j] for p in simplex[:-1]) / n for j in range(n)]
        worst = simplex[-1]
        reflect = [centroid[j] + (centroid[j] - worst[j]) for j in range(n)]
        fr = f(reflect)
        if fr < values[0]:
            expand = [centroid[j] + 2.0 * (centroid[j] - worst[j]) for j in range(n)]
            fe = f(expand)
            simplex[-1], values[-1] = (expand, fe) if fe < fr else (reflect, fr)
        elif fr < values[-2]:
            simplex[-1], values[-1] = reflect, fr
        else:
            contract = [centroid[j] + 0.5 * (worst[j] - centroid[j]) for j in range(n)]
            fc = f(contract)
            if fc < values[-1]:
                simplex[-1], values[-1] = contract, fc
            else:
                best = simplex[0]
                simplex = [best] + [[best[j] + 0.5 * (p[j] - best[j]) for j in range(n)]
                                    for p in simplex[1:]]
                values = [values[0]] + [f(p) for p in simplex[1:]]
    return simplex[int(min(range(n + 1), key=lambda k: values[k]))]


__all__ = ["add_hidden_lineages", "fit_birth_death_rates"]
