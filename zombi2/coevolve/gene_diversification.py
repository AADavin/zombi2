"""Gene-content-dependent diversification — the ``genes:species`` edge (key innovations + HGT).

Here a lineage's **gene content** sets its speciation and extinction rates, so the tree's shape
depends on the genes it carries — and, because a gene is gained horizontally at a rate that grows
with how prevalent it already is (transfer is frequency-dependent), the tree and the gene content
must grow **together** in one forward Gillespie. This is the ``genes:species`` arrow of the
:doc:`coevolution model <models/coevolution>`, and the one edge that genuinely merges the
species-tree and gene processes rather than layering them.

To keep it tractable (and exact under ZOMBI's independent-family assumption) only a small panel of
binary **driver** ("key innovation") families rides in the forward loop — those are the ones that
touch diversification. The rest of the genome does not affect the tree, so it is simulated
**afterward** on the finished tree with the ordinary ``zombi2 genomes`` / :func:`simulate_genomes`.

Each live lineage carries a set of present drivers ``S`` and competes to:

* **speciate** at ``λ(S) = λ₀ · exp(Σ_{d∈S} βλ_d)`` — both daughters inherit ``S``;
* **go extinct** at ``μ(S) = μ₀ · exp(Σ_{d∈S} βμ_d)`` — the lineage dies;
* **lose** a present driver ``d`` at rate ``loss`` — deletion;
* **originate** an absent driver ``d`` at rate ``origination`` — de-novo gain;
* **receive** an absent driver ``d`` by **transfer** at rate ``t·carriers[d]/(n−1)`` — a driver in
  more of the ``n`` live genomes is donated more often (rich-get-richer, and self-saturating).

    import zombi2 as z
    m = z.GeneDiversification(2, driver_speciation=1.5, transfer=0.6, root_drivers=1)
    res = z.simulate_gene_diversification(m, age=4.0, seed=1)
    res.tree                 # complete tree (its shape was driven by gene content)
    res.tip_prevalence()     # fraction of extant tips carrying each driver
"""

from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np

from zombi2.tree import Tree, TreeNode
from zombi2.species.forward import _name


# --------------------------------------------------------------------------- model
def _as_vector(x, k, name):
    """A scalar (broadcast to length ``k``) or a length-``k`` sequence, as a float array."""
    arr = np.atleast_1d(np.asarray(x, dtype=float))
    if arr.size == 1:
        arr = np.full(k, float(arr[0]))
    if arr.shape != (k,):
        raise ValueError(f"{name} must be a scalar or a length-{k} vector, got shape {arr.shape}")
    return arr


class GeneDiversification:
    """Gene-content-dependent diversification: ``K`` binary driver families shape the tree.

    Parameters
    ----------
    n_drivers:
        Number of binary driver ("key innovation") families, ``K ≥ 1``.
    lambda0, mu0:
        Base speciation / extinction rates of a lineage carrying **no** drivers (``≥ 0``).
    driver_speciation, driver_extinction:
        Per-driver effect on the **log** speciation / extinction rate — a scalar (shared by all
        ``K``) or a length-``K`` vector. A present driver ``d`` multiplies ``λ`` by ``exp(βλ_d)``
        and ``μ`` by ``exp(βμ_d)``; ``driver_speciation > 0`` is a key innovation (a radiation),
        ``driver_extinction`` lets a driver instead (or also) buffer/raise extinction.
    loss, origination, transfer:
        Per-driver rates: ``loss`` deletes a present driver; ``origination`` gains an absent one de
        novo; ``transfer`` is the per-donor HGT rate (a driver in more live genomes spreads faster).
    root_drivers:
        Drivers present at the root — a count ``m`` (the first ``m`` drivers) or an explicit list of
        indices (default ``0`` = none; drivers then enter only by origination).
    """

    def __init__(self, n_drivers, *, lambda0=1.0, mu0=0.2,
                 driver_speciation=1.0, driver_extinction=0.0,
                 loss=0.1, origination=0.05, transfer=0.5, root_drivers=0,
                 cladogenetic_loss=0.0, cladogenetic_gain=0.0):
        k = int(n_drivers)
        if k < 1:
            raise ValueError(f"n_drivers must be >= 1, got {n_drivers}")
        if lambda0 < 0 or mu0 < 0:
            raise ValueError("lambda0 and mu0 must be >= 0")
        if loss < 0 or origination < 0 or transfer < 0:
            raise ValueError("loss, origination and transfer must be >= 0")
        if not (0.0 <= cladogenetic_loss <= 1.0) or not (0.0 <= cladogenetic_gain <= 1.0):
            raise ValueError("cladogenetic_loss and cladogenetic_gain must be probabilities in [0, 1]")
        self.n_drivers = k
        self.lambda0 = float(lambda0)
        self.mu0 = float(mu0)
        self.beta_lambda = _as_vector(driver_speciation, k, "driver_speciation")
        self.beta_mu = _as_vector(driver_extinction, k, "driver_extinction")
        self.loss = float(loss)
        self.origination = float(origination)
        self.transfer = float(transfer)
        # cladogenetic burst on the drivers at each speciation (the species:genes arrow). With both
        # > 0 this is the species<->genes JOINT model ("co-diversification"): the same drivers set
        # the rates AND are reshuffled at every split, so a burst can hand one daughter a key
        # innovation and not its sister. Default 0 = pure genes:species (drivers change only along
        # branches).
        self.cladogenetic_loss = float(cladogenetic_loss)
        self.cladogenetic_gain = float(cladogenetic_gain)
        if np.isscalar(root_drivers):
            m = int(root_drivers)
            if not (0 <= m <= k):
                raise ValueError(f"root_drivers count must be in [0, {k}], got {m}")
            self.root_set = frozenset(range(m))
        else:
            s = frozenset(int(i) for i in root_drivers)
            if any(i < 0 or i >= k for i in s):
                raise ValueError(f"root_drivers indices must be in [0, {k - 1}]")
            self.root_set = s

    def rates(self, present):
        """``(λ, μ)`` for a lineage carrying the driver-index set ``present``."""
        if present:
            idx = list(present)
            lam = self.lambda0 * float(np.exp(self.beta_lambda[idx].sum()))
            mu = self.mu0 * float(np.exp(self.beta_mu[idx].sum()))
        else:
            lam, mu = self.lambda0, self.mu0
        return lam, mu

    def null(self, kind="neutral", **kwargs):
        """Decoupled **null** for the ``genes:species`` arrow (gene content → diversification).
        See :doc:`the null-models guide </guide/coevolution_nulls>`.

        ``"neutral"`` zeroes every driver's effect on the rates (β = 0): the drivers still spread
        and vary, but no longer set λ/μ, so the tree is constant-rate. The character-independent
        (``"cid"``) null for this edge is a **workflow**, not a model transform — grow the
        driver-shaped tree, then analyse a *neutral overlay genome* while withholding the drivers;
        the CLI ``--null cid`` builds it for you.
        """
        kind = kind.lower()
        if kind == "neutral":
            m = copy.copy(self)
            m.beta_lambda = np.zeros_like(self.beta_lambda)
            m.beta_mu = np.zeros_like(self.beta_mu)
            return m
        if kind == "cid":
            raise TypeError(
                "the genes:species CID null is a workflow (grow the driver-shaped tree, then "
                "observe a neutral overlay genome and withhold the drivers), not a model "
                "transform; use the CLI `--null cid` — see docs/guide/coevolution_nulls.md")
        if kind == "timing":
            raise ValueError("genes:species has no 'timing' null (its driver is a gene state, not "
                             "an at-speciation event); use kind='neutral' (or 'cid' via the CLI)")
        raise ValueError(f"unknown null kind {kind!r}; expected 'neutral'")

    @property
    def is_co_diversification(self) -> bool:
        """True when a cladogenetic burst is active — the species<->genes joint model."""
        return self.cladogenetic_loss > 0.0 or self.cladogenetic_gain > 0.0

    def __repr__(self):
        return f"GeneDiversification(n_drivers={self.n_drivers}, transfer={self.transfer:g})"


@dataclass
class GeneDiversificationResult:
    """The outcome of :func:`simulate_gene_diversification`.

    ``node_drivers`` maps every node to the (frozen) set of driver indices present there — the tips
    are the observable data; internal nodes are the true ancestral gene content. ``.tree`` is the
    **complete** tree whose shape the drivers produced (extinct leaves carry ``is_extant=False``).
    """

    tree: Tree
    model: GeneDiversification
    node_drivers: dict

    def driver_names(self):
        return [f"D{i}" for i in range(self.model.n_drivers)]

    def tip_prevalence(self):
        """Fraction of the extant tips carrying each driver (length ``K``)."""
        tips = self.tree.extant_leaves()
        k = self.model.n_drivers
        if not tips:
            return [0.0] * k
        counts = [0] * k
        for tip in tips:
            for d in self.node_drivers[tip]:
                counts[d] += 1
        return [c / len(tips) for c in counts]

    def to_tsv(self, nodes: str = "all") -> str:
        """A ``node`` × driver 0/1 presence table (one column per driver)."""
        if nodes == "extant":
            selected = self.tree.extant_leaves()
        elif nodes == "leaves":
            selected = self.tree.leaves()
        elif nodes == "all":
            selected = self.tree.nodes()
        else:
            raise ValueError("nodes must be 'extant', 'leaves', or 'all'")
        names = self.driver_names()
        lines = ["node\t" + "\t".join(names)]
        for node in selected:
            present = self.node_drivers[node]
            lines.append(node.name + "\t" + "\t".join(
                "1" if d in present else "0" for d in range(self.model.n_drivers)))
        return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- engine
class _Lineage:
    """A live lineage: its growing node and the set of driver indices it currently carries."""

    __slots__ = ("node", "drivers")

    def __init__(self, node, drivers):
        self.node = node
        self.drivers = set(drivers)


def _burst_drivers(present, k, clado_loss, clado_gain, rng):
    """A daughter's driver set after a cladogenetic burst at speciation (the ``species:genes``
    arrow over a fixed binary panel): drop each *present* driver with probability ``clado_loss``,
    and gain each *absent* driver with probability ``clado_gain``."""
    out = set()
    for d in range(k):
        if d in present:
            if clado_loss <= 0.0 or rng.random() >= clado_loss:
                out.add(d)                            # retained through the burst
        elif clado_gain > 0.0 and rng.random() < clado_gain:
            out.add(d)                                # gained at the burst
    return out


def _simulate_once(model, age, n_tips, rng, max_lineages):
    """One forward trial from a crown of two lineages sharing the root driver set.

    Exact-Gillespie (rates recomputed each event, since transfer depends on the live population).
    Returns ``(root, end, node_drivers)`` or ``None`` to reject (extinct / fewer than 2 survivors).
    """
    k = model.n_drivers
    lam0, mu0 = model.lambda0, model.mu0
    beta_l, beta_m = model.beta_lambda, model.beta_mu
    loss, orig, tr = model.loss, model.origination, model.transfer
    clado_loss, clado_gain = model.cladogenetic_loss, model.cladogenetic_gain
    bursts = clado_loss > 0.0 or clado_gain > 0.0            # species:genes arrow active = co-div
    root_set = model.root_set

    carriers = [0] * k                            # how many live lineages carry each driver
    node_drivers = {}

    root = TreeNode(name="", time=0.0)
    node_drivers[root] = frozenset(root_set)
    live: list[_Lineage] = []

    def add(node, drivers):
        lin = _Lineage(node, drivers)
        live.append(lin)
        for d in lin.drivers:
            carriers[d] += 1
        return lin

    def drop(lin):                                # lineage leaves the live set
        for d in lin.drivers:
            carriers[d] -= 1

    for _ in range(2):                            # crown: two lineages carrying the root drivers
        child = TreeNode(name="", time=0.0)
        root.add_child(child)
        add(child, root_set)

    t = 0.0
    end = None
    while True:
        n = len(live)
        if n == 0:
            return None
        if n > max_lineages:
            raise RuntimeError(
                f"gene-diversification tree exceeded max_lineages={max_lineages}; a driver has "
                "likely fixed and run away — lower driver_speciation/transfer or raise max_lineages")

        total_copies = sum(carriers)
        t_unit = tr / (n - 1) if (n > 1 and tr > 0.0) else 0.0      # per-copy transfer contribution
        # per-lineage total rates (present drivers iterated directly; absent handled via globals)
        rates = []
        rtot = 0.0
        for lin in live:
            S = lin.drivers
            s = len(S)
            if s:
                lam = lam0 * float(np.exp(beta_l[list(S)].sum()))
                mu = mu0 * float(np.exp(beta_m[list(S)].sum()))
                present_carriers = sum(carriers[d] for d in S)
            else:
                lam, mu, present_carriers = lam0, mu0, 0
            loss_i = loss * s
            gain_i = orig * (k - s) + t_unit * (total_copies - present_carriers)   # absent drivers
            ri = lam + mu + loss_i + gain_i
            rates.append((ri, lam, mu, loss_i))
            rtot += ri

        if rtot <= 0.0:                            # nothing can happen
            if age is not None:
                end = age
                break
            return None                            # n_tips mode can never be reached

        if n_tips is not None and n == n_tips:
            # present strictly after the N-th birth: last event + Exp(total rate). See forward._grow.
            end = t + rng.exponential(1.0 / rtot)
            break

        dt = rng.exponential(1.0 / rtot)
        if age is not None and t + dt >= age:
            end = age
            break
        t += dt

        # pick the lineage (weighted by its total rate), then the event within it
        u = rng.random() * rtot
        idx = 0
        acc = 0.0
        for j, (ri, *_rest) in enumerate(rates):
            acc += ri
            if u < acc:
                idx = j
                break
        lin = live[idx]
        ri, lam, mu, loss_i = rates[idx]
        S = lin.drivers
        e = rng.random() * ri
        if e < lam:                                # speciation
            lin.node.time = t
            node_drivers[lin.node] = frozenset(S)
            drop(lin)
            live[idx] = live[-1]
            live.pop()
            for _ in range(2):                     # each daughter: burst the drivers, then live on
                child = TreeNode(name="", time=t)
                lin.node.add_child(child)
                child_S = _burst_drivers(S, k, clado_loss, clado_gain, rng) if bursts else set(S)
                add(child, child_S)
        elif e < lam + mu:                         # extinction
            lin.node.time = t
            lin.node.is_extant = False
            node_drivers[lin.node] = frozenset(S)
            drop(lin)
            live[idx] = live[-1]
            live.pop()
        elif e < lam + mu + loss_i:                # lose a present driver (uniform among present)
            d = int(rng.choice(list(S)))
            S.discard(d)
            carriers[d] -= 1
        else:                                      # gain an absent driver (origination or transfer)
            absent = [d for d in range(k) if d not in S]
            weights = [orig + t_unit * carriers[d] for d in absent]
            wsum = sum(weights)
            d = absent[int(rng.choice(len(absent), p=[w / wsum for w in weights]))] if wsum > 0 \
                else absent[int(rng.integers(len(absent)))]
            S.add(d)
            carriers[d] += 1

    for lin in live:                               # survivors reach the present
        lin.node.time = end
        lin.node.is_extant = True
        lin.node.sampled = True
        node_drivers[lin.node] = frozenset(lin.drivers)
    if len(live) < 2:
        return None
    return root, end, node_drivers


def simulate_gene_diversification(
    model: GeneDiversification,
    *,
    age: float | None = None,
    n_tips: int | None = None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_attempts: int = 10_000,
    max_lineages: int = 1_000_000,
) -> GeneDiversificationResult:
    """Grow a tree jointly with a set of driver gene families that shape its diversification.

    Provide exactly one stopping condition — ``age`` (grow for this crown age; the tip count is
    random) or ``n_tips`` (grow until this many extant lineages coexist; the age is random). The
    run starts from a crown of two lineages carrying the model's ``root_drivers`` and is conditioned
    on at least two extant survivors.

    Returns a :class:`GeneDiversificationResult` whose ``.tree`` is the **complete** tree (extinct
    lineages kept) and ``.node_drivers`` the exact per-node driver presence. Overlay the neutral
    genome afterward with :func:`~zombi2.simulate_genomes` (or ``zombi2 genomes``) on ``.tree``.
    """
    if (age is None) == (n_tips is None):
        raise ValueError("provide exactly one of `age` or `n_tips`")
    if age is not None and age <= 0:
        raise ValueError(f"age must be > 0, got {age}")
    if n_tips is not None and n_tips < 2:
        raise ValueError(f"n_tips must be >= 2, got {n_tips}")
    if rng is None:
        rng = np.random.default_rng(seed)

    for _ in range(max_attempts):
        result = _simulate_once(model, age, n_tips, rng, max_lineages)
        if result is not None:
            root, end, node_drivers = result
            tree = Tree(root, end)
            _name(tree)
            return GeneDiversificationResult(tree=tree, model=model, node_drivers=node_drivers)

    raise RuntimeError(
        f"gene-diversification produced no surviving tree in {max_attempts} attempts "
        "(the clade kept going extinct); raise max_attempts or lower the extinction rates")


def simulate_co_diversification(
    model: GeneDiversification,
    *,
    age: float | None = None,
    n_tips: int | None = None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    max_attempts: int = 10_000,
    max_lineages: int = 1_000_000,
) -> GeneDiversificationResult:
    """The **species<->genes joint model** — ``genes:species`` *and* ``species:genes`` at once.

    The same panel of driver families both **sets the diversification rates** (``genes:species``:
    a lineage carrying key innovations speciates faster) *and* is **reshuffled by a cladogenetic
    burst at every speciation** (``species:genes``: at each split a daughter drops each driver it
    carries with probability ``cladogenetic_loss`` and gains each absent one with probability
    ``cladogenetic_gain``). Because a burst can hand one daughter a key innovation and not its
    sister, speciation *itself* seeds diversification-rate heterogeneity — the genomic analogue of
    ClaSSE. One arrow points into S, so the tree is an **output** (grown jointly; take ``age`` or
    ``n_tips``, not a ``-t`` tree).

    ``model`` must have ``cladogenetic_loss > 0`` or ``cladogenetic_gain > 0`` (otherwise the
    ``species:genes`` arrow is off and this is plain :func:`simulate_gene_diversification`). Returns
    the same :class:`GeneDiversificationResult`; overlay the neutral bulk genome afterward with
    :func:`~zombi2.simulate_genomes` on ``.tree`` as usual.
    """
    if not model.is_co_diversification:
        raise ValueError(
            "simulate_co_diversification is the species<->genes joint model — set "
            "cladogenetic_loss and/or cladogenetic_gain > 0 (the species:genes burst); with both 0 "
            "this reduces to genes:species, so call simulate_gene_diversification instead")
    return simulate_gene_diversification(model, age=age, n_tips=n_tips, seed=seed, rng=rng,
                                         max_attempts=max_attempts, max_lineages=max_lineages)
