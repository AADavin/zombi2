"""The ALEml_undated / GeneRax ``UndatedDTL`` reconciliation likelihood.

This is the *undated* ALE model: the species tree carries no dates, each branch has the same
event *odds*, and a transfer may land on any other branch. It is the model people actually
run on real data (ALEml_undated, ALErax), so it is the natural first engine — and, unlike the
dated model, it has closed-form limits we can check by hand. It is **not** a faithful match to
ZOMBI2's dated, contemporaneous-transfer simulator; that is what the forthcoming ``DatedDTL``
engine is for.

Parameterisation (per branch, as in GeneRax). Given duplication/transfer/loss *rates*
``d, t, l`` (here dimensionless odds relative to a speciation), one traversal step of a branch
consumes exactly one of four slots — duplicate, transfer, lose, or speciate/sample — with

    denom = 1 + d + t + l,   pD = d/denom,  pT = t/denom,  pL = l/denom,  pS = 1/denom.

Extinction ``E[e]`` — the probability a single gene copy at the top of branch ``e`` leaves no
sampled descendant — is the coupled fixed point (``Ē_e`` = mean extinction over the *other*
branches, since a transfer cannot target its own donor)::

    E[e] = pL + pD·E[e]²  + pT·E[e]·Ē_e  + pS·(E[left]·E[right]   if e internal
                                               0                  if e a sampled leaf)

The gene-tree DP fills, per gene node ``u`` (post-order) and species branch ``e``, ``P[u,e]`` =
probability of the subtree at ``u`` given ``u``'s lineage sits at the top of ``e``. Writing
``a = P[left(u)]``, ``b = P[right(u)]`` and ``P̄_e(x) = mean over branches ≠ e of x``::

    A[e] =  pS·(a[f]·b[g] + a[g]·b[f])          # S  (speciation into daughters f,g of e)
          + 2·pD·a[e]·b[e]                       # D  (duplication on e; ordered mother/daughter)
          + pT·(a[e]·P̄_e(b) + b[e]·P̄_e(a))       # T  (one child stays on e, one is transferred out)
       (+ pS   if u is a leaf sampled from species e)

    P[u,e] = ( A[e]
             + pS·(P[u,f]·E[g] + P[u,g]·E[f])    # SL  (speciate, follow one daughter, lose the other)
             + pT·E[e]·P̄_e(P[u]) )               # TL  (donor copy lost, transferred copy carries u)
             / ( 1 − 2·pD·E[e] − pT·Ē_e )        # DL + TL  (the *other* copy is lost, u stays on e)

The SL term references child branches (already filled) and the TL term references the mean of
``P[u]`` over all branches, so the ``e``-recursion is itself a small fixed point, solved by a
few Gauss-Seidel sweeps.

Two exact limits pin the S / SL / loss / extinction machinery (see ``tests/test_alelite.py``):

* d = t = 0, gene tree perfectly matching a species subtree with ``k`` tips ⇒
  ``P = pS^(2k−1)`` (one slot per speciation and per tip sampling; nothing to lose).
* d = t = 0, species tree ``((A,B))`` but the gene present only in ``A`` ⇒
  ``P = pS²·pL = l/(1+l)³`` (root speciates, A sampled, the B copy is lost).

Origination (where the family enters the tree): ``"root"`` conditions on the family being
present on the root branch — exact for ZOMBI2's root-seeded ``initial_size`` families;
``"uniform"`` averages the root gene node over all branches.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .genetree import GeneTree
from .species import SpeciesTree

_MAX_ITERS = 500
_TOL = 1e-13


@dataclass(frozen=True, slots=True)
class UndatedDTL:
    """Undated ALE rates (per-branch odds). See the module docstring for the parameterisation."""

    dup: float = 0.0
    transfer: float = 0.0
    loss: float = 0.0

    def probs(self) -> tuple[float, float, float, float]:
        """Return ``(pD, pT, pL, pS)`` — the normalised per-slot probabilities."""
        d, t, lo = self.dup, self.transfer, self.loss
        if min(d, t, lo) < 0:
            raise ValueError("dup/transfer/loss rates must be >= 0")
        denom = 1.0 + d + t + lo
        return d / denom, t / denom, lo / denom, 1.0 / denom


def _mean_over(vec: list[float], total: float, e: int, nb, n: int) -> float:
    """Mean of ``vec`` over the transfer recipients of branch ``e``. ``nb is None`` = every
    other branch (plain undated, ``(total - vec[e])/(n-1)``); otherwise the explicit
    time-overlapping list (reldated)."""
    if nb is None:
        return (total - vec[e]) / (n - 1) if n > 1 else 0.0
    lst = nb[e]
    return math.fsum(vec[f] for f in lst) / len(lst) if lst else 0.0


def _transfer_neighbors(sp: SpeciesTree, transfers: str):
    """Per-branch transfer-recipient sets. ``"global"`` → ``None`` (any branch, plain undated);
    ``"dated"`` → for each branch the branches whose time interval **overlaps** it (reldated:
    transfers stay time-consistent, using the dates but no sub-branch integration)."""
    if transfers == "global":
        return None
    if transfers != "dated":
        raise ValueError(f"transfers must be 'global' or 'dated', got {transfers!r}")
    n = sp.n
    tol = 1e-9 * max((b.time for b in sp.branches), default=1.0)
    nb: list[list[int]] = [[] for _ in range(n)]
    for e in range(n):
        be = sp.branches[e]
        for f in range(n):
            if f == e:
                continue
            bf = sp.branches[f]
            overlap = min(be.time, bf.time) - max(be.parent_time, bf.parent_time)
            if overlap > tol:
                nb[e].append(f)
    return nb


def extinction(sp: SpeciesTree, model: UndatedDTL, nb=None) -> list[float]:
    """Solve the coupled extinction fixed point ``E[e]`` for every species branch.

    ``nb`` is the transfer neighborhood (``None`` = plain undated; a per-branch list =
    reldated), controlling which branches a transfer's other copy can land on."""
    pD, pT, pL, pS = model.probs()
    n = sp.n
    E = [0.0] * n
    for _ in range(_MAX_ITERS):
        total = math.fsum(E)
        delta = 0.0
        for e in range(n):  # post-order: children before parents
            b = sp.branches[e]
            ebar = _mean_over(E, total, e, nb, n)
            child = 0.0 if b.is_leaf else E[b.left] * E[b.right]
            new = pL + pD * E[e] * E[e] + pT * E[e] * ebar + pS * child
            delta = max(delta, abs(new - E[e]))
            total += new - E[e]
            E[e] = new
        if delta < _TOL:
            break
    return E


def _propagate(Pu: list[float], A: list[float], sp: SpeciesTree,
               E: list[float], pD: float, pT: float, pS: float, nb=None) -> None:
    """In-place solve of the per-gene-node ``e``-recursion (SL/DL/TL fixed point) given the
    event-independent accumulator ``A``. Overwrites ``Pu``. ``nb`` = transfer neighborhood."""
    n = sp.n
    total_E = math.fsum(E)
    for _ in range(_MAX_ITERS):
        total = math.fsum(Pu)
        delta = 0.0
        for e in range(n):
            b = sp.branches[e]
            ebar = _mean_over(E, total_E, e, nb, n)
            pbar = _mean_over(Pu, total, e, nb, n)
            sl = 0.0 if b.is_leaf else pS * (Pu[b.left] * E[b.right] + Pu[b.right] * E[b.left])
            denom = 1.0 - 2.0 * pD * E[e] - pT * ebar
            new = (A[e] + sl + pT * E[e] * pbar) / denom
            delta = max(delta, abs(new - Pu[e]))
            total += new - Pu[e]
            Pu[e] = new
        if delta < _TOL:
            break


def undated_loglik(gene_tree: GeneTree, sp: SpeciesTree, model: UndatedDTL,
                   *, origination: str = "root", transfers: str = "global") -> float:
    """Log of the marginal reconciliation likelihood ``P(gene_tree | sp, model)``.

    Sums over every reconciliation of the (fixed) gene tree against the species tree under the
    undated DTL model. ``origination`` is ``"root"`` (family present on the root branch) or
    ``"uniform"`` (root gene node averaged over all branches). ``transfers`` is ``"global"``
    (plain undated — a transfer may land on any branch) or ``"dated"`` (**reldated** — a
    transfer may only land on a branch that overlaps the donor in time, using the species-tree
    dates to stay time-consistent without full slicing; see :func:`reldated_loglik`).
    """
    pD, pT, pL, pS = model.probs()
    n = sp.n
    nb = _transfer_neighbors(sp, transfers)
    E = extinction(sp, model, nb)

    P: list[list[float]] = [None] * gene_tree.n  # type: ignore[list-item]
    for u in range(gene_tree.n):  # gene tree post-order: children before parents
        g = gene_tree.nodes[u]
        A = [0.0] * n
        if g.is_leaf:
            leaf = sp.leaf_index.get(g.species)
            if leaf is None:
                raise KeyError(f"gene tip species {g.species!r} is not a species-tree leaf")
            A[leaf] = pS
        else:
            a, b = P[g.left], P[g.right]
            atot, btot = math.fsum(a), math.fsum(b)
            for e in range(n):
                br = sp.branches[e]
                term = 2.0 * pD * a[e] * b[e]  # D
                if not br.is_leaf:  # S
                    f, gg = br.left, br.right
                    term += pS * (a[f] * b[gg] + a[gg] * b[f])
                if pT > 0.0:  # T: one child transferred to a recipient branch
                    abar = _mean_over(a, atot, e, nb, n)
                    bbar = _mean_over(b, btot, e, nb, n)
                    term += pT * (a[e] * bbar + b[e] * abar)
                A[e] = term
        Pu = [0.0] * n
        _propagate(Pu, A, sp, E, pD, pT, pS, nb)
        P[u] = Pu

    root_P = P[gene_tree.root]
    if origination == "root":
        like = root_P[sp.root]
    elif origination == "uniform":
        like = math.fsum(root_P) / n
    else:
        raise ValueError(f"origination must be 'root' or 'uniform', got {origination!r}")
    if like <= 0.0:
        return -math.inf
    return math.log(like)


def reldated_loglik(gene_tree: GeneTree, sp: SpeciesTree, model: UndatedDTL,
                    *, origination: str = "root") -> float:
    """The **reldated** likelihood: the undated model with transfers restricted to branches
    that overlap the donor in time.

    This is the middle ground between :func:`undated_loglik` (transfers to any branch — cheap
    but time-inconsistent) and the fully time-sliced :func:`~zombi2.alelite.dated_loglik`
    (exact but costly). It reuses the undated per-branch recursion but forbids transfers to
    non-contemporaneous branches, using the species-tree dates — the trick recent ALE versions
    use to keep undated transfers consistent without sub-branch integration. Same
    :class:`UndatedDTL` rates (per-branch odds).
    """
    return undated_loglik(gene_tree, sp, model, origination=origination, transfers="dated")
