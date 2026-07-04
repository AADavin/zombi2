"""The dated DTL reconciliation likelihood — faithful to ZOMBI2's simulator.

Unlike the undated model, this engine honours the species-tree **dates** and restricts
transfers to lineages alive at the **same instant** — exactly ZOMBI2's forward process
(``GenomeSimulator._choose_recipient`` draws the recipient uniformly among contemporaneous
branches). It is the Szöllősi-2013 dated model: the tree is cut into time slices at the node
times, and within a slice all alive branches coexist and exchange transfers.

Rates are **per unit time** (ZOMBI2's native δ/τ/λ), so this likelihood is directly comparable
to the values a ZOMBI2 simulation was generated under — that is what makes inject→recover work.

Time convention (from :mod:`zombi2.tree`): ``t`` increases from the root (``0``) to the present
(``T``); a branch ``e`` spans ``(parent_time_e, node_time_e]``. Everything integrates
**backward**, from the present toward the root.

Extinction ``E_e(t)`` — probability a copy on branch ``e`` at time ``t`` leaves no sampled
descendant — obeys the coupled ODE (``Ē_e`` = mean extinction over the *other* branches alive
at ``t``, since a transfer cannot target its own donor)::

    dE_e/dt = (δ+τ+λ)·E_e − λ − δ·E_e² − τ·E_e·Ē_e

with ``E_e(T)=0`` at sampled leaves, and ``E_e(node_time) = E_left·E_right`` where a branch
splits (a copy present at a speciation is inherited by both daughters).

The gene-tree DP carries ``f_{u,e}(t)`` = probability of the observed subtree at gene node
``u`` given ``u``'s lineage sits on branch ``e`` at time ``t``. It is a linear coupled ODE in
``f`` (``E`` known) with a source ``S`` that injects ``u``'s defining event::

    df_{u,e}/dt = f_{u,e}·[(δ+τ+λ) − 2δE_e − τĒ_e] − τ·E_e·f̄_e − S_{u,e}

    S_{u,e}(t) = δ·f_{v,e}·f_{w,e}                                  # duplication of u on e
               + (τ/(n−1))·Σ_{h≠e}(f_{v,e}·f_{w,h} + f_{w,e}·f_{v,h})  # transfer of one child off e

where ``v,w`` are ``u``'s children and ``f̄_e`` the mean of ``f_u`` over the other alive
branches (a transferred copy arriving on ``e`` while the donor copy is lost). At a species
node the speciation source ``f_{v,left}·f_{w,right}+f_{v,right}·f_{w,left}`` and the
speciation-loss term ``f_{u,left}·E_right+f_{u,right}·E_left`` enter as boundary jumps. Gene
leaves seed ``f_{u,s}(T)=1`` at their species' leaf branch.

Origination: ``"root"`` places the family as one lineage on the root branch at ``t=0`` (exact
for ZOMBI2's root-seeded ``initial_size`` families); ``"uniform"`` averages the root gene node
over the tops of all branches.

The integration is an explicit backward Euler-style sweep on a per-slice sub-grid; ``n_steps``
(sub-steps per slice, auto-raised to keep ``dt·(δ+τ+λ)`` small) controls resolution. The pure
loops here are the part a Rust core will replace.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .genetree import GeneTree
from .species import SpeciesTree

_EPS = 1e-9


@dataclass(frozen=True, slots=True)
class DatedDTL:
    """Per-unit-time DTL rates for the dated model (ZOMBI2's native δ/τ/λ)."""

    dup: float = 0.0
    transfer: float = 0.0
    loss: float = 0.0

    def __post_init__(self):
        if min(self.dup, self.transfer, self.loss) < 0:
            raise ValueError("dup/transfer/loss rates must be >= 0")


class _DatedEngine:
    """Slices the species tree, integrates extinction, then runs the gene-tree DP.

    One engine is built per (species tree, rates, resolution); it can score many gene trees
    against the same background (extinction is computed once).
    """

    def __init__(self, sp: SpeciesTree, model: DatedDTL, n_steps: int = 100):
        self.sp = sp
        self.d, self.t, self.lam = model.dup, model.transfer, model.loss
        self.tot = self.d + self.t + self.lam
        self._build_slices(n_steps)
        self.E_bottom: dict[int, float] = {}   # branch -> E at its root-ward end (parent_time)
        self._compute_extinction()

    # --- slice structure ---------------------------------------------------
    def _build_slices(self, n_steps: int) -> None:
        sp = self.sp
        T = max(b.time for b in sp.branches)
        self.T = T
        # boundaries = every node time (leaves at T, internal nodes, root at 0)
        times = sorted({round(b.time, 12) for b in sp.branches} | {0.0})
        self.grids: list[dict] = []
        for k in range(len(times) - 1):
            lo, hi = times[k], times[k + 1]
            alive = [i for i, b in enumerate(sp.branches)
                     if b.parent_time <= lo + _EPS and b.time >= hi - _EPS and b.length > _EPS]
            # sub-steps: at least n_steps, more if needed to keep dt*(total rate) small
            span = hi - lo
            D = max(n_steps, math.ceil(span * self.tot * 4.0)) if span > 0 else 1
            dt = span / D
            col = {e: j for j, e in enumerate(alive)}
            self.grids.append({"lo": lo, "hi": hi, "alive": alive, "col": col,
                               "D": D, "dt": dt, "E": None})

    def _split_children(self, e: int) -> tuple[int, int] | None:
        b = self.sp.branches[e]
        return None if b.is_leaf else (b.left, b.right)

    # --- extinction --------------------------------------------------------
    def _compute_extinction(self) -> None:
        d, t, lam, tot = self.d, self.t, self.lam, self.tot
        # process slices present-first (reverse: last grid is nearest present)
        for g in reversed(self.grids):
            alive, col, D, dt, hi = g["alive"], g["col"], g["D"], g["dt"], g["hi"]
            n = len(alive)
            # boundary E at t_hi for each alive branch
            E = [0.0] * n
            for e in alive:
                b = self.sp.branches[e]
                if abs(b.time - hi) < _EPS:          # e ends at t_hi
                    ch = self._split_children(e)
                    E[col[e]] = 0.0 if ch is None else self.E_bottom[ch[0]] * self.E_bottom[ch[1]]
                else:                                 # e continues from above
                    E[col[e]] = self.E_bottom[e]
            rows = [E[:]]
            for _ in range(D):                        # backward Euler, decreasing t
                tot_E = math.fsum(E)
                new = [0.0] * n
                for j in range(n):
                    Ej = E[j]
                    ebar = (tot_E - Ej) / (n - 1) if n > 1 else 0.0
                    dE_dt = tot * Ej - lam - d * Ej * Ej - t * Ej * ebar
                    new[j] = Ej - dt * dE_dt
                E = new
                rows.append(E[:])
            g["E"] = rows                             # rows[0]=t_hi ... rows[D]=t_lo
            for e in alive:                           # record value at root-ward end (t_lo)
                self.E_bottom[e] = E[col[e]]

    # --- gene-tree DP ------------------------------------------------------
    def gene_loglik(self, gt: GeneTree, origination: str = "root") -> float:
        d, t, tot = self.d, self.t, self.tot
        sp = self.sp
        # f_bottom[u][branch] = f_{u}(branch, branch.parent_time); f_slices[u][k] = rows
        f_bottom: list[dict[int, float]] = [dict() for _ in range(gt.n)]
        f_rows: list[list] = [None] * gt.n            # per node: list over slices of E-shaped rows

        for u in range(gt.n):                         # gene tree post-order
            node = gt.nodes[u]
            internal = not node.is_leaf
            v, w = (node.left, node.right) if internal else (None, None)
            per_slice: list = [None] * len(self.grids)
            fb = f_bottom[u]

            for gi in range(len(self.grids) - 1, -1, -1):   # present-first
                g = self.grids[gi]
                alive, col, D, dt, hi = g["alive"], g["col"], g["D"], g["dt"], g["hi"]
                Erows = g["E"]
                n = len(alive)
                # boundary f at t_hi
                f = [0.0] * n
                for e in alive:
                    b = sp.branches[e]
                    j = col[e]
                    if abs(b.time - hi) < _EPS:       # e ends at t_hi
                        ch = self._split_children(e)
                        if ch is None:                # leaf branch at present: observation
                            f[j] = 1.0 if (node.is_leaf and sp.leaf_index.get(node.species) == e) else 0.0
                        else:
                            a, c = ch
                            Ea, Ec = self.E_bottom[a], self.E_bottom[c]
                            f[j] = fb.get(a, 0.0) * Ec + fb.get(c, 0.0) * Ea   # SL
                            if internal:              # speciation of u maps here
                                fva, fvc = f_bottom[v].get(a, 0.0), f_bottom[v].get(c, 0.0)
                                fwa, fwc = f_bottom[w].get(a, 0.0), f_bottom[w].get(c, 0.0)
                                f[j] += fva * fwc + fvc * fwa
                    else:                             # continues from above
                        f[j] = fb.get(e, 0.0)
                rows = [f[:]]
                fv_slice = f_rows[v][gi] if internal else None
                fw_slice = f_rows[w][gi] if internal else None
                for s in range(D):
                    Erow = Erows[s]
                    tot_f = math.fsum(f)
                    tot_E = math.fsum(Erow)
                    fv_row = fv_slice[s] if internal else None
                    fw_row = fw_slice[s] if internal else None
                    if internal and n > 1:
                        sum_fv, sum_fw = math.fsum(fv_row), math.fsum(fw_row)
                    new = [0.0] * n
                    for j in range(n):
                        fj = f[j]
                        Ej = Erow[j]
                        ebar = (tot_E - Ej) / (n - 1) if n > 1 else 0.0
                        pbar = (tot_f - fj) / (n - 1) if n > 1 else 0.0
                        homog = fj * (tot - 2 * d * Ej - t * ebar) - t * Ej * pbar
                        src = 0.0
                        if internal:
                            src += d * fv_row[j] * fw_row[j]
                            if n > 1:
                                src += t / (n - 1) * (fv_row[j] * (sum_fw - fw_row[j])
                                                      + fw_row[j] * (sum_fv - fv_row[j]))
                        new[j] = fj - dt * homog + dt * src
                    f = new
                    rows.append(f[:])
                per_slice[gi] = rows
                for e in alive:
                    fb[e] = f[col[e]]
            f_rows[u] = per_slice

        return self._originate(gt, f_bottom, origination)

    def _originate(self, gt: GeneTree, f_bottom, origination: str) -> float:
        sp = self.sp
        root_u = gt.root
        if origination == "uniform":
            reals = [i for i, b in enumerate(sp.branches) if b.length > _EPS]
            like = math.fsum(f_bottom[root_u].get(e, 0.0) for e in reals) / len(reals)
        elif origination == "root":
            a, c = sp.branches[sp.root].left, sp.branches[sp.root].right
            Ea, Ec = self.E_bottom[a], self.E_bottom[c]
            fr_a, fr_c = f_bottom[root_u].get(a, 0.0), f_bottom[root_u].get(c, 0.0)
            like = fr_a * Ec + fr_c * Ea                       # SL at the root split
            node = gt.nodes[root_u]
            if not node.is_leaf:                               # gene root == root speciation
                v, w = node.left, node.right
                fva, fvc = f_bottom[v].get(a, 0.0), f_bottom[v].get(c, 0.0)
                fwa, fwc = f_bottom[w].get(a, 0.0), f_bottom[w].get(c, 0.0)
                like += fva * fwc + fvc * fwa
        else:
            raise ValueError(f"origination must be 'root' or 'uniform', got {origination!r}")
        if like <= 0.0:
            return -math.inf
        return math.log(like)


def dated_loglik(gene_tree: GeneTree, sp: SpeciesTree, model: DatedDTL,
                 *, origination: str = "root", n_steps: int = 100) -> float:
    """Log marginal reconciliation likelihood ``P(gene_tree | sp, model)`` under the dated model.

    ``n_steps`` sets the per-slice integration resolution (increase to check convergence);
    ``origination`` is ``"root"`` or ``"uniform"``. See the module docstring for the model.
    """
    return _DatedEngine(sp, model, n_steps).gene_loglik(gene_tree, origination)


def dated_extinction(sp: SpeciesTree, model: DatedDTL, *, n_steps: int = 100) -> dict[int, float]:
    """Extinction probability at the **root-ward end** of every branch — ``{branch_index: E}``.
    Exposed for validation against closed forms (e.g. the birth-death limit τ=0)."""
    return dict(_DatedEngine(sp, model, n_steps).E_bottom)


def dated_joint_loglik(gene_trees, sp: SpeciesTree, model: DatedDTL, *,
                       origination: str = "root", n_extinct: int = 0, n_steps: int = 100) -> float:
    """Joint dated log-lik of many gene trees sharing one species tree and rates.

    Builds the background (extinction) once and reuses it across every tree — the efficient
    path for rate inference / inject-recover. ``n_extinct`` adds ``k·log P(no survivor)`` for
    ``k`` families that were seeded but left no extant copy; include it (with the count of
    fully-extinct families) so the likelihood of a *set* of families seeded together is
    unbiased. ``P(no survivor)`` is taken under the same ``origination`` model.
    """
    eng = _DatedEngine(sp, model, n_steps)
    ll = 0.0
    if n_extinct:
        if origination == "root":
            a, c = sp.branches[sp.root].left, sp.branches[sp.root].right
            pe = eng.E_bottom[a] * eng.E_bottom[c]
        elif origination == "uniform":
            reals = [i for i, b in enumerate(sp.branches) if b.length > _EPS]
            pe = math.fsum(eng.E_bottom[e] for e in reals) / len(reals)
        else:
            raise ValueError(f"origination must be 'root' or 'uniform', got {origination!r}")
        ll += n_extinct * (math.log(pe) if pe > 0 else -math.inf)
    for gt in gene_trees:
        ll += eng.gene_loglik(gt, origination)
    return ll
