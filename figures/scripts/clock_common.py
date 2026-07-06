"""Shared helpers for the molecular-clock / sequence-evolution figures (Ch15, Ch16).

The clock *classes* (StrictClock, UncorrelatedLogNormalClock, ...) live in
``zombi2.clocks`` on ``main``; this branch predates them. So these figures
re-implement each clock's rate draw here, following the exact formulas of
Chapter 16, seeded with ``numpy``'s default_rng for reproducibility. Every clock
returns two dicts keyed by ete3-node name:

  * ``rate[name]``      -- the (branch-average) rate multiplier on the branch
                           leading to that node; used to *paint* the branch.
  * ``subst_len[name]`` -- the substitution length of that branch, i.e. the time
                           integral of the rate: for the memoryless clocks just
                           ``rate * dist``; for CIR the sub-stepped integral.

A phylogram is then the tree with branch lengths ``subst_len`` instead of time.

All figures share ONE canonical species tree (``build_tree``) so they read as a
family, exactly as the BM / OU trait figures do.

Run nothing directly; imported by the ``fig_clock_*`` / ``fig_seq_*`` scripts.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from zombi2.species import BirthDeath, simulate_species_tree

from model_common import zombi_to_ete3
# viridis colormap shared with the trait figures (the house continuous ramp)
from fig_trait_bm import VIRIDIS, hexc, viridis  # noqa: F401  (re-exported)

# --- the one canonical tree, reused across every clock figure ----------------
CANON = dict(n_tips=14, age=1.0, seed=3)


def build_tree(n_tips: int | None = None, age: float | None = None, seed: int | None = None):
    """The shared time-calibrated species tree (ete3), branch lengths in TIME."""
    p = dict(CANON)
    if n_tips is not None:
        p["n_tips"] = n_tips
    if age is not None:
        p["age"] = age
    if seed is not None:
        p["seed"] = seed
    zt = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=p["n_tips"], age=p["age"],
                               direction="backward", seed=p["seed"])
    return zombi_to_ete3(zt)


# --- layout: time-from-root and a leaf-ordered y for every node --------------
def node_times(tree) -> tuple[dict, float]:
    """{name: time from root} for every node, plus the present (max tip time)."""
    tfo = {}
    for n in tree.traverse("preorder"):
        tfo[n.name] = 0.0 if n.is_root() else tfo[n.up.name] + n.dist
    return tfo, max(tfo.values())


def leaf_ys(tree) -> tuple[dict, int]:
    """{name: y-index}; leaves get 0..nleaf-1 (draw order), internals the mean of
    their children -- the standard ladderized vertical layout."""
    ys, i = {}, 0
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            ys[n.name] = float(i)
            i += 1
        else:
            ys[n.name] = sum(ys[c.name] for c in n.children) / len(n.children)
    return ys, i


def subst_dist_to_root(tree, subst_len: dict) -> dict:
    """Cumulative substitution distance from the root, per node -- the x-axis of a
    phylogram. Root = 0; a child = parent + its own branch's substitution length."""
    dr = {}
    for n in tree.traverse("preorder"):
        dr[n.name] = 0.0 if n.is_root() else dr[n.up.name] + subst_len[n.name]
    return dr


# --- rate -> colour: a shared LOG scale so a colour means the same rate in every
#     panel (rate 1 sits mid-viridis; faster warmer, slower cooler) -----------
RATE_LO, RATE_HI = 1.0 / 3.0, 3.0   # colour clamps (3x slow ... 3x fast)


def rate_rgb(r: float, lo: float = RATE_LO, hi: float = RATE_HI):
    """viridis colour for a rate on a symmetric log scale in [lo, hi]."""
    r = max(lo, min(hi, r))
    t = (math.log(r) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return viridis(t)


def rate_hex(r: float, lo: float = RATE_LO, hi: float = RATE_HI) -> str:
    return hexc(rate_rgb(r, lo, hi))


# =============================================================================
# The clock family.  Each returns (rate, subst_len), dicts keyed by node name.
# The root's own (stub) branch is given rate ``root_rate`` where meaningful.
# =============================================================================
def strict(tree, rate: float = 1.0):
    """One rate everywhere: the phylogram is the chronogram uniformly stretched."""
    R, S = {}, {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            R[n.name] = rate
            continue
        R[n.name] = rate
        S[n.name] = rate * n.dist
    return R, S


def uncorrelated_lognormal(tree, sigma: float = 0.5, seed: int = 0):
    """Each branch an i.i.d. lognormal(mean 1) multiplier: exp(N(-s^2/2, s))."""
    rng = np.random.default_rng(seed)
    R, S = {}, {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            R[n.name] = 1.0
            continue
        r = float(math.exp(rng.normal(-sigma * sigma / 2.0, sigma)))
        R[n.name] = r
        S[n.name] = r * n.dist
    return R, S


def uncorrelated_gamma(tree, shape: float = 3.0, seed: int = 0):
    """Each branch an i.i.d. Gamma(shape, 1/shape): mean 1, variance 1/shape."""
    rng = np.random.default_rng(seed)
    R, S = {}, {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            R[n.name] = 1.0
            continue
        r = float(rng.gamma(shape, 1.0 / shape))
        R[n.name] = r
        S[n.name] = r * n.dist
    return R, S


def white_noise(tree, sigma: float = 0.5, seed: int = 0):
    """Branch rate ~ Gamma with mean 1 and variance sigma^2 / dt: long branches
    average the noise away (rate near 1), short branches are highly variable."""
    rng = np.random.default_rng(seed)
    R, S = {}, {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            R[n.name] = 1.0
            continue
        dt = max(n.dist, 1e-9)
        var = sigma * sigma / dt
        shape = 1.0 / var                # mean 1 => shape*scale = 1, scale = var
        r = float(rng.gamma(shape, var))
        R[n.name] = r
        S[n.name] = r * n.dist
    return R, S


def autocorrelated_lognormal(tree, sigma: float = 0.3, seed: int = 0, root_rate: float = 1.0):
    """Geometric random walk down the tree: R_child = R_parent * exp(N(0, s*sqrt(l)))."""
    rng = np.random.default_rng(seed)
    R, S = {}, {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            R[n.name] = root_rate
            continue
        step = math.exp(rng.normal(0.0, sigma * math.sqrt(max(n.dist, 0.0))))
        r = R[n.up.name] * step
        R[n.name] = r
        S[n.name] = r * n.dist
    return R, S


def cir(tree, theta: float = 3.0, sigma: float = 0.5, mean: float = 1.0, seed: int = 0,
        substeps_per_unit: int = 240):
    """Cox-Ingersoll-Ross mean-reverting rate, autocorrelated AND varying WITHIN a
    branch. dr = theta*(mean - r) dt + sigma*sqrt(r) dW, Euler sub-stepped; the
    rate is continuous down the tree (child starts where its parent ended).

    Returns (rate, subst_len, segments) where ``segments[name]`` is a list of
    (rate, dt) pieces for the branch -- used by the CIR figure to show the
    within-branch path; ``rate[name]`` is the time-average over the branch."""
    rng = np.random.default_rng(seed)
    end_rate = {}          # rate value at the END of each node's branch
    R, S, SEG = {}, {}, {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            end_rate[n.name] = mean
            R[n.name] = mean
            continue
        r = end_rate[n.up.name]
        nsteps = max(1, int(round(n.dist * substeps_per_unit)))
        dt = n.dist / nsteps
        segs, integ = [], 0.0
        for _ in range(nsteps):
            r_mid = r  # left-point rate over the step
            segs.append((r_mid, dt))
            integ += r_mid * dt
            r = r + theta * (mean - r) * dt + sigma * math.sqrt(max(r, 0.0)) * math.sqrt(dt) * rng.normal()
            r = max(r, 1e-6)
        end_rate[n.name] = r
        S[n.name] = integ
        R[n.name] = integ / n.dist if n.dist > 0 else mean
        SEG[n.name] = segs
    return R, S, SEG


def discrete_bin(tree, bins=(0.25, 0.5, 1.0, 2.0, 4.0), switch_rate: float = 1.0, seed: int = 0,
                 substeps_per_unit: int = 240, start_bin: int | None = None):
    """GTDB-style ordered rate bins with a nearest-neighbour Markov walk along the
    tree (autocorrelated, within-branch segments). Same (rate, subst_len, segments)
    shape as :func:`cir`, but rate takes only the discrete bin values."""
    rng = np.random.default_rng(seed)
    bins = list(bins)
    nb = len(bins)
    start = nb // 2 if start_bin is None else start_bin
    end_bin = {}
    R, S, SEG = {}, {}, {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            end_bin[n.name] = start
            R[n.name] = bins[start]
            continue
        b = end_bin[n.up.name]
        nsteps = max(1, int(round(n.dist * substeps_per_unit)))
        dt = n.dist / nsteps
        segs, integ = [], 0.0
        for _ in range(nsteps):
            segs.append((bins[b], dt))
            integ += bins[b] * dt
            # nearest-neighbour jump with prob switch_rate*dt, direction 50/50 (reflected at ends)
            if rng.random() < switch_rate * dt:
                if b == 0:
                    b = 1
                elif b == nb - 1:
                    b = nb - 2
                else:
                    b += 1 if rng.random() < 0.5 else -1
        end_bin[n.name] = b
        S[n.name] = integ
        R[n.name] = integ / n.dist if n.dist > 0 else bins[start]
        SEG[n.name] = segs
    return R, S, SEG


# --- continuous 1-D sample paths (for the CIR figure, off-tree) --------------
def cir_paths(T: float, n_paths: int, theta: float, sigma: float, mean: float,
              seed: int = 0, r0: float | None = None, nstep: int = 600):
    """n_paths CIR trajectories r(t) on [0, T]; returns (t, list_of_paths)."""
    rng = np.random.default_rng(seed)
    dt = T / nstep
    t = np.linspace(0.0, T, nstep + 1)
    out = []
    for _ in range(n_paths):
        r = mean if r0 is None else r0
        path = [r]
        for _ in range(nstep):
            r = r + theta * (mean - r) * dt + sigma * math.sqrt(max(r, 0.0)) * math.sqrt(dt) * rng.normal()
            r = max(r, 1e-6)
            path.append(r)
        out.append(np.array(path))
    return t, out


def gbm_paths(T: float, n_paths: int, sigma: float, seed: int = 0, r0: float = 1.0,
              nstep: int = 600):
    """n_paths geometric-random-walk (autocorrelated-lognormal) rate trajectories
    r(t) = r0 * exp(sigma * W(t)); returns (t, list_of_paths)."""
    rng = np.random.default_rng(seed)
    dt = T / nstep
    t = np.linspace(0.0, T, nstep + 1)
    out = []
    for _ in range(n_paths):
        logr = math.log(r0)
        path = [r0]
        for _ in range(nstep):
            logr += sigma * math.sqrt(dt) * rng.normal()
            path.append(math.exp(logr))
        out.append(np.array(path))
    return t, out
