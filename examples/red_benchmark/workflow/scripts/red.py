"""Relative Evolutionary Divergence (RED) and the benchmark's pure building blocks.

This module holds the *side-effect-free* core shared by the Snakemake scripts and the unit
tests: the RED estimator (Parks et al. 2018), a clock factory that turns a config spec into a
ZOMBI2 relaxed clock, a time-tree simulator, and the per-run metrics. Nothing here imports
Snakemake, so it is importable and testable on its own.

The RED question (Rinke et al. 2021, GTDB): after a tree's branch rates are perturbed across
lineages (time -> substitutions), does RED recover the true *relative* node ages? ZOMBI2 knows
the true ages, so the ground truth is built in.
"""

from __future__ import annotations

import zlib

import numpy as np

import zombi2 as z
from zombi2.species.forward import simulate_forward
from zombi2.tools import relative_evolutionary_divergence
from zombi2.tree import prune, read_newick


# --------------------------------------------------------------------------- RED
def compute_red(tree, bl: dict) -> dict:
    """Relative Evolutionary Divergence of every node, using the branch lengths ``bl`` (node->len).

    RED is the whole point of the benchmark, so this calls the **shipped** tool,
    :func:`zombi2.tools.relative_evolutionary_divergence` (root = 0, leaves = 1), rather than a
    private copy. The perturbed branch lengths are fed in through the tool's ``branch_length`` hook,
    so RED is computed on the (unchanged) tree topology with the substitution branch lengths — the
    example thus exercises exactly the code a user gets from ``zombi2 tools red``.
    """
    return relative_evolutionary_divergence(tree, branch_length=lambda n: bl[n])


# --------------------------------------------------------------------------- clocks
def _log_bins(spread: float, n_bins: int) -> list:
    """``n_bins`` log-spaced rate multipliers symmetric about 1.0, spanning ``spread`` each way.

    RED is scale-invariant, so only the *ratio* (max/min ~ spread**2) matters, not the level.
    """
    return list(np.exp(np.linspace(-np.log(spread), np.log(spread), n_bins)))


def build_clock(spec: dict):
    """Turn a perturbation spec (a ``config['perturbations']`` entry) into a ZOMBI2 clock.

    ``clock`` selects the model; the remaining keys are its parameters:

    - ``strict`` -> :class:`StrictClock` (rate 1: branch length == time, so RED is exact).
    - ``ratevar`` -> :class:`RateVariation` (the GTDB discrete-bin clock): ``switch_rate``,
      ``spread`` and ``n_bins`` (bins are log-spaced), optional ``up_bias``.
    - ``aln`` -> :class:`AutocorrelatedLogNormalClock` (``sigma``).
    - ``cir`` -> :class:`CIRClock` (``theta``, ``sigma``).
    - ``ucln`` -> :class:`UncorrelatedLogNormalClock` (``sigma``).
    - ``whitenoise`` -> :class:`WhiteNoiseClock` (``sigma``).
    """
    kind = spec["clock"]
    if kind == "strict":
        return z.StrictClock(rate=float(spec.get("rate", 1.0)))
    if kind == "ratevar":
        bins = _log_bins(float(spec.get("spread", 5.5)), int(spec.get("n_bins", 15)))
        return z.RateVariation(bins=bins, switch_rate=float(spec["switch_rate"]),
                               up_bias=float(spec.get("up_bias", 0.5)))
    if kind == "aln":
        return z.AutocorrelatedLogNormalClock(float(spec["sigma"]))
    if kind == "cir":
        return z.CIRClock(theta=float(spec["theta"]), sigma=float(spec["sigma"]))
    if kind == "ucln":
        from zombi2.sequences.clocks import UncorrelatedLogNormalClock
        return UncorrelatedLogNormalClock(float(spec["sigma"]))
    if kind == "whitenoise":
        return z.WhiteNoiseClock(float(spec["sigma"]))
    raise ValueError(f"unknown clock kind {kind!r} in perturbation spec {spec!r}")


def strength_of(spec: dict) -> float:
    """The scalar 'perturbation strength' of a spec, for labelling/x-axes.

    ``switch_rate`` for the discrete-bin clock, ``sigma`` for the continuous ones, ``0`` for a
    strict clock. This is the *nominal* knob; the *realized* across-lineage spread is measured
    per run as ``fold_range`` (a better common x-axis, since clocks differ in units)."""
    if spec["clock"] == "strict":
        return 0.0
    return float(spec.get("switch_rate", spec.get("sigma", 0.0)))


# --------------------------------------------------------------------------- trees
def simulate_time_tree(model: str, n_tips: int, birth: float, death: float, seed: int):
    """An ultrametric time tree with known node ages.

    ``model='yule'`` (death 0) grows a pure-birth tree forward — ultrametric by construction, as
    in the prototype. ``model='bd'`` grows a birth-death tree (extinct lineages included) and
    prunes to the extant, reconstructed tree, so RED still sees all tips at the present but the
    branch-length / node-age distribution now carries an extinction signature. Raises if the
    result is not ultrametric (the RED benchmark assumes tips == present).
    """
    tree = simulate_forward(z.BirthDeath(birth=birth, death=death), n_tips=n_tips, seed=seed)
    if model == "bd" and death > 0:
        tree = prune(tree, keep="extant")
        if tree is None:
            raise RuntimeError(f"birth-death tree went extinct (seed={seed}); raise n_tips or lower death")
    # prune() keeps the ORIGINAL crown time on the reconstructed root, so rebase the tree to root=0
    # and recompute total_age as the true root-to-tip depth. Without this the returned object
    # violates the "root at time 0, total_age == depth" contract and the (1-RED)*total_age
    # calibration would be wrong on it (the workflow's Newick round-trip masks this; a direct
    # caller — e.g. the unit tests — would not). No-op for an already-rooted-at-0 yule tree.
    shift = tree.root.time
    if shift:
        for n in tree.nodes_preorder():
            n.time -= shift
    tree.total_age = max(lf.time for lf in tree.leaves())
    tip_times = [lf.time for lf in tree.leaves()]
    spread = max(tip_times) - min(tip_times)
    if spread > 1e-6 * (tree.total_age or 1.0):
        raise RuntimeError(f"time tree is not ultrametric (tip-time spread={spread:.3g}); "
                           f"model={model!r} death={death}")
    return tree


# --------------------------------------------------------------------------- metrics
def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation via numpy only (no scipy dependency)."""
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def red_recovery(tree, scaled, total_age: float):
    """Compare RED-recovered node ages to the truth for one perturbed realization.

    ``scaled`` is a :class:`RateScaledTree` (from ``clock.scale``). Returns
    ``(points, metrics)`` where ``points`` is a list of ``(name, true_age, recovered_age)`` over
    the internal (aged) nodes, and ``metrics`` carries Pearson r, Spearman r, nRMSE (of the root
    age) and the realized ``fold_range`` (max/min branch rate across lineages).
    """
    internal = tree.internal_nodes()
    red = compute_red(tree, scaled.branch_lengths)
    true_age = np.array([total_age - n.time for n in internal])
    recovered = np.array([(1.0 - red[n]) * total_age for n in internal])   # calibrate to root age
    rates = np.array([scaled.branch_rate[n] for n in tree.nodes_preorder()
                      if n.parent is not None and n in scaled.branch_rate])
    # Across-lineage rate heterogeneity, two ways: fold_range (max/min) is intuitive but
    # outlier-dominated for heavy-tailed clocks; rate_p95_p5 (95th/5th-percentile ratio) is the
    # robust measure and the sensible common x-axis when comparing clocks.
    fold = float(rates.max() / rates.min()) if rates.size and rates.min() > 0 else float("nan")
    p5 = float(np.percentile(rates, 5)) if rates.size else float("nan")
    p95_p5 = float(np.percentile(rates, 95) / p5) if rates.size and p5 > 0 else float("nan")
    logsd = float(np.std(np.log(rates))) if rates.size and rates.min() > 0 else float("nan")
    metrics = dict(
        n_internal=len(internal),
        total_age=float(total_age),
        pearson_r=float(np.corrcoef(true_age, recovered)[0, 1]),
        spearman_r=_spearman(true_age, recovered),
        nrmse=float(np.sqrt(np.mean((true_age - recovered) ** 2)) / total_age),
        fold_range=fold,
        rate_p95_p5=p95_p5,
        rate_logsd=logsd,
    )
    points = [(n.name, float(t), float(r)) for n, t, r in zip(internal, true_age, recovered)]
    return points, metrics


# --------------------------------------------------------------------------- misc
def derive_seed(*parts) -> int:
    """A stable non-negative 31-bit seed from arbitrary wildcard parts (reproducible re-runs)."""
    return zlib.crc32("|".join(str(p) for p in parts).encode()) & 0x7FFFFFFF


def load_tree(path: str):
    """Read a Newick time tree back into a ZOMBI2 ``Tree`` (node times reconstructed)."""
    with open(path) as fh:
        return read_newick(fh.read())
