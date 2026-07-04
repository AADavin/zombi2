"""Simulation analysis for ABC inference of the coupled (Potts) gene-family model.

Generates every figure in `report/report.tex` and prints the numbers the report cites.
All results are deterministic (fixed seeds). Run from anywhere:

    python analyses/coupling_abc/run_analysis.py

Compute: ~2-3 min (two reference sets of coupled simulations, reused across figures).
"""
from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import zombi2 as z
from zombi2.coupling import pathway_blocks

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

# --- house style (ink on white, muted accents, clean spines) ----------------------
INK, MUTED, BLUE, ORANGE = "#1A1A1A", "#6A6A6A", "#2B6CB0", "#C05621"
plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 200, "font.size": 10.5, "font.family": "serif",
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlelocation": "left",
    "axes.titleweight": "bold", "axes.titlepad": 8, "legend.frameon": False,
})


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


# --- model + simulation -----------------------------------------------------------
SIZE, BLOCKS = 4, 15
N = SIZE * BLOCKS
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.5), n_tips=60, age=1.0, seed=1)
S = len(tree.extant_leaves())
H_NEUTRAL = lambda J: 0.33 - 2.35 * J           # holds prevalence ~0.5 (calibrated)


def spec(J, h):
    return pathway_blocks([SIZE] * BLOCKS, within=J, between=0.0, h=h,
                          base_loss=1.0, transfer=0.4, beta=1.0)


def presence(J, h, seed):
    return (z.simulate_coupled(tree, spec(J, h), seed=seed).profiles.matrix > 0).astype(float)


# --- summary statistics on a presence matrix P (families x species) ---------------
def s_freq(P):                                   # marginal: gene frequency spectrum
    return np.bincount(P.sum(1).astype(int), minlength=S + 1)[1:S + 1].astype(float)


def _corr(P):
    keep = P.std(1) > 0
    return np.corrcoef(P[keep]) if keep.sum() >= 5 else None


def s_hist(P):                                   # co-occurrence: correlation histogram
    C = _corr(P)
    if C is None:
        return np.zeros(12)
    return np.histogram(C[np.triu_indices(C.shape[0], 1)], bins=12, range=(-1, 1))[0].astype(float)


def s_eig(P, k=8):                               # co-occurrence: top eigenvalues
    C = _corr(P)
    if C is None:
        return np.zeros(k)
    ev = np.sort(np.linalg.eigvalsh(C))[::-1]
    return np.pad(ev[:k], (0, max(0, k - len(ev))))


def s_clust(P, thr=0.35):                         # co-occurrence: module/clustering features
    C = _corr(P)
    if C is None:
        return np.zeros(3)
    A = (C > thr).astype(float)
    np.fill_diagonal(A, 0.0)
    tri = float(np.trace(A @ A @ A)) / 6.0
    deg = A.sum(1)
    triples = float(np.sum(deg * (deg - 1))) / 2.0
    return np.array([A.sum() / 2.0, tri, 3 * tri / triples if triples > 0 else 0.0])


def within_between_gap(P):
    """Diagnostic (uses known blocks): mean within- vs between-pathway presence correlation."""
    keep = P.std(1) > 0
    idx = np.where(keep)[0]
    if len(idx) < 3:
        return np.nan
    C = np.corrcoef(P[keep])
    iu = np.triu_indices(len(idx), k=1)
    same = (idx[iu[0]] // SIZE) == (idx[iu[1]] // SIZE)
    corr = C[iu]
    return corr[same].mean() - corr[~same].mean()


def abc_posterior(ref, target, Js, accept=0.05):
    """Scale-free rejection ABC: return the accepted J values."""
    sd = ref.std(0); sd[sd == 0] = 1.0
    d = np.linalg.norm((ref - target) / sd, axis=1)
    k = max(1, round(accept * len(Js)))
    return Js[np.argsort(d)[:k]]


def med_ci(x):
    return float(np.median(x)), float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))


TARGET_SEEDS = range(700, 708)      # average RMSE over several target datasets per J (stable)


def robust_rmse(ref, summ_fn, hfun, Js, trueJs):
    """Mean recovery RMSE over TARGET_SEEDS target datasets per true J (a stable estimate)."""
    err = []
    for tJ in trueJs:
        for sd in TARGET_SEEDS:
            Pt = presence(tJ, hfun(tJ), sd)
            err.append((np.median(abc_posterior(ref, summ_fn(Pt), Js)) - tJ) ** 2)
    return float(np.sqrt(np.mean(err)))


# --- reference simulation sets (generated once, reused across figures) -------------
N_SIMS = 500
rng = np.random.default_rng(0)
Js = rng.uniform(0.0, 1.5, N_SIMS)
TRUE = [0.0, 0.3, 0.6, 0.9, 1.2]
results = {"model": dict(n_species=S, n_families=N, block_size=SIZE, blocks=BLOCKS,
                         n_sims=N_SIMS, prior="U(0, 1.5)")}

t0 = time.perf_counter()
P_real = [presence(J, -0.5, 1000 + i) for i, J in enumerate(Js)]           # realistic regime
P_neut = [presence(J, H_NEUTRAL(J), 5000 + i) for i, J in enumerate(Js)]   # prevalence-neutral
print(f"generated 2x{N_SIMS} reference simulations in {time.perf_counter()-t0:.0f}s "
      f"(S={S} species, N={N} families)")


# ==================================================================================
# Figure 1 — the co-occurrence signal: correlation heatmaps (matched prevalence)
# ==================================================================================
def fig_heatmaps():
    Pc = presence(1.2, -2.4, 1)          # coupled, prevalence ~0.5
    Pi = presence(0.0, 0.2, 1)           # independent, prevalence ~0.5
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.5))
    for ax, P, ttl in [(axes[0], Pc, "a  Coupled (block J > 0)"),
                       (axes[1], Pi, "b  Independent (J = 0)")]:
        keep = P.std(1) > 0
        C = np.corrcoef(P[keep])
        im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
        ax.set_title(ttl)
        ax.set_xlabel("gene family"); ax.set_ylabel("gene family")
        ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label("presence correlation", fontsize=9)
    fig.suptitle("Gene-family co-occurrence (both at prevalence $\\approx$ 0.5)",
                 x=0.09, ha="left", fontweight="bold", fontsize=12)
    save(fig, "fig1_heatmaps")
    results["heatmap"] = dict(coupled_prev=float(Pc.mean()), indep_prev=float(Pi.mean()),
                              coupled_triangles=float(s_clust(Pc)[1]),
                              indep_triangles=float(s_clust(Pi)[1]))


# ==================================================================================
# Figure 2 — ABC recovers coupling strength J (realistic regime, marginal+clustering)
# ==================================================================================
def fig_recover_realistic():
    ref = np.array([np.concatenate([s_freq(P), s_clust(P)]) for P in P_real])
    fitted, los, his = [], [], []
    for tJ in TRUE:
        Pt = presence(tJ, -0.5, 777)
        m, lo, hi = med_ci(abc_posterior(ref, np.concatenate([s_freq(Pt), s_clust(Pt)]), Js))
        fitted.append(m); los.append(lo); his.append(hi)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.plot([0, 1.5], [0, 1.5], color=MUTED, ls=(0, (4, 3)), lw=1.1, label="truth")
    ax.errorbar(TRUE, fitted, yerr=[np.array(fitted) - los, np.array(his) - fitted],
                fmt="o", color=BLUE, ms=7, capsize=3, lw=1.4, label="ABC posterior (median, 95% CI)")
    ax.set_xlabel("true coupling  $J$"); ax.set_ylabel("inferred coupling  $\\hat J$")
    ax.set_title("ABC recovers the coupling strength")
    ax.set_xlim(-0.1, 1.5); ax.set_ylim(-0.1, 1.6); ax.legend(loc="upper left")
    rmse = robust_rmse(ref, lambda P: np.concatenate([s_freq(P), s_clust(P)]),
                       lambda J: -0.5, Js, TRUE)         # stable RMSE over many targets
    ax.text(0.98, 0.05, f"RMSE = {rmse:.3f}", transform=ax.transAxes, ha="right", color=MUTED)
    save(fig, "fig2_recover_realistic")
    results["recover_realistic"] = dict(true=TRUE, fitted=fitted, lo=los, hi=his, rmse=rmse)


# ==================================================================================
# Figure 3 — why clustering: module signal vs J, and RMSE by summary (fixed prevalence)
# ==================================================================================
def fig_summary_comparison():
    # (a) module signal vs J at fixed prevalence (within-block gap + label-blind triangles)
    gaps, tris = [], []
    for J in TRUE:
        g = np.mean([within_between_gap(presence(J, H_NEUTRAL(J), s)) for s in range(1, 6)])
        tr = np.mean([s_clust(presence(J, H_NEUTRAL(J), s))[1] for s in range(1, 6)])
        gaps.append(g); tris.append(tr)

    # (b) recover-J RMSE per summary at fixed prevalence (co-occurrence only; marginal is blind).
    # Averaged over several target datasets per J so the ranking is stable, not a single-draw fluke.
    summ_fns = {"marginal\n(blind)": s_freq, "histogram": s_hist,
                "eigenvalues": s_eig, "clustering": s_clust}
    refs = {k: np.array([f(P) for P in P_neut]) for k, f in summ_fns.items()}
    rmse = {k: robust_rmse(refs[k], f, H_NEUTRAL, Js, TRUE) for k, f in summ_fns.items()}

    fig, (a, b) = plt.subplots(1, 2, figsize=(11.8, 4.7))
    fig.subplots_adjust(wspace=0.62)
    a.plot(TRUE, gaps, "o-", color=BLUE, lw=1.6, ms=6)
    a2 = a.twinx(); a2.spines["top"].set_visible(False)
    a2.plot(TRUE, tris, "s--", color=ORANGE, lw=1.4, ms=5)
    a.set_xlabel("coupling  $J$"); a.set_ylabel("within$-$between correlation", color=BLUE)
    a2.set_ylabel("triangle count", color=ORANGE)
    a.set_title("a  Module signal vs $J$", fontsize=11)
    a.tick_params(axis="y", colors=BLUE); a2.tick_params(axis="y", colors=ORANGE)

    names = list(rmse); vals = [rmse[n] for n in names]
    cols = [MUTED, BLUE, BLUE, BLUE]                  # blind marginal vs the co-occurrence group
    b.bar(range(len(names)), vals, color=cols, width=0.62)
    b.set_xticks(range(len(names))); b.set_xticklabels(names, fontsize=9)
    b.set_ylabel("recovery RMSE")
    b.set_title("b  Co-occurrence beats the blind marginal", fontsize=11)
    for i, v in enumerate(vals):
        b.text(i, v + 0.007, f"{v:.2f}", ha="center", fontsize=9)
    b.set_ylim(0, max(vals) * 1.15)
    save(fig, "fig3_summary_comparison")
    results["fixed_prevalence"] = dict(J=TRUE, within_between_gap=gaps, triangles=tris, rmse=rmse)


# ==================================================================================
# Figure 4 — the payoff: clustering recovers J where the marginal goes blind
# ==================================================================================
def fig_recover_neutral():
    ref_m = np.array([s_freq(P) for P in P_neut])
    ref_c = np.array([s_clust(P) for P in P_neut])
    fm, fc = [], []
    for tJ in TRUE:
        Pt = presence(tJ, H_NEUTRAL(tJ), 777)
        fm.append(med_ci(abc_posterior(ref_m, s_freq(Pt), Js)))
        fc.append(med_ci(abc_posterior(ref_c, s_clust(Pt), Js)))
    fm, fc = np.array(fm), np.array(fc)
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    ax.plot([0, 1.5], [0, 1.5], color=MUTED, ls=(0, (4, 3)), lw=1.1, label="truth")
    ax.errorbar(np.array(TRUE) - 0.015, fm[:, 0], yerr=[fm[:, 0] - fm[:, 1], fm[:, 2] - fm[:, 0]],
                fmt="s", color=MUTED, ms=6, capsize=3, lw=1.2, label="marginal only (blind)")
    ax.errorbar(np.array(TRUE) + 0.015, fc[:, 0], yerr=[fc[:, 0] - fc[:, 1], fc[:, 2] - fc[:, 0]],
                fmt="o", color=BLUE, ms=7, capsize=3, lw=1.4, label="clustering (co-occurrence)")
    ax.set_xlabel("true coupling  $J$"); ax.set_ylabel("inferred coupling  $\\hat J$")
    ax.set_title("Prevalence-neutral: clustering sees what the marginal can't")
    ax.set_xlim(-0.1, 1.5); ax.set_ylim(-0.1, 1.7); ax.legend(loc="upper left")
    save(fig, "fig4_recover_neutral")
    results["recover_neutral"] = dict(
        true=TRUE, marginal=fm.tolist(), clustering=fc.tolist(),
        rmse_marginal=robust_rmse(ref_m, s_freq, H_NEUTRAL, Js, TRUE),
        rmse_clustering=robust_rmse(ref_c, s_clust, H_NEUTRAL, Js, TRUE))


for f in (fig_heatmaps, fig_recover_realistic, fig_summary_comparison, fig_recover_neutral):
    t = time.perf_counter()
    f()
    print(f"  {f.__name__}: {time.perf_counter()-t:.0f}s")

(HERE / "results.json").write_text(json.dumps(results, indent=2))
print("\nwrote figures to", FIG, "and results.json")
print(json.dumps({k: results[k] for k in ("heatmap", "recover_realistic", "fixed_prevalence",
                                          "recover_neutral")}, indent=2, default=float))
