"""Forward map / ABC over (inversion rate × mean inversion length), inversions only.

At each grid cell we simulate the genome down the dated tree (several seeds), measure the synteny
observables on a density-matched core, and compare to the real data:

  * gene-order conservation (uns) + conserved block size (blk)  -> the fit targets
  * orientation flips (flip), reversed-block size (rev)          -> checks / size probe

Per clade it writes an ``.npz`` (grid axes, the distance surface, and every observable surface) and
a ``.json`` (best cell + pseudo-posterior marginals for the rate and the length). Where the distance
is minimal the simulation matches the data; the length runs unconstrained (the rate–length ridge)
while the rate is pinned. Ported from ``synteny/scripts/forward_map.py`` + ``fit_mc.py``.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import pathlib

import numpy as np

import sim
import stats

HERE = pathlib.Path(__file__).parent
REFERENCE_L = 4.0        # inversion length (genes) at which the rate is read off the ridge — a few genes
                         # (Keogh et al. 2000); the length itself is not identifiable (see the report)


def markers(n_total: int, n_core: int) -> set[int]:
    """The ``n_core`` observed markers, spread evenly across all ``n_total`` genes so their density
    matches the real single-copy core (~n_core/n_total) — the block-size statistic is density-sensitive."""
    n_core = min(n_core, n_total)
    idx = (np.arange(n_core) * (n_total / n_core)).astype(int) + 1
    return {int(i) for i in idx}


def _sim_summary(job):
    tree_path, (inv, L), G, n_core, K, bins, divergence, seed = job
    tree, namemap = sim.load_dated_tree(tree_path)
    g = sim.simulate_signed_order(tree, namemap, inversion=float(inv), mean_length=float(L),
                                  n_total=G, n_chrom=K, seed=int(seed))
    mk = markers(G, n_core)
    marked = {sp: [(c, f, o) for (c, f, o) in order if f in mk] for sp, order in g.items()}
    return stats.summary_vector(marked, divergence, bins=bins, core=mk)


def build(tree_path, thetas, G, n_core, K, bins, divergence, seed0, workers):
    jobs = [(tree_path, tuple(t), G, n_core, K, bins, divergence, seed0 + i)
            for i, t in enumerate(thetas)]
    workers = workers or max(1, mp.cpu_count() - 1)
    if workers == 1:
        return [_sim_summary(j) for j in jobs]
    # fork (not macOS-default spawn): workers inherit the imports and never re-run the caller
    with mp.get_context("fork").Pool(workers) as pool:
        return pool.map(_sim_summary, jobs)


def adaptive_scales(sums, keys):
    """Per-key scale = s.d. across the grid, so conservation (0–1) and block size (1–hundreds)
    weigh comparably in the distance."""
    M = np.array([[s.get(k, np.nan) for k in keys] for s in sums], float)
    out = {}
    for j, k in enumerate(keys):
        col = M[:, j][~np.isnan(M[:, j])]
        out[k] = float(np.std(col)) if col.size and np.std(col) > 0 else 1.0
    return out


def sweep(name, tree_path, obs_path, K, G, bins, invs, Ls, n_seed=3, seed0=90001,
          workers=None, out=None):
    real = stats.load_signed(obs_path)
    core = stats.core_single_copy(real, 1.0)
    n_core = len(core)
    divergence = stats.pairwise_divergence(tree_path)
    target = stats.summary_vector(real, divergence, bins=bins, core=core)
    fit_keys = stats.summary_keys(bins)
    all_keys = fit_keys + stats.flip_keys(bins) + stats.rev_keys(bins)

    grid = [(float(inv), float(L)) for L in Ls for inv in invs]          # row-major: L outer, inv inner
    thetas, cell_of = [], []
    for gi, (inv, L) in enumerate(grid):
        for _ in range(n_seed):
            thetas.append((inv, L)); cell_of.append(gi)
    sums = build(tree_path, thetas, G, n_core, K, bins, divergence, seed0, workers)

    nI, nL = len(invs), len(Ls)
    cell_sum = [{} for _ in grid]
    for gi in range(len(grid)):
        cs = [sums[j] for j in range(len(sums)) if cell_of[j] == gi]
        for k in all_keys:
            v = [s.get(k, np.nan) for s in cs]
            v = [x for x in v if not np.isnan(x)]
            cell_sum[gi][k] = float(np.mean(v)) if v else np.nan

    surf = {k: np.full((nL, nI), np.nan) for k in all_keys}
    for gi in range(len(grid)):
        li, ii = gi // nI, gi % nI
        for k in all_keys:
            surf[k][li, ii] = cell_sum[gi][k]

    scales = adaptive_scales(cell_sum + [target], fit_keys)
    w = {k: 1.0 / scales[k] ** 2 for k in fit_keys}
    dist = np.full((nL, nI), np.nan)
    for gi in range(len(grid)):
        li, ii = gi // nI, gi % nI
        dist[li, ii] = stats.distance(cell_sum[gi], target, fit_keys, w)

    save = {"invs": np.asarray(invs), "Ls": np.asarray(Ls), "dist": dist,
            "n_core": n_core, "bins": np.asarray(bins)}
    for k in all_keys:
        save[f"surf_{k}"] = surf[k]
        save[f"tgt_{k}"] = float(target.get(k, np.nan))
    np.savez(out, **save)

    bi = np.unravel_index(np.nanargmin(dist), dist.shape)
    # the headline reading: fix the inversion length at a few genes (literature) and read the rate off
    # the ridge — the length is unidentifiable, so we do not fit it (see the report).
    ref_li = int(np.argmin(np.abs(np.asarray(Ls) - REFERENCE_L)))
    rate_at_ref = float(invs[np.nanargmin(dist[ref_li])])
    eps = np.nanpercentile(dist, 10)
    wgrid = np.exp(-(dist / eps) ** 2); wgrid[np.isnan(wgrid)] = 0
    pinv = wgrid.sum(axis=0); pinv /= pinv.sum()
    cdf = np.cumsum(pinv)
    lo = float(invs[np.searchsorted(cdf, 0.05)]); hi = float(invs[np.searchsorted(cdf, 0.95)])
    med = float(invs[np.searchsorted(cdf, 0.50)])
    pL = wgrid.sum(axis=1); pL /= pL.sum()
    cdfL = np.cumsum(pL)
    Llo = float(Ls[np.searchsorted(cdfL, 0.05)]); Lhi = float(Ls[np.searchsorted(cdfL, 0.95)])
    Lmed = float(Ls[np.searchsorted(cdfL, 0.50)])
    est = {"name": name, "rate_at_reference_length": rate_at_ref, "reference_length": REFERENCE_L,
           "inversion": {"median": med, "q05": lo, "q95": hi},
           "mean_length": {"median": Lmed, "q05": Llo, "q95": Lhi},
           "best_cell": {"inversion": float(invs[bi[1]]), "mean_length": float(Ls[bi[0]]),
                         "dist": float(dist[bi])},
           "n_core": n_core, "K": K, "G": G, "bins": list(bins)}
    json.dump(est, open(str(out).replace(".npz", ".json"), "w"), indent=2)
    print(f"[{name}] n_core={n_core} K={K}  rate@L={REFERENCE_L:.0f}genes={rate_at_ref:.2e}/gene·Myr  "
          f"| marginal median={med:.2e} 90%CI=[{lo:.2e},{hi:.2e}]  L median={Lmed:.0f} (90%CI [{Llo:.0f},{Lhi:.0f}] = unconstrained)  "
          f"dist_min={dist[bi]:.3f}", flush=True)
    return est


CLADES = {
    "lachancea": dict(name="Lachancea", tree="data/lachancea/tree_myr.nwk",
                      obs="data/lachancea/observed_signed_order.json", K=8, bins=(0., 70., 130., 200.)),
    "kluyveromyces": dict(name="Kluyveromyces+Eremothecium", tree="data/kluyveromyces/tree_myr.nwk",
                          obs="data/kluyveromyces/observed_signed_order.json", K=6, bins=(0., 40., 120., 200.)),
}


def main(G=5000, n_seed=3, workers=None):
    invs = np.round(np.logspace(np.log10(3e-5), np.log10(1.5e-3), 24), 10)
    Ls = np.array([1., 2., 4., 8., 16., 32., 64., 128., 256.])          # mean inversion length (genes)
    (HERE / "results").mkdir(exist_ok=True)
    for key, c in CLADES.items():
        sweep(c["name"], str(HERE / c["tree"]), str(HERE / c["obs"]), c["K"], G, c["bins"],
              invs, Ls, n_seed=n_seed, workers=workers,
              out=str(HERE / "results" / f"{key}_forward_map.npz"))


if __name__ == "__main__":
    main()
