"""Snakemake step 3: aggregate every per-run metrics.tsv into one summary + figures.

Writes ``summary.tsv`` (one row per run) and two figures:

- ``red_curve.png`` — RED accuracy (Pearson r and nRMSE) vs the *realized* across-lineage rate
  spread (``fold_range``), one line per (tree config × clock family). This is the generalization
  the single-tree prototype cannot show: where RED degrades, across tree sizes and clock models.
- ``red_scatter.png`` — the prototype-style true-vs-recovered scatter for a showcase tree across
  a few perturbation levels (pooled over replicates), for direct visual reproduction.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

FLOAT_COLS = {"strength", "total_age", "pearson_r", "spearman_r", "nrmse", "fold_range",
              "rate_p95_p5", "rate_logsd"}
INT_COLS = {"n_tips", "tseed", "cseed", "seed", "n_internal"}
XKEY = "rate_p95_p5"      # robust across-lineage spread (outlier-safe), the figures' x-axis
XLABEL = "across-lineage rate spread (95th/5th-percentile ratio, ×)"


def _read_metrics(path: str) -> dict:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    row = rows[0]
    for k in list(row):
        if k in FLOAT_COLS:
            row[k] = float(row[k])
        elif k in INT_COLS:
            row[k] = int(row[k])
    row["_dir"] = os.path.dirname(path)
    return row


def write_summary(rows: list, path: str) -> None:
    cols = [c for c in rows[0] if c != "_dir"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])


# --------------------------------------------------------------------------- figures
_STYLE = {"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
          "axes.edgecolor": "#444", "axes.linewidth": 0.8, "figure.dpi": 130}
_COLORS = ["#2f6f8f", "#3f8f5a", "#c98a3a", "#a4443f", "#7a5ea8", "#4c4c4c"]
_LINESTYLE = {"ratevar": "-", "aln": "--", "cir": ":", "whitenoise": "-.", "ucln": (0, (3, 1, 1, 1)),
              "strict": "-"}


def figure_curve(rows: list, path: str) -> None:
    """Pearson r and nRMSE vs across-lineage rate spread — one marker per perturbation level.

    Colour = tree config, line style = clock model; each marker aggregates the replicates of one
    ``(treeset, clock, pert)`` at its mean realized spread (x) with ±sd error bars, and markers of
    the same ``(treeset, clock)`` are joined into a degradation curve.
    """
    from matplotlib.lines import Line2D

    treesets = sorted({r["treeset"] for r in rows}, key=lambda t: rows_lookup(rows, t))
    clocks = sorted({r["clock"] for r in rows})
    color = {t: _COLORS[i % len(_COLORS)] for i, t in enumerate(treesets)}

    agg = defaultdict(list)                       # (treeset, clock, pert) -> runs
    for r in rows:
        agg[(r["treeset"], r["clock"], r["pert"])].append(r)
    lines = defaultdict(list)                     # (treeset, clock) -> [(x, r, r_sd, e, e_sd)]
    for (treeset, clock, _pert), rs in agg.items():
        lines[(treeset, clock)].append((
            float(np.mean([b[XKEY] for b in rs])),
            float(np.mean([b["pearson_r"] for b in rs])), float(np.std([b["pearson_r"] for b in rs])),
            100 * float(np.mean([b["nrmse"] for b in rs])), 100 * float(np.std([b["nrmse"] for b in rs])),
        ))

    with plt.rc_context(_STYLE):
        fig, (axr, axe) = plt.subplots(1, 2, figsize=(12, 5))
        for (treeset, clock), pts in sorted(lines.items()):
            pts.sort()
            xs = [p[0] for p in pts]
            ls = _LINESTYLE.get(clock, "-")
            axr.errorbar(xs, [p[1] for p in pts], yerr=[p[2] for p in pts], color=color[treeset],
                         ls=ls, marker="o", ms=4, lw=1.4, capsize=2)
            axe.errorbar(xs, [p[3] for p in pts], yerr=[p[4] for p in pts], color=color[treeset],
                         ls=ls, marker="o", ms=4, lw=1.4, capsize=2)
        for ax in (axr, axe):
            ax.set_xscale("log")
            ax.set_xlabel(XLABEL)
        axr.set_ylabel("Pearson r (true vs RED-recovered age)")
        axe.set_ylabel("nRMSE (% of root age)")
        axr.set_title("RED age recovery vs rate heterogeneity", loc="left")
        axe.set_title("RED error vs rate heterogeneity", loc="left")
        # two compact legends: colour = tree config, line style = clock model
        tleg = [Line2D([], [], color=color[t], lw=2, label=t) for t in treesets]
        cleg = [Line2D([], [], color="#555", lw=1.6, ls=_LINESTYLE.get(c, "-"), label=c) for c in clocks]
        leg1 = axr.legend(handles=tleg, title="tree", fontsize=8, title_fontsize=8,
                          frameon=False, loc="lower left")
        axr.add_artist(leg1)
        axr.legend(handles=cleg, title="clock", fontsize=8, title_fontsize=8,
                   frameon=False, loc="upper right")
        fig.tight_layout()
        fig.savefig(path, bbox_inches="tight", facecolor="white")
        plt.close(fig)


def rows_lookup(rows, treeset):
    # sort treesets by tip count for a stable legend/colour order
    for r in rows:
        if r["treeset"] == treeset:
            return r["n_tips"]
    return 0


def figure_scatter(rows: list, showcase: dict, path: str) -> None:
    """Prototype-style true-vs-recovered scatter for one tree across perturbation levels."""
    tree = showcase.get("tree")
    perts = showcase.get("perturbations", [])
    pooled = {p: {"true": [], "rec": [], "r": [], "spread": []} for p in perts}
    total_age = None
    for r in rows:
        if r["treeset"] != tree or r["pert"] not in pooled:
            continue
        total_age = r["total_age"]
        pts = _load_points(os.path.join(r["_dir"], "points.csv"))
        pooled[r["pert"]]["true"].extend(pts[0])
        pooled[r["pert"]]["rec"].extend(pts[1])
        pooled[r["pert"]]["r"].append(r["pearson_r"])
        pooled[r["pert"]]["spread"].append(r[XKEY])
    perts = [p for p in perts if pooled[p]["true"]]
    if not perts:
        return
    n = len(perts)
    ncol = 2 if n > 1 else 1
    nrow = (n + ncol - 1) // ncol
    with plt.rc_context(_STYLE):
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.4 * nrow), squeeze=False)
        for ax, p in zip(axes.flat, perts):
            d = pooled[p]
            ax.plot([0, total_age], [0, total_age], ls="--", lw=1, color="#888", zorder=1)
            ax.scatter(d["true"], d["rec"], s=6, alpha=0.12, color="#2f6f8f",
                       edgecolors="none", zorder=2)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(0, total_age * 1.02)
            ax.set_ylim(0, total_age * 1.02)
            ax.set_title(f"{p}  (~{np.mean(d['spread']):.0f}× spread)", loc="left", fontsize=11)
            ax.text(0.04, 0.95, f"r = {np.mean(d['r']):.3f} ± {np.std(d['r']):.3f}",
                    transform=ax.transAxes, va="top", fontsize=9.5,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#ddd", alpha=0.9))
            ax.set_xlabel("true node age")
            ax.set_ylabel("RED-recovered age")
        for ax in axes.flat[n:]:
            ax.axis("off")
        fig.suptitle(f"RED node-age recovery — showcase tree '{tree}' (pooled replicates)", y=1.0)
        fig.tight_layout()
        fig.savefig(path, bbox_inches="tight", facecolor="white")
        plt.close(fig)


def _load_points(path: str):
    true, rec = [], []
    with open(path) as fh:
        r = csv.reader(fh)
        next(r, None)
        for row in r:
            true.append(float(row[1]))
            rec.append(float(row[2]))
    return true, rec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--fig-dir", required=True)
    ap.add_argument("--showcase", default="{}", help="JSON: {tree, perturbations}")
    ap.add_argument("--runs-dir", help="discover metrics.tsv under this dir (avoids an argv-length "
                    "limit on large sweeps); alternative to passing files positionally")
    ap.add_argument("metrics", nargs="*", help="per-run metrics.tsv files (if --runs-dir is unset)")
    args = ap.parse_args()

    import glob
    paths = list(args.metrics)
    if args.runs_dir:
        paths += sorted(glob.glob(os.path.join(args.runs_dir, "**", "metrics.tsv"), recursive=True))
    if not paths:
        ap.error("no metrics files (pass --runs-dir or metrics paths)")
    rows = [_read_metrics(p) for p in paths]
    write_summary(rows, args.out_summary)
    os.makedirs(args.fig_dir, exist_ok=True)
    figure_curve(rows, os.path.join(args.fig_dir, "red_curve.png"))
    figure_scatter(rows, json.loads(args.showcase), os.path.join(args.fig_dir, "red_scatter.png"))

    # console summary + the two invariant self-checks
    print(f"[summarize] {len(rows)} runs -> {args.out_summary}")
    by = defaultdict(list)
    for r in rows:
        by[(r["treeset"], r["pert"])].append(r["pearson_r"])
    for (treeset, pert), rs in sorted(by.items()):
        print(f"  {treeset:14} {pert:14} r={np.mean(rs):.4f} ± {np.std(rs):.4f}  (n={len(rs)})")
    strict = [r for r in rows if r["clock"] == "strict"]
    if strict:
        # Pearson r is affine-invariant, so it alone does not prove exact age *recovery*; also
        # assert calibration (nRMSE ~ 0). Together these are the "RED recovers time exactly" check.
        min_r = min(r["pearson_r"] for r in strict)
        max_e = max(r["nrmse"] for r in strict)
        assert min_r > 0.9999, f"strict-clock r should be ~1.0, got min {min_r}"
        assert max_e < 1e-6, f"strict-clock nRMSE should be ~0 (calibration), got max {max_e}"
        print(f"  [self-check] strict-clock r = {np.mean([r['pearson_r'] for r in strict]):.6f}, "
              f"nRMSE ≤ {max_e:.1e} (≡ RED recovers time exactly)")


if __name__ == "__main__":
    main()
