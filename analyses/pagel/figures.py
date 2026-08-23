#!/usr/bin/env python3
"""The one figure: rejection rates by arm, the cox pairing against the control, with
Wilson intervals. Reads results.json beside it, writes figures/pagel.{png,pdf}.

    python figures.py
"""
from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"

INK, MUTED, FAINT = "#1a1a1a", "#8a8a8a", "#c9c9c9"
ARMS = [("feedback", "feedback\n(both connections)"),
        ("gen2trait", "the cox family\ndrives the switch rate"),
        ("trait2gen", "the habitat\ndrives the loss rate"),
        ("null", "null\n(no connection)")]

plt.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.bbox": "tight",
})


def main() -> int:
    d = json.loads((HERE / "results.json").read_text())
    FIG.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    x = np.arange(len(ARMS)) * 1.15
    for off, pair, color, label in ((-0.19, "cox", INK, "habitat and the cox family"),
                                    (0.19, "ctrl", FAINT, "habitat and the control")):
        rates = [d[a][pair]["rate"] for a, _ in ARMS]
        los = [d[a][pair]["rate"] - d[a][pair]["wilson95"][0] for a, _ in ARMS]
        his = [d[a][pair]["wilson95"][1] - d[a][pair]["rate"] for a, _ in ARMS]
        ax.bar(x + off, rates, width=0.36, color=color, label=label,
               yerr=np.array([los, his]), error_kw={"ecolor": MUTED, "capsize": 3, "lw": 1.1})
    ax.axhline(0.05, color=MUTED, lw=0.9, ls=(0, (3, 3)))
    ax.annotate("the nominal 5%", (x[-1] - 0.87, 0.065), color=MUTED, fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in ARMS], fontsize=9.5)
    ax.set_ylabel("share of replicates rejecting independence")
    ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, loc="upper right")
    for ext in ("png", "pdf"):
        out = FIG / f"pagel.{ext}"
        fig.savefig(out, dpi=300)
        print(" ", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
