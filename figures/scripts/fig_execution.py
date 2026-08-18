"""Figure: the three kinds of run — independent, conditioned, joint.

Chapter 8's execution triptych: level 1 and level 2 as ellipses, one panel per kind of run.
A plain arrow is the hierarchy's order; the black circle on the conditioned arrow's tail is
the finished driver the second run reads; the joint panel has one arrow with two heads,
because neither end can be finished first.

Run:  python figures/scripts/fig_execution.py
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent

PANELS = [
    ("independent", "-|>", False, ["level 1 first, level 2 second,", "independent"]),
    ("conditioned", "-|>", True, ["level 1 first, level 2 second,", "level 2 depends on level 1"]),
    ("joint", "<|-|>", False, ["both levels simulated", "simultaneously"]),
]
PW, GAP = 350.0, 40.0
W = 3 * PW + 2 * GAP + 2 * 21
H = 500.0
RX, RY = 100.0, 40.0
FRAME_TOP, FRAME_BOT = 54.0, 384.0


def main() -> None:
    out = ROOT / "manual" / "book" / "figures" / "execution.svg"
    mid = (FRAME_TOP + FRAME_BOT) / 2
    y_top, y_bot = mid - 70, mid + 70
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=190)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_axis_off(); ax.set_aspect("equal")
    for i, (name, arrow, dot, sub) in enumerate(PANELS):
        x0 = 21 + i * (PW + GAP)
        cx = x0 + PW / 2
        ax.add_patch(FancyBboxPatch((x0, FRAME_TOP), PW, FRAME_BOT - FRAME_TOP,
                                    boxstyle="round,pad=0", facecolor="none",
                                    edgecolor="#8a8a8a", linewidth=1.3))
        ax.text(x0 + 26, FRAME_TOP, f" {name} ", fontsize=15, style="italic", color="#555",
                ha="left", va="center", backgroundcolor="white")
        for yc, level in ((y_top, "level 1"), (y_bot, "level 2")):
            ax.add_patch(Ellipse((cx, yc), 2 * RX, 2 * RY, facecolor="#f2f2f0",
                                 edgecolor="#2b2b2b", linewidth=1.6))
            ax.text(cx, yc, level, ha="center", va="center", fontsize=19)
        a0 = y_top + RY + 4
        ax.add_patch(FancyArrowPatch((cx, a0), (cx, y_bot - RY), arrowstyle=arrow,
                                     mutation_scale=24, linewidth=1.8, color="#2b2b2b",
                                     shrinkA=0, shrinkB=4))
        if dot:
            ax.add_patch(Circle((cx, a0 + 2), 7.5, facecolor="#2b2b2b", edgecolor="none"))
        for j, line in enumerate(sub):
            ax.text(cx, 420 + j * 26, line, ha="center", va="center", fontsize=14.5,
                    color="#444")
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
