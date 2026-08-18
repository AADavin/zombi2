"""Figure: a connection in the abstract — Driver, Link, Target, no example attached.

Chapter 8's first figure: the two boxes carry only the part's name and the level it sits
at; the arrow between them is the link.

Run:  python figures/scripts/fig_connection_abstract.py
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent

W, H = 1150.0, 300.0
BW, BH = 300.0, 130.0


def main() -> None:
    out = ROOT / "manual" / "book" / "figures" / "connection_abstract.svg"
    y = 150.0
    xl, xr = 60.0, W - 60.0 - BW
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=190)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_axis_off(); ax.set_aspect("equal")
    for x, level, word in ((xl, "level 1", "Driver"), (xr, "level 2", "Target")):
        ax.add_patch(FancyBboxPatch((x, y - BH / 2), BW, BH, boxstyle="round,pad=12",
                                    facecolor="#f2f2f0", edgecolor="#2b2b2b", linewidth=1.6))
        ax.text(x + BW / 2, y - 30, level, ha="center", va="center", fontsize=14,
                style="italic", color="#666")
        ax.text(x + BW / 2, y + 12, word, ha="center", va="center", fontsize=27)
    a0, a1 = xl + BW + 24, xr - 24
    ax.add_patch(FancyArrowPatch((a0, y), (a1, y), arrowstyle="-|>", mutation_scale=28,
                                 linewidth=2.0, color="#2b2b2b"))
    ax.text((a0 + a1) / 2, y - 26, "Link", ha="center", va="center", fontsize=22)
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
