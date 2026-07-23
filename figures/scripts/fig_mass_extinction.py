"""Figure: the mass-extinction model — a tree-wide survival pulse.

Chapter 3 introduces ``mass_extinctions=[(time, fraction), ...]`` as the one model that does not
fit the modifier framework: at a chosen instant a whole cohort of standing lineages dies at once,
a pulse rather than a steady rate. This figure grows the COMPLETE tree with the current engine —
survivors solid, extinct lineages dashed, exactly as ``fig_species_tree_extinct`` draws them — and
marks the pulse as a vertical wall where a cohort of lineages terminates together. An aligned
lineages-through-time curve below, sharing the tree's time axis, shows the diversity crash at the
wall and the recovery after it.

The tree is **simulated here**, by the engine the chapter documents, so the figure cannot drift
from what the code does. The pulse victims are read off the run: the extinct leaves whose death
time lands on the pulse. Bounding by ``total_time`` (not tip count) keeps the pulse at a fixed
place on the axis.

House style: B&W (the LTT curve is black, not an accent), ASCII text, no title inside the figure —
the manual captions it.

Run:  python figures/scripts/fig_mass_extinction.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

from zombi2 import species
from zombi_style import save, FONT, INK, MUTED, FS_LABEL, FS_ANNOT, FS_TICK

# One birth-death run with a single hard pulse. Low death keeps the dashed extinct lineages a
# garnish before the wall; the seed is one where enough survive the 75%-kill to leave a full tree.
BIRTH, DEATH = 1.1, 0.15
PULSE_TIME, KILL = 2.5, 0.75      # kill three-quarters of the standing lineages at time 2.5
TOTAL_TIME, SEED = 5.0, 2

# Canvas: a wide landscape, tree on top and the LTT strip below it, sharing one time axis.
W, H = 1400, 800
XL, XR = 108, 1330                # the shared plotting band, left to right
TREE_TOP, TREE_H = 150, 396
LTT_TOP, LTT_H = 588, 132
DASH = "6,5"                      # extinct-lineage dashes, matching fig_species_tree_extinct
BAND = "#e9e9e9"                  # faint fill marking the instant of the pulse


def _survivors(nodes: dict) -> dict:
    """Map node id -> True if the node or any descendant is a present-day survivor (fate 'extant')."""
    surv: dict[int, bool] = {}

    def rec(i: int) -> bool:
        n = nodes[i]
        # list-then-any, never `any(rec(c) ...)`: the generator short-circuits and would leave
        # a node's later children unvisited (so their id is missing from `surv`).
        s = (n.fate == "extant") if not n.children else any([rec(c) for c in n.children])
        surv[i] = s
        return s

    root = next(i for i, n in nodes.items() if n.parent is None)
    rec(root)
    return surv


def _order(nodes: dict) -> dict:
    """Ladderised leaf order: each node gets a y-rank, internal nodes the mean of their children."""
    y: dict[int, float] = {}
    counter = [0]

    def size(i: int) -> int:
        n = nodes[i]
        return 1 if not n.children else sum(size(c) for c in n.children)

    def rec(i: int) -> None:
        n = nodes[i]
        if not n.children:
            y[i] = counter[0]
            counter[0] += 1
        else:
            for c in sorted(n.children, key=size):
                rec(c)
            y[i] = sum(y[c] for c in n.children) / len(n.children)

    root = next(i for i, n in nodes.items() if n.parent is None)
    rec(root)
    return y, counter[0]


def render() -> None:
    result = species.simulate_species_tree(birth=BIRTH, death=DEATH,
                                            mass_extinctions=[(PULSE_TIME, KILL)],
                                            total_time=TOTAL_TIME, seed=SEED)
    nodes = result.complete_tree.nodes
    present = max(n.end_time for n in nodes.values())
    surv = _survivors(nodes)
    yrank, nleaf = _order(nodes)

    dy = min(12.0, TREE_H / max(1, nleaf - 1))
    top = TREE_TOP + (TREE_H - dy * (nleaf - 1)) / 2

    def X(t: float) -> float:
        return XL + t / present * (XR - XL)

    def Y(i: int) -> float:
        return top + yrank[i] * dy

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))

    # present reference and the pulse wall, both spanning the tree and the LTT strip
    d.append(draw.Line(X(present), TREE_TOP - 18, X(present), LTT_TOP + LTT_H,
                       stroke="#cccccc", stroke_width=1.0, stroke_dasharray="2,4"))
    pw = X(PULSE_TIME)
    d.append(draw.Rectangle(pw - 7, TREE_TOP - 18, 14, (LTT_TOP + LTT_H) - (TREE_TOP - 18), fill=BAND))
    d.append(draw.Line(pw, TREE_TOP - 18, pw, LTT_TOP + LTT_H, stroke=INK, stroke_width=1.6,
                       stroke_dasharray="5,4"))

    # the tree, drawn by hand so each segment is solid (survivor) or dashed (no survivor)
    def seg(x1, y1, x2, y2, dashed):
        if dashed:
            d.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=2.0,
                               stroke_dasharray=DASH, stroke_linecap="butt"))
        else:
            d.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=2.6, stroke_linecap="round"))

    victims = []
    for i, n in nodes.items():
        px = X(0) - 14 if n.parent is None else X(nodes[n.parent].end_time)
        seg(px, Y(i), X(n.end_time), Y(i), not surv[i])
        for c in (n.children or ()):
            seg(X(n.end_time), Y(i), X(n.end_time), Y(c), not surv[c])
        if not n.children and n.fate == "extinct" and abs(n.end_time - PULSE_TIME) < 1e-6:
            victims.append(i)
    for i in victims:                                    # the cohort the pulse killed
        d.append(draw.Circle(X(PULSE_TIME), Y(i), 3.2, fill=INK))

    # legend: solid = surviving lineage, dashed = extinct lineage (top-left, clear of the crown)
    lx, ly = 44, 96
    d.append(draw.Line(lx, ly, lx + 34, ly, stroke=INK, stroke_width=2.6, stroke_linecap="round"))
    d.append(draw.Text("surviving lineage", FS_LABEL, lx + 46, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx + 260, ly, lx + 294, ly, stroke=INK, stroke_width=2.0,
                       stroke_dasharray=DASH, stroke_linecap="butt"))
    d.append(draw.Text("extinct lineage", FS_LABEL, lx + 306, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    # --- lineages-through-time, on the same X axis ---
    grid = [present * k / 700 for k in range(701)]

    def alive(t):
        return sum(1 for n in nodes.values() if n.birth_time < t <= n.end_time)

    counts = [alive(t) for t in grid]
    cmax = max(counts)

    def LY(c):
        return LTT_TOP + LTT_H - c / cmax * LTT_H

    d.append(draw.Line(XL, LTT_TOP + LTT_H, XR, LTT_TOP + LTT_H, stroke="#bdbdbd", stroke_width=1.2))
    d.append(draw.Line(XL, LTT_TOP, XL, LTT_TOP + LTT_H, stroke="#bdbdbd", stroke_width=1.2))
    pts = []
    for t, c in zip(grid, counts):
        pts += [X(t), LY(c)]
    d.append(draw.Lines(*pts, close=False, fill="none", stroke=INK, stroke_width=2.8,
                        stroke_linejoin="round"))
    for c in (0, cmax):
        d.append(draw.Text(str(c), FS_TICK, XL - 12, LY(c), font_family=FONT, text_anchor="end",
                           dominant_baseline="central", fill=MUTED))
    d.append(draw.Text("lineages", FS_LABEL, XL - 46, LTT_TOP + LTT_H / 2, font_family=FONT,
                       text_anchor="middle", fill=MUTED,
                       transform=f"rotate(-90 {XL - 46} {LTT_TOP + LTT_H / 2})"))
    d.append(draw.Text("mass extinction", FS_ANNOT, pw + 16, LTT_TOP + 24, font_family=FONT,
                       text_anchor="start", fill=MUTED, font_style="italic"))

    # shared time axis, under the LTT strip
    ya = LTT_TOP + LTT_H
    for k in range(6):
        t = present * k / 5
        d.append(draw.Line(X(t), ya, X(t), ya + 6, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(f"{t:.0f}" if t == int(t) else f"{t:.1f}", FS_TICK, X(t), ya + 26,
                           font_family=FONT, text_anchor="middle", fill=MUTED))
    d.append(draw.Text("time (origin to present)", FS_LABEL, (XL + XR) / 2, ya + 56, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    save(d, "mass_extinction")
    n_ext = sum(1 for n in nodes.values() if n.fate == "extant")
    print(f"  ({nleaf} leaves, {len(victims)} pulse victims, {n_ext} survivors)")


if __name__ == "__main__":
    render()
