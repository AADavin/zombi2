"""Figure: the stem — the founding lineage's own stretch of time, before the first split.

A forward run starts from **one** lineage at time 0 and that lineage lives for a while before it
speciates. The figure names the three points that stretch sits between: the **origin** at time 0, the
**crown** at the first split, and the present at the tips. Everything between origin and crown is the
**stem**, and it is ordinary simulated time — genes are gained and lost along it, traits drift along
it — which is why every tree ZOMBI2 writes gives its root a branch length.

Grown out of the old `fig_age_crown.py`, which compared two settings of an `age_type=` argument the
API no longer has. The drawing is the same; what it names is not.

House style: B&W, ASCII text. No title inside the figure — the manual captions it.

Run:  python figures/scripts/fig_stem.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

from zombi_style import save, FONT, INK, MUTED, FS_LABEL, FS_ANNOT, FS_TICK


W, H = 760, 355

# a small fixed 4-tip ultrametric tree, node times normalised 0 (crown) .. 1 (present)
TIP_T = 1.0
NODE_U_T, NODE_L_T = 0.46, 0.58        # the two internal nodes below the crown
STEM_FRACTION = 0.30                   # how much of the run is over before the crown


def draw_tree(d, xc, xp, ytop, ybot):
    """Draw the 4-tip ultrametric tree from crown x=xc to present x=xp, spanning ytop..ybot."""
    x = lambda t: xc + t * (xp - xc)                       # noqa: E731
    ys = [ytop + (ybot - ytop) * k / 3 for k in range(4)]  # 4 tip rows
    yU, yL = (ys[0] + ys[1]) / 2, (ys[2] + ys[3]) / 2
    yC = (yU + yL) / 2
    lw = 2.6

    def branch(x0, x1, y):
        d.append(draw.Line(x0, y, x1, y, stroke=INK, stroke_width=lw))

    def connector(xx, ya, yb):
        d.append(draw.Line(xx, ya, xx, yb, stroke=INK, stroke_width=lw))

    connector(x(0.0), yU, yL)                              # crown -> two internal nodes
    branch(x(0.0), x(NODE_U_T), yU)
    branch(x(0.0), x(NODE_L_T), yL)
    connector(x(NODE_U_T), ys[0], ys[1])                   # internal U -> tips 0,1
    connector(x(NODE_L_T), ys[2], ys[3])                   # internal L -> tips 2,3
    for k in (0, 1):
        branch(x(NODE_U_T), x(TIP_T), ys[k])
    for k in (2, 3):
        branch(x(NODE_L_T), x(TIP_T), ys[k])
    return yC, ys


def span(d, x0, x1, y, label, *, italic=False):
    """A measured span with end ticks, labelled below."""
    d.append(draw.Line(x0, y, x1, y, stroke=INK, stroke_width=1.8))
    for xx in (x0, x1):
        d.append(draw.Line(xx, y - 7, xx, y + 7, stroke=INK, stroke_width=1.8))
    d.append(draw.Text(label, FS_LABEL, (x0 + x1) / 2, y + 26, font_family=FONT,
                       text_anchor="middle", fill=INK,
                       font_weight="normal" if italic else "bold",
                       font_style="italic" if italic else "normal"))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))

    ytop, ybot = 48, 208
    x0, x1 = 90, 680                                       # origin ... present
    xcrown = x0 + STEM_FRACTION * (x1 - x0)

    yC, _ = draw_tree(d, xcrown, x1, ytop, ybot)
    d.append(draw.Line(x0, yC, xcrown, yC, stroke=INK, stroke_width=2.6))   # the stem branch

    # origin / crown share one baseline above the branch, so the letters sit at one height
    yrow = yC - 16
    d.append(draw.Circle(x0, yC, 5.5, fill=INK))
    d.append(draw.Text("origin", FS_TICK, x0, yrow, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    d.append(draw.Circle(xcrown, yC, 5.5, fill=INK))
    d.append(draw.Text("crown", FS_TICK, xcrown + 12, yrow, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("first split", FS_ANNOT, xcrown + 12, yC + 30, font_family=FONT,
                       text_anchor="start", fill=MUTED, font_style="italic"))

    span(d, x0, xcrown, 250, "stem")
    span(d, x0, x1, 306, "time (origin to present)")

    save(d, "stem")


if __name__ == "__main__":
    render()
