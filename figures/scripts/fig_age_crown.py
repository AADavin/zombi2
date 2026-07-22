"""Figure: what `age` measures — crown age vs stem (origin) age.

`simulate_species_tree(..., age=A, age_type=...)` fixes the tree's timescale, but the reference
point differs:
  * `age_type="crown"` — A is the depth from the CROWN (the first speciation, i.e. the root of the
    reconstructed tree) to the present.
  * `age_type="stem"` — A is the depth from the ORIGIN (the start of the lineage) to the present, so
    a stem branch precedes the crown and the crown subtree is correspondingly shorter.

Two side-by-side panels show the SAME total age A (same left-to-right span) read against the two
reference points; in the stem panel part of A is spent on the stem before the crown.

House style: B&W, one centered title, ASCII text.

Run:  python figures/scripts/fig_age_crown.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

from zombi_style import save, FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK


W, H = 1160, 380

# a small fixed 4-tip ultrametric tree, node times normalized 0 (crown) .. 1 (present)
TIP_T = 1.0
NODE_U_T, NODE_L_T = 0.46, 0.58        # the two internal nodes below the crown


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

    # crown -> two internal nodes
    connector(x(0.0), yU, yL)
    branch(x(0.0), x(NODE_U_T), yU)
    branch(x(0.0), x(NODE_L_T), yL)
    # internal U -> tips 0,1 ; internal L -> tips 2,3
    connector(x(NODE_U_T), ys[0], ys[1])
    connector(x(NODE_L_T), ys[2], ys[3])
    for k in (0, 1):
        branch(x(NODE_U_T), x(TIP_T), ys[k])
    for k in (2, 3):
        branch(x(NODE_L_T), x(TIP_T), ys[k])
    return yC, ys


def node_dot(d, x, y, label, below=False):
    d.append(draw.Circle(x, y, 5.5, fill=INK))
    d.append(draw.Text(label, FS_TICK, x, y + (26 if below else -14), font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))


def age_bracket(d, x0, x1, y, label):
    d.append(draw.Line(x0, y, x1, y, stroke=INK, stroke_width=1.8))
    for xx in (x0, x1):
        d.append(draw.Line(xx, y - 7, xx, y + 7, stroke=INK, stroke_width=1.8))
    d.append(draw.Text(label, FS_LABEL, (x0 + x1) / 2, y + 26, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))

    ytop, ybot = 108, 300
    ybase = 340                                            # age-bracket baseline
    ylab = 58                                              # panel header baseline
    # --- Panel A: crown ---
    ax0, ax1 = 90, 520                                     # crown ... present (full age)
    d.append(draw.Text("age_type = 'crown'", FS_LABEL, (ax0 + ax1) / 2, ylab, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    yC, _ = draw_tree(d, ax0, ax1, ytop, ybot)
    node_dot(d, ax0, yC, "crown")
    age_bracket(d, ax0, ax1, ybase, "age")

    # --- Panel B: stem ---
    bx0, bx1 = 640, 1070                                   # origin ... present (same total age)
    stem = 0.30 * (bx1 - bx0)                              # part of the age is the stem
    xcrown = bx0 + stem
    d.append(draw.Text("age_type = 'stem'", FS_LABEL, (bx0 + bx1) / 2, ylab, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    yC2, _ = draw_tree(d, xcrown, bx1, ytop, ybot)
    d.append(draw.Line(bx0, yC2, xcrown, yC2, stroke=INK, stroke_width=2.6))   # stem branch
    # origin / stem / crown share one baseline above the branch (consistent letter size)
    yrow = yC2 - 14
    d.append(draw.Circle(bx0, yC2, 5.5, fill=INK))
    d.append(draw.Text("origin", FS_TICK, bx0, yrow, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    d.append(draw.Text("stem", FS_TICK, (bx0 + xcrown) / 2, yrow, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    d.append(draw.Circle(xcrown, yC2, 5.5, fill=INK))
    d.append(draw.Text("crown", FS_TICK, xcrown + 10, yrow, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    age_bracket(d, bx0, bx1, ybase, "age")

    name = "age_crown"
    save(d, name)


if __name__ == "__main__":
    render()
