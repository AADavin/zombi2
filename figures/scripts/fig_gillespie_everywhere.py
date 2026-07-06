"""Figure: it's all Gillespie -- the same engine at every level of ZOMBI2.

Each level of a ZOMBI2 simulation supplies a different *bag of events* with their own
rates, but they are all realised by the one loop on the right. Species trees race
speciation against extinction; gene families race duplication, transfer, loss and
origination; a discrete trait races its between-state transitions.

Event markers follow the manual convention: filled square = duplication, black arc =
transfer, cross = loss / extinction, ring = origination. The speciation glyph is a lineage
splitting into two forward in time (one ancestor above, two descendants below).

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_gillespie_everywhere.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, ACCENT, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1240, 516


def text(d, s, x, y, size, *, anchor="middle", fill=INK, weight="normal"):
    d.append(draw.Text(s, size, x, y, font_family=FONT, text_anchor=anchor,
                       dominant_baseline="central", fill=fill, font_weight=weight))


def arrowhead(d, hx, hy, ax, ay, col=INK, ah=13.0):
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy,
                        hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


# ---- event glyphs (small, centred on (x, y)) --------------------------------
def g_speciation(d, x, y):
    # one ancestor above splitting into two descendants below (forward in time)
    c = ACCENT["speciation"]
    d.append(draw.Line(x, y - 13, x, y, stroke=c, stroke_width=3.0))
    d.append(draw.Line(x, y, x - 11, y + 12, stroke=c, stroke_width=3.0))
    d.append(draw.Line(x, y, x + 11, y + 12, stroke=c, stroke_width=3.0))


def g_cross(d, x, y, col=INK):
    r = 9
    d.append(draw.Line(x - r, y - r, x + r, y + r, stroke=col, stroke_width=3.0))
    d.append(draw.Line(x - r, y + r, x + r, y - r, stroke=col, stroke_width=3.0))


def g_square(d, x, y):
    s = 16
    d.append(draw.Rectangle(x - s / 2, y - s / 2, s, s, fill=INK))


def g_transfer(d, x, y):
    p = draw.Path(fill="none", stroke=INK, stroke_width=3.0)
    p.M(x - 12, y + 8).Q(x, y - 20, x + 12, y + 8)
    d.append(p)
    arrowhead(d, x + 12, y + 8, x + 4, y - 4, INK, ah=9.0)


def g_ring(d, x, y):
    d.append(draw.Circle(x, y, 8, fill="none", stroke=INK, stroke_width=3.0))


def g_state(d, x, y):
    d.append(draw.Rectangle(x - 22, y - 8, 15, 15, fill="#dddddd", stroke=INK, stroke_width=1.6))
    d.append(draw.Rectangle(x + 8, y - 8, 15, 15, fill=INK))
    d.append(draw.Line(x - 5, y, x + 6, y, stroke=INK, stroke_width=2.4))
    arrowhead(d, x + 6, y, x - 3, y, INK, ah=7.0)


# (title, rate line, [ (glyph, label), ... ])
LEVELS = [
    ("Species tree", "rates: lambda, mu",
     [(g_speciation, "speciation"), (g_cross, "extinction")]),
    ("Gene families", "rates: one per event",
     [(g_square, "duplication"), (g_transfer, "transfer"), (g_cross, "loss"),
      (g_ring, "origination")]),
    ("Discrete trait (Mk)", "rates: Q-matrix entries",
     [(g_state, "state change")]),
]

CARD_X, CARD_W, CARD_H = 50, 700, 104
YS = [143, 283, 423]
GLYPH_CENTER = 526          # centre of the event-glyph group (right of the text block)
GSTEP = 100                 # spacing between glyphs; the 4-wide row must stay inside the card


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    text(d, "One engine, many events: it's all Gillespie", W / 2, 42, FS_TITLE, weight="bold")

    for (title, rate, events), cy in zip(LEVELS, YS):
        d.append(draw.Rectangle(CARD_X, cy - CARD_H / 2, CARD_W, CARD_H, rx=16, ry=16,
                                fill="white", stroke=INK, stroke_width=2.4))
        text(d, title, CARD_X + 28, cy - 22, FS_LABEL, anchor="start", weight="bold")
        text(d, rate, CARD_X + 28, cy + 8, FS_TICK, anchor="start", fill=MUTED)
        # event glyphs, centred as a group and kept clear of the text block and the border
        n = len(events)
        start_x = GLYPH_CENTER - (n - 1) * GSTEP / 2
        for i, (glyph, label) in enumerate(events):
            gx = start_x + i * GSTEP
            glyph(d, gx, cy - 4)
            text(d, label, gx, cy + 22, FS_TICK)

    # ---- the engine box on the right ----------------------------------------
    ex0, ew = CARD_X + CARD_W + 64, 374
    ey0 = YS[0] - CARD_H / 2
    eh = (YS[-1] + CARD_H / 2) - ey0
    d.append(draw.Rectangle(ex0, ey0, ew, eh, rx=18, ry=18, fill="#f2f2f2", stroke=INK,
                            stroke_width=3.0))
    text(d, "one Gillespie loop", ex0 + ew / 2, ey0 + 40, FS_LABEL, weight="bold")
    step_rows = [
        "1.  total rate R = sum",
        "     of every event's rate",
        "",
        "2.  when:  draw dt ~",
        "     Exponential(R)",
        "",
        "3.  what:  pick event i",
        "     with prob. (rate)/R",
        "",
        "4.  apply it, then repeat",
    ]
    ry = ey0 + 80
    for row in step_rows:
        bold = row[:2] in ("1.", "2.", "3.", "4.")
        text(d, row, ex0 + 32, ry, FS_ANNOT, anchor="start", fill=(INK if bold else MUTED))
        ry += 28

    # arrows from each card into the engine box
    for cy in YS:
        d.append(draw.Line(CARD_X + CARD_W, cy, ex0 - 3, ey0 + eh / 2, stroke=MUTED,
                           stroke_width=2.4))
        arrowhead(d, ex0 - 3, ey0 + eh / 2, CARD_X + CARD_W, cy, MUTED, ah=12.0)

    name = "gillespie_everywhere"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
