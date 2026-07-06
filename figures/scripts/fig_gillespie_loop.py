"""Figure: the Gillespie algorithm as a loop.

A flowchart of the exact ("direct method") algorithm ZOMBI2 runs: from the current
state, compute every event's rate and their total R; draw a waiting time from an
Exponential(R); stop if the clock has passed the target age, otherwise choose one
event proportional to its rate, apply it, and go round again.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_gillespie_loop.py
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

W, H = 1000, 900
CX = 430                 # centre x of the main column
BW, BH = 380, 74         # box width / height


def text(d, s, x, y, size, *, anchor="middle", fill=INK, weight="normal"):
    d.append(draw.Text(s, size, x, y, font_family=FONT, text_anchor=anchor,
                       dominant_baseline="central", fill=fill, font_weight=weight))


def lines(d, rows, x, y, size, **kw):
    n = len(rows)
    for i, row in enumerate(rows):
        text(d, row, x, y + (i - (n - 1) / 2) * (size + 8), size, **kw)


def arrowhead(d, hx, hy, ax, ay, col=INK, ah=13.0):
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy,
                        hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


def varrow(d, x, y1, y2, col=INK, lw=2.6):
    d.append(draw.Line(x, y1, x, y2 - 3, stroke=col, stroke_width=lw))
    arrowhead(d, x, y2, x, y1, col)


def box(d, cy, rows, *, fill="white", stroke=INK, rx=14, w=BW, h=BH):
    d.append(draw.Rectangle(CX - w / 2, cy - h / 2, w, h, rx=rx, ry=rx, fill=fill,
                            stroke=stroke, stroke_width=2.6))
    lines(d, rows, CX, cy, FS_ANNOT)


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    text(d, "The Gillespie loop", CX, 40, FS_TITLE, weight="bold")

    y_start = 118
    y_rate = 250
    y_wait = 382
    y_dec = 520          # decision diamond centre
    y_choose = 662
    y_apply = 786

    # 1. start (stadium)
    box(d, y_start, ["start:  t = 0,  initial state"], fill="#f2f2f2", rx=37, h=66)
    varrow(d, CX, y_start + 33, y_rate - BH / 2)

    # 2. compute rates
    box(d, y_rate, ["compute every event's rate", "total rate  R = sum of the rates"])
    varrow(d, CX, y_rate + BH / 2, y_wait - BH / 2)

    # 3. draw waiting time
    box(d, y_wait, ["draw waiting time  dt ~ Exponential(R)", "advance the clock  t = t + dt"])
    varrow(d, CX, y_wait + BH / 2, y_dec - 66)

    # 4. decision diamond
    dw, dh = 320, 130
    d.append(draw.Lines(CX, y_dec - dh / 2, CX + dw / 2, y_dec, CX, y_dec + dh / 2,
                        CX - dw / 2, y_dec, close=True, fill="white", stroke=INK,
                        stroke_width=2.6))
    lines(d, ["clock past the target age?", "(or nothing left alive)"], CX, y_dec, 16)
    # yes -> stop box on the right
    stop_x = CX + dw / 2 + 250
    d.append(draw.Line(CX + dw / 2, y_dec, stop_x - 118, y_dec, stroke=INK, stroke_width=2.6))
    arrowhead(d, stop_x - 118, y_dec, CX + dw / 2, y_dec, INK)
    text(d, "yes", CX + dw / 2 + 66, y_dec - 18, FS_TICK, fill=MUTED, weight="bold")
    d.append(draw.Rectangle(stop_x - 112, y_dec - 40, 224, 80, rx=14, ry=14, fill="#f2f2f2",
                            stroke=ACCENT["loss"], stroke_width=2.8))
    lines(d, ["stop:", "return the history"], stop_x, y_dec, FS_ANNOT, weight="bold")
    # no -> continue down
    varrow(d, CX, y_dec + dh / 2, y_choose - BH / 2)
    text(d, "no", CX + 24, y_dec + dh / 2 + 22, FS_TICK, fill=MUTED, weight="bold")

    # 5. choose the event
    box(d, y_choose, ["choose one event, i,",
                      "with probability  (rate of i) / R"])
    varrow(d, CX, y_choose + BH / 2, y_apply - BH / 2)

    # 6. apply the event
    box(d, y_apply, ["apply event i to the state,", "and record it in the history"])

    # loop-back arrow: from the apply box, down the left margin, up to the rate box
    lx = CX - BW / 2 - 70
    d.append(draw.Line(CX - BW / 2, y_apply, lx, y_apply, stroke=INK, stroke_width=2.6))
    d.append(draw.Line(lx, y_apply, lx, y_rate, stroke=INK, stroke_width=2.6))
    d.append(draw.Line(lx, y_rate, CX - BW / 2 - 3, y_rate, stroke=INK, stroke_width=2.6))
    arrowhead(d, CX - BW / 2, y_rate, lx, y_rate, INK)
    text(d, "repeat", lx - 16, (y_rate + y_apply) / 2, FS_ANNOT, anchor="end",
         fill=MUTED, weight="bold")

    name = "gillespie_loop"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
