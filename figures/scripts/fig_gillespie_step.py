"""Figure: the anatomy of one Gillespie step.

Three panels, left to right, walk through a single iteration of the algorithm:

  1. list every event the current state can undergo, each with its rate; sum them to
     the total rate R;
  2. draw the waiting time to the next event -- an Exponential(R) draw, mean 1/R;
  3. choose which event fires, each with probability (its rate) / R.

House style: B&W plus the ZOMBI Set1 event accents (duplication blue, transfer green,
loss red) so the same events keep the same colour across the manual.

Run:  python figures/scripts/fig_gillespie_step.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

from zombi_style import save, FONT, INK, MUTED, ACCENT, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK


W, H = 1260, 512

# The three illustrative events and their (aggregate) rates. Round numbers so the
# probabilities are clean: R = 6, shares 1/2, 1/3, 1/6.
EVENTS = [
    ("duplication", 3.0, ACCENT["duplication"]),
    ("transfer",    2.0, ACCENT["transfer"]),
    ("loss",        1.0, ACCENT["loss"]),
]
R = sum(r for _, r, _ in EVENTS)


def text(d, s, x, y, size, *, anchor="middle", fill=INK, weight="normal", italic=False):
    d.append(draw.Text(s, size, x, y, font_family=FONT, text_anchor=anchor,
                       dominant_baseline="central", fill=fill, font_weight=weight,
                       font_style=("italic" if italic else "normal")))


def arrowhead(d, hx, hy, ax, ay, col, ah=12.0):
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy,
                        hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


def flow_arrow(d, x1, y, x2, col=MUTED, lw=2.6):
    d.append(draw.Line(x1, y, x2 - 3, y, stroke=col, stroke_width=lw, stroke_linecap="round"))
    arrowhead(d, x2, y, x1, y, col, ah=13.0)


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    text(d, "The anatomy of one Gillespie step", W / 2, 40, FS_TITLE, weight="bold")

    top = 120                      # top of the panel content
    band_y = (top + H - 60) / 2    # vertical centre line for the connecting arrows

    # ---- Panel 1: events and their rates -------------------------------------
    ax0 = 60
    text(d, "1.  the events and their rates", ax0, top, FS_LABEL, anchor="start",
         weight="bold")
    bar_x = ax0 + 150              # where the rate bars start
    scale = 34.0                   # px per unit rate
    row_y = top + 58
    for name, r, col in EVENTS:
        d.append(draw.Rectangle(ax0, row_y - 9, 18, 18, fill=col))
        text(d, name, ax0 + 30, row_y, FS_ANNOT, anchor="start")
        d.append(draw.Rectangle(bar_x, row_y - 11, r * scale, 22, fill=col, opacity=0.85))
        text(d, f"{r:.0f}", bar_x + r * scale + 16, row_y, FS_ANNOT, anchor="start",
             fill=MUTED)
        row_y += 56
    # total rate
    ty = row_y + 8
    d.append(draw.Line(ax0, ty - 20, ax0 + 340, ty - 20, stroke=INK, stroke_width=1.4))
    text(d, "total rate  R = 3 + 2 + 1 = 6", ax0, ty + 8, FS_LABEL, anchor="start",
         weight="bold")

    flow_arrow(d, ax0 + 360, band_y, ax0 + 410, col=MUTED)

    # ---- Panel 2: the waiting time (exponential) -----------------------------
    bx0 = 500
    text(d, "2.  when?  draw a waiting time", bx0, top, FS_LABEL, anchor="start",
         weight="bold")
    px, py = bx0 + 8, top + 60
    pw, ph = 300, 236
    d.append(draw.Line(px, py + ph, px + pw, py + ph, stroke=INK, stroke_width=2.0))   # x axis
    d.append(draw.Line(px, py, px, py + ph, stroke=INK, stroke_width=2.0))             # y axis
    mean = 1.0 / R
    xmax = 4.5 * mean
    # exponential density f(x) = R exp(-R x), normalised so the peak (R, at x=0) is the box top
    pts = []
    N = 90
    for i in range(N + 1):
        x = xmax * i / N
        f = R * math.exp(-R * x)
        sx = px + (x / xmax) * pw
        sy = py + ph - (f / R) * ph
        pts.append((sx, sy))
    path = draw.Path(fill="none", stroke=INK, stroke_width=3.0, stroke_linejoin="round")
    path.M(*pts[0])
    for sx, sy in pts[1:]:
        path.L(sx, sy)
    d.append(path)
    # mean marker (dashed vertical at x = 1/R)
    mx = px + (mean / xmax) * pw
    d.append(draw.Line(mx, py + 8, mx, py + ph, stroke=MUTED, stroke_width=2.0,
                       stroke_dasharray="5,5"))
    text(d, "mean = 1/R", mx + 8, py + 24, FS_TICK, anchor="start", fill=MUTED)
    # a sampled draw
    xs = 0.62 * mean
    sxs = px + (xs / xmax) * pw
    fys = py + ph - (R * math.exp(-R * xs) / R) * ph
    d.append(draw.Line(sxs, py + ph, sxs, fys, stroke=ACCENT["origination"], stroke_width=3.0))
    d.append(draw.Circle(sxs, fys, 6, fill=ACCENT["origination"]))
    text(d, "dt", sxs - 4, py + ph + 26, FS_ANNOT, anchor="middle",
         fill=ACCENT["origination"], weight="bold")
    text(d, "waiting time", px + pw / 2, py + ph + 52, FS_TICK, fill=MUTED)
    text(d, "dt ~ Exponential(R)", px + pw / 2, py - 22, FS_ANNOT, weight="bold")

    flow_arrow(d, bx0 + 330, band_y, bx0 + 380, col=MUTED)

    # ---- Panel 3: which event fires ------------------------------------------
    cx0 = 900
    text(d, "3.  which?  pick an event", cx0, top, FS_LABEL, anchor="start", weight="bold")
    # a single bar of length R, split proportionally: the roulette wheel, laid flat
    bar_left = cx0 + 8
    bar_w = 300
    bar_top = top + 78
    bar_h = 54
    x = bar_left
    seg_centers = {}
    for name, r, col in EVENTS:
        w = (r / R) * bar_w
        d.append(draw.Rectangle(x, bar_top, w, bar_h, fill=col, opacity=0.9,
                                stroke="white", stroke_width=2.0))
        frac = {"duplication": "1/2", "transfer": "1/3", "loss": "1/6"}[name]
        text(d, frac, x + w / 2, bar_top + bar_h / 2, FS_ANNOT, fill="white", weight="bold")
        seg_centers[name] = x + w / 2
        x += w
    text(d, "0", bar_left, bar_top + bar_h + 22, FS_TICK, fill=MUTED)
    text(d, "R", bar_left + bar_w, bar_top + bar_h + 22, FS_TICK, fill=MUTED)
    # spinner: a dart pointing into the transfer segment
    spin_x = seg_centers["transfer"]
    d.append(draw.Line(spin_x, bar_top - 34, spin_x, bar_top - 2, stroke=INK, stroke_width=2.6))
    arrowhead(d, spin_x, bar_top - 2, spin_x, bar_top - 34, INK, ah=14.0)
    text(d, "u x R", spin_x, bar_top - 48, FS_TICK, fill=MUTED)
    text(d, "chosen: transfer", bar_left + bar_w / 2, bar_top + bar_h + 62, FS_ANNOT,
         weight="bold", fill=ACCENT["transfer"])
    text(d, "each event fires with", bar_left + bar_w / 2, bar_top + bar_h + 102, FS_TICK,
         fill=MUTED)
    text(d, "probability (its rate) / R", bar_left + bar_w / 2, bar_top + bar_h + 126,
         FS_TICK, fill=MUTED)

    name = "gillespie_step"
    save(d, name)


if __name__ == "__main__":
    render()
