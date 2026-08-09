"""Figure: a rate is a probability, so the count in a fixed interval is Poisson.

A rate does not deliver a fixed number of events per unit of time. Over any given unit the
number that actually fire is random: with rate lambda, the count follows a Poisson
distribution with mean lambda. Two panels (a low and a higher rate) show how the count is
scattered around the mean, never exactly equal to it -- often zero or one at a low rate.

Run:  python figures/scripts/fig_gillespie_poisson.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

from zombi_style import save, FONT, INK, MUTED, ACCENT, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK


_LIFT = 56          # the band the in-figure title used to take
W, H = 1200, 444
KMAX = 10


def text(d, s, x, y, size, *, anchor="middle", fill=INK, weight="normal", italic=False):
    d.append(draw.Text(s, size, x, y, font_family=FONT, text_anchor=anchor,
                       dominant_baseline="central", fill=fill, font_weight=weight,
                       font_style=("italic" if italic else "normal")))


def poisson_pmf(lam, k):
    return math.exp(-lam) * lam ** k / math.factorial(k)


def panel(d, x0, lam, label):
    px, py = x0 + 40, 120          # plot origin (top-left of plot area)
    pw, ph = 440, 250
    base = py + ph
    step = pw / (KMAX + 1)
    bw = step * 0.64
    pmax = max(poisson_pmf(lam, k) for k in range(KMAX + 1))

    # axes
    d.append(draw.Line(px, base, px + pw, base, stroke=INK, stroke_width=2.0))
    d.append(draw.Line(px, py, px, base, stroke=INK, stroke_width=2.0))

    # bars
    for k in range(KMAX + 1):
        cx = px + (k + 0.5) * step
        h = ph * poisson_pmf(lam, k) / pmax
        d.append(draw.Rectangle(cx - bw / 2, base - h, bw, h, fill=INK, opacity=0.82))
        if k <= 8 or lam > 2:
            text(d, str(k), cx, base + 22, FS_TICK, fill=MUTED)

    # mean marker
    mx = px + (lam + 0.5) * step
    d.append(draw.Line(mx, py - 6, mx, base, stroke=ACCENT["loss"], stroke_width=2.4,
                       stroke_dasharray="5,5"))
    text(d, f"mean = {lam:g}", mx, py - 22, FS_ANNOT, fill=ACCENT["loss"], weight="bold")

    text(d, label, px + pw / 2, py + ph + 58, FS_LABEL, fill=INK)
    text(d, "number of events", px + pw / 2, py + ph + 90, FS_TICK, fill=MUTED)


def render():
    # No title inside the figure (figures/STYLE.md): the manual captions every one, and a title
    # repeats the caption in a second voice a few millimetres above it. `_LIFT` is the band the
    # title used to occupy, taken back off the top so the drawing keeps its own coordinates.
    d = draw.Drawing(W, H + _LIFT, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H + _LIFT, fill="white"))

    panel(d, 40, 1, "rate = 1 per unit of time")
    panel(d, 620, 4, "rate = 4 per unit of time")

    name = "gillespie_poisson"
    lifted = draw.Group(transform="translate(0,-56)")
    for element in d.elements:
        lifted.append(element)
    out = draw.Drawing(W, H, origin=(0, 0))
    out.append(draw.Rectangle(0, 0, W, H, fill="white"))
    out.append(lifted)
    save(out, name)


if __name__ == "__main__":
    render()
