"""Figure: a rate is a probability, so the count in a fixed interval is Poisson.

A rate does not deliver a fixed number of events per unit of time. Over any given unit the
number that actually fire is random: with rate lambda, the count follows a Poisson
distribution with mean lambda. Two panels (a low and a higher rate) show how the count is
scattered around the mean, never exactly equal to it -- often zero or one at a low rate.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_gillespie_poisson.py
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

W, H = 1200, 500
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
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    text(d, "How many events happen in one unit of time?", W / 2, 42, FS_TITLE, weight="bold")

    panel(d, 40, 1, "rate = 1 per unit of time")
    panel(d, 620, 4, "rate = 4 per unit of time")

    name = "gillespie_poisson"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
