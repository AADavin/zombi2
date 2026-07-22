"""Figure (Ch16): the per-branch rate distributions of the uncorrelated clocks.

Where the family figure shows the clocks *on a tree*, this one shows the distribution a
single branch draws its rate from -- and how each clock's one knob controls the spread.
All three are centred on mean rate 1, so on average the tree's length is preserved; they
differ in how far individual branches stray from 1.

  * Lognormal: `sigma` is the spread. sigma -> 0 collapses to a spike at 1 (the strict clock).
  * Gamma: the spread is set by `shape`, but INVERSELY -- a *large* shape concentrates
    rates near 1, a small shape spreads them out.
  * White noise: same gamma shape, but the variance is sigma^2 / (branch length). A short
    branch draws a wildly variable rate; a long branch averages the noise away toward 1.
    This branch-length dependence is what sets white noise apart.

House style: B&W line plots, curves labelled in place (no legend box), one centered title.

Run:  python figures/scripts/fig_clock_distributions.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1300, 560
XMAX = 3.0
GRID = np.linspace(1e-3, XMAX, 700)


def lognormal_pdf(x, sigma):
    mu = -sigma * sigma / 2.0
    return np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma * sigma)) / (x * sigma * math.sqrt(2 * math.pi))


def gamma_pdf(x, k, theta):
    return x ** (k - 1) * np.exp(-x / theta) / (theta ** k * math.gamma(k))


# curve shade ramp: light grey -> black, keyed to the "amount of spread" (light = tight)
SHADES = ["#111111", "#6f6f6f", "#b3b3b3"]


def panel(d, ox, oy, pw, ph, title, subtitle, curves, ymax=None):
    """curves: list of (y_array, label, shade, label_x). Draws framed axes + curves."""
    x_at = lambda v: ox + (v / XMAX) * pw               # noqa: E731
    ymax = ymax or max(float(np.nanmax(y)) for y, *_ in curves) * 1.08
    y_at = lambda v: oy + ph - (v / ymax) * ph          # noqa: E731

    d.append(draw.Text(title, FS_LABEL, ox + pw / 2, oy - 40, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text(subtitle, FS_TICK, ox + pw / 2, oy - 18, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    # axes
    d.append(draw.Line(ox, oy + ph, ox + pw, oy + ph, stroke=INK, stroke_width=1.6))
    d.append(draw.Line(ox, oy, ox, oy + ph, stroke=INK, stroke_width=1.6))
    # rate = 1 guide
    d.append(draw.Line(x_at(1.0), oy, x_at(1.0), oy + ph, stroke="#dcdcdc", stroke_width=1.2))
    d.append(draw.Text("1", FS_TICK, x_at(1.0), oy + ph + 20, font_family=FONT,
                       text_anchor="middle", fill=MUTED))
    for v in (0, XMAX):
        d.append(draw.Line(x_at(v), oy + ph, x_at(v), oy + ph + 5, stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"{v:g}", FS_TICK, x_at(v), oy + ph + 20, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
    d.append(draw.Text("branch rate", FS_TICK, ox + pw / 2, oy + ph + 42, font_family=FONT,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("density", FS_TICK, ox - 14, oy + ph / 2, font_family=FONT,
                       text_anchor="middle", fill=INK, transform=f"rotate(-90,{ox - 14},{oy + ph / 2})"))
    # curves
    for y, label, shade, lx in curves:
        yy = np.clip(y, 0, ymax)
        pts = []
        for xv, yv in zip(GRID, yy):
            pts += [x_at(xv), y_at(yv)]
        d.append(draw.Lines(*pts, close=False, fill="none", stroke=shade, stroke_width=3.0,
                            stroke_linejoin="round"))
        # label sits ON the curve at x = lx, so the three labels cascade apart
        j = int(np.argmin(np.abs(GRID - lx)))
        ylab = min(y_at(yy[j]) - 8, oy + ph - 8)
        d.append(draw.Text(label, FS_TICK, x_at(lx), ylab, font_family=FONT,
                           text_anchor="start", fill=shade, font_weight="bold"))


def main():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("What rate does one branch draw? The uncorrelated clocks", FS_TITLE,
                       W / 2, 44, font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    pw, ph = 340, 320
    oy = 150
    xs = [70, 70 + pw + 40, 70 + 2 * (pw + 40)]

    # A -- lognormal: sigma controls spread
    panel(d, xs[0], oy, pw, ph, "Lognormal", "knob: sigma  (0 = strict clock)",
          [(lognormal_pdf(GRID, 0.25), "sigma 0.25", SHADES[2], 1.2),
           (lognormal_pdf(GRID, 0.5), "sigma 0.5", SHADES[1], 1.75),
           (lognormal_pdf(GRID, 1.0), "sigma 1.0", SHADES[0], 2.35)],
          ymax=1.9)

    # B -- gamma: shape controls spread INVERSELY
    panel(d, xs[1], oy, pw, ph, "Gamma", "knob: shape  (large = tight)",
          [(gamma_pdf(GRID, 10.0, 0.1), "shape 10", SHADES[2], 1.15),
           (gamma_pdf(GRID, 3.0, 1 / 3), "shape 3", SHADES[1], 1.7),
           (gamma_pdf(GRID, 1.0, 1.0), "shape 1", SHADES[0], 2.25)],
          ymax=1.5)

    # C -- white noise: variance set by BRANCH LENGTH
    sig = 0.5
    def wn(dt):
        var = sig * sig / dt
        return gamma_pdf(GRID, 1.0 / var, var)
    panel(d, xs[2], oy, pw, ph, "White noise", "variance = sigma^2 / branch length",
          [(wn(2.5), "long branch", SHADES[2], 1.15),
           (wn(0.6), "medium", SHADES[1], 1.7),
           (wn(0.15), "short branch", SHADES[0], 2.3)],
          ymax=1.9)

    name = "clock_distributions"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    main()
