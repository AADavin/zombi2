"""Figure (Ch16): the autocorrelated clocks as rate *processes* over time.

Both autocorrelated clocks let the rate evolve continuously, but they differ in one
decisive way -- is there anything pulling the rate back? Each panel plots many sample
paths of the instantaneous rate against elapsed time.

  * Left (autocorrelated lognormal): a geometric random walk. Nothing restores it, so the
    paths fan out without bound -- the longer you wait, the further a lineage's rate can
    stray from 1, and the tree's total length wanders with it.
  * Right (Cox-Ingersoll-Ross): a mean-reverting diffusion, dr = theta*(mean - r) dt +
    sigma*sqrt(r) dW. The drift term pulls every path back toward the long-run mean, so the
    spread stabilises and the total tree length stays close to mean x time. CIR also varies
    the rate *within* a branch (the path wiggles along a single lineage), unlike the
    lognormal walk, which only jumps at nodes.

House style: B&W paths, the long-run mean drawn as a dashed reference, ASCII text.

Run:  python figures/scripts/fig_clock_cir.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

import clock_common as C
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1300, 620
T = 3.0
YMAX = 3.4
N_PATHS = 16
MEAN = 1.0


def panel(d, ox, oy, pw, ph, title, subtitle, t, paths, note):
    x_at = lambda v: ox + (v / T) * pw                    # noqa: E731
    y_at = lambda v: oy + ph - (v / YMAX) * ph            # noqa: E731 (unclamped; clipped below)

    d.append(draw.Text(title, FS_LABEL, ox + pw / 2, oy - 40, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text(subtitle, FS_TICK, ox + pw / 2, oy - 18, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    # frame
    d.append(draw.Line(ox, oy + ph, ox + pw, oy + ph, stroke=INK, stroke_width=1.6))
    d.append(draw.Line(ox, oy, ox, oy + ph, stroke=INK, stroke_width=1.6))
    # sample paths (thin, translucent ink), clipped to the plot box so paths that
    # run off the top simply exit the frame rather than piling up on the ceiling
    clip = draw.ClipPath()
    clip.append(draw.Rectangle(ox, oy, pw, ph))
    d.append(clip)
    g = draw.Group(clip_path=clip)
    for p in paths:
        pts = []
        for tv, rv in zip(t, p):
            pts += [x_at(tv), y_at(rv)]
        g.append(draw.Lines(*pts, close=False, fill="none", stroke=INK, stroke_width=1.3,
                            stroke_opacity=0.34, stroke_linejoin="round"))
    d.append(g)
    # long-run mean reference
    d.append(draw.Line(ox, y_at(MEAN), ox + pw, y_at(MEAN), stroke=INK, stroke_width=1.8,
                       stroke_dasharray="7,5"))
    # label sits over the densest band of paths, so give it a white plate (halo) to stay legible;
    # lift it a little further off the dashed line as well.
    _lab = "long-run mean = 1"
    _lx, _ly = ox + pw - 6, y_at(MEAN) - 16
    _lw, _lh = 0.56 * FS_TICK * len(_lab), FS_TICK + 8
    d.append(draw.Rectangle(_lx - _lw, _ly - _lh / 2 - 2, _lw + 8, _lh, rx=4, ry=4,
                            fill="white", fill_opacity=0.86, stroke="none"))
    d.append(draw.Text(_lab, FS_TICK, _lx, _ly, font_family=FONT,
                       text_anchor="end", dominant_baseline="central", fill=INK, font_weight="bold"))
    # y ticks
    for v in (0, 1, 2, 3):
        d.append(draw.Line(ox - 5, y_at(v), ox, y_at(v), stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"{v:g}", FS_TICK, ox - 10, y_at(v) + 4, font_family=FONT,
                           text_anchor="end", fill=MUTED))
    # x ticks
    for v in (0, 1, 2, 3):
        d.append(draw.Line(x_at(v), oy + ph, x_at(v), oy + ph + 5, stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"{v:g}", FS_TICK, x_at(v), oy + ph + 20, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
    d.append(draw.Text("elapsed time along the lineage", FS_TICK, ox + pw / 2, oy + ph + 42,
                       font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("rate", FS_TICK, ox - 34, oy + ph / 2, font_family=FONT,
                       text_anchor="middle", fill=INK, transform=f"rotate(-90,{ox - 34},{oy + ph / 2})"))
    # takeaway note, bottom-left inside the panel
    d.append(draw.Text(note, FS_TICK, ox + 12, oy + 22, font_family=FONT, text_anchor="start",
                       fill=MUTED))


def main():
    t_g, gbm = C.gbm_paths(T, N_PATHS, sigma=0.6, seed=11)
    t_c, cir = C.cir_paths(T, N_PATHS, theta=3.5, sigma=0.6, mean=MEAN, seed=11)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Does anything pull the rate back? Random walk vs mean reversion",
                       FS_TITLE, W / 2, 44, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    pw, ph = 520, 400
    oy = 150
    panel(d, 80, oy, pw, ph, "Autocorrelated lognormal", "geometric random walk -- no restoring force",
          t_g, gbm, "paths fan out; length wanders")
    panel(d, 80 + pw + 100, oy, pw, ph, "Cox-Ingersoll-Ross", "mean-reverting diffusion",
          t_c, cir, "paths pulled back; length stabilises")

    name = "clock_cir"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    main()
