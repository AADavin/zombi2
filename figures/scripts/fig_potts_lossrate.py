"""Figure: coupling enters the dynamics through the loss rate (Potts, idea 3).

The couplings do not act directly on gain; they act on how fast a present family is LOST. Each
present family i feels a local field f_i = h_i + sum_j J_ij sigma_j (a sum over its present
partners), and is lost at rate  loss_i = base_loss * exp(-beta * f_i). So a family surrounded by
present partners has a large field, an exponentially small loss rate, and is retained; a family
with no partners present is lost fast. That is what makes a coupled module stay together.

House style: colour (didactic), one centered bold title, ASCII text, legend clear of the curve.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_potts_lossrate.py
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))   # zombi_style

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, COOCCUR, AVOID, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent
NAME = "potts_lossrate"

W, H = 1020, 720
POS, NEG = COOCCUR, AVOID                         # protected / fast-loss regimes


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Coupling enters through the loss rate", FS_TITLE, W / 2, 52,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # the two equations, centered under the title (ASCII only)
    d.append(draw.Text("local field   f_i  =  h_i  +  sum_j J_ij sigma_j     (over present partners)",
                       FS_ANNOT, W / 2, 100, font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("loss rate   loss_i  =  base_loss * exp(-beta * f_i)",
                       FS_ANNOT, W / 2, 134, font_family=FONT, text_anchor="middle", fill=INK))

    # loss-vs-field curve, generously sized. Range chosen so the whole visible curve is the
    # smooth exponential (exp(1.4) = 4.06 < lmax) with no clipping plateau at the top-left.
    bx, by, bw, bh = 190, 190, 640, 380
    fmin, fmax, lmax = -1.4, 3.0, 4.2

    def X(f):
        return bx + (f - fmin) / (fmax - fmin) * bw

    def Y(v):
        return by + bh - min(v, lmax) / lmax * bh

    # axes
    d.append(draw.Line(bx, by, bx, by + bh, stroke="#bdbdbd", stroke_width=1.6))
    d.append(draw.Line(bx, by + bh, bx + bw, by + bh, stroke="#bdbdbd", stroke_width=1.6))
    # field = 0 guide line
    d.append(draw.Line(X(0), by, X(0), by + bh, stroke="#dcdcdc", stroke_width=1.4,
                       stroke_dasharray="4,4"))

    # the exponential curve
    pts = []
    f = fmin
    while f <= fmax + 1e-9:
        pts += [X(f), Y(math.exp(-f))]
        f += 0.02
    d.append(draw.Lines(*pts, close=False, fill="none", stroke=INK, stroke_width=3.4,
                        stroke_linejoin="round"))

    # axis titles
    d.append(draw.Text("loss rate", FS_LABEL, bx - 44, by + bh / 2, font_family=FONT,
                       text_anchor="middle", fill=MUTED,
                       transform=f"rotate(-90 {bx - 44} {by + bh / 2})"))
    d.append(draw.Text("local field  f", FS_LABEL, bx + bw / 2, by + bh + 46, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    # base_loss reference point at f = 0
    d.append(draw.Circle(X(0), Y(1), 5, fill="#999"))
    d.append(draw.Text("base_loss", FS_TICK, X(0) + 12, Y(1) - 14, font_family=FONT,
                       text_anchor="start", fill=MUTED))

    # regime 1: incomplete module (partners absent) -> fast loss
    d.append(draw.Circle(X(-1.0), Y(math.exp(1)), 8, fill=NEG, stroke=INK, stroke_width=1.6))
    d.append(draw.Text("module incomplete", FS_TICK, X(-1.0) + 20, Y(math.exp(1)) - 12,
                       font_family=FONT, text_anchor="start", fill=NEG, font_weight="bold"))
    d.append(draw.Text("(partners absent): fast loss", FS_TICK, X(-1.0) + 20, Y(math.exp(1)) + 12,
                       font_family=FONT, text_anchor="start", fill=NEG, font_weight="bold"))

    # regime 2: complete module (partners present) -> protected
    d.append(draw.Circle(X(2.0), Y(math.exp(-2)), 8, fill=POS, stroke=INK, stroke_width=1.6))
    d.append(draw.Text("module complete", FS_TICK, X(2.0), Y(math.exp(-2)) - 44,
                       font_family=FONT, text_anchor="middle", fill=POS, font_weight="bold"))
    d.append(draw.Text("(partners present): protected", FS_TICK, X(2.0), Y(math.exp(-2)) - 20,
                       font_family=FONT, text_anchor="middle", fill=POS, font_weight="bold"))

    # caption
    d.append(draw.Text("more present partners  =>  lower loss, so a complete module stays together",
                       FS_TICK, W / 2, by + bh + 96, font_family=FONT, text_anchor="middle",
                       fill=MUTED))

    out = OUT_DIR / NAME
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{NAME}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{NAME}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{NAME}.svg / .png")


if __name__ == "__main__":
    render()
