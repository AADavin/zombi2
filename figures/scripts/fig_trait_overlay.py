"""Figure: the two trait-model families — how a branch is traversed.

Every trait in ZOMBI2 runs on one overlay engine: a value is dropped at the root and evolved
branch by branch. The two families differ only in how a single branch is traversed, and that is
what this opener contrasts:

  * CONTINUOUS (left): the whole branch is one draw. Along a branch of length t the change is a
    single Normal draw whose mean and variance the model sets; nothing between the endpoints is
    simulated because nothing needs to be.
  * DISCRETE (right): the Markov jumps are simulated exactly. The chain steps through states with
    random waiting times, so the realized history along the branch — a stochastic character map —
    comes for free.

House style: B&W, one centered title, ASCII text.

Run:  python figures/scripts/fig_trait_overlay.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 560
GREY = "#9a9a9a"


def axes(d, x0, y0, w, h, xlabel, ylabel):
    d.append(draw.Line(x0, y0, x0, y0 - h, stroke=INK, stroke_width=1.6))          # y
    d.append(draw.Line(x0, y0, x0 + w, y0, stroke=INK, stroke_width=1.6))          # x
    d.append(draw.Text(xlabel, FS_TICK, x0 + w / 2, y0 + 34, font_family=FONT,
                       text_anchor="middle", fill=MUTED))
    d.append(draw.Text(ylabel, FS_TICK, x0 - 34, y0 - h / 2, font_family=FONT, text_anchor="middle",
                       fill=MUTED, transform=f"rotate(-90, {x0 - 34}, {y0 - h / 2})"))


# --------------------------------------------------------------------------- panel A: continuous
def panel_continuous(d, x0, ytop, w, h):
    ybase = ytop + h
    d.append(draw.Text("continuous", FS_LABEL, x0 + w / 2, ytop - 26,
                       font_family=FONT, text_anchor="middle", fill=INK, font_weight="bold"))
    axes(d, x0, ybase, w, h, "time (one branch)", "trait value")
    ymid = ytop + h / 2
    root = (x0, ymid)
    d.append(draw.Circle(*root, 5, fill=INK))
    d.append(draw.Text("root value", FS_TICK, x0 + 10, ymid - 14, font_family=FONT,
                       text_anchor="start", fill=MUTED))

    # annotation at the top, clear of the data
    d.append(draw.Text("endpoint = one Normal draw", FS_ANNOT, x0 + w / 2, ytop + 8, font_family=FONT,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("(the path is never simulated)", FS_TICK, x0 + w / 2, ytop + 30,
                       font_family=FONT, text_anchor="middle", fill=MUTED))

    # normal endpoint distribution at the far end of the branch (kept inside the plot box).
    # sd is deliberately modest so the curve's peak stays clear of the subtitle above.
    xend = x0 + 0.66 * w
    sd = h * 0.105
    ys = np.linspace(ymid - 3 * sd, ymid + 3 * sd, 80)
    dens = np.exp(-0.5 * ((ys - ymid) / sd) ** 2)
    bell = draw.Path(fill="none", stroke=INK, stroke_width=2.4)
    bell.M(xend, ys[0])
    for yy, dd in zip(ys, dens):
        bell.L(xend + dd * 0.26 * w, yy)
    d.append(bell)
    d.append(draw.Line(xend, ys[0], xend, ys[-1], stroke=GREY, stroke_width=1.2))

    # a handful of sampled endpoints on the branch-end line; one is "the" draw
    draws = [ymid + sd * z for z in (-1.6, -0.7, 0.15, 0.9, 1.7)]
    chosen = draws[2]
    d.append(draw.Line(x0, ymid, xend, chosen, stroke=GREY, stroke_width=1.6, stroke_dasharray="7,6"))
    for yy in draws:
        d.append(draw.Circle(xend, yy, 4.5, fill="white", stroke=INK, stroke_width=1.6))
    d.append(draw.Circle(xend, chosen, 6.0, fill=INK, stroke=INK, stroke_width=1.6))


# --------------------------------------------------------------------------- panel B: discrete
def panel_discrete(d, x0, ytop, w, h):
    ybase = ytop + h
    d.append(draw.Text("discrete", FS_LABEL, x0 + w / 2, ytop - 26,
                       font_family=FONT, text_anchor="middle", fill=INK, font_weight="bold"))
    axes(d, x0, ybase, w, h, "time (one branch)", "state")

    levels = [ybase - h * f for f in (0.2, 0.5, 0.8)]      # states 0,1,2 (bottom..top)
    labels = ["0", "1", "2"]
    for yy, lab in zip(levels, labels):
        d.append(draw.Line(x0, yy, x0 + w, yy, stroke="#e4e4e4", stroke_width=1.0))
        d.append(draw.Text(lab, FS_TICK, x0 - 12, yy, font_family=FONT, text_anchor="end",
                           dominant_baseline="central", fill=MUTED))

    # a hand-drawn CTMC path: state jumps at random waiting times
    rng = np.random.default_rng(3)
    seq = [1]                                              # start in state 1
    times = [0.0]
    t = 0.0
    while t < 1.0:
        t += rng.exponential(0.28)
        if t >= 1.0:
            break
        cur = seq[-1]
        nxt = rng.choice([s for s in (0, 1, 2) if s != cur])
        times.append(t)
        seq.append(nxt)
    times.append(1.0)
    xt = lambda u: x0 + u * w                              # noqa: E731
    for i, s in enumerate(seq):
        xa, xb = xt(times[i]), xt(times[i + 1])
        d.append(draw.Line(xa, levels[s], xb, levels[s], stroke=INK, stroke_width=4.0,
                           stroke_linecap="round"))
        if i > 0:                                          # vertical riser at the jump
            d.append(draw.Line(xa, levels[seq[i - 1]], xa, levels[s], stroke=INK, stroke_width=2.0))
            d.append(draw.Circle(xa, levels[s], 4.0, fill=INK))
    d.append(draw.Text("every jump simulated exactly", FS_ANNOT, x0 + w / 2, ytop + 6,
                       font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("(a stochastic character map)", FS_TICK, x0 + w / 2, ytop + 28,
                       font_family=FONT, text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- render
def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Two ways a branch is traversed", FS_TITLE, W / 2, 46, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    panel_continuous(d, 110, 150, 360, 300)
    panel_discrete(d, 720, 150, 360, 300)

    name = "trait_overlay"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
