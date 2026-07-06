"""Figure: QuaSSE -- continuous-trait state-dependent diversification.

The trait is a real number that diffuses by Brownian motion along the branches, and the speciation
rate is a *function* of the current value: lambda(x) rises with x here, so high-trait lineages branch
faster and the tree fills with high (yellow) values.

  * Panel A (the model): lambda(x) (rising sigmoid) and mu(x) (flat) as curves over the trait axis;
    the trait axis is tinted with the same viridis ramp used to paint the tree.
  * Panel B (a realization): one simulate_sse(QuaSSE) tree, each branch painted by its trait value
    (viridis); the extant tips carry colored chips. The high-value (yellow) lineages proliferate.

House style: viridis for the continuous trait (as in the BM/OU figures), one centered title.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_sse_quasse.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2.coevolve import QuaSSE, simulate_sse

from fig_trait_bm import viridis, hexc, VIRIDIS
from fig_trait_pagel import _layout
from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 700
N_TIPS = 14

# lambda(x) rises sigmoidally with the trait; mu flat; slow diffusion
SPEC = lambda x: 0.6 + 2.6 / (1.0 + np.exp(-1.6 * x))     # noqa: E731
EXT = lambda x: 0.3                                        # noqa: E731
SIGMA2, RATE_BOUND, X0 = 0.5, 6.0, 0.0
XLO, XHI = -3.2, 3.2


# --------------------------------------------------------------------------- panel A: the model
def panel_model(d, ox, oy, pw, ph):
    d.append(draw.Text("A   the model", FS_LABEL, ox, oy - 18, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    x_at = lambda x: ox + (x - XLO) / (XHI - XLO) * pw      # noqa: E731
    rmax = SPEC(XHI) * 1.08
    y_at = lambda r: oy + ph - (r / rmax) * ph              # noqa: E731

    # viridis strip along the trait axis
    grad = draw.LinearGradient(ox, 0, ox + pw, 0)
    for t, c in VIRIDIS:
        grad.add_stop(t, hexc(c))
    d.append(grad)
    d.append(draw.Rectangle(ox, oy + ph + 10, pw, 12, fill=grad, stroke=INK, stroke_width=0.8))
    d.append(draw.Text("trait value  x", FS_TICK, ox + pw / 2, oy + ph + 46, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    # axes
    d.append(draw.Line(ox, oy, ox, oy + ph, stroke=INK, stroke_width=1.4))
    d.append(draw.Line(ox, oy + ph, ox + pw, oy + ph, stroke=INK, stroke_width=1.4))
    d.append(draw.Text("rate", FS_TICK, ox - 14, oy + ph / 2, font_family=FONT,
                       text_anchor="middle", dominant_baseline="central", fill=MUTED,
                       transform=f"rotate(-90 {ox - 14} {oy + ph / 2})"))

    xs = np.linspace(XLO, XHI, 160)
    lam = draw.Path(fill="none", stroke=INK, stroke_width=3.4)
    for i, x in enumerate(xs):
        (lam.M if i == 0 else lam.L)(x_at(x), y_at(SPEC(x)))
    d.append(lam)
    mu = draw.Path(fill="none", stroke=MUTED, stroke_width=3.0, stroke_dasharray="7,5")
    for i, x in enumerate(xs):
        (mu.M if i == 0 else mu.L)(x_at(x), y_at(EXT(x)))
    d.append(mu)
    # place the lambda(x) label in the open space above the curve's shoulder (not touching it):
    # sit it at the plateau height but shifted left, where the rising curve is still well below.
    d.append(draw.Text("lambda(x)", FS_ANNOT, x_at(0.35), y_at(SPEC(2.9)) - 20, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    d.append(draw.Text("mu(x)", FS_ANNOT, x_at(-2.1), y_at(EXT(-2.1)) - 14, font_family=FONT,
                       text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- panel B: a realization
def _pick(seed_range):
    model = QuaSSE(speciation=SPEC, extinction=EXT, sigma2=SIGMA2, rate_bound=RATE_BOUND, x0=X0)
    best = None
    for s in seed_range:
        res = simulate_sse(model, n_tips=N_TIPS, seed=s)
        ete = zombi_to_ete3(res.tree)
        if len(ete.get_leaves()) > 26:
            continue
        vals = {n.name: float(v) for n, v in res.node_values.items()}
        extant = [vals[n.name] for n in ete.get_leaves() if n.is_extant]
        spread = max(vals.values()) - min(vals.values())
        score = np.mean(extant) + 0.3 * spread                # high-value tips + visible spread
        if best is None or score > best[0]:
            best = (score, s, res, ete, vals)
    return best


def panel_realization(d, ox, oy, pw, ph):
    _, seed, res, ete, vals = _pick(range(1, 160))
    lo, hi = min(vals.values()), max(vals.values())
    norm = lambda v: (v - lo) / (hi - lo) if hi > lo else 0.5   # noqa: E731
    col = lambda name: hexc(viridis(norm(vals[name])))          # noqa: E731

    tfo, present, ys, nleaf = _layout(ete)
    x_at = lambda t: ox + 40 + (t / present) * (pw - 150)       # noqa: E731
    y_at = lambda k: oy + 40 + (k / max(1, nleaf - 1)) * (ph - 96)   # noqa: E731

    d.append(draw.Text("B   a simulated realization", FS_LABEL, ox, oy - 18, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        # paint the branch by the mean of its endpoints' values (a smooth read of the diffusion)
        mid = (vals[n.name] + vals[n.up.name]) / 2
        c = hexc(viridis(norm(mid)))
        d.append(draw.Line(x_at(tfo[n.up.name]), y, x_at(tfo[n.name]), y, stroke=c,
                           stroke_width=4.4, stroke_linecap="butt"))
        if n.is_leaf() and not n.is_extant:
            xx, ah = x_at(tfo[n.name]), 5.5
            d.append(draw.Line(xx - ah, y - ah, xx + ah, y + ah, stroke=INK, stroke_width=2.0))
            d.append(draw.Line(xx - ah, y + ah, xx + ah, y - ah, stroke=INK, stroke_width=2.0))
    for n in ete.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=col(n.name), stroke_width=3.2))

    colc = ox + pw - 34
    for n in ete.get_leaves():
        if n.is_extant:
            y = y_at(ys[n.name])
            d.append(draw.Rectangle(colc - 8, y - 8, 16, 16, fill=col(n.name),
                                    stroke=INK, stroke_width=1.0))

    # colour bar
    grad = draw.LinearGradient(ox, 0, ox + 200, 0)
    for t, c in VIRIDIS:
        grad.add_stop(t, hexc(c))
    d.append(grad)
    by = oy + ph + 6
    d.append(draw.Rectangle(ox, by, 200, 14, fill=grad, stroke=INK, stroke_width=0.8))
    d.append(draw.Text(f"{lo:+.1f}", FS_TICK, ox, by + 30, font_family=FONT, text_anchor="start", fill="#555"))
    d.append(draw.Text(f"{hi:+.1f}", FS_TICK, ox + 200, by + 30, font_family=FONT, text_anchor="end", fill="#555"))
    d.append(draw.Text("trait value", FS_TICK, ox + 100, by + 30, font_family=FONT,
                       text_anchor="middle", fill=MUTED))
    return seed


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Continuous-trait diversification (QuaSSE)", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    panel_model(d, 110, 190, 360, 300)
    seed = panel_realization(d, 620, 150, 520, 400)

    name = "sse_quasse"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png  (tree seed {seed})")


if __name__ == "__main__":
    render()
