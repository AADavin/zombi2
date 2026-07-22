"""Figure: ZOMBI2's continuous-trait models as diffusion processes.

The continuous trait models — Brownian motion, Ornstein-Uhlenbeck and Early-Burst — are
Gaussian diffusions: a trait value wanders through time, and the model sets how its variance
accumulates. This is the continuous companion to the discrete Markov-chain gallery. Each
panel shows real sample paths, stepped through the model's own ``evolve()`` (so they are
exactly ZOMBI2's dynamics), with the theoretical variance envelope.

Rendered in colour and B&W.  Run:  python figures/scripts/fig_trait_diffusion.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2 import BrownianMotion, EarlyBurst, OrnsteinUhlenbeck

from zombi_style import FONT, INK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1260, 580
T, STEPS, NPATH = 1.0, 260, 7
PY0, PH, PW = 172, 278, 330
COLS = [90, 475, 860]                                     # plot-box left edges
BOX_COLOR, REF = "#bdbdbd", "#8a8a8a"


def path(model, x0, rng):
    dt = T / STEPS
    xs, x, t = [x0], x0, 0.0
    for _ in range(STEPS):
        x, _ = model.evolve(x, dt, t, rng)
        t += dt
        xs.append(x)
    return np.linspace(0, T, STEPS + 1), np.array(xs)


# --------------------------------------------------------------------------- plotting
def _x(box, dr, t):
    px0, _, pw, _ = box
    return px0 + (t - dr[0]) / (dr[1] - dr[0]) * pw


def _y(box, dr, y):
    _, py0, _, ph = box
    return py0 + ph - (y - dr[2]) / (dr[3] - dr[2]) * ph


def polyline(d, box, dr, ts, ys, color, lw=1.5, opacity=1.0, dash=None):
    pts = []
    for t, y in zip(ts, ys):
        pts += [_x(box, dr, t), _y(box, dr, max(dr[2], min(dr[3], y)))]
    kw = {"stroke": color, "stroke_width": lw, "fill": "none", "stroke_opacity": opacity,
          "stroke_linejoin": "round", "stroke_linecap": "round"}
    if dash:
        kw["stroke_dasharray"] = dash
    d.append(draw.Lines(*pts, close=False, **kw))


def axes(d, box, dr, cx):
    px0, py0, pw, ph = box
    d.append(draw.Rectangle(px0, py0, pw, ph, fill="none", stroke=BOX_COLOR, stroke_width=1.3))
    d.append(draw.Lines(px0 + pw + 4, py0 + ph, px0 + pw - 4, py0 + ph - 4,
                        px0 + pw - 4, py0 + ph + 4, close=True, fill=INK))   # time arrowhead
    d.append(draw.Text("0", 11, px0, py0 + ph + 18, font_family=FONT, text_anchor="middle", fill="#777"))
    d.append(draw.Text("present", 11, px0 + pw, py0 + ph + 18, font_family=FONT, text_anchor="middle", fill="#777"))
    d.append(draw.Text("time", 12.5, cx, py0 + ph + 34, font_family=FONT, text_anchor="middle", fill="#777"))
    d.append(draw.Text("trait value", 12.5, px0 - 26, py0 + ph / 2, font_family=FONT, text_anchor="middle",
                       fill="#777", transform=f"rotate(-90 {px0 - 26} {py0 + ph / 2})"))


def refline(d, box, dr, y, color, label):
    yy = _y(box, dr, y)
    d.append(draw.Line(box[0], yy, box[0] + box[2], yy, stroke=color, stroke_width=1.4, stroke_dasharray="5,4"))
    d.append(draw.Text(label, 12, box[0] + box[2] - 4, yy - 7, font_family=FONT, text_anchor="end",
                       dominant_baseline="central", fill=color, font_weight="bold"))


def panel(d, i, title, formula, caption, model_paths, dr, envelope, refl, color):
    box = (COLS[i], PY0, PW, PH)
    cx = COLS[i] + PW / 2
    axes(d, box, dr, cx)
    d.append(draw.Text(title, 16, cx, 118, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    d.append(draw.Text(formula, 13.5, cx, 140, font_family=FONT, text_anchor="middle", fill="#555"))
    if refl:
        refline(d, box, dr, *refl)
    if envelope:
        te, mid, env = envelope
        polyline(d, box, dr, te, mid + env, REF, lw=1.3, dash="4,3")
        polyline(d, box, dr, te, mid - env, REF, lw=1.3, dash="4,3")
    for ts, ys in model_paths:
        polyline(d, box, dr, ts, ys, color, lw=1.5, opacity=0.72)
    d.append(draw.Text(caption, 12.5, cx, PY0 + PH + 52, font_family=FONT, text_anchor="middle", fill="#666"))


# --------------------------------------------------------------------------- render
def render(mode):
    rng = np.random.default_rng(7)
    te = np.linspace(0, T, 60)

    bm = BrownianMotion(sigma2=1.0)
    bm_paths = [path(bm, 0.0, rng) for _ in range(NPATH)]
    bm_env = (te, np.zeros_like(te), 2 * np.sqrt(1.0 * te))

    ou = OrnsteinUhlenbeck(sigma2=0.5, alpha=3.0, theta=2.0)
    ou_paths = [path(ou, x0, rng) for x0 in np.linspace(-1.0, 5.0, NPATH)]
    ssd = np.sqrt(0.5 / (2 * 3.0))
    ou_env = (te, np.full_like(te, 2.0), np.full_like(te, 2 * ssd))

    eb = EarlyBurst(sigma2=9.0, rate=-4.5)
    eb_paths = [path(eb, 0.0, rng) for _ in range(NPATH)]
    eb_env = (te, np.zeros_like(te), 2 * np.sqrt(9.0 * (1 - np.exp(-4.5 * te)) / 4.5))

    if mode == "bw":
        cbm = cou = ceb = "#555555"
    else:
        cbm, cou, ceb = "#4477AA", "#AA3377", "#228833"

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Continuous-trait models as diffusion processes", 20, 40, 40,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("sample paths of the trait through time (stepped through each model's evolve()); "
                       "dashed = 95% variance envelope", 13.5, 40, 64, font_family=FONT,
                       text_anchor="start", fill="#777"))
    panel(d, 0, "Brownian motion (BM)", "dx = σ·dW",
          "no optimum — variance grows without bound", bm_paths, (0, T, -3.2, 3.2), bm_env,
          (0.0, REF, "x0"), cbm)
    panel(d, 1, "Ornstein-Uhlenbeck (OU)", "dx = α·(θ - x)·dt + σ·dW",
          "pulled toward the optimum θ — variance saturates", ou_paths, (0, T, -1.8, 5.8), ou_env,
          (2.0, cou, "θ"), cou)
    panel(d, 2, "Early burst (EB)", "σ²(t) = σ²·exp(r·t),  r < 0",
          "rate decays through time — disparity front-loaded", eb_paths, (0, T, -3.2, 3.2), eb_env,
          (0.0, REF, "x0"), ceb)

    name = "trait_diffusion" if mode == "color" else "trait_diffusion_bw"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)


def main():
    for mode in ("color", "bw"):
        render(mode)
    print("wrote trait_diffusion (+_bw)")


if __name__ == "__main__":
    main()
