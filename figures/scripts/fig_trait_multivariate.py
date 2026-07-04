"""Figure: multivariate Brownian motion — two traits evolving jointly.

``MultivariateBrownian`` diffuses a vector-valued trait with a rate (covariance) matrix R;
the off-diagonal of R makes the traits evolve in a correlated way. This figure plots the
joint (X, Y) plane: a cloud of endpoints (each stepped through the model's own ``evolve()``),
the 95% covariance ellipse, and a few sample trajectories — for a correlated R (tilted
ellipse) versus an independent R (upright circle).

Rendered in colour and B&W.  Run:  python figures/scripts/fig_trait_multivariate.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2 import MultivariateBrownian

from zombi_style import FONT, INK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1120, 600
T, NPTS, NTRAJ, TSTEPS = 1.0, 280, 3, 140
DR = (-3.3, 3.3)
BS = 384
BOXES = [(110, 150), (628, 150)]                             # (x, y) top-left of each square box
MODE = "color"


def _mx(box, x):
    return box[0] + (x - DR[0]) / (DR[1] - DR[0]) * BS


def _my(box, y):
    return box[1] + BS - (y - DR[0]) / (DR[1] - DR[0]) * BS


def endpoints(model, rng, n):
    return np.array([model.evolve(np.zeros(2), T, 0.0, rng)[0] for _ in range(n)])


def trajectory(model, rng):
    dt, x, pts = T / TSTEPS, np.zeros(2), [np.zeros(2)]
    for _ in range(TSTEPS):
        x, _ = model.evolve(x, dt, 0.0, rng)
        pts.append(x.copy())
    return np.array(pts)


def ellipse_pts(R):
    C = T * R
    vals, vecs = np.linalg.eigh(C)
    th = np.linspace(0, 2 * np.pi, 90)
    return np.array([2 * np.sqrt(vals[0]) * np.cos(t) * vecs[:, 0]
                     + 2 * np.sqrt(vals[1]) * np.sin(t) * vecs[:, 1] for t in th])


def panel(d, box, R, title, pt_c, ell_c, traj_c, rng):
    bx, by = box
    cx0, cy0 = _mx(box, 0), _my(box, 0)
    d.append(draw.Rectangle(bx, by, BS, BS, fill="none", stroke="#cccccc", stroke_width=1.2))
    d.append(draw.Line(bx, cy0, bx + BS, cy0, stroke="#d9d9d9", stroke_width=1))     # X axis
    d.append(draw.Line(cx0, by, cx0, by + BS, stroke="#d9d9d9", stroke_width=1))     # Y axis
    d.append(draw.Text("trait X", 12, bx + BS - 6, cy0 - 8, font_family=FONT, text_anchor="end", fill="#999"))
    d.append(draw.Text("trait Y", 12, cx0 + 8, by + 14, font_family=FONT, text_anchor="start", fill="#999"))

    for (x, y) in endpoints(R_model(R), rng, NPTS):
        if DR[0] < x < DR[1] and DR[0] < y < DR[1]:
            d.append(draw.Circle(_mx(box, x), _my(box, y), 2.6, fill=pt_c, fill_opacity=0.5))

    ell = ellipse_pts(R)
    pts = []
    for (x, y) in ell:
        pts += [_mx(box, x), _my(box, y)]
    d.append(draw.Lines(*pts, close=True, fill="none", stroke=ell_c, stroke_width=2.4))

    for _ in range(NTRAJ):
        tr = trajectory(R_model(R), rng)
        pp = []
        for (x, y) in tr:
            pp += [_mx(box, max(DR[0], min(DR[1], x))), _my(box, max(DR[0], min(DR[1], y)))]
        d.append(draw.Lines(*pp, close=False, fill="none", stroke=traj_c, stroke_width=1.5,
                            stroke_opacity=0.85, stroke_linejoin="round"))
        d.append(draw.Circle(_mx(box, tr[-1][0]), _my(box, tr[-1][1]), 3.4, fill=traj_c))
    d.append(draw.Circle(cx0, cy0, 3.6, fill="white", stroke=INK, stroke_width=1.6))   # root

    d.append(draw.Text(title, 16, bx + BS / 2, by - 16, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))


def R_model(R):
    return MultivariateBrownian(R)


def render(mode):
    global MODE
    MODE = mode
    if mode == "bw":
        pt_c, ell_c, traj_c = "#b4b4b4", INK, "#555555"
    else:
        pt_c, ell_c, traj_c = "#8aa9c8", "#c0436a", "#2f5c86"

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Multivariate Brownian motion — traits evolve jointly", 20, 40, 40,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("the rate matrix R sets each trait's variance and their covariance;  "
                       "dots = endpoints, ring = 95% ellipse, lines = sample paths", 13.5, 40, 64,
                       font_family=FONT, text_anchor="start", fill="#777"))

    panel(d, BOXES[0], np.array([[1.0, 0.75], [0.75, 1.0]]), "Correlated  (ρ = 0.75)",
          pt_c, ell_c, traj_c, np.random.default_rng(3))
    panel(d, BOXES[1], np.array([[1.0, 0.0], [0.0, 1.0]]), "Independent  (ρ = 0)",
          pt_c, ell_c, traj_c, np.random.default_rng(5))

    d.append(draw.Text("R = [[1, 0.75], [0.75, 1]]", 13, BOXES[0][0] + BS / 2, 570, font_family=FONT,
                       text_anchor="middle", fill="#555"))
    d.append(draw.Text("R = [[1, 0], [0, 1]]", 13, BOXES[1][0] + BS / 2, 570, font_family=FONT,
                       text_anchor="middle", fill="#555"))

    name = "trait_multivariate" if mode == "color" else "trait_multivariate_bw"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)


def main():
    for mode in ("color", "bw"):
        render(mode)
    print("wrote trait_multivariate (+_bw)")


if __name__ == "__main__":
    main()
