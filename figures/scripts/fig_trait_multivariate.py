"""Figure: multivariate Brownian motion — two traits evolving jointly.

``MultivariateBrownian`` diffuses a vector-valued trait with a rate (covariance) matrix R; the
off-diagonal of R makes the traits evolve in a correlated way. This figure plots the joint (X, Y)
plane: a cloud of endpoints (each stepped through the model's own ``evolve()``), the 95% covariance
ellipse, and a few sample trajectories — for a correlated R (tilted ellipse) versus an independent
R (upright circle).

House style: B&W, one centered title, ASCII text.

Run:  python figures/scripts/fig_trait_multivariate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2.traits import MultivariateBrownian

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 640
T, NPTS, NTRAJ, TSTEPS = 1.0, 280, 3, 140
DR = (-3.3, 3.3)
BS = 430                                                     # box side
BOXES = [(120, 150), (650, 150)]
PT_C, ELL_C, TRAJ_C = "#b4b4b4", INK, "#555555"              # grey dots, black ellipse, dark paths


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


def panel(d, box, R, title, rng):
    bx, by = box
    cx0, cy0 = _mx(box, 0), _my(box, 0)
    d.append(draw.Rectangle(bx, by, BS, BS, fill="none", stroke="#cccccc", stroke_width=1.2))
    d.append(draw.Line(bx, cy0, bx + BS, cy0, stroke="#d9d9d9", stroke_width=1))
    d.append(draw.Line(cx0, by, cx0, by + BS, stroke="#d9d9d9", stroke_width=1))
    d.append(draw.Text("trait X", FS_TICK, bx + BS - 6, cy0 - 10, font_family=FONT,
                       text_anchor="end", fill=MUTED))
    d.append(draw.Text("trait Y", FS_TICK, cx0 + 10, by + 18, font_family=FONT,
                       text_anchor="start", fill=MUTED))

    for (x, y) in endpoints(MultivariateBrownian(R), rng, NPTS):
        if DR[0] < x < DR[1] and DR[0] < y < DR[1]:
            d.append(draw.Circle(_mx(box, x), _my(box, y), 3.0, fill=PT_C, fill_opacity=0.55))

    pts = []
    for (x, y) in ellipse_pts(R):
        pts += [_mx(box, x), _my(box, y)]
    d.append(draw.Lines(*pts, close=True, fill="none", stroke=ELL_C, stroke_width=2.6))

    for _ in range(NTRAJ):
        tr = trajectory(MultivariateBrownian(R), rng)
        pp = []
        for (x, y) in tr:
            pp += [_mx(box, max(DR[0], min(DR[1], x))), _my(box, max(DR[0], min(DR[1], y)))]
        d.append(draw.Lines(*pp, close=False, fill="none", stroke=TRAJ_C, stroke_width=1.7,
                            stroke_opacity=0.85, stroke_linejoin="round"))
        d.append(draw.Circle(_mx(box, tr[-1][0]), _my(box, tr[-1][1]), 4.0, fill=TRAJ_C))
    d.append(draw.Circle(cx0, cy0, 4.0, fill="white", stroke=INK, stroke_width=1.8))

    d.append(draw.Text(title, FS_LABEL, bx + BS / 2, by - 18, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Multivariate Brownian motion", FS_TITLE, W / 2, 46, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    panel(d, BOXES[0], np.array([[1.0, 0.75], [0.75, 1.0]]), "Correlated  (rho = 0.75)",
          np.random.default_rng(3))
    panel(d, BOXES[1], np.array([[1.0, 0.0], [0.0, 1.0]]), "Independent  (rho = 0)",
          np.random.default_rng(5))
    d.append(draw.Text("R = [[1, 0.75], [0.75, 1]]", FS_TICK, BOXES[0][0] + BS / 2, 602,
                       font_family=FONT, text_anchor="middle", fill=MUTED))
    d.append(draw.Text("R = [[1, 0], [0, 1]]", FS_TICK, BOXES[1][0] + BS / 2, 602,
                       font_family=FONT, text_anchor="middle", fill=MUTED))
    d.append(draw.Text("dots = endpoints,  ring = 95% ellipse,  lines = sample paths", FS_TICK,
                       W / 2, 632, font_family=FONT, text_anchor="middle", fill=MUTED))

    name = "trait_multivariate"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
