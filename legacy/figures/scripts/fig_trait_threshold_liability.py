"""Figure: the threshold model — a continuous liability read out as a discrete state.

Felsenstein's threshold model bridges the continuous and discrete trait families: an
unobserved **liability** diffuses by Brownian motion, and the observed discrete state is
simply which side of a threshold the liability is on. This schematic shows one liability
path (stepped through the model's own ``evolve()``), the threshold, and the discrete state
track it produces — flipping exactly at each crossing.

Rendered in colour and B&W.  Run:  python figures/scripts/fig_trait_threshold_liability.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2 import ThresholdModel

from zombi_style import FONT, INK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 940, 520
T, STEPS = 1.0, 320
PX0, PY0, PW, PH = 96, 132, 760, 250
THR = 0.0
YR = (-2.6, 2.6)

COLOR = {0: "#4477AA", 1: "#EE6677"}          # below / above
BAND_C = {0: "#e6eef6", 1: "#fbe4e8"}
BW = {0: "#4a4a4a", 1: "#9a9a9a"}
BAND_BW = {0: "#ededed", 1: "#f7f7f7"}
MODE = "color"


def scol(s):
    return (BW if MODE == "bw" else COLOR)[s]


def band(s):
    return (BAND_BW if MODE == "bw" else BAND_C)[s]


def sim(model, rng):
    dt, x, t, xs = T / STEPS, model.x0, 0.0, [model.x0]
    for _ in range(STEPS):
        x, _ = model.evolve(x, dt, t, rng)
        t += dt
        xs.append(x)
    return np.linspace(0, T, STEPS + 1), np.array(xs)


def _x(t):
    return PX0 + t / T * PW


def _y(y):
    return PY0 + PH - (y - YR[0]) / (YR[1] - YR[0]) * PH


def crossings(ts, xs):
    """Insert threshold-crossing points; return densified (t, x, state) with state per point."""
    out_t, out_x = [ts[0]], [xs[0]]
    for i in range(1, len(ts)):
        x0, x1 = xs[i - 1], xs[i]
        if (x0 - THR) * (x1 - THR) < 0:                        # crossed
            f = (THR - x0) / (x1 - x0)
            out_t.append(ts[i - 1] + f * (ts[i] - ts[i - 1]))
            out_x.append(THR)
        out_t.append(ts[i])
        out_x.append(x1)
    return np.array(out_t), np.array(out_x)


def render(mode):
    global MODE
    MODE = mode
    model = ThresholdModel(thresholds=[THR], sigma2=1.3, x0=-0.9)
    for seed in range(600):                                    # a path with a few well-spread crossings
        rng = np.random.default_rng(seed)
        ts, xs = sim(model, rng)
        cti = ts[np.where((xs[:-1] - THR) * (xs[1:] - THR) < 0)[0]]
        if (2 <= len(cti) <= 3 and (len(cti) < 2 or np.min(np.diff(cti)) > 0.14)
                and xs.min() > YR[0] + 0.25 and xs.max() < YR[1] - 0.25):
            break
    ct, cx = crossings(ts, xs)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The threshold model — a continuous liability read out as a discrete state",
                       19, 40, 40, font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("the liability diffuses by Brownian motion; the observed state = which side "
                       "of the threshold it is on", 13.5, 40, 64, font_family=FONT,
                       text_anchor="start", fill="#777"))

    # tinted bands + threshold
    d.append(draw.Rectangle(PX0, _y(YR[1]), PW, _y(THR) - _y(YR[1]), fill=band(1)))
    d.append(draw.Rectangle(PX0, _y(THR), PW, _y(YR[0]) - _y(THR), fill=band(0)))
    d.append(draw.Rectangle(PX0, PY0, PW, PH, fill="none", stroke="#bdbdbd", stroke_width=1.3))
    d.append(draw.Line(PX0, _y(THR), PX0 + PW, _y(THR), stroke=INK, stroke_width=1.6, stroke_dasharray="6,4"))
    d.append(draw.Text("threshold", 12, PX0 + PW - 4, _y(THR) - 8, font_family=FONT, text_anchor="end",
                       font_weight="bold", fill=INK))

    # liability path, coloured by side, split at crossings
    seg_t, seg_x = [ct[0]], [cx[0]]
    for i in range(1, len(ct)):
        seg_t.append(ct[i])
        seg_x.append(cx[i])
        if cx[i] == THR or i == len(ct) - 1:
            s = 1 if (max(seg_x) + min(seg_x)) / 2 > THR else 0
            mid = 0.5 * (seg_x[0] + seg_x[-1])
            s = 1 if mid > THR else 0
            pts = []
            for t, y in zip(seg_t, seg_x):
                pts += [_x(t), _y(y)]
            d.append(draw.Lines(*pts, close=False, fill="none", stroke=scol(s),
                                stroke_width=2.4, stroke_linejoin="round"))
            seg_t, seg_x = [ct[i]], [cx[i]]

    # crossing dots
    for t, y in zip(ct, cx):
        if y == THR:
            d.append(draw.Circle(_x(t), _y(y), 3.2, fill="white", stroke=INK, stroke_width=1.4))

    d.append(draw.Text("liability", 12.5, PX0 - 26, PY0 + PH / 2, font_family=FONT, text_anchor="middle",
                       fill="#777", transform=f"rotate(-90 {PX0 - 26} {PY0 + PH / 2})"))
    d.append(draw.Lines(PX0 + PW + 4, PY0 + PH, PX0 + PW - 4, PY0 + PH - 4, PX0 + PW - 4, PY0 + PH + 4,
                        close=True, fill=INK))
    d.append(draw.Text("time", 12.5, PX0 + PW / 2, PY0 + PH + 20, font_family=FONT, text_anchor="middle", fill="#777"))

    # discrete-state track below
    ty, th = PY0 + PH + 44, 26
    d.append(draw.Text("observed state", 12.5, PX0 - 10, ty + th / 2, font_family=FONT, text_anchor="end",
                       dominant_baseline="central", fill=INK))
    marks = [0.0] + list(ct[cx == THR]) + [T]
    cur = 0 if xs[0] < THR else 1
    for a, b in zip(marks, marks[1:]):
        xa, xb = _x(a), _x(b)
        fill = scol(cur)
        d.append(draw.Rectangle(xa, ty, xb - xa, th, fill=fill, stroke="white", stroke_width=1.2))
        if xb - xa > 46:
            tc = "white" if (0.299 * int(fill[1:3], 16) + 0.587 * int(fill[3:5], 16)
                             + 0.114 * int(fill[5:7], 16)) < 150 else INK
            d.append(draw.Text("above" if cur else "below", 11, (xa + xb) / 2, ty + th / 2,
                               font_family=FONT, text_anchor="middle", dominant_baseline="central", fill=tc))
        cur = 1 - cur

    # legend
    ly = ty + th + 24
    d.append(draw.Rectangle(PX0, ly - 8, 15, 15, fill=scol(0)))
    d.append(draw.Text("below threshold", 12, PX0 + 20, ly, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Rectangle(PX0 + 150, ly - 8, 15, 15, fill=scol(1)))
    d.append(draw.Text("above threshold", 12, PX0 + 170, ly, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Text("two thresholds give three ordered states, and so on", 12, PX0 + 360, ly,
                       font_family=FONT, text_anchor="start", fill="#999", font_style="italic"))

    name = "trait_threshold_liability" if mode == "color" else "trait_threshold_liability_bw"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)


def main():
    for mode in ("color", "bw"):
        render(mode)
    print("wrote trait_threshold_liability (+_bw)")


if __name__ == "__main__":
    main()
