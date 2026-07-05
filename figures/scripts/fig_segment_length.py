"""Figure: an event's segment length (number of genes) is geometric in `extension`.

Three panels for extension = 0.3, 0.6, 0.8. Each shows the geometric PMF
P(L=k) = (1-ext)*ext**(k-1) with the mean 1/(1-ext) marked. House style.

    python figures/scripts/fig_segment_length.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT = Path(__file__).resolve().parent.parent / "segment_length"

W, H = 1180, 420
EXTS = [0.3, 0.6, 0.8]
KMAX = 12
YMAX = 0.72
BARGREY = "#c9c9c9"

PLOT_TOP, PLOT_H = 116, 224
BASE = PLOT_TOP + PLOT_H
PW, GAP, LEFT = 320, 55, 96


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("How many genes does an event affect?", FS_TITLE, W / 2, 44,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    for p, ext in enumerate(EXTS):
        x0 = LEFT + p * (PW + GAP)
        mean = 1.0 / (1.0 - ext)
        d.append(draw.Text(f"extension = {ext}", FS_LABEL, x0 + PW / 2, 84,
                           font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
        d.append(draw.Line(x0, BASE, x0 + PW, BASE, stroke="#bdbdbd", stroke_width=1.2))
        d.append(draw.Line(x0, PLOT_TOP, x0, BASE, stroke="#bdbdbd", stroke_width=1.2))

        for k in range(1, KMAX + 1):
            prob = (1 - ext) * ext ** (k - 1)
            cx = x0 + (k - 0.5) / KMAX * PW
            bh = prob / YMAX * PLOT_H
            d.append(draw.Rectangle(cx - 9, BASE - bh, 18, bh, fill=BARGREY, stroke=INK, stroke_width=1.0))

        mx = x0 + (mean - 0.5) / KMAX * PW
        d.append(draw.Line(mx, PLOT_TOP - 4, mx, BASE, stroke=INK, stroke_width=1.6, stroke_dasharray="6,4"))
        d.append(draw.Text(f"mean = {mean:.1f}", FS_ANNOT, mx + 7, PLOT_TOP + 10,
                           font_family=FONT, text_anchor="start", fill=INK))

        for k in (1, 4, 8, 12):
            tx = x0 + (k - 0.5) / KMAX * PW
            d.append(draw.Line(tx, BASE, tx, BASE + 5, stroke="#999999", stroke_width=1))
            d.append(draw.Text(str(k), FS_TICK, tx, BASE + 22, font_family=FONT,
                               text_anchor="middle", fill="#777777"))
        d.append(draw.Text("genes in the segment", FS_LABEL, x0 + PW / 2, BASE + 48,
                           font_family=FONT, text_anchor="middle", fill="#555555"))

        if p == 0:
            for yv in (0.0, 0.35, 0.70):
                yy = BASE - yv / YMAX * PLOT_H
                d.append(draw.Text(f"{yv:.2f}", FS_TICK, x0 - 10, yy, font_family=FONT,
                                   text_anchor="end", dominant_baseline="central", fill="#777777"))
            d.append(draw.Text("probability", FS_LABEL, x0 - 56, (PLOT_TOP + BASE) / 2,
                               font_family=FONT, text_anchor="middle", fill="#555555",
                               transform=f"rotate(-90 {x0 - 56} {(PLOT_TOP + BASE) / 2})"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "segment_length.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "segment_length.png"), scale=300 / 72.0)
    print("wrote segment_length")


if __name__ == "__main__":
    render()
