"""Figure: pathway modules, and a genome that loses an incomplete module (Potts, idea 1).

In ZOMBI2's coupled ("Potts") gene-family model a genome is a present/absent vector over a fixed
panel of gene families. The families are organised into pathway MODULES that are internally
coupled (J > 0), so members protect one another from loss.

  * Panel A - the panel: 20 families in 4 modules of 5, numbered 1..5 within each module.
  * Panel B - one genome, drawn UNORDERED (a bag of present families). Three modules are complete
    and one (module 4) is nearly gone: its single surviving gene has no present partners, so its
    local field is small, its loss rate is high, and it is lost -- the module empties out.

The "incomplete module is lost" claim is exactly the model's loss law
loss_i = base_loss * exp(-beta * f_i) with f_i = h_i + sum over present partners J_ij (see
zombi2/coupling.py): fewer present partners -> smaller field -> faster loss.

House style: categorical module colour (STYLE.md exception), one centered bold title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_potts_genome.py
"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))   # zombi_style

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, MODULE_COLORS, AVOID, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent
NAME = "potts_genome"

W, H = 1200, 660
R = 20                                              # family dot radius


def dot(d, x, y, color, label=None, r=R, faded=False):
    kw = dict(fill=color, stroke=INK, stroke_width=2.0)
    if faded:
        kw = dict(fill=color, fill_opacity=0.30, stroke="#bfc6cf", stroke_width=2.0)
    d.append(draw.Circle(x, y, r, **kw))
    if label is not None:
        d.append(draw.Text(label, FS_TICK, x, y, font_family=FONT, text_anchor="middle",
                           dominant_baseline="central", fill="white", font_weight="bold"))


def module_group(d, gx, gy, mi):
    """Panel A: one module = faint container + 5 numbered dots + a colour label."""
    color = MODULE_COLORS[mi]
    d.append(draw.Rectangle(gx - 108, gy - 72, 216, 150, rx=14, fill="#f7f7f7",
                            stroke="#e4e4e4", stroke_width=1.2))
    d.append(draw.Text(f"module {mi + 1}", FS_TICK, gx, gy - 50, font_family=FONT,
                       text_anchor="middle", fill=color, font_weight="bold"))
    top = [(gx - 58, gy - 2), (gx, gy - 2), (gx + 58, gy - 2)]
    bot = [(gx - 29, gy + 46), (gx + 29, gy + 46)]
    for k, (x, y) in enumerate(top + bot):
        dot(d, x, y, color, label=str(k + 1))


def cross(d, x, y, color, r=15):
    for a, b in (((-r, -r), (r, r)), ((-r, r), (r, -r))):
        d.append(draw.Line(x + a[0], y + a[1], x + b[0], y + b[1], stroke=color,
                           stroke_width=5, stroke_linecap="round"))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Modules are kept or lost together", FS_TITLE, W / 2, 48,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # divider between the two panels
    d.append(draw.Line(600, 92, 600, 590, stroke="#e0e0e0", stroke_width=1.4))

    # --- Panel A: the family panel, 4 modules of 5 ---
    d.append(draw.Text("A", FS_LABEL, 44, 108, font_family=FONT, text_anchor="start",
                       font_weight="bold", fill=INK))
    d.append(draw.Text("the panel: 20 families in 4 modules", FS_ANNOT, 78, 108, font_family=FONT,
                       text_anchor="start", fill=MUTED))
    for mi, (gx, gy) in enumerate([(175, 250), (445, 250), (175, 470), (445, 470)]):
        module_group(d, gx, gy, mi)

    # --- Panel B: one genome, drawn unordered ---
    d.append(draw.Text("B", FS_LABEL, 636, 108, font_family=FONT, text_anchor="start",
                       font_weight="bold", fill=INK))
    d.append(draw.Text("one genome, drawn unordered", FS_ANNOT, 670, 108, font_family=FONT,
                       text_anchor="start", fill=MUTED))

    gx0, gy0, gw, gh = 650, 150, 474, 356
    d.append(draw.Rectangle(gx0, gy0, gw, gh, rx=16, fill="#f7f7f7", stroke="#d7d7d7",
                            stroke_width=1.6))

    # present families: modules 1-3 complete (5 each), module 4 nearly gone (1 of 5).
    cols = [732, 826, 920, 1014]
    rows = [242, 322, 402, 472]
    grid = [[0, 1, 2, 0],
            [1, 2, 0, 1],
            [2, 0, 1, 3],       # the single teal (module 4) survivor at row 2, col 3
            [2, 1, 0, 2]]
    jit = [[(-6, 4), (7, -5), (-4, 6), (5, 5)],
           [(6, -6), (-7, 4), (4, 7), (-5, -4)],
           [(-6, -5), (5, 6), (-4, -6), (6, 4)],
           [(4, -5), (-6, 5), (7, -4), (-5, 6)]]
    lone = None
    for r, row in enumerate(grid):
        for c, mi in enumerate(row):
            x, y = cols[c] + jit[r][c][0], rows[r] + jit[r][c][1]
            dot(d, x, y, MODULE_COLORS[mi], faded=(mi == 3))
            if mi == 3:
                lone = (x, y)

    # flag the lone survivor of the incomplete module as lost
    cross(d, lone[0], lone[1], AVOID, r=14)
    d.append(draw.Line(lone[0] + 18, lone[1] + 10, lone[0] + 44, lone[1] + 40,
                       stroke=MUTED, stroke_width=1.4))
    d.append(draw.Text("lost", FS_TICK, lone[0] + 48, lone[1] + 50, font_family=FONT,
                       text_anchor="start", fill=AVOID, font_weight="bold"))

    d.append(draw.Text("module 4 is incomplete (1 of 5): its lone gene has no partners,",
                       FS_TICK, gx0 + gw / 2, gy0 + gh + 34, font_family=FONT,
                       text_anchor="middle", fill=MUTED))
    d.append(draw.Text("so a small field, a fast loss -- the module empties out",
                       FS_TICK, gx0 + gw / 2, gy0 + gh + 58, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    # shared module colour key, centered along the bottom
    ky = 628
    labels = [f"module {i + 1}" for i in range(4)]
    xs = [230, 430, 630, 830]
    for x, lab, col in zip(xs, labels, MODULE_COLORS):
        d.append(draw.Circle(x, ky, 11, fill=col, stroke=INK, stroke_width=1.4))
        d.append(draw.Text(lab, FS_TICK, x + 20, ky, font_family=FONT, text_anchor="start",
                           dominant_baseline="central", fill=INK))

    out = OUT_DIR / NAME
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{NAME}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{NAME}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{NAME}.svg / .png")


if __name__ == "__main__":
    render()
