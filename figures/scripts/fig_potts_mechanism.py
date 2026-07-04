"""Figure: how the Potts gene-family coupling model works (didactic, 4 steps).

ZOMBI2's ``PottsRates`` makes gene families NON-independent. Each genome is a presence/absence
vector over a fixed panel; pairwise couplings J make families gain/lose together. The coupling
enters through the LOSS rate — a present family i is lost at ``base_loss·exp(-β·f_i)`` with the
local field ``f_i = h_i + Σ_j J_ij·σ_j`` — while gain is field-blind HGT whose gene is then
selectively retained. This plate walks through it in four steps.

Colour, didactic.  Run:  python figures/scripts/fig_potts_mechanism.py
"""

from __future__ import annotations

import math
from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK

OUT = Path(__file__).resolve().parent.parent / "potts_mechanism"

W, H = 1240, 860
MODA, MODB = "#4477AA", "#E08A3C"                 # two pathway modules
POS, NEG = "#2f8f4e", "#cc4b3c"                   # J>0 co-occur, J<0 avoid
FAINT = "#dfe4ea"
COL = {0: MODA, 1: MODA, 2: MODA, 3: MODB, 4: MODB, 5: MODB}


def header(d, x, y, num, text):
    d.append(draw.Text(f"{num}", 19, x, y, font_family=FONT, text_anchor="start",
                       font_weight="bold", fill=INK))
    d.append(draw.Text(text, 16, x + 26, y, font_family=FONT, text_anchor="start",
                       font_weight="bold", fill=INK))


def cell(d, x, y, present, color, s=38):
    if present:
        d.append(draw.Rectangle(x, y, s, s, fill=color, stroke=INK, stroke_width=1.4))
    else:
        d.append(draw.Rectangle(x, y, s, s, fill="white", stroke="#bfc6cf", stroke_width=1.4))
    return s


def check(d, x, y, color):
    d.append(draw.Lines(x - 6, y + 1, x - 1, y + 7, x + 8, y - 6, close=False, fill="none",
                        stroke=color, stroke_width=3, stroke_linecap="round", stroke_linejoin="round"))


def cross(d, x, y, color):
    for a, b in (((-5, -5), (5, 5)), ((-5, 5), (5, -5))):
        d.append(draw.Line(x + a[0], y + a[1], x + b[0], y + b[1], stroke=color, stroke_width=3,
                           stroke_linecap="round"))


# --------------------------------------------------------------------------- panel 1
def panel_genome(d):
    header(d, 56, 130, "1", "A genome = which families it has  (present / absent vector σ)")
    sigma = [1, 1, 0, 1, 0, 0]
    x0, y = 96, 168
    for i, on in enumerate(sigma):
        x = x0 + i * 50
        cell(d, x, y, on, COL[i])
        d.append(draw.Text(f"F{i}", 13, x + 19, y + 58, font_family=FONT, text_anchor="middle", fill=INK))
        d.append(draw.Text(str(on), 13, x + 19, y + 78, font_family=FONT, text_anchor="middle",
                           font_weight="bold", fill="#888"))
    d.append(draw.Text("present", 12.5, x0 + 300 + 26, y + 12, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    cell(d, x0 + 300 + 8, y, 1, "#9aa7b4", s=16)
    d.append(draw.Rectangle(x0 + 300 + 8, y + 28, 16, 16, fill="white", stroke="#bfc6cf", stroke_width=1.4))
    d.append(draw.Text("absent", 12.5, x0 + 300 + 26, y + 36, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Text("a fixed panel of N gene families; colour = pathway module", 12.5, 96, 290,
                       font_family=FONT, text_anchor="start", fill="#777"))


# --------------------------------------------------------------------------- panel 2
def panel_coupling(d):
    header(d, 640, 130, "2", "Families are coupled by J")
    # horizontal legend (top, clear of the nodes)
    ly = 166
    d.append(draw.Line(660, ly, 686, ly, stroke=POS, stroke_width=3))
    d.append(draw.Text("J > 0  (partners co-occur)", 12.5, 692, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(940, ly, 966, ly, stroke=NEG, stroke_width=2.6, stroke_dasharray="6,4"))
    d.append(draw.Text("J < 0  (partners avoid)", 12.5, 972, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    pos = {0: (760, 216), 1: (720, 272), 2: (810, 270),
           3: (1030, 216), 4: (990, 272), 5: (1080, 270)}
    for a, b in [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)]:
        d.append(draw.Line(*pos[a], *pos[b], stroke=POS, stroke_width=3))
    d.append(draw.Line(*pos[2], *pos[4], stroke=NEG, stroke_width=2.6, stroke_dasharray="6,4"))
    for i, (x, y) in pos.items():
        d.append(draw.Circle(x, y, 17, fill=COL[i], stroke=INK, stroke_width=1.5))
        d.append(draw.Text(f"F{i}", 12, x, y, font_family=FONT, text_anchor="middle",
                           dominant_baseline="central", fill="white", font_weight="bold"))
    d.append(draw.Text("module A", 12.5, 765, 306, font_family=FONT, text_anchor="middle", fill=MODA, font_weight="bold"))
    d.append(draw.Text("module B", 12.5, 1035, 306, font_family=FONT, text_anchor="middle", fill=MODB, font_weight="bold"))


# --------------------------------------------------------------------------- panel 3
def panel_loss(d):
    header(d, 56, 372, "3", "Coupling enters through the LOSS rate")
    d.append(draw.Text("local field   f_i  =  h_i  +  Σ_j  J_ij · σ_j        (sum over present partners)",
                       14, 96, 404, font_family=FONT, text_anchor="start", fill=INK))
    d.append(draw.Text("loss rate   loss_i  =  base_loss · e^(−β · f_i)", 14, 96, 430, font_family=FONT,
                       text_anchor="start", fill=INK))

    # loss vs field curve
    bx, by, bw, bh = 132, 462, 356, 232
    fmin, fmax, lmax = -2.0, 3.0, 4.0

    def X(f):
        return bx + (f - fmin) / (fmax - fmin) * bw

    def Y(v):
        return by + bh - min(v, lmax) / lmax * bh
    d.append(draw.Line(bx, by, bx, by + bh, stroke="#bdbdbd", stroke_width=1.2))
    d.append(draw.Line(bx, by + bh, bx + bw, by + bh, stroke="#bdbdbd", stroke_width=1.2))
    d.append(draw.Line(X(0), by, X(0), by + bh, stroke="#dcdcdc", stroke_width=1, stroke_dasharray="3,3"))
    pts = []
    f = fmin
    while f <= fmax + 1e-9:
        pts += [X(f), Y(math.exp(-f))]
        f += 0.05
    d.append(draw.Lines(*pts, close=False, fill="none", stroke=INK, stroke_width=2.4, stroke_linejoin="round"))
    d.append(draw.Text("loss rate", 12.5, bx - 28, by + bh / 2, font_family=FONT, text_anchor="middle",
                       fill="#777", transform=f"rotate(-90 {bx - 28} {by + bh / 2})"))
    d.append(draw.Text("local field  f", 12.5, bx + bw / 2, by + bh + 24, font_family=FONT,
                       text_anchor="middle", fill="#777"))
    d.append(draw.Circle(X(0), Y(1), 3, fill="#999"))
    d.append(draw.Text("base_loss", 11, X(0) + 6, Y(1) - 7, font_family=FONT, text_anchor="start", fill="#999"))
    # two regimes
    d.append(draw.Circle(X(-1.0), Y(math.exp(1)), 5, fill=NEG, stroke=INK, stroke_width=1.2))
    d.append(draw.Text("no partners: fast loss", 12, X(-1.0) + 12, Y(math.exp(1)), font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=NEG, font_weight="bold"))
    d.append(draw.Circle(X(2.0), Y(math.exp(-2)), 5, fill=POS, stroke=INK, stroke_width=1.2))
    d.append(draw.Text("partners present: protected", 12, X(2.0), Y(math.exp(-2)) - 16, font_family=FONT,
                       text_anchor="middle", fill=POS, font_weight="bold"))
    d.append(draw.Text("more present + partners  ->  lower loss, so coupled families stay together",
                       12.5, 96, 776, font_family=FONT, text_anchor="start", fill="#777"))


# --------------------------------------------------------------------------- panel 4
def panel_hgt(d):
    header(d, 640, 372, "4", "Gain = HGT (blind); retention = selective")

    def scene(y0, present, label, kept):
        # genome F0 [i=F1 arrives] F2
        x0 = 700
        cell(d, x0, y0, present, MODA); d.append(draw.Text("F0", 12, x0 + 19, y0 + 54, font_family=FONT, text_anchor="middle", fill=INK))
        xm = x0 + 60
        if kept:
            cell(d, xm, y0, True, MODA)
        else:
            d.append(draw.Rectangle(xm, y0, 38, 38, fill=MODA, fill_opacity=0.28, stroke="#bfc6cf",
                                    stroke_width=1.4, stroke_dasharray="4,3"))
        d.append(draw.Text("F1", 12, xm + 19, y0 + 54, font_family=FONT, text_anchor="middle", fill=INK))
        x2 = x0 + 120
        cell(d, x2, y0, present, MODA); d.append(draw.Text("F2", 12, x2 + 19, y0 + 54, font_family=FONT, text_anchor="middle", fill=INK))
        # HGT arrow into F1 from above
        d.append(draw.Line(xm + 19, y0 - 34, xm + 19, y0 - 4, stroke=INK, stroke_width=2))
        d.append(draw.Lines(xm + 19, y0 - 2, xm + 13, y0 - 12, xm + 25, y0 - 12, close=True, fill=INK))
        d.append(draw.Text("HGT", 11, xm + 42, y0 - 22, font_family=FONT, text_anchor="start", fill="#777"))
        # outcome
        if kept:
            check(d, x2 + 66, y0 + 19, POS)
            d.append(draw.Text("kept", 15, x2 + 82, y0 + 19, font_family=FONT, text_anchor="start",
                               dominant_baseline="central", font_weight="bold", fill=POS))
        else:
            cross(d, x2 + 66, y0 + 19, NEG)
            d.append(draw.Text("purged", 15, x2 + 82, y0 + 19, font_family=FONT, text_anchor="start",
                               dominant_baseline="central", font_weight="bold", fill=NEG))
        d.append(draw.Text(label, 12, 700, y0 + 74, font_family=FONT, text_anchor="start", fill="#555"))

    scene(500, True, "partners F0, F2 present:  the transferred F1 is retained", True)
    scene(650, False, "partners F0, F2 absent:  F1 has no protection, lost fast", False)
    d.append(draw.Text("differential retention of transferred genes writes J into the profiles",
                       12.5, 640, 760, font_family=FONT, text_anchor="start", fill="#777"))


def main():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The Potts model — gene families that gain and lose together", 21, 40, 44,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("non-independent gene families: pairwise couplings J make partners co-occur "
                       "(+) or avoid each other (−).  ZOMBI2: PottsRates / simulate_coupled", 13.5, 40, 70,
                       font_family=FONT, text_anchor="start", fill="#777"))
    d.append(draw.Line(620, 100, 620, 812, stroke="#e8e8e8", stroke_width=1.2))
    d.append(draw.Line(48, 330, 1192, 330, stroke="#e8e8e8", stroke_width=1.2))
    panel_genome(d)
    panel_coupling(d)
    panel_loss(d)
    panel_hgt(d)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "potts_mechanism.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "potts_mechanism.png"), scale=300 / 72.0)
    print("wrote potts_mechanism")


if __name__ == "__main__":
    main()
