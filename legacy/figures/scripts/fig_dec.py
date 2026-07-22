"""Figure: the DEC model of geographic-range evolution (Ree & Smith 2008).

A species' trait is its geographic RANGE — a non-empty subset of the areas. DEC evolves it
with two processes, both shown here:

* ANAGENETIC (along a branch): a continuous-time Markov chain over ranges — dispersal gains
  an area, local extinction loses one. Drawn as the range lattice (single areas at the
  bottom, unions above); moving up = dispersal, down = extinction.
* CLADOGENETIC (at a speciation): the ancestral range is inherited by the two daughters by
  narrow sympatry, subset sympatry or vicariance.

Colour. Run:  python figures/scripts/fig_dec.py
"""

from __future__ import annotations

import math
from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1210, 720
# Areas are categorical identities, so they are the one place this figure uses colour
# (Adrian's call): a distinct, colour-blind-safe hue per area (Paul Tol 'bright').
AREA_COLOR = {"A": "#4477AA", "B": "#EE6677", "C": "#228833"}
CELL = 26
FS_AREA = 15                                             # small area letter inside a range cell
DISP = "#2f8f4e"                                          # dispersal (gain) arrowhead — green
EXT = "#cc4b3c"                                           # local-extinction (loss) arrowhead — red


def _acol(a):
    return AREA_COLOR[a]


def _lum(h):
    return 0.299 * int(h[1:3], 16) + 0.587 * int(h[3:5], 16) + 0.114 * int(h[5:7], 16)


def range_chip(d, cx, cy, areas):
    """A geographic range = a row of coloured (or grey) area cells."""
    areas = sorted(areas)
    total = len(areas) * CELL
    x0 = cx - total / 2
    for i, a in enumerate(areas):
        x, fill = x0 + i * CELL, _acol(a)
        d.append(draw.Rectangle(x, cy - CELL / 2, CELL, CELL, fill=fill, stroke=INK, stroke_width=1.3))
        # smaller area letter, centred in its cell (explicit y offset = reliable vertical centring)
        d.append(draw.Text(a, FS_AREA, x + CELL / 2, cy + 0.34 * FS_AREA, font_family=FONT,
                           text_anchor="middle", fill="white" if _lum(fill) < 150 else INK,
                           font_weight="bold"))
    d.append(draw.Rectangle(x0, cy - CELL / 2, total, CELL, fill="none", stroke=INK, stroke_width=1.3))
    return total


def _head(d, x, y, ang, color, size=8.5):
    d.append(draw.Lines(x, y, x - size * math.cos(ang - 0.42), y - size * math.sin(ang - 0.42),
                        x - size * math.cos(ang + 0.42), y - size * math.sin(ang + 0.42),
                        close=True, fill=color))


def lattice_edge(d, lower, upper):
    """A rung of the range lattice: dispersal (green, up) and extinction (red, down)."""
    lx, ly = lower
    ux, uy = upper
    a = math.atan2(uy - ly, ux - lx)
    sh = 20
    sx, sy = lx + sh * math.cos(a), ly + sh * math.sin(a)
    ex, ey = ux - sh * math.cos(a), uy - sh * math.sin(a)
    d.append(draw.Line(sx, sy, ex, ey, stroke="#b7b7b7", stroke_width=1.8))
    _head(d, ex, ey, a, DISP)                            # into the larger range = dispersal
    _head(d, sx, sy, a + math.pi, EXT)                   # into the smaller range = extinction


# --------------------------------------------------------------------------- panel A: anagenesis
def panel_a(d):
    d.append(draw.Text("A   along a branch: dispersal & extinction", FS_LABEL, 60, 150,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    yl0, yl1, yl2 = 580, 420, 260
    pos = {
        "A": (170, yl0), "B": (330, yl0), "C": (490, yl0),
        "AB": (170, yl1), "AC": (330, yl1), "BC": (490, yl1),
        "ABC": (330, yl2),
    }
    edges = [("A", "AB"), ("A", "AC"), ("B", "AB"), ("B", "BC"), ("C", "AC"), ("C", "BC"),
             ("AB", "ABC"), ("AC", "ABC"), ("BC", "ABC")]
    for lo, up in edges:
        lattice_edge(d, pos[lo], pos[up])
    for name, (x, y) in pos.items():
        range_chip(d, x, y, name)

    # legend (two stacked rows, kept within panel A)
    lx, ly = 60, 648
    d.append(draw.Line(lx, ly, lx + 26, ly, stroke="#b7b7b7", stroke_width=1.8))
    _head(d, lx + 26, ly, 0.0, DISP)
    d.append(draw.Text("dispersal: gain an area (rate d)", FS_TICK, lx + 36, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    ly2 = ly + 30
    d.append(draw.Line(lx, ly2, lx + 26, ly2, stroke="#b7b7b7", stroke_width=1.8))
    _head(d, lx, ly2, math.pi, EXT)
    d.append(draw.Text("local extinction: lose an area (rate e)", FS_TICK, lx + 36, ly2,
                       font_family=FONT, text_anchor="start", dominant_baseline="central", fill=INK))


# --------------------------------------------------------------------------- panel B: cladogenesis
def clado(d, cx, cy, anc, d1, d2, name, desc):
    range_chip(d, cx, cy, anc)
    ybar, yd = cy + 30, cy + 70
    x1, x2 = cx - 66, cx + 66
    d.append(draw.Line(cx, cy + CELL / 2, cx, ybar, stroke=INK, stroke_width=2))
    d.append(draw.Line(x1, ybar, x2, ybar, stroke=INK, stroke_width=2))
    d.append(draw.Line(x1, ybar, x1, yd - CELL / 2, stroke=INK, stroke_width=2))
    d.append(draw.Line(x2, ybar, x2, yd - CELL / 2, stroke=INK, stroke_width=2))
    range_chip(d, x1, yd, d1)
    range_chip(d, x2, yd, d2)
    d.append(draw.Text(name, FS_ANNOT, cx + 160, cy - 8, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", font_weight="bold", fill=INK))
    d.append(draw.Text(desc, FS_TICK, cx + 160, cy + 16, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=MUTED))


def panel_b(d):
    x0 = 700
    d.append(draw.Line(x0 - 30, 130, x0 - 30, 690, stroke="#e2e2e2", stroke_width=1.2))
    d.append(draw.Text("B   at a speciation: cladogenesis", FS_LABEL, x0, 150, font_family=FONT,
                       text_anchor="start", font_weight="bold", fill=INK))
    cx = x0 + 70
    clado(d, cx, 250, "A", "A", "A", "Narrow sympatry",
          "both daughters keep the range")
    clado(d, cx, 430, "AB", "A", "AB", "Subset sympatry",
          "one area, and the full range")
    clado(d, cx, 600, "AB", "A", "B", "Vicariance",
          "range splits into complements")


# --------------------------------------------------------------------------- render
def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The DEC model of geographic-range evolution", FS_TITLE, W / 2, 52,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    panel_a(d)
    panel_b(d)
    name = "dec"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)


def main():
    render()
    print("wrote dec")


if __name__ == "__main__":
    main()
