"""Figure: it's all Gillespie -- the same engine at every level of ZOMBI2.

Each level of a ZOMBI2 simulation supplies a different *bag of events* with their own
rates, but they are all realised by the one loop on the right. Species trees race
speciation against extinction; gene families race duplication, transfer, loss and
origination; a discrete trait races its between-state transitions. Only sequence
substitution along a fixed branch steps outside the loop (a single matrix operation).

Event markers follow the manual convention: filled square = duplication, black arc =
transfer, cross = loss / extinction, ring = origination.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_gillespie_everywhere.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, ACCENT, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1300, 760


def text(d, s, x, y, size, *, anchor="middle", fill=INK, weight="normal", italic=False):
    d.append(draw.Text(s, size, x, y, font_family=FONT, text_anchor=anchor,
                       dominant_baseline="central", fill=fill, font_weight=weight,
                       font_style=("italic" if italic else "normal")))


def arrowhead(d, hx, hy, ax, ay, col=INK, ah=13.0):
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy,
                        hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


# ---- event glyphs (small, centred on (x, y)) --------------------------------
def g_speciation(d, x, y):
    d.append(draw.Line(x, y + 11, x, y - 2, stroke=ACCENT["speciation"], stroke_width=3.0))
    d.append(draw.Line(x, y - 2, x - 11, y - 13, stroke=ACCENT["speciation"], stroke_width=3.0))
    d.append(draw.Line(x, y - 2, x + 11, y - 13, stroke=ACCENT["speciation"], stroke_width=3.0))


def g_cross(d, x, y, col=INK):
    r = 9
    d.append(draw.Line(x - r, y - r, x + r, y + r, stroke=col, stroke_width=3.0))
    d.append(draw.Line(x - r, y + r, x + r, y - r, stroke=col, stroke_width=3.0))


def g_square(d, x, y):
    s = 16
    d.append(draw.Rectangle(x - s / 2, y - s / 2, s, s, fill=INK))


def g_transfer(d, x, y):
    p = draw.Path(fill="none", stroke=INK, stroke_width=3.0)
    p.M(x - 12, y + 8).Q(x, y - 20, x + 12, y + 8)
    d.append(p)
    arrowhead(d, x + 12, y + 8, x + 4, y - 4, INK, ah=9.0)


def g_ring(d, x, y):
    d.append(draw.Circle(x, y, 8, fill="none", stroke=INK, stroke_width=3.0))


def g_state(d, x, y):
    d.append(draw.Rectangle(x - 22, y - 8, 15, 15, fill="#dddddd", stroke=INK, stroke_width=1.6))
    d.append(draw.Rectangle(x + 8, y - 8, 15, 15, fill=INK))
    d.append(draw.Line(x - 5, y, x + 6, y, stroke=INK, stroke_width=2.4))
    arrowhead(d, x + 6, y, x - 3, y, INK, ah=7.0)


LEVELS = [
    ("Species tree", "rates: lambda, mu  per lineage",
     [(g_speciation, "speciation"), (g_cross, "extinction")]),
    ("Gene families", "rates: duplication, transfer, loss, origination",
     [(g_square, "duplication"), (g_transfer, "transfer"), (g_cross, "loss"),
      (g_ring, "origination")]),
    ("Discrete trait (Mk)", "rates: the entries of the Q matrix",
     [(g_state, "state change")]),
]


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    text(d, "One engine, many events: it's all Gillespie", W / 2, 42, FS_TITLE,
         weight="bold")

    card_x, card_w = 50, 700
    card_h = 150
    ys = [150, 340, 530]
    ys = [y + 30 for y in ys]

    for (title, rates, events), cy in zip(LEVELS, ys):
        d.append(draw.Rectangle(card_x, cy - card_h / 2, card_w, card_h, rx=16, ry=16,
                                fill="white", stroke=INK, stroke_width=2.4))
        text(d, title, card_x + 26, cy - card_h / 2 + 34, FS_LABEL, anchor="start",
             weight="bold")
        text(d, rates, card_x + 26, cy - card_h / 2 + 66, FS_TICK, anchor="start",
             fill=MUTED)
        # event glyphs in a row across the lower half of the card
        gx = card_x + 70
        spacing = min(180, (card_w - 120) / max(1, len(events)))
        gx = card_x + 90
        for glyph, label in events:
            glyph(d, gx, cy + 34)
            text(d, label, gx, cy + 62, FS_TICK)
            gx += spacing

    # ---- the engine box on the right ----------------------------------------
    ex0, ew = 880, 380
    ey0, eh = 150, 460
    d.append(draw.Rectangle(ex0, ey0, ew, eh, rx=18, ry=18, fill="#f2f2f2", stroke=INK,
                            stroke_width=3.0))
    ecx = ex0 + ew / 2
    text(d, "one Gillespie loop", ecx, ey0 + 46, FS_LABEL, weight="bold")
    step_rows = [
        "1.  total rate R = sum",
        "     of every event's rate",
        "",
        "2.  when:  draw",
        "     dt ~ Exponential(R)",
        "",
        "3.  what:  pick event i",
        "     with prob. (rate)/R",
        "",
        "4.  apply it, then repeat",
    ]
    ry = ey0 + 108
    for row in step_rows:
        weight = "bold" if row[:2] in ("1.", "2.", "3.", "4.") else "normal"
        text(d, row, ex0 + 34, ry, FS_ANNOT, anchor="start",
             fill=(INK if weight == "bold" else MUTED))
        ry += 34

    # arrows from each card into the engine box
    for cy in ys:
        x1 = card_x + card_w
        d.append(draw.Line(x1, cy, ex0 - 3, ey0 + eh / 2, stroke=MUTED, stroke_width=2.4))
        arrowhead(d, ex0 - 3, ey0 + eh / 2, x1, cy, MUTED, ah=12.0)

    # ---- footnote: the one exception ----------------------------------------
    text(d, "The one exception: substituting a sequence along a fixed branch needs only the "
            "branch endpoints,", W / 2, H - 54, FS_TICK, fill=MUTED, italic=True)
    text(d, "so ZOMBI2 does it in a single matrix step, P(t) = exp(Q t), instead of "
            "event-by-event.", W / 2, H - 28, FS_TICK, fill=MUTED, italic=True)

    name = "gillespie_everywhere"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
