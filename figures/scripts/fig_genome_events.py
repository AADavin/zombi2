"""Figures: the rearrangement events of the ordered-genome model, as before -> after
ring pairs. One figure per event (duplication, loss, inversion, transposition,
transfer), each acting on a 2-3 gene segment, rendered in colour and in B&W.

Run:  python figures/scripts/fig_genome_events.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw

from fig_genome_circular import Gene, color_map, draw_ring, initial_genes
from fig_genome_circular_bw import make_patterns
from zombi_style import FONT, INK

OUT = Path(__file__).resolve().parent.parent / "genome_events"

CW, CH = 1000, 520
LC, RC, RR = (235, 250), (765, 250), 150         # ring centres + radius

BASE = initial_genes(12)


def dup():
    seg = BASE[3:6]                                # families 4,5,6
    return BASE[:6] + seg + BASE[6:], [3, 4, 5], [3, 4, 5, 6, 7, 8]


def loss():
    return BASE[:3] + BASE[6:], [3, 4, 5], []


def inversion():
    inv = [Gene(g.family, -g.orientation) for g in reversed(BASE[3:6])]   # reverse + flip
    return BASE[:3] + inv + BASE[6:], [3, 4, 5], [3, 4, 5]


def transposition():
    seg, rest = BASE[3:6], BASE[:3] + BASE[6:]
    return rest[:6] + seg + rest[6:], [3, 4, 5], [6, 7, 8]


def transfer():
    seg = [Gene("13", 1), Gene("14", -1)]         # arrives from a donor genome
    return BASE[:6] + seg + BASE[6:], [], [6, 7]


EVENTS = [("Duplication", dup, "a segment is copied in tandem"),
          ("Loss", loss, "a segment is deleted"),
          ("Inversion", inversion, "a segment is reversed; its strands flip"),
          ("Transposition", transposition, "a segment moves elsewhere"),
          ("Transfer", transfer, "a segment arrives from a donor genome")]

BW_FILLS_ORDER = ["diag", "#606060", "cross", "#bcbcbc", "dots", "diag2", "#8c8c8c",
                  "horiz", "#e0e0e0", "vert", "white", "#3a3a3a", "black", "check"]


def render(name, after, hlb, hla, desc, mono):
    d = draw.Drawing(CW, CH, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, CW, CH, fill="white"))

    allgenes = BASE + after
    if mono:
        pats = make_patterns(d)
        fills = [pats[x] if x in pats else x for x in BW_FILLS_ORDER]
        order = list(dict.fromkeys(g.family for g in allgenes))
        fam = {f: fills[i % len(fills)] for i, f in enumerate(order)}
        fill_of, hl = (lambda f: fam[f]), "#e6e6e6"
    else:
        fam = color_map(allgenes)
        fill_of, hl = (lambda f: fam[f]), "#f6e7bf"

    draw_ring(d, BASE, fill_of, cx=LC[0], cy=LC[1], r=RR, center_text=None, label_off=30,
              highlight=set(hlb) or None, hl_fill=hl, fs_label=12)
    draw_ring(d, after, fill_of, cx=RC[0], cy=RC[1], r=RR, center_text=None, label_off=30,
              highlight=set(hla) or None, hl_fill=hl, fs_label=12)

    d.append(draw.Text("before", 15, LC[0], LC[1] + RR + 52, font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("after", 15, RC[0], RC[1] + RR + 52, font_family=FONT, text_anchor="middle", fill=INK))

    ax0, ax1, ay = 432, 568, 250                   # central arrow
    d.append(draw.Line(ax0, ay, ax1 - 12, ay, stroke=INK, stroke_width=3))
    d.append(draw.Lines(ax1, ay, ax1 - 15, ay - 8, ax1 - 15, ay + 8, close=True, fill=INK))
    d.append(draw.Text(name, 20, (ax0 + ax1) / 2, ay - 18, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    d.append(draw.Text(desc, 14, CW / 2, CH - 22, font_family=FONT, text_anchor="middle", fill="#555"))

    sub = OUT / ("bw" if mono else "color")
    sub.mkdir(parents=True, exist_ok=True)
    stem = sub / name.lower()
    stem.with_suffix(".svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(stem.with_suffix(".png")), scale=300 / 72.0)


def main():
    for name, fn, desc in EVENTS:
        after, hlb, hla = fn()
        for mono in (False, True):
            render(name, after, hlb, hla, desc, mono)
    print(f"wrote {len(EVENTS)} events × color+bw to {OUT}/")


if __name__ == "__main__":
    main()
