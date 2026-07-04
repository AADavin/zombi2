"""Monochrome variant of the circular-genome figure: family fills are hatch patterns
and grey shades instead of colours (no colour). 12 genes, each its own family.

Run:  python figures/scripts/fig_genome_circular_bw.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw

from fig_genome_circular import H, W, draw_ring, initial_genes

OUT_STEM = Path(__file__).resolve().parent.parent / "genome_circular_bw" / "genome_circular_bw"


def make_patterns(d, s=6, sw=1.3):
    def hatch(pid, lines):
        p = draw.Pattern(s, s, id=pid, patternUnits="userSpaceOnUse")
        p.append(draw.Rectangle(0, 0, s, s, fill="white"))
        for x1, y1, x2, y2 in lines:
            p.append(draw.Line(x1, y1, x2, y2, stroke="black", stroke_width=sw))
        d.append(p)
        return p

    pats = {
        "diag":  hatch("g_diag",  [(o, s, o + s, 0) for o in (-s, 0, s)]),      # /
        "diag2": hatch("g_diag2", [(o, 0, o + s, s) for o in (-s, 0, s)]),      # \
        "horiz": hatch("g_horiz", [(0, s * 0.3, s, s * 0.3), (0, s * 0.7, s, s * 0.7)]),
        "vert":  hatch("g_vert",  [(s * 0.3, 0, s * 0.3, s), (s * 0.7, 0, s * 0.7, s)]),
        "cross": hatch("g_cross", [(0, s * 0.5, s, s * 0.5), (s * 0.5, 0, s * 0.5, s)]),
    }
    pd = draw.Pattern(s, s, id="g_dots", patternUnits="userSpaceOnUse")
    pd.append(draw.Rectangle(0, 0, s, s, fill="white"))
    pd.append(draw.Circle(s / 2, s / 2, 1.35, fill="black"))
    d.append(pd)
    pats["dots"] = pd

    pc = draw.Pattern(s, s, id="g_check", patternUnits="userSpaceOnUse")
    pc.append(draw.Rectangle(0, 0, s, s, fill="white"))
    pc.append(draw.Rectangle(0, 0, s / 2, s / 2, fill="black"))
    pc.append(draw.Rectangle(s / 2, s / 2, s / 2, s / 2, fill="black"))
    d.append(pc)
    pats["check"] = pc
    return pats


def main():
    genes = initial_genes(12)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    pats = make_patterns(d)
    # 12 distinct black-and-white fills (patterns + grey shades)
    fills = [pats["diag"], "#555555", pats["cross"], "#bfbfbf", pats["dots"], pats["diag2"],
             "#8f8f8f", pats["horiz"], "#e2e2e2", pats["vert"], "white", "#333333"]
    fam_fill = {g.family: fills[i % len(fills)] for i, g in enumerate(genes)}
    draw_ring(d, genes, lambda f: fam_fill[f])

    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    OUT_STEM.with_suffix(".svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT_STEM.with_suffix(".png")),
                     scale=300 / 72.0)
    print(f"wrote {OUT_STEM}.svg / .png  ({len(genes)} genes, each its own family)")


if __name__ == "__main__":
    main()
