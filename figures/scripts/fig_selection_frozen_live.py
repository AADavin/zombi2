"""Figure (experimental / language-model selection): the two modes, frozen vs live.

The critic (ESM2) can be read in two ways as a gene family evolves down its tree:

  * FROZEN -- read the critic ONCE, on the root protein. Each site's amino-acid preference is baked in
    and the sites then evolve independently (no epistasis). One critic call per family: cheap, exact,
    embarrassingly parallel.
  * LIVE   -- re-read the critic on the CURRENT sequence every `refresh` substitutions/site along each
    lineage, so each site feels the others' current states (within-gene epistasis). Many critic calls;
    refresh -> infinity recovers frozen.

House style: one centered bold title, ASCII text, a single accent for the 'critic read' marker.
Run:  python figures/scripts/fig_selection_frozen_live.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_TICK

OUT = Path(__file__).resolve().parent.parent.parent / "docs" / "img"
W, H = 1120, 500
READ = "#2f7d84"        # a "critic read" marker (muted teal)
BRANCH_W = 3.0


def _tree_segments(ox, oy, w, h):
    """Horizontal 3-tip cladogram; returns (h_segments, v_segments) as ((x0,y0,x1,y1),...)."""
    cy = oy + h / 2
    s1 = ox + 0.34 * w
    s2 = ox + 0.63 * w
    tipx = ox + w - 8
    yA, yB, yC = oy + 0.10 * h, oy + 0.56 * h, oy + 0.92 * h
    ymid = (yB + yC) / 2
    hseg = [(ox, cy, s1, cy), (s1, yA, tipx, yA), (s1, ymid, s2, ymid),
            (s2, yB, tipx, yB), (s2, yC, tipx, yC)]
    vseg = [(s1, yA, s1, yC), (s2, yB, s2, yC)]
    return hseg, vseg, [(tipx, yA), (tipx, yB), (tipx, yC)]


def _draw_tree(d, ox, oy, w, h):
    hseg, vseg, tips = _tree_segments(ox, oy, w, h)
    for x0, y0, x1, y1 in hseg + vseg:
        d.append(draw.Line(x0, y0, x1, y1, stroke=INK, stroke_width=BRANCH_W, stroke_linecap="round"))
    return hseg, tips


def _read_marker(d, x, y):
    d.append(draw.Circle(x, y, 7.5, fill=READ, stroke="white", stroke_width=1.4))


def _panel(d, ox, oy, w, h, title, mode):
    d.append(draw.Text(title, FS_LABEL, ox + w / 2, oy - 18, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    hseg, tips = _draw_tree(d, ox, oy, w, h)
    if mode == "frozen":
        _read_marker(d, hseg[0][0] + 2, hseg[0][1])                 # one read, at the root
    else:
        for (x0, y0, x1, y1) in hseg:                               # reads spaced along every branch
            n = max(1, int((x1 - x0) // 42))
            for k in range(n):
                _read_marker(d, x0 + (k + 0.5) * (x1 - x0) / n, y0)
    for (tx, ty) in tips:
        d.append(draw.Rectangle(tx, ty - 6, 12, 12, fill=INK, stroke="white", stroke_width=1))


def _caption(d, ox, oy, w, lines):
    for k, (txt, bold) in enumerate(lines):
        d.append(draw.Text(txt, FS_TICK, ox, oy + k * 22, font_family=FONT, text_anchor="start",
                           fill=INK if bold else MUTED, font_weight="bold" if bold else "normal"))


def main():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Two ways to read the critic: frozen vs live", FS_TITLE, W / 2, 44,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # legend (top-left, clear of the panels)
    _read_marker(d, 60, 84)
    d.append(draw.Text("= one ESM2 read", FS_TICK, 76, 90, font_family=FONT, text_anchor="start", fill=INK))

    pw, ph = 420, 210
    _panel(d, 70, 150, pw, ph, "FROZEN (default)", "frozen")
    _panel(d, 630, 150, pw, ph, "LIVE (epistatic)", "live")

    _caption(d, 70, 410, pw, [
        ("read the critic ONCE, on the root protein", True),
        ("sites evolve independently (no epistasis)", False),
        ("one call per gene: cheap, exact, parallel", False),
    ])
    _caption(d, 630, 410, pw, [
        ("re-read the critic every 'refresh' subs/site", True),
        ("each site feels the others (within-gene epistasis)", False),
        ("many calls; refresh = infinity recovers frozen", False),
    ])

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "selection_frozen_live.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(),
                     write_to=str(OUT / "selection_frozen_live.png"), scale=200 / 72.0)
    print(f"wrote {OUT}/selection_frozen_live.svg / .png")


if __name__ == "__main__":
    main()
