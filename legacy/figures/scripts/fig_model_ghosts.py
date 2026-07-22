"""Model figure: backward simulation + grafted ghost lineages, as three panels.

Backward simulation yields only the reconstructed tree of survivors (panel A).
``add_ghost_lineages`` then grafts back the extinct 'ghost' lineages the same
birth-death process would have produced — a *random* overlay, so different seeds
give different hidden histories on the SAME observed tree (panels B and C).

  * reconstructed (observed) lineage -> solid, reaches the present
  * ghost lineage (grafted, extinct) -> dashed dead-end

Run:  python figures/scripts/fig_model_ghosts.py
"""

from __future__ import annotations

import re
from pathlib import Path

import cairosvg
import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, add_ghost_lineages, simulate_species_tree

from model_common import annotate_depths, draw_skeleton, mark_observed, zombi_to_ete3
from zombi_style import INK, species_style

OUT_STEM = Path(__file__).resolve().parent.parent / "model_ghosts" / "model_ghosts"

MODEL = BirthDeath(1.0, 0.7)
N_TIPS, AGE, RECON_SEED = 12, 2.0, 11
GHOST_SEEDS = (11, 13)          # two realizations for panels B and C
PW, PH = 820, 680               # per-panel size (landscape: wide, height sized to content)
# The three panels compose into one very wide (~4:1) figure, so on the page every
# label is scaled down hard. Fonts are therefore set MUCH larger than a single-panel
# figure would use, so the axis, panel headers, and legend stay readable in the manual.
PRESENT_LINE = "#c9c9c9"


def reconstructed():
    return simulate_species_tree(MODEL, n_tips=N_TIPS, age=AGE, direction="backward", seed=RECON_SEED)


def panel(label, subtitle, ghost_seed=None):
    ztree = reconstructed()
    if ghost_seed is not None:
        ztree = add_ghost_lineages(ztree, MODEL, seed=ghost_seed)
    tree = zombi_to_ete3(ztree)
    present = annotate_depths(tree)
    mark_observed(tree)

    style = species_style(width=PW, height=PH, margin=118, font_size=34)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    ys = [l.y_coord for l in tree.get_leaves()]
    present_x = d.root_x + present * d.sf
    d.drawing.append(draw.Line(present_x, min(ys) - 12, present_x, max(ys) + 12,
                               stroke=PRESENT_LINE, stroke_width=1.0, stroke_dasharray="2,4"))
    draw_skeleton(d, tree)

    ticks = [round(present * i / 2, 6) for i in range(3)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.1f}" for t in ticks],
                    label="Time", tick_size=8.0, padding=16.0, stroke_width=1.6)

    d.add_text(label, x=-PW / 2 + 24, y=-PH / 2 + 60, font_size=54, color=INK, weight="bold")
    d.add_text(subtitle, x=-PW / 2 + 82, y=-PH / 2 + 60, font_size=38, color=INK)
    return d


FOOTER = 110                    # strip below the panels for a shared legend


def _footer_legend(cx, y):
    """Raw-SVG one-line solid/dashed key centred at x = cx."""
    items = [("reconstructed lineage", False), ("ghost lineage", True)]
    out = []
    x = cx - 470
    for lab, dash in items:
        dash_attr = ' stroke-dasharray="8,6" stroke-linecap="butt"' if dash else ' stroke-linecap="round"'
        out.append(f'<line x1="{x}" y1="{y}" x2="{x + 54}" y2="{y}" stroke="{INK}" '
                   f'stroke-width="4.2"{dash_attr} />')
        out.append(f'<text x="{x + 66}" y="{y}" font-size="38" font-family="Helvetica" '
                   f'text-anchor="start" dominant-baseline="middle" fill="{INK}">{lab}</text>')
        x += 560
    return "\n".join(out)


def compose(drawers, out_stem):
    """Lay the panel drawings out left-to-right into one SVG (+ PNG), shared footer legend."""
    n = len(drawers)
    W, H = PW * n, PH + FOOTER
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
             f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
             f'<rect x="0" y="0" width="{W}" height="{H}" fill="white" />']
    for i, d in enumerate(drawers):
        svg = re.sub(r"<\?xml[^>]*\?>\s*", "", d.drawing.as_svg())
        svg = svg.replace("<svg ", f'<svg x="{PW * i}" y="0" ', 1)
        parts.append(svg)
    parts.append(_footer_legend(W / 2, PH + FOOTER / 2))
    parts.append("</svg>")
    parent = "\n".join(parts)
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    (out_stem.with_suffix(".svg")).write_text(parent, encoding="utf-8")
    cairosvg.svg2png(bytestring=parent.encode("utf-8"), write_to=str(out_stem.with_suffix(".png")),
                     scale=300 / 72.0)


def main():
    a = panel("A", "Reconstructed tree")
    b = panel("B", "Ghost lineages (seed 1)", ghost_seed=GHOST_SEEDS[0])
    c = panel("C", "Ghost lineages (seed 2)", ghost_seed=GHOST_SEEDS[1])
    compose([a, b, c], OUT_STEM)
    ng = [sum(1 for l in zombi_to_ete3(add_ghost_lineages(reconstructed(), MODEL, seed=s)).get_leaves()
              if not l.is_extant) for s in GHOST_SEEDS]
    print(f"wrote {OUT_STEM}.svg / .png  (12 reconstructed; ghosts B={ng[0]}, C={ng[1]})")


if __name__ == "__main__":
    main()
