"""Model figure: the Fossilized Birth-Death (FBD) process.

Forward simulation with a fossil (serial) sampling rate psi and removal r < 1, so
lineages leave dated fossil samples through time and some continue past sampling
(sampled ancestors). We observe the extant tips (reaching the present) plus the
fossils; the unsampled extinct lineages are the dashed 'dark' part of the tree.

  * extant sampled tip   -> solid lineage reaching the present
  * fossil sample        -> solid black diamond (a dated fossil tip, or a
                            sampled ancestor sitting on a continuing branch)
  * extinct / unsampled  -> dashed lineage (never observed)

Run:  python figures/scripts/fig_model_fbd.py
"""

from __future__ import annotations

from pathlib import Path

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, simulate_species_tree

from model_common import (annotate_depths, draw_fossils, draw_skeleton, fossil_nodes,
                          mark_observed, zombi_to_ete3)
from zombi_style import INK, PANEL, species_style, FS_TITLE, FS_LABEL, FS_TICK

OUT_STEM = Path(__file__).resolve().parent.parent / "model_fbd" / "model_fbd"

# forward FBD: psi = fossil sampling rate, removal < 1 keeps sampled ancestors
MODEL = BirthDeath(0.9, 0.35, fossilization=0.6, removal=0.5)
AGE, SEED = 2.8, 31
PRESENT_LINE = "#c9c9c9"


def add_legend(d, x, y):
    """Single-column key (shared font scale): solid extant / fossil diamond / dashed unsampled."""
    fam = d.style.font_family
    sw = d.style.branch_stroke_width
    L, cy, row = 34, y, FS_LABEL * 1.9
    # extant (solid)
    d.drawing.append(draw.Line(x, cy, x + L, cy, stroke=INK, stroke_width=sw, stroke_linecap="round"))
    d.drawing.append(draw.Text("extant lineage", FS_LABEL, x + L + 14, cy, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))
    cy += row
    # fossil (diamond)
    d._draw_shape_at(x + L / 2, cy, "square", INK, r=7.0, stroke=PANEL, stroke_width=1.1, rotation=45.0)
    d.drawing.append(draw.Text("fossil sample", FS_LABEL, x + L + 14, cy, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))
    cy += row
    # unsampled/extinct (dashed)
    d.drawing.append(draw.Line(x, cy, x + L, cy, stroke=INK, stroke_width=sw,
                               stroke_dasharray="6,5", stroke_linecap="butt"))
    d.drawing.append(draw.Text("extinct / unsampled", FS_LABEL, x + L + 14, cy, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))


def main():
    ztree = simulate_species_tree(MODEL, age=AGE, direction="forward", seed=SEED)
    tree = zombi_to_ete3(ztree)
    present = annotate_depths(tree)
    mark_observed(tree)

    n_leaves = len(tree.get_leaves())
    # Wide/short landscape canvas: fixed width, height sized to the tip count with
    # a cap so the tree stays wider than it is tall (extra top headroom leaves a
    # clean band for the centered title above the legend).
    style = species_style(width=1300, height=min(760, 44 * n_leaves + 300), margin=118)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    ys = [l.y_coord for l in tree.get_leaves()]
    present_x = d.root_x + present * d.sf
    d.drawing.append(draw.Line(present_x, min(ys) - 14, present_x, max(ys) + 14,
                               stroke=PRESENT_LINE, stroke_width=1.0, stroke_dasharray="2,4"))

    draw_skeleton(d, tree)
    draw_fossils(d, fossil_nodes(tree))

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)

    # title: one short bold line, horizontally centered at the top
    d.drawing.append(draw.Text("The fossilized birth-death model", FS_TITLE, 0,
                               -style.height / 2 + 44, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    # legend: single top-left column, below the title and clear of the tree
    add_legend(d, x=-style.width / 2 + 34, y=-style.height / 2 + 96)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    n_foss = len(fossil_nodes(tree))
    n_ext = sum(1 for l in tree.get_leaves() if l.is_extant)
    print(f"wrote {OUT_STEM}.svg / .png  ({n_ext} extant, {n_foss} fossils, {n_leaves} leaves)")


if __name__ == "__main__":
    main()
