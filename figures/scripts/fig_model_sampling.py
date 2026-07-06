"""Model figure: incomplete (extant) sampling, rho < 1.

Forward birth-death where only a fraction rho of the species alive at the present
are sampled. Unlike extinction, these missing lineages are *alive today* — they
just aren't in your data. The sample (solid) is therefore missing two different
things: lineages that went extinct (dashed, stop short) AND extant lineages that
were not sampled (dashed, reach the present, marked with an open circle).

  * sampled species       -> solid, reaches the present
  * unsampled extant      -> dashed, reaches the present, open circle at the tip
  * extinct               -> dashed, stops short of the present

Run:  python figures/scripts/fig_model_sampling.py
"""

from __future__ import annotations

from pathlib import Path

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, simulate_species_tree

from model_common import annotate_depths, draw_skeleton, mark_observed, zombi_to_ete3
from zombi_style import INK, PANEL, species_style, FS_TITLE, FS_LABEL, FS_TICK

OUT_STEM = Path(__file__).resolve().parent.parent / "model_sampling" / "model_sampling"

MODEL = BirthDeath(1.0, 0.35, sampling_fraction=0.6)     # rho = 0.6
AGE, SEED = 2.8, 5
PRESENT_LINE = "#c9c9c9"
CIRCLE_R = 6.0


def unsampled_extant(tree, present, tol):
    """Leaves alive at the present but not sampled (is_extant False, sampled False, at present)."""
    return [l for l in tree.get_leaves()
            if not l.is_extant and not l.sampled and abs(l.depth - present) <= tol]


def add_legend(d, x, y):
    """Single-column key (shared font scale): sampled / unsampled-extant (open circle) / extinct."""
    fam, sw = d.style.font_family, d.style.branch_stroke_width
    L, row, cy = 34, FS_LABEL * 1.9, y

    def item(label, dashed, circle=False):
        nonlocal cy
        kw = dict(stroke=INK, stroke_width=sw, stroke_linecap="round")
        if dashed:
            kw = dict(stroke=INK, stroke_width=sw, stroke_dasharray="6,5", stroke_linecap="butt")
        d.drawing.append(draw.Line(x, cy, x + L, cy, **kw))
        tx = x + L + 14
        if circle:
            d.drawing.append(draw.Circle(x + L + 8, cy, 6.0, fill=PANEL, stroke=INK, stroke_width=2))
            tx = x + L + 22
        d.drawing.append(draw.Text(label, FS_LABEL, tx, cy, font_family=fam,
                                   text_anchor="start", dominant_baseline="central", fill=INK))
        cy += row

    item("sampled", False)
    item("unsampled extant", True, circle=True)
    item("extinct", True)


def main():
    ztree = simulate_species_tree(MODEL, age=AGE, direction="forward", seed=SEED)
    tree = zombi_to_ete3(ztree)
    present = annotate_depths(tree)
    mark_observed(tree)

    n_leaves = len(tree.get_leaves())
    # extra-wide landscape canvas so the tree fills the PDF page width (height capped
    # so it stays wider than it is tall); extra top headroom for the title + legend
    style = species_style(width=1440, height=min(860, 30 * n_leaves + 260), margin=118)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    ys = [l.y_coord for l in tree.get_leaves()]
    present_x = d.root_x + present * d.sf
    d.drawing.append(draw.Line(present_x, min(ys) - 14, present_x, max(ys) + 14,
                               stroke=PRESENT_LINE, stroke_width=1.0, stroke_dasharray="2,4"))

    draw_skeleton(d, tree)

    # open circle at each unsampled-but-extant tip (alive today, not sampled)
    for lf in unsampled_extant(tree, present, tol=1e-6 * present):
        x, y = lf.coordinates
        d.drawing.append(draw.Circle(x, y, CIRCLE_R, fill=PANEL, stroke=INK, stroke_width=2))

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)

    # title: one short bold line, horizontally centered at the top (ASCII "rho")
    d.drawing.append(draw.Text("Incomplete sampling (rho = 0.6)", FS_TITLE, 0,
                               -style.height / 2 + 44, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    # legend: single column in the top-left, just under the title and clear of the
    # crown (which radiates to the right)
    add_legend(d, x=-style.width / 2 + 34, y=-style.height / 2 + 120)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    n_s = sum(1 for l in tree.get_leaves() if l.is_extant)
    n_u = len(unsampled_extant(tree, present, tol=1e-6 * present))
    print(f"wrote {OUT_STEM}.svg / .png  ({n_s} sampled, {n_u} unsampled-extant, {n_leaves} leaves)")


if __name__ == "__main__":
    main()
