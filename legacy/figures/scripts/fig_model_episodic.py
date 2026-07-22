"""Model figure: the episodic (skyline) birth-death process.

Speciation/extinction rates change in time epochs (rate shifts). Here a recent
epoch of fast diversification (high λ, low μ) sits on top of an older, quieter
epoch — visible as a burst of branching once you cross the shift line.

The rate epochs are shown as shaded time bands labelled with their λ (speciation)
and μ (extinction); a few lineages still go extinct (dashed).

Run:  python figures/scripts/fig_model_episodic.py
"""

from __future__ import annotations

from pathlib import Path

import drawsvg as draw

import phylustrator as ph
from zombi2 import EpisodicBirthDeath, simulate_species_tree

from model_common import annotate_depths, draw_skeleton, mark_observed, zombi_to_ete3
from zombi_style import INK, MUTED, species_style, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_STEM = Path(__file__).resolve().parent.parent / "model_episodic" / "model_episodic"

# four rate epochs (present -> past), alternating fast / slow diversification
BIRTH = [1.5, 0.7, 1.5, 0.6]
DEATH = [0.15, 0.35, 0.15, 0.3]
SHIFTS = [1.0, 2.0, 3.0]
AGE, SEED = 4.0, 22            # a densely populated tree (~200 tips)
BAND = "#e9e9e9"          # grey for the alternating shaded epochs


def epoch_bounds(present):
    """Tree-time boundaries [t0, t1] and (λ, μ) for each epoch, present→past."""
    ages = [0.0] + list(SHIFTS) + [present]
    out = []
    for i in range(len(BIRTH)):
        t_lo, t_hi = present - ages[i + 1], present - ages[i]   # older edge, younger edge
        out.append((max(0.0, t_lo), t_hi, BIRTH[i], DEATH[i]))
    return out


def main():
    model = EpisodicBirthDeath(birth=BIRTH, death=DEATH, shifts=SHIFTS)
    ztree = simulate_species_tree(model, age=AGE, direction="forward", seed=SEED)
    tree = zombi_to_ete3(ztree)
    present = annotate_depths(tree)
    mark_observed(tree)

    n_leaves = len(tree.get_leaves())
    # Wide/short landscape canvas; a fixed height keeps the ~200-tip tree from
    # exploding vertically, so the leaves pack tight and the branches thin out.
    style = species_style(width=1320, height=840, margin=118,
                          branch_stroke_width=0.9)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    ys = [l.y_coord for l in tree.get_leaves()]
    y0, y1 = min(ys) - 10, max(ys) + 10

    # alternating white / grey rate-epoch bands + a one-line birth/death label on each
    # ("b" = birth/speciation rate, "d" = death/extinction rate; ASCII only)
    for i, (t_lo, t_hi, lam, mu) in enumerate(epoch_bounds(present)):
        x_lo = d.root_x + t_lo * d.sf
        x_hi = d.root_x + t_hi * d.sf
        if i % 2 == 1:                                      # shade every other epoch
            d.drawing.append(draw.Rectangle(x_lo, y0, x_hi - x_lo, y1 - y0, fill=BAND))
        d.drawing.append(draw.Text(f"b = {lam:g},  d = {mu:g}", FS_ANNOT, (x_lo + x_hi) / 2, y0 - 14,
                                   font_family=style.font_family, text_anchor="middle", fill=INK))
    # shift line(s)
    for s in SHIFTS:
        xs = d.root_x + (present - s) * d.sf
        d.drawing.append(draw.Line(xs, y0, xs, y1, stroke=MUTED, stroke_width=1.2,
                                   stroke_dasharray="3,4"))

    draw_skeleton(d, tree)

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)

    # title: one short bold line, horizontally centered at the top
    d.drawing.append(draw.Text("Episodic (skyline) birth-death", FS_TITLE, 0,
                               -style.height / 2 + 44, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    # legend, single column in the open upper-left corner (the oldest epoch's
    # first band is empty up top before the crown fans out): a gloss for the
    # per-epoch b/d labels + the dashed = extinct convention
    lx, ly = -style.width / 2 + 34, y0 + 30
    d.drawing.append(draw.Text("b = birth rate,  d = death rate", FS_LABEL, lx, ly,
                               font_family=style.font_family, text_anchor="start",
                               dominant_baseline="central", fill=INK))
    d.drawing.append(draw.Line(lx, ly + 32, lx + 34, ly + 32, stroke=INK,
                               stroke_width=2.2,
                               stroke_dasharray="6,5", stroke_linecap="butt"))
    d.drawing.append(draw.Text("extinct lineage", FS_LABEL, lx + 46, ly + 32,
                               font_family=style.font_family, text_anchor="start",
                               dominant_baseline="central", fill=INK))

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    n_ext = sum(1 for l in tree.get_leaves() if l.is_extant)
    print(f"wrote {OUT_STEM}.svg / .png  ({n_ext} extant, {n_leaves} leaves)")


if __name__ == "__main__":
    main()
