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
from zombi_style import INK, MUTED, species_style

OUT_STEM = Path(__file__).resolve().parent.parent / "model_episodic" / "model_episodic"

# four rate epochs (present -> past), alternating fast / slow diversification
BIRTH = [1.5, 0.7, 1.5, 0.6]
DEATH = [0.15, 0.35, 0.15, 0.3]
SHIFTS = [0.75, 1.5, 2.25]
AGE, SEED = 3.0, 16
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
    style = species_style(width=920, height=max(680, 42 * n_leaves + 200))
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    ys = [l.y_coord for l in tree.get_leaves()]
    y0, y1 = min(ys) - 22, max(ys) + 22

    # alternating white / grey rate-epoch bands + a one-line λ/μ label on each
    for i, (t_lo, t_hi, lam, mu) in enumerate(epoch_bounds(present)):
        x_lo = d.root_x + t_lo * d.sf
        x_hi = d.root_x + t_hi * d.sf
        if i % 2 == 1:                                      # shade every other epoch
            d.drawing.append(draw.Rectangle(x_lo, y0, x_hi - x_lo, y1 - y0, fill=BAND))
        d.drawing.append(draw.Text(f"λ = {lam:g},  μ = {mu:g}", 14, (x_lo + x_hi) / 2, y0 - 12,
                                   font_family=style.font_family, text_anchor="middle", fill=INK))
    # shift line(s)
    for s in SHIFTS:
        xs = d.root_x + (present - s) * d.sf
        d.drawing.append(draw.Line(xs, y0, xs, y1, stroke=MUTED, stroke_width=1.2,
                                   stroke_dasharray="3,4"))

    d.add_text("Episodic (skyline) birth–death", x=-style.width / 2 + 34, y=y0 - 36,
               font_size=16, color=INK, weight="bold")

    draw_skeleton(d, tree)

    ticks = [round(present * i / 4, 6) for i in range(5)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.2f}" for t in ticks],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    n_ext = sum(1 for l in tree.get_leaves() if l.is_extant)
    print(f"wrote {OUT_STEM}.svg / .png  ({n_ext} extant, {n_leaves} leaves)")


if __name__ == "__main__":
    main()
