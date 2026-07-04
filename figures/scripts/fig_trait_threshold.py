"""Trait figure: Felsenstein's threshold model, drawn as a phenogram.

An unobserved continuous *liability* evolves by Brownian motion; the observed
discrete state is simply which side of a threshold the liability is on. A
phenogram (liability on y, time on x) shows the mechanism directly: lineages
diffuse, and the discrete state flips exactly when a lineage crosses the
threshold line.

  * y = liability, x = time; dashed line = threshold
  * each lineage coloured by its current discrete state (below / above)
  * a segment that crosses the threshold changes colour at the crossing

Run:  python figures/scripts/fig_trait_threshold.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/aadria/Desktop/CLAUDE/ZOMBI2-traits")   # traits-enabled zombi2

import cairosvg
import drawsvg as draw

from zombi2 import BirthDeath, ThresholdModel, simulate_species_tree, simulate_traits

from zombi_style import FONT, INK

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_threshold/trait_threshold")

N_TIPS, AGE, TREE_SEED = 12, 1.0, 2
SIGMA2, THRESH, TRAIT_SEED = 1.6, 0.0, 3
STATES = ["below", "above"]
PALETTE = {"below": "#4477AA", "above": "#EE6677"}      # below / above threshold

W, H = 940, 660
ML, MR, MT, MB = 84, 190, 44, 68                        # margins (wide right for legend)


def state_of(v):
    return STATES[1] if v >= THRESH else STATES[0]


def nice_ticks(lo, hi, n=5):
    return [round(lo + (hi - lo) * i / n, 2) for i in range(n + 1)]


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    res = simulate_traits(ztree, ThresholdModel(thresholds=[THRESH], sigma2=SIGMA2, x0=0.0,
                                                states=STATES), seed=TRAIT_SEED)
    liab = {n: float(v) for n, v in res.node_values.items()}
    present = ztree.total_age
    vlo, vhi = min(liab.values()), max(liab.values())
    padv = (vhi - vlo) * 0.10 or 1.0
    vlo, vhi = vlo - padv, vhi + padv

    def X(t):
        return ML + (t / present) * (W - ML - MR)

    def Y(v):
        return H - MB - (v - vlo) / (vhi - vlo) * (H - MT - MB)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    fam = FONT

    # threshold line (labelled in the legend, to avoid overlapping the lineages)
    yt = Y(THRESH)
    d.append(draw.Line(ML, yt, W - MR, yt, stroke="#999999", stroke_width=1.4, stroke_dasharray="6,5"))

    # axes
    d.append(draw.Line(ML, H - MB, W - MR, H - MB, stroke=INK, stroke_width=1.6))   # x
    d.append(draw.Line(ML, MT, ML, H - MB, stroke=INK, stroke_width=1.6))           # y
    for t in [round(AGE * i / 4, 6) for i in range(5)]:
        x = X(t)
        d.append(draw.Line(x, H - MB, x, H - MB + 6, stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"{t:.2f}", 13, x, H - MB + 22, font_family=fam, text_anchor="middle", fill=INK))
    for v in nice_ticks(vlo, vhi, 5):
        y = Y(v)
        d.append(draw.Line(ML - 6, y, ML, y, stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"{v:g}", 13, ML - 10, y, font_family=fam, text_anchor="end",
                           dominant_baseline="middle", fill=INK))
    d.append(draw.Text("Time (root to present)", 15, (ML + W - MR) / 2, H - 16, font_family=fam,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("Liability", 15, 22, (MT + H - MB) / 2, font_family=fam, text_anchor="middle",
                       fill=INK, transform=f"rotate(-90, 22, {(MT + H - MB) / 2})"))

    # lineages (colour by discrete state; split a segment at a threshold crossing)
    for n in ztree.nodes_preorder():
        if n.parent is None:
            continue
        t0, v0, t1, v1 = n.parent.time, liab[n.parent], n.time, liab[n]
        s0, s1 = state_of(v0), state_of(v1)
        if s0 == s1:
            d.append(draw.Line(X(t0), Y(v0), X(t1), Y(v1), stroke=PALETTE[s0], stroke_width=3.0,
                               stroke_linecap="round"))
        else:
            f = (THRESH - v0) / (v1 - v0)
            tc = t0 + (t1 - t0) * f
            d.append(draw.Line(X(t0), Y(v0), X(tc), yt, stroke=PALETTE[s0], stroke_width=3.0, stroke_linecap="round"))
            d.append(draw.Line(X(tc), yt, X(t1), Y(v1), stroke=PALETTE[s1], stroke_width=3.0, stroke_linecap="round"))

    # speciation dots + tip markers
    for n in ztree.nodes_preorder():
        if not n.is_leaf():
            d.append(draw.Circle(X(n.time), Y(liab[n]), 2.6, fill=INK))
        else:
            d.append(draw.Circle(X(n.time), Y(liab[n]), 5.0, fill=PALETTE[state_of(liab[n])],
                                 stroke="white", stroke_width=1.2))

    # legend (right margin)
    lx, ly = W - MR + 20, MT + 20
    d.append(draw.Text("Threshold model", 16, lx, ly, font_weight="bold", font_family=fam, text_anchor="start"))
    for i, s in enumerate(STATES):
        cy = ly + 26 + i * 24
        d.append(draw.Line(lx, cy, lx + 26, cy, stroke=PALETTE[s], stroke_width=3.0, stroke_linecap="round"))
        d.append(draw.Text(f"state: {s}", 14, lx + 34, cy, font_family=fam, text_anchor="start",
                           dominant_baseline="middle", fill=INK))
    cy = ly + 26 + 2 * 24
    d.append(draw.Line(lx, cy, lx + 26, cy, stroke="#999999", stroke_width=1.4, stroke_dasharray="6,5"))
    d.append(draw.Text("threshold", 14, lx + 34, cy, font_family=fam, text_anchor="start",
                       dominant_baseline="middle", fill="#666"))

    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    OUT_STEM.with_suffix(".svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT_STEM.with_suffix(".png")),
                     scale=300 / 72.0)
    n_cross = sum(1 for n in ztree.nodes_preorder() if n.parent is not None
                  and state_of(liab[n]) != state_of(liab[n.parent]))
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, {n_cross} threshold crossings)")


if __name__ == "__main__":
    main()
