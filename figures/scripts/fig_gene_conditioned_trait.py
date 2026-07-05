"""Figure: genes:traits -- gene-conditioned trait.

A binary *modifier* gene is gained and lost along the tree, and its presence sets a continuous
trait's OU optimum: a lineage that carries the gene is pulled toward theta_present, one that lacks it
toward theta_absent. A discrete genomic event thus reads out as a shift in a continuous phenotype.

  * Panel A (the mechanism): trait value vs time for one lineage -- it sits near theta_absent while
    the gene is absent (light), then the gene is gained (marker) and the trait climbs to
    theta_present (heavy). Two dashed lines mark the two optima.
  * Panel B (a realization): the tree drawn heavy where the modifier is present and light where it
    is absent; to the right, each tip's trait value as a dot on a shared axis with the two optima.
    Carriers sit near theta_present, non-carriers near theta_absent.

House style: B&W, one centered title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_gene_conditioned_trait.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.coevolve import GeneConditionedTrait, simulate_gene_conditioned_trait

from fig_trait_pagel import _layout
from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 700
GREY = "#9a9a9a"
N_TIPS, AGE, TREE_SEED = 12, 1.0, 3
TH_ABS, TH_PRES, ALPHA, SIGMA2 = 0.0, 5.0, 3.0, 0.35
GENE_GAIN, GENE_LOSS = 0.55, 0.55


# --------------------------------------------------------------------------- panel A: the mechanism
def panel_model(d, ox, oy, pw, ph):
    d.append(draw.Text("A   the mechanism", FS_LABEL, ox, oy - 16, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    vlo, vhi = -1.0, 6.0
    x_at = lambda t: ox + t * pw                              # noqa: E731  t in [0,1]
    y_at = lambda v: oy + ph - (v - vlo) / (vhi - vlo) * ph   # noqa: E731

    d.append(draw.Line(ox, oy, ox, oy + ph, stroke=INK, stroke_width=1.4))
    d.append(draw.Line(ox, oy + ph, ox + pw, oy + ph, stroke=INK, stroke_width=1.4))
    for v, lab in ((TH_ABS, "theta_absent"), (TH_PRES, "theta_present")):
        d.append(draw.Line(ox, y_at(v), ox + pw, y_at(v), stroke=MUTED, stroke_width=2.0,
                           stroke_dasharray="6,5"))
        d.append(draw.Text(lab, FS_TICK, ox + pw, y_at(v) - 8, font_family=FONT,
                           text_anchor="end", fill=MUTED))
    d.append(draw.Text("trait", FS_TICK, ox - 12, oy + ph / 2, font_family=FONT,
                       text_anchor="middle", dominant_baseline="central", fill=MUTED,
                       transform=f"rotate(-90 {ox - 12} {oy + ph / 2})"))
    d.append(draw.Text("time", FS_TICK, ox + pw / 2, oy + ph + 26, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    # a trajectory: near theta_absent (light), gene gained at t=0.45 (marker), climbs to present (heavy)
    tg = 0.45
    rng = np.random.default_rng(3)
    absent = [(t, TH_ABS + rng.normal(0, 0.28)) for t in np.linspace(0.02, tg, 22)]
    pres = []
    v = absent[-1][1]
    for t in np.linspace(tg, 0.98, 34):
        v += ALPHA * (TH_PRES - v) * 0.03 + rng.normal(0, 0.22)
        pres.append((t, v))
    pa = draw.Path(fill="none", stroke=GREY, stroke_width=3.0)
    for i, (t, v) in enumerate(absent):
        (pa.M if i == 0 else pa.L)(x_at(t), y_at(v))
    d.append(pa)
    pp = draw.Path(fill="none", stroke=INK, stroke_width=4.4)
    (pp.M)(x_at(absent[-1][0]), y_at(absent[-1][1]))
    for t, v in pres:
        pp.L(x_at(t), y_at(v))
    d.append(pp)
    gx, gy = x_at(tg), y_at(absent[-1][1])
    d.append(draw.Circle(gx, gy, 8.5, fill="white", stroke=INK, stroke_width=1.8))
    d.append(draw.Text("+", FS_TICK, gx, gy + 1, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=INK, font_weight="bold"))
    d.append(draw.Text("gene gained", FS_TICK, gx, gy + 24, font_family=FONT, text_anchor="middle",
                       fill=INK, font_style="italic"))


# --------------------------------------------------------------------------- panel B: a realization
def _pick():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                 direction="backward", seed=TREE_SEED)
    best = None
    for s in range(1, 120):
        res = simulate_gene_conditioned_trait(
            tree, GeneConditionedTrait(gene_gain=GENE_GAIN, gene_loss=GENE_LOSS,
                                       theta_absent=TH_ABS, theta_present=TH_PRES,
                                       alpha=ALPHA, sigma2=SIGMA2), seed=s)
        pres = {n.name: int(v) for n, v in res.gene_presence().items()}
        tv = {n.name: float(v) for n, v in res.trait_values().items()}
        carr = [tv[k] for k in pres if pres[k] == 1]
        noncarr = [tv[k] for k in pres if pres[k] == 0]
        if not carr or not noncarr:
            continue
        sep = np.mean(carr) - np.mean(noncarr)               # clean separation, both groups present
        balance = -abs(len(carr) - len(noncarr))
        score = sep + 0.3 * balance
        if best is None or score > best[0]:
            best = (score, tree, res)
    return best


def panel_realization(d, ox, oy, pw, ph):
    _, tree, res = _pick()
    gmap = {n.name: int(v) for n, v in res.gene.node_values.items()}
    gene_chg = {}
    for node, t, frm, to in res.gene.changes():
        gene_chg.setdefault(node.name, []).append((t, to))
    tv = {n.name: float(v) for n, v in res.trait_values().items()}

    ete = zombi_to_ete3(tree)
    tfo, present, ys, nleaf = _layout(ete)
    tw = 250
    x_at = lambda t: ox + (t / present) * tw                 # noqa: E731
    y_at = lambda k: oy + 34 + (k / max(1, nleaf - 1)) * (ph - 80)   # noqa: E731

    d.append(draw.Text("B   a simulated realization", FS_LABEL, ox, oy - 16, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    def seg(x1, x2, y, on):
        d.append(draw.Line(x1, y, x2, y, stroke=INK if on else GREY,
                           stroke_width=5.2 if on else 2.4, stroke_linecap="butt"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        t0, cur = tfo[n.up.name], gmap[n.up.name]
        for tt, to in sorted(gene_chg.get(n.name, [])):
            seg(x_at(t0), x_at(tt), y, cur == 1)
            t0, cur = tt, to
        seg(x_at(t0), x_at(tfo[n.name]), y, cur == 1)
    for n in ete.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))

    # trait-value strip aligned to the tips
    vlo, vhi = -1.0, 6.5
    ax0 = ox + tw + 60
    axw = pw - tw - 90
    tx = lambda v: ax0 + (v - vlo) / (vhi - vlo) * axw        # noqa: E731
    for v, lab in ((TH_ABS, "theta_absent"), (TH_PRES, "theta_present")):
        d.append(draw.Line(tx(v), y_at(0) - 16, tx(v), y_at(nleaf - 1) + 12, stroke=MUTED,
                           stroke_width=2.0, stroke_dasharray="6,5"))
        d.append(draw.Text(lab, FS_TICK, tx(v), y_at(0) - 24, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
    for n in ete.get_leaves():
        y = y_at(ys[n.name])
        on = gmap[n.name] == 1
        d.append(draw.Circle(tx(tv[n.name]), y, 6.5, fill=INK if on else "white",
                             stroke=INK, stroke_width=1.8))
    d.append(draw.Text("trait value at the tip", FS_TICK, ax0 + axw / 2, y_at(nleaf - 1) + 34,
                       font_family=FONT, text_anchor="middle", fill=MUTED))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Gene-conditioned trait (genes to traits)", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    ly = 82
    d.append(draw.Line(W / 2 - 300, ly, W / 2 - 268, ly, stroke=INK, stroke_width=5.2))
    d.append(draw.Text("modifier present", FS_TICK, W / 2 - 260, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(W / 2 - 110, ly, W / 2 - 78, ly, stroke=GREY, stroke_width=2.4))
    d.append(draw.Text("absent", FS_TICK, W / 2 - 70, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Circle(W / 2 + 40, ly, 6.5, fill=INK, stroke=INK, stroke_width=1.8))
    d.append(draw.Text("tip trait value", FS_TICK, W / 2 + 54, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    panel_model(d, 90, 220, 300, 300)
    panel_realization(d, 500, 200, 660, 400)

    name = "gene_conditioned_trait"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
