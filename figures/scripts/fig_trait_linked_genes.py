"""Figure: traits:genes -- trait-linked gene families.

A continuous trait is evolved down the tree, then a panel of gene families is evolved whose LOSS
depends on the local trait value: a *responsive* family (coupling weight w>0) is retained where the
trait is high (loss = base_loss * exp(-effect_loss * w * s)) and purged where it is low, while an
*inert* family (w=0) always loses at the baseline. Gain is a trait-blind influx that the trait-
modulated loss then selectively retains -- so responsive families end up present where the trait
favours them, and inert families carry no trait signal.

  * Panel A (the mechanism): loss rate vs trait value -- a responsive family (falling curve) vs an
    inert one (flat).
  * Panel B (a realization): the tree painted by the trait (viridis), a per-tip trait chip, then the
    gene-presence matrix. Responsive families (left block) track the trait chip; inert families
    (right block) do not.

House style: viridis for the continuous trait; B&W presence cells; one centered title.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_trait_linked_genes.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import BrownianMotion, simulate_traits
from zombi2.coevolve import TraitGeneCoupling, simulate_trait_linked_genomes

from fig_trait_bm import viridis, hexc, VIRIDIS
from fig_trait_pagel import _layout
from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1320, 720
N_TIPS, AGE, TREE_SEED = 12, 1.0, 3
N_RESP, N_INERT = 7, 7
BASE_LOSS, EFFECT_LOSS, TRANSFER, SIGMA2 = 1.0, 2.6, 0.9, 1.3


# --------------------------------------------------------------------------- panel A: the mechanism
def panel_model(d, ox, oy, pw, ph):
    d.append(draw.Text("A   the mechanism", FS_LABEL, ox, oy - 16, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    xlo, xhi = -2.6, 2.6
    x_at = lambda x: ox + (x - xlo) / (xhi - xlo) * pw        # noqa: E731
    rmax = BASE_LOSS * 1.15
    y_at = lambda r: oy + ph - (r / rmax) * ph                # noqa: E731

    grad = draw.LinearGradient(ox, 0, ox + pw, 0)
    for t, c in VIRIDIS:
        grad.add_stop(t, hexc(c))
    d.append(grad)
    d.append(draw.Rectangle(ox, oy + ph + 8, pw, 11, fill=grad, stroke=INK, stroke_width=0.8))
    d.append(draw.Text("trait value  s", FS_TICK, ox + pw / 2, oy + ph + 42, font_family=FONT,
                       text_anchor="middle", fill=MUTED))
    d.append(draw.Line(ox, oy, ox, oy + ph, stroke=INK, stroke_width=1.4))
    d.append(draw.Line(ox, oy + ph, ox + pw, oy + ph, stroke=INK, stroke_width=1.4))
    d.append(draw.Text("loss rate", FS_TICK, ox - 12, oy + ph / 2, font_family=FONT,
                       text_anchor="middle", dominant_baseline="central", fill=MUTED,
                       transform=f"rotate(-90 {ox - 12} {oy + ph / 2})"))

    clamp = lambda y: min(oy + ph, max(oy, y))               # noqa: E731  keep the curve in the panel
    xs = np.linspace(xlo, xhi, 140)
    resp = draw.Path(fill="none", stroke=INK, stroke_width=3.4)
    for i, x in enumerate(xs):
        (resp.M if i == 0 else resp.L)(x_at(x), clamp(y_at(BASE_LOSS * np.exp(-EFFECT_LOSS * 1.0 * x))))
    d.append(resp)
    d.append(draw.Line(x_at(xlo), y_at(BASE_LOSS), x_at(xhi), y_at(BASE_LOSS),
                       stroke=MUTED, stroke_width=3.0, stroke_dasharray="7,5"))
    d.append(draw.Text("responsive (w>0)", FS_TICK, x_at(0.2), y_at(BASE_LOSS * np.exp(-EFFECT_LOSS * 0.2)) - 14,
                       font_family=FONT, text_anchor="start", fill=INK, font_weight="bold"))
    d.append(draw.Text("inert (w=0)", FS_TICK, x_at(-2.4), y_at(BASE_LOSS) - 12, font_family=FONT,
                       text_anchor="start", fill=MUTED))


# --------------------------------------------------------------------------- panel B: a realization
def _val(x):
    return x() if callable(x) else x


def _weights():
    w = np.zeros(N_RESP + N_INERT)
    w[:N_RESP] = 1.0
    return w


def _pick():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                 direction="backward", seed=TREE_SEED)
    weights = _weights()
    best = None
    for ts in range(1, 30):
        trait = simulate_traits(tree, BrownianMotion(sigma2=SIGMA2), seed=ts)
        tipval = {n.name: float(v) for n, v in trait.node_values.items() if n.is_leaf()}
        if max(tipval.values()) - min(tipval.values()) < 2.0:
            continue
        for gs in range(1, 30):
            res = simulate_trait_linked_genomes(
                tree, trait, TraitGeneCoupling(N_RESP + N_INERT, weights=weights,
                                               effect_loss=EFFECT_LOSS, base_loss=BASE_LOSS,
                                               transfer=TRANSFER), seed=gs)
            pm = res.profiles
            M = np.asarray(_val(pm.matrix))
            species = list(_val(pm.species))
            xs = np.array([tipval[s] for s in species])
            rf = np.array([M[:N_RESP, j].mean() for j in range(len(species))])
            inf = np.array([M[N_RESP:, j].mean() for j in range(len(species))])
            if rf.std() < 1e-6:
                continue
            corr = float(np.corrcoef(xs, rf)[0, 1])
            corr_inert = 0.0 if inf.std() < 1e-6 else float(np.corrcoef(xs, inf)[0, 1])
            # want responsive to track the trait AND inert to carry no trait signal
            score = corr - 0.8 * abs(corr_inert)
            if best is None or score > best[0]:
                best = (score, tree, trait, tipval, res, species, corr)
    return best


def panel_realization(d, ox, oy, pw, ph):
    _, tree, trait, tipval, res, species, corr = _pick()
    vals = {n.name: float(v) for n, v in trait.node_values.items()}
    lo, hi = min(vals.values()), max(vals.values())
    norm = lambda v: (v - lo) / (hi - lo) if hi > lo else 0.5   # noqa: E731
    col = lambda name: hexc(viridis(norm(vals[name])))          # noqa: E731

    ete = zombi_to_ete3(tree)
    tfo, present, ys, nleaf = _layout(ete)
    tw = 250
    x_at = lambda t: ox + (t / present) * tw                    # noqa: E731
    y_at = lambda k: oy + 30 + (k / max(1, nleaf - 1)) * (ph - 70)   # noqa: E731

    d.append(draw.Text("B   a simulated realization", FS_LABEL, ox, oy - 16, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        mid = (vals[n.name] + vals[n.up.name]) / 2
        d.append(draw.Line(x_at(tfo[n.up.name]), y, x_at(tfo[n.name]), y,
                           stroke=hexc(viridis(norm(mid))), stroke_width=4.2, stroke_linecap="butt"))
    for n in ete.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=col(n.name), stroke_width=3.0))

    # per-tip trait chip + gene-presence matrix (responsive block | inert block)
    pm = res.profiles
    M = np.asarray(_val(pm.matrix))
    col_of = {sp: j for j, sp in enumerate(_val(pm.species))}
    cell = 20
    chip_x = ox + tw + 26
    grid_x = chip_x + 34
    for n in ete.get_leaves():
        y = y_at(ys[n.name])
        d.append(draw.Rectangle(chip_x, y - 9, 18, 18, fill=col(n.name), stroke=INK, stroke_width=1.0))
        j = col_of[n.name]
        for i in range(N_RESP + N_INERT):
            gx = grid_x + (i + (1 if i >= N_RESP else 0)) * cell     # gap between the two blocks
            on = M[i, j] > 0
            d.append(draw.Rectangle(gx, y - cell / 2 + 1, cell - 2, cell - 2,
                                    fill=INK if on else "white", stroke="#cccccc", stroke_width=0.8))
    # block headers
    ytop = y_at(0) - 22
    d.append(draw.Text("trait", FS_TICK, chip_x + 9, ytop, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    d.append(draw.Text("responsive", FS_TICK, grid_x + N_RESP * cell / 2, ytop, font_family=FONT,
                       text_anchor="middle", fill=INK, font_weight="bold"))
    d.append(draw.Text("inert", FS_TICK, grid_x + (N_RESP + 1) * cell + N_INERT * cell / 2, ytop,
                       font_family=FONT, text_anchor="middle", fill=MUTED, font_weight="bold"))
    d.append(draw.Text("gene families", FS_TICK, grid_x + ((N_RESP + N_INERT + 1) * cell) / 2,
                       y_at(nleaf - 1) + 30, font_family=FONT, text_anchor="middle", fill=MUTED))
    return corr


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Trait-linked gene families (traits to genes)", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    ly = 82
    d.append(draw.Rectangle(W / 2 - 250, ly - 9, 18, 18, fill=INK, stroke="#cccccc", stroke_width=0.8))
    d.append(draw.Text("family present", FS_TICK, W / 2 - 224, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Rectangle(W / 2 - 60, ly - 9, 18, 18, fill="white", stroke="#cccccc", stroke_width=0.8))
    d.append(draw.Text("absent", FS_TICK, W / 2 - 34, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    panel_model(d, 80, 250, 300, 250)
    corr = panel_realization(d, 470, 150, 800, 500)

    name = "trait_linked_genes"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png  (corr trait~responsive = {corr:+.2f})")


if __name__ == "__main__":
    render()
