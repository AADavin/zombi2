"""Figure: species:genes -- gradual vs punctuational genome evolution.

The genomic twin of the anagenetic/cladogenetic trait figure. A genome (a set of gene families) is
evolved down ONE shared tree two ways:

  * Panel A (gradual / anagenetic): families are lost and gained ALONG the branches; the amount of
    gene-content turnover on a branch scales with its length, so sister tips have similar genomes.
  * Panel B (punctuational / cladogenetic): gene content changes only in a BURST at each speciation
    (a founder-effect upheaval); sister tips can differ sharply in gene content.

Each branch carries a mark whose area is the amount of gene-content change on it -- a circle at the
branch midpoint when the change is gradual, a diamond at the node when it is punctuational. The tip
bars show each extant genome's size.

House style: B&W, one centered title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_punctuational_genome.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.coevolve import CladogeneticGenome, simulate_cladogenetic_genome

from fig_trait_pagel import _layout
from model_common import zombi_to_ete3
from zombi_style import (FONT, INK, MUTED, STATE_ON,
                         FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK)

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1220, 664
GREY = "#9a9a9a"
N_TIPS, AGE, TREE_SEED = 11, 1.0, 3
INIT = 24

# This figure is not branch-state-coded (both panels have black trees; shape distinguishes
# gradual circles from punctuational diamonds). The tasteful colour touch is a pale-teal tint on
# the genome-size bars; the tree and event markers stay INK per house style. Default output is the
# colour version (-> punctuational_genome.svg, embedded); a preserved B&W copy is *_bw.svg.
BAR_COL = STATE_ON


def _by_name(d):
    if callable(d):
        d = d()
    return {(k.name if hasattr(k, "name") else k): v for k, v in d.items()}


def _panel(d, ox, oy, pw, ph, tree, res, mode, header):
    genomes = _by_name(res.node_genomes)
    sizes = _by_name(res.genome_sizes)
    ete = zombi_to_ete3(tree)
    tfo, present, ys, nleaf = _layout(ete)
    x_at = lambda t: ox + 34 + (t / present) * (pw - 190)      # noqa: E731
    y_at = lambda k: oy + 44 + (k / max(1, nleaf - 1)) * (ph - 104)   # noqa: E731

    d.append(draw.Text(header, FS_LABEL, ox, oy - 6, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        d.append(draw.Line(x_at(tfo[n.up.name]), y, x_at(tfo[n.name]), y, stroke=INK,
                           stroke_width=2.6, stroke_linecap="butt"))
        change = len(genomes[n.name] ^ genomes[n.up.name])    # families gained or lost on this branch
        if change:
            r = 3.0 + 1.5 * math.sqrt(change)
            if mode == "grad":                                # gradual: circle at branch midpoint
                mx = (x_at(tfo[n.up.name]) + x_at(tfo[n.name])) / 2
                d.append(draw.Circle(mx, y, r, fill="white", stroke=INK, stroke_width=2.0))
            else:                                             # punctuational: diamond at the node
                nx = x_at(tfo[n.up.name])
                d.append(draw.Lines(nx, y - r, nx + r, y, nx, y + r, nx - r, y,
                                    close=True, fill=INK, stroke=INK, stroke_width=1.5))
    for n in ete.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.6))

    # tip genome-size bars
    smax = max(sizes[n.name] for n in ete.get_leaves())
    bx = ox + pw - 138
    bw = 120
    for n in ete.get_leaves():
        y = y_at(ys[n.name])
        w = bw * sizes[n.name] / smax
        d.append(draw.Rectangle(bx, y - 6, w, 12, fill=BAR_COL, stroke="none"))
        d.append(draw.Text(str(sizes[n.name]), FS_TICK, bx + w + 6, y, font_family=FONT,
                           text_anchor="start", dominant_baseline="central", fill=MUTED))
    d.append(draw.Text("genome size", FS_TICK, bx, oy + 20, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    base = oy + ph - 8
    d.append(draw.Line(x_at(0), base, x_at(present), base, stroke=INK, stroke_width=1.4))
    for k in range(5):
        tv = present * k / 4
        d.append(draw.Line(x_at(tv), base, x_at(tv), base + 6, stroke=INK, stroke_width=1.4))
    d.append(draw.Text("time (root to present)", FS_TICK, (x_at(0) + x_at(present)) / 2, base + 24,
                       font_family=FONT, text_anchor="middle", fill=MUTED))


def render(bw=False):
    global BAR_COL
    BAR_COL = INK if bw else STATE_ON

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Gradual vs punctuational genome (species to genes)",
                       FS_TITLE, W / 2, 46, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))

    ly = 82
    lty = ly + 0.34 * FS_TICK   # L#1: text baseline vertically centred on the legend marker
    d.append(draw.Circle(W / 2 - 250, ly, 6.5, fill="white", stroke=INK, stroke_width=2.0))
    d.append(draw.Text("gene turnover along a branch", FS_TICK, W / 2 - 236, lty, font_family=FONT,
                       text_anchor="start", fill=INK))
    dx = W / 2 + 70
    r = 6.5
    d.append(draw.Lines(dx, ly - r, dx + r, ly, dx, ly + r, dx - r, ly, close=True, fill=INK))
    d.append(draw.Text("turnover at a speciation", FS_TICK, dx + 14, lty, font_family=FONT,
                       text_anchor="start", fill=INK))

    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                 direction="backward", seed=TREE_SEED)
    grad = simulate_cladogenetic_genome(
        tree, CladogeneticGenome(INIT, loss=0.6, origination=0.6,
                                 cladogenetic_loss=0.0, cladogenetic_gain=0.0), seed=4)
    punc = simulate_cladogenetic_genome(
        tree, CladogeneticGenome(INIT, loss=0.0, origination=0.0,
                                 cladogenetic_loss=0.16, cladogenetic_gain=3.0), seed=4)

    _panel(d, 70, 158, 520, 470, tree, grad, "grad",
           "A   gradual: gene content drifts along branches")
    _panel(d, 660, 158, 520, 470, tree, punc, "punc",
           "B   punctuational: gene content bursts at each split")

    name = "punctuational_genome"
    suffix = "_bw" if bw else ""
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}{suffix}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}{suffix}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}{suffix}.svg / .png")


if __name__ == "__main__":
    render(bw=False)   # colour -> punctuational_genome.svg (embedded)
    render(bw=True)    # preserved B&W -> punctuational_genome_bw.svg
