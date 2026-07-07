"""Figure (Ch15): the nucleotide substitution models, and an alignment they produce.

Once a gene tree is a phylogram, sequence evolution draws an actual alignment down it under
a substitution model -- a rate matrix Q over the four bases. The figure has two parts.

Top -- the model progression as exchange graphs on the four bases (purines A, G on top;
pyrimidines C, T below). Edge width = exchange rate; node area = stationary base frequency.
  * JC69   : all exchanges equal, all bases equal -- one parameter-free model.
  * K80    : transitions (A<->G, C<->T) faster than transversions by kappa; bases still equal.
  * HKY85  : K80's transition bias PLUS unequal base frequencies (bigger nodes = commoner base).
  * GTR    : all six exchange rates free, plus unequal frequencies -- the general model.

Bottom -- one HKY85 alignment simulated down a small phylogram: a base per site at each tip,
substitutions accumulating along the branches. Longer branches carry more changes.

House style: B&W graphs; the four bases are a categorical set, so they (and only they) get
four colour-blind-safe colours, as the DEC areas do.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_seq_subst_models.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

import clock_common as C
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1340, 900

# four bases -- categorical colours (Paul Tol 'bright'), the sanctioned exception
BASES = "ACGT"
BASE_COL = {"A": "#4477AA", "C": "#EE6677", "G": "#228833", "T": "#CCBB44"}
PURINES, PYRIM = "AG", "CT"


def graph_panel(d, ox, oy, s, title, subtitle, rates, freqs):
    """One exchange graph. rates: dict of frozenset({x,y}) -> rate; freqs: dict base->freq."""
    pos = {"A": (ox, oy), "G": (ox + s, oy), "C": (ox, oy + s), "T": (ox + s, oy + s)}
    rmax = max(rates.values())
    ew = lambda r: 1.2 + 8.8 * (r / rmax)                 # noqa: E731  edge width by rate
    nr = lambda b: 11 + 41 * math.sqrt(freqs[b])          # noqa: E731  node radius by freq

    d.append(draw.Text(title, FS_LABEL, ox + s / 2, oy - 66, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text(subtitle, FS_TICK, ox + s / 2, oy - 44, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))

    # edges first (behind nodes); transitions solid dark, transversions lighter
    order = [("A", "G"), ("C", "T"), ("A", "C"), ("G", "T"), ("A", "T"), ("G", "C")]
    for x, y in order:
        r = rates[frozenset((x, y))]
        is_ts = (x in PURINES and y in PURINES) or (x in PYRIM and y in PYRIM)
        col = INK if is_ts else "#9a9a9a"
        (x1, y1), (x2, y2) = pos[x], pos[y]
        d.append(draw.Line(x1, y1, x2, y2, stroke=col, stroke_width=ew(r), stroke_linecap="round"))
    # nodes
    for b, (x, y) in pos.items():
        d.append(draw.Circle(x, y, nr(b), fill=BASE_COL[b], stroke=INK, stroke_width=1.4))
        d.append(draw.Text(b, FS_LABEL, x, y + 0.34 * FS_LABEL, font_family=FONT,
                           text_anchor="middle", fill="white", font_weight="bold"))


def expm_reversible(Q):
    """exp(Q) via eigendecomposition (Q is a valid 4x4 generator over one unit of time)."""
    w, V = np.linalg.eig(Q)
    return (V @ np.diag(np.exp(w)) @ np.linalg.inv(V)).real


def hky_Q(kappa, pi):
    pi = np.array(pi, float)
    idx = {b: i for i, b in enumerate(BASES)}
    Q = np.zeros((4, 4))
    for i, bi in enumerate(BASES):
        for j, bj in enumerate(BASES):
            if i == j:
                continue
            ts = (bi in PURINES and bj in PURINES) or (bi in PYRIM and bj in PYRIM)
            Q[i, j] = (kappa if ts else 1.0) * pi[j]
        Q[i, i] = -Q[i].sum()
    # normalise to one expected substitution per unit time
    mu = -sum(pi[i] * Q[i, i] for i in range(4))
    return Q / mu, pi, idx


def simulate_alignment(tree, dist, ncol, kappa, pi, seed):
    """Evolve an ncol-site sequence down the phylogram under HKY85; return {node: seq}."""
    rng = np.random.default_rng(seed)
    Q, pi, idx = hky_Q(kappa, pi)
    inv = {i: b for b, i in idx.items()}
    seqs = {}
    root = tree.get_tree_root()
    seqs[root.name] = rng.choice(4, size=ncol, p=pi)
    for n in tree.traverse("preorder"):
        if n.is_root():
            continue
        t = dist[n.name] - dist[n.up.name]
        P = expm_reversible(Q * max(t, 0.0))
        P = np.clip(P, 0, None)
        P /= P.sum(axis=1, keepdims=True)
        parent = seqs[n.up.name]
        child = np.empty(ncol, dtype=int)
        for c in range(ncol):
            child[c] = rng.choice(4, p=P[parent[c]])
        seqs[n.name] = child
    return {k: "".join(inv[i] for i in v) for k, v in seqs.items()}


def alignment_panel(d, ox, oy, tree, seqs, dist):
    ys, nleaf = C.leaf_ys(tree)
    leaves = tree.get_leaves()
    maxd = max(dist[l.name] for l in leaves)
    tw = 300                                              # tree width
    x_at = lambda v: ox + (v / maxd) * tw                 # noqa: E731
    rh = 30
    y_at = lambda k: oy + 30 + k * rh                     # noqa: E731 (one row per leaf, by draw order)
    # map each leaf to a compact row index
    leaf_rows = {lf.name: i for i, lf in enumerate(leaves)}
    ynode = {}
    for n in tree.traverse("postorder"):
        ynode[n.name] = (y_at(leaf_rows[n.name]) if n.is_leaf()
                         else sum(ynode[c.name] for c in n.children) / len(n.children))

    d.append(draw.Text("An HKY85 alignment simulated down the phylogram", FS_LABEL, ox, oy,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    # tree
    for n in tree.traverse():
        if n.is_root():
            continue
        d.append(draw.Line(x_at(dist[n.up.name]), ynode[n.name], x_at(dist[n.name]), ynode[n.name],
                           stroke=INK, stroke_width=2.4, stroke_linecap="round"))
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(dist[n.name])
            yy = [ynode[c.name] for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))
    # aligned sequences at the tips
    ncol = len(next(iter(seqs.values())))
    cw = 15
    ax = ox + tw + 40
    # column shading guides
    for lf in leaves:
        r = leaf_rows[lf.name]
        y = y_at(r)
        seq = seqs[lf.name]
        for c, base in enumerate(seq):
            d.append(draw.Rectangle(ax + c * cw, y - rh / 2 + 3, cw - 1.4, rh - 6,
                                    fill=BASE_COL[base], stroke="none"))
            d.append(draw.Text(base, 13, ax + c * cw + cw / 2, y + 1, font_family="Courier",
                               text_anchor="middle", dominant_baseline="central", fill="white",
                               font_weight="bold"))
    # (base-colour legend removed -- each alignment cell already carries its letter)
    return ax + ncol * cw


def main():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Substitution models over the four bases", FS_TITLE, W / 2, 44,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # --- top: four model graphs -------------------------------------------
    def rates_dict(ag, ct, ac, gt, at, gc):
        return {frozenset(("A", "G")): ag, frozenset(("C", "T")): ct, frozenset(("A", "C")): ac,
                frozenset(("G", "T")): gt, frozenset(("A", "T")): at, frozenset(("G", "C")): gc}
    eqf = {b: 0.25 for b in BASES}
    unf = {"A": 0.15, "C": 0.34, "G": 0.36, "T": 0.15}
    s = 150
    gy = 205
    xs = [90, 90 + 300, 90 + 600, 90 + 900]
    graph_panel(d, xs[0], gy, s, "JC69", "equal rates, equal freqs",
                rates_dict(1, 1, 1, 1, 1, 1), eqf)
    graph_panel(d, xs[1], gy, s, "K80", "transition bias (kappa)",
                rates_dict(3, 3, 1, 1, 1, 1), eqf)
    graph_panel(d, xs[2], gy, s, "HKY85", "kappa + unequal freqs",
                rates_dict(3, 3, 1, 1, 1, 1), unf)
    graph_panel(d, xs[3], gy, s, "GTR", "six free rates + freqs",
                rates_dict(3.1, 2.4, 0.7, 1.3, 0.5, 1.8), unf)

    # shared legend under the graphs -- laid out left-to-right, then centred on the page so the
    # three blocks sit as one balanced row (previously it ran off to the right and the transition
    # label collided with the transversion swatch).
    lgy = gy + s + 78
    ty = lgy + 0.34 * FS_TICK
    cwf = 0.55 * FS_TICK                                   # approx character width at FS_TICK
    seg, g1, gap = 30, 12, 50                              # swatch line, line->text gap, block gap
    t1 = "transition (A<->G, C<->T)"
    t2 = "transversion"
    t3 = "edge width = rate;   node area = base frequency"
    wA = seg + g1 + cwf * len(t1)
    wB = seg + g1 + cwf * len(t2)
    wC = cwf * len(t3)
    x = W / 2 - (wA + gap + wB + gap + wC) / 2
    d.append(draw.Line(x, lgy, x + seg, lgy, stroke=INK, stroke_width=6))
    d.append(draw.Text(t1, FS_TICK, x + seg + g1, ty, font_family=FONT, text_anchor="start", fill=INK))
    x += wA + gap
    d.append(draw.Line(x, lgy, x + seg, lgy, stroke="#9a9a9a", stroke_width=3))
    d.append(draw.Text(t2, FS_TICK, x + seg + g1, ty, font_family=FONT, text_anchor="start", fill=INK))
    x += wB + gap
    d.append(draw.Text(t3, FS_TICK, x, ty, font_family=FONT, text_anchor="start", fill=MUTED))

    # --- bottom: a simulated alignment ------------------------------------
    tree = C.build_tree(n_tips=8, seed=11)
    _, present = C.node_times(tree)
    rate, sub = C.autocorrelated_lognormal(tree, 0.4, seed=2)
    dist = C.subst_dist_to_root(tree, sub)
    seqs = simulate_alignment(tree, dist, ncol=44, kappa=4.0,
                              pi=[unf["A"], unf["C"], unf["G"], unf["T"]], seed=7)
    alignment_panel(d, 90, gy + s + 170, tree, seqs, dist)

    name = "seq_subst_models"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    main()
