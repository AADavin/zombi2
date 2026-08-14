"""Shared rendering utilities for the ZOMBI2 examples-gallery preview.

Each example is an :class:`Example` — id, title, caption, tag, and a ``render(out_png)`` function that
simulates with ZOMBI2 and plots with Phylustrator (optionally compositing a matplotlib companion
panel). Level modules (``traits.py``, ``species.py``, …) each expose an ``EXAMPLES`` list; ``build.py``
renders them all and regenerates ``gallery.html``.
"""

from __future__ import annotations

import hashlib
import inspect
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Callable

# Run against the zombi2 in this repo (gallery/ sits at the repo root), installed or not.
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import cm, colors

import phylustrator as ph

# Phylustrator's line-safe ramp, mirrored into matplotlib. The gallery draws trees with Phylustrator
# and the diagram beside them with matplotlib, so a colormap named in one has to exist in the other or
# the two halves of a figure disagree. `viridis_dark` is viridis stopped before it goes pale, because
# a branch a few pixels wide loses the light end against white (see the Phylustrator colour module).
_LINE_SAFE = {
    "viridis_dark": [(68, 1, 84), (72, 26, 108), (71, 47, 125), (65, 68, 135),
                     (57, 86, 140), (49, 104, 142), (42, 120, 142), (35, 136, 142),
                     (31, 152, 139), (34, 168, 132), (53, 183, 121), (84, 197, 104)],
    "magma_dark": [(0, 0, 4), (11, 9, 36), (32, 17, 75), (59, 15, 112),
                   (87, 21, 126), (114, 31, 129), (140, 41, 129), (168, 50, 125),
                   (196, 60, 117), (222, 73, 104), (241, 96, 93), (250, 127, 94)],
}
for _name, _anchors in _LINE_SAFE.items():
    if _name not in matplotlib.colormaps:
        matplotlib.colormaps.register(colors.LinearSegmentedColormap.from_list(
            _name, [tuple(v / 255 for v in rgb) for rgb in _anchors]))

# --- one dial for the whole gallery ---------------------------------------
BW = 4.0                # branch width
MARGIN = 45
TREE_W, TREE_H = 1000, 760      # tree-only figures
COMP_W, COMP_H = 1000, 560      # a tree destined for a panel below it


@dataclass
class Example:
    id: str
    title: str
    caption: str          # HTML — may contain &nbsp;, entities
    tag: str
    render: Callable[[str], None]
    code: str | None = None   # clean snippet shown on the detail view (falls back to render source)


def style(w: int = TREE_W, h: int = TREE_H) -> ph.Style:
    return ph.Style(width=w, height=h, margin=MARGIN, branch_width=BW)


def node_values(result) -> dict:
    """{'n<id>': value} keyed to match Phylustrator's node names."""
    return {f"n{i}": v for i, v in result.node_values.items()}


def dashed_extinct(tree, ct) -> set:
    """The set of node names whose whole subtree is extinct (their branches draw dashed).

    ``ct`` is the ZOMBI2 complete tree, and the names come from its own `labels()` — a lineage that
    went extinct is ``e<id>``, not ``n<id>``. Spelling the prefix here rather than at each call site
    is deliberate: when the extinct marking was introduced every caller kept building ``n<id>``, the
    set matched nothing, and four figures silently lost their dashing."""
    extinct_leaf_names = {ct.labels()[n.id] for n in (ct.nodes[_i] for _i in ct.extinct_leaves())}
    dashed = set()
    for node in tree.walk("postorder"):
        if node.is_leaf:
            if node.name in extinct_leaf_names:
                dashed.add(node.name)
        elif node.children and all(c.name in dashed for c in node.children):
            dashed.add(node.name)
    return dashed


# --- matplotlib companion panels (precise time-aligned compositing) --------

def _extent(present: float, canvas_w: int = COMP_W, margin: int = MARGIN):
    """Map the tree image's pixels to time so px(t) == t: the tree spans ``canvas_w - 2*margin``
    px for ``present`` time units, so each margin is ``margin * present / (canvas_w - 2*margin)``
    time units of overhang."""
    off = margin * present / (canvas_w - 2 * margin)
    return -off, present + off


def composite_below(tree_png: str, present: float, out: str, panel, ylabel: str, *,
                    tree_w: int = COMP_W, margin: int = MARGIN, figsize=(10, 8.2),
                    height_ratios=(3.1, 1.0), axis_fontsize: float | None = None) -> None:
    """Tree on top, a time-indexed matplotlib panel below, sharing the exact time axis. ``tree_w`` /
    ``margin`` must match the pixel size the tree PNG was rendered at (so the axes line up)."""
    img = mpimg.imread(tree_png)
    fig, (axt, axp) = plt.subplots(2, 1, figsize=figsize, height_ratios=list(height_ratios),
                                   sharex=True, gridspec_kw={"hspace": 0.06})
    left, right = _extent(present, tree_w, margin)
    axt.imshow(img, extent=[left, right, 0, 1], aspect="auto")
    axt.set_axis_off()
    panel(axp)
    axp.set_xlim(left, right)
    axp.set_xlabel("time")
    axp.set_ylabel(ylabel)
    if axis_fontsize:
        axp.xaxis.label.set_size(axis_fontsize)
        axp.yaxis.label.set_size(axis_fontsize)
        axp.tick_params(labelsize=axis_fontsize * 0.82)
    for spine in ("top", "right"):
        axp.spines[spine].set_visible(False)
    fig.savefig(out, dpi=125, bbox_inches="tight")
    plt.close(fig)


def composite_beside(tree_png: str, out: str, panel, figsize=(12, 6), ratios=(3, 1.15),
                     *, geometry=None, wspace: float = 0.15, inset=None) -> None:
    """Tree on the left, a standalone matplotlib panel (its own axes) on the right.

    ``geometry`` is the rendered figure's `Figure.geometry()` — pass it for a **row-aligned** panel
    (one bar per tip): the panel's y axis is set to the tree's own canvas height, so a bar drawn at
    ``y = tip.y`` sits at the height of that tip rather than merely in the same order. It must be the
    canvas height and not the PNG's, which Phylustrator renders at 2x.

    ``inset`` is ``(rect, draw_fn)`` — a small schematic drawn over the tree, as `composite_markov`
    does for a Markov chain, for a figure whose model is worth stating as a picture."""
    img = mpimg.imread(tree_png)
    fig, (axt, axp) = plt.subplots(1, 2, figsize=figsize,
                                   gridspec_kw={"width_ratios": list(ratios), "wspace": wspace})
    axt.imshow(img)
    axt.set_axis_off()
    if inset is not None:
        rect, draw_fn = inset
        ax_in = axt.inset_axes(rect)
        draw_fn(ax_in)
        ax_in.set_axis_off()
    panel(axp)
    if geometry is not None:
        axp.set_ylim(geometry.size[1], 0)      # the tree's own coordinate space, top-down
        axp.set_yticks([])
    for spine in ("top", "right"):
        axp.spines[spine].set_visible(False)
    fig.savefig(out, dpi=125, bbox_inches="tight")
    plt.close(fig)


def render_tree_for_composite(tree, out_png: str, values: dict | None = None, *,
                              width: int = COMP_W, height: int = COMP_H,
                              branch_width: float = BW) -> None:
    fig = ph.trees.plot(tree, style=ph.Style(width=width, height=height, margin=MARGIN,
                                       branch_width=branch_width))
    if values is not None:
        fig = fig + ph.trees.color_branches(values)
    fig.save(out_png)


def draw_markov(ax, states, palette, rates, *, symbol="λ") -> None:
    """Draw a continuous-time Markov chain: a coloured circle per state (labelled with the rate it
    drives) and **bidirectional** transition arrows between states."""
    import itertools
    import matplotlib.patches as mpatches

    n = len(states)
    if n == 2:
        rad, pos = 0.205, {states[0]: (0.26, 0.55), states[1]: (0.74, 0.55)}
    else:
        rad, (cx, cy, R) = 0.185, (0.5, 0.47, 0.28)
        pos = {s: (cx + R * math.cos(math.pi / 2 - 2 * math.pi * i / n),
                   cy + R * math.sin(math.pi / 2 - 2 * math.pi * i / n)) for i, s in enumerate(states)}
    mx = sum(p[0] for p in pos.values()) / n
    my = sum(p[1] for p in pos.values()) / n
    for a, b in itertools.combinations(states, 2):          # one arrow each way, on opposite sides
        p, q = pos[a], pos[b]
        dx, dy = q[0] - p[0], q[1] - p[1]
        L = math.hypot(dx, dy) or 1.0
        ux, uy, ox, oy = dx / L, dy / L, -dy / L, dx / L
        for sgn, off in ((1, 0.058), (-1, -0.058)):
            (sx, sy), (ex, ey) = (p, q) if sgn > 0 else (q, p)
            a0 = (sx + ux * rad * sgn + ox * off, sy + uy * rad * sgn + oy * off)
            b0 = (ex - ux * rad * sgn + ox * off, ey - uy * rad * sgn + oy * off)
            ax.add_patch(mpatches.FancyArrowPatch(a0, b0, arrowstyle="-|>", mutation_scale=17,
                                                  lw=1.9, color="#555", zorder=2))
    for s in states:
        x, y = pos[s]
        ax.add_patch(mpatches.Circle((x, y), rad, facecolor=palette[s], edgecolor="white",
                                     lw=2.5, zorder=3))
        ax.text(x, y, s, ha="center", va="center", color="white",
                fontsize=19 if len(s) <= 4 else 16, fontweight="bold", zorder=4)
        if n == 2:                                          # rate label below; else radially outward
            lx, ly, va = x, y - rad - 0.07, "top"
        else:
            ddx, ddy = x - mx, y - my
            dl = math.hypot(ddx, ddy) or 1.0
            lx, ly = x + ddx / dl * (rad + 0.11), y + ddy / dl * (rad + 0.11)
            va = "bottom" if ddy >= 0 else "top"
        ax.text(lx, ly, f"{symbol} = {rates[s]}", ha="center", va=va, fontsize=22, zorder=4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")


def draw_grid_markov(ax, rates, xpal, ypal, *, labels=("00", "01", "10", "11")) -> None:
    """A 2×2 CTMC over two binary characters (compound states 00/01/10/11): only single-bit-flip
    arrows, their width scaling with the rate. Each directed transition is a curved arrow; the two
    directions of an edge bow to opposite sides so they never overlap. Each circle is split
    X-half | Y-half using the two binary palettes, so the diagram is also the colour key for the tree's
    two lanes."""
    import matplotlib.patches as mpatches
    import matplotlib.patheffects as pe

    # a touch below centre: the rate on the top arrow sits above it, and at 0.72 that label was
    # pinched between the arrow and the circles. The caption at 0.015 leaves room underneath.
    pos = {labels[0]: (0.27, 0.68), labels[1]: (0.73, 0.68),
           labels[2]: (0.27, 0.24), labels[3]: (0.73, 0.24)}
    rad = 0.108
    mx = max(rates.values()) or 1.0
    curve = 0.32

    def lw(r):
        return 1.2 + 2.3 * (r / mx)          # arrow width ~ rate (slimmer than before)

    for key, r in rates.items():                     # every directed transition, bowed to its left
        a, b = key.split("->")
        p, q = pos[a], pos[b]
        dx, dy = q[0] - p[0], q[1] - p[1]
        L = math.hypot(dx, dy) or 1.0
        ux, uy, ox, oy = dx / L, dy / L, -dy / L, dx / L        # (ox, oy) = left normal of a->b
        a0 = (p[0] + ux * rad, p[1] + uy * rad)
        b0 = (q[0] - ux * rad, q[1] - uy * rad)
        ax.add_patch(mpatches.FancyArrowPatch(a0, b0, connectionstyle=f"arc3,rad={curve}",
                                              arrowstyle="-|>", mutation_scale=13, lw=lw(r),
                                              color="#555", zorder=2))
        ap = ((a0[0] + b0[0]) / 2 + ox * curve * L * 0.5, (a0[1] + b0[1]) / 2 + oy * curve * L * 0.5)
        ax.text(ap[0] + ox * 0.05, ap[1] + oy * 0.05, f"{r:g}", ha="center", va="center",
                fontsize=13.5, color="#8a8a8a", zorder=4)
    for s, (x, y) in pos.items():
        ax.add_patch(mpatches.Wedge((x, y), rad, 90, 270, facecolor=xpal[s[0]],
                                    edgecolor="white", lw=2.2, zorder=3))     # left half = X
        ax.add_patch(mpatches.Wedge((x, y), rad, -90, 90, facecolor=ypal[s[1]],
                                    edgecolor="white", lw=2.2, zorder=3))      # right half = Y
        ax.text(x, y, s, ha="center", va="center", color="#1a1a1a", fontsize=17, fontweight="bold",
                zorder=4, path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
    ax.text(0.5, 0.015, "left half = X,   right half = Y", ha="center",
            va="bottom", fontsize=11.5, color="#8a8a8a")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")


def composite_model_realization(realization_png: str, out: str, draw_model, *,
                                figsize=(13, 6.4), ratios=(1.0, 1.6)) -> None:
    """Two panels side by side, like the reference: the **model** drawn by ``draw_model(ax)`` on the
    left, a rendered **realization** PNG on the right."""
    img = mpimg.imread(realization_png)
    fig, (axm, axr) = plt.subplots(1, 2, figsize=figsize,
                                   gridspec_kw={"width_ratios": list(ratios), "wspace": 0.02})
    draw_model(axm)
    axm.set_axis_off()
    axr.imshow(img)
    axr.set_axis_off()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def composite_markov(tree_png: str, out: str, draw_fn, *, loc=(0.02, 0.09, 0.34, 0.36)) -> None:
    """Place a tree PNG as the full background and draw a Markov-chain inset in the bottom-left."""
    img = mpimg.imread(tree_png)
    fig = plt.figure(figsize=(img.shape[1] / 150, img.shape[0] / 150))
    bg = fig.add_axes([0, 0, 1, 1])
    bg.imshow(img)
    bg.set_axis_off()
    inset = fig.add_axes(loc)
    draw_fn(inset)
    inset.set_axis_off()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def _draw_key(ax, key) -> None:
    """One row's colour key, drawn into a thin axes: what the colours on the tree above it mean.

    Two shapes, because the trees come in two kinds. A **dict** ``{state: colour}`` is a discrete
    trait or a gene's presence, and each entry gets a thick line the weight of a branch. A **tuple**
    ``(cmap, low, [middle,] high)`` is a continuous trait, and gets the ramp itself with its ends
    named — a swatch per value would be meaningless there.
    """
    ax.set_axis_off()
    if isinstance(key, dict):
        handles = [plt.Line2D([], [], color=colour, lw=7, solid_capstyle="butt")
                   for colour in key.values()]
        ax.legend(handles, list(key), loc="center left", ncol=len(key), frameon=False,
                  fontsize=13, handlelength=1.5, handletextpad=0.6, columnspacing=2.2,
                  borderpad=0.0, borderaxespad=0.0)
        return
    cmap, *names = key
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    left, right = 0.055, 0.235                          # the ramp, in axes fractions
    ax.imshow([[i / 255 for i in range(256)]], aspect="auto", cmap=cmap,
              extent=(left, right, 0.32, 0.68), zorder=2)
    ax.text(left - 0.012, 0.5, names[0], ha="right", va="center", fontsize=13)
    ax.text(right + 0.012, 0.5, names[-1], ha="left", va="center", fontsize=13)
    if len(names) == 3:                                 # a diverging scale needs its middle named
        ax.text((left + right) / 2, 0.14, names[1], ha="center", va="top", fontsize=12,
                color="#555555")


def conditioned_figure(out, ct, layers, values, tipcol, diagram, *, label="genome size (genes)"):
    """The conditioning-figure layout: the tree painted by the driver, beside one bar per tip, with
    the driver·mapping·target diagram small on top.

    ``layers`` are the Phylustrator layers that colour the tree, ``values`` is ``{tip name: number}``
    for the bars, ``diagram`` the kwargs for :func:`conditioning_png`, and ``label`` names what the
    bars measure. The bar quantity is whatever the target *did* — genome size where a genome rate was
    driven, root-to-tip substitutions where the substitution rate was, inversions where the inversion
    rate was — so the driver is on the tree and its consequence beside it."""
    fig = ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
                        style=ph.Style(width=900, height=900, margin=92, branch_width=3.0))
    for layer in layers:                      # no legend on the tree — the diagram is the key
        fig = fig + layer
    fig = fig + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)
    real = out.replace(".png", "_real.png")
    ph.beside(fig, ph.genomes.bars(values, colors=tipcol, label=label,
                                   tick_size=20, label_size=26),
              width=1150, tree_fraction=0.58, footer=36).save(real)
    diag = conditioning_png(out.replace(".png", "_diag.png"), **diagram)
    fig2 = plt.figure(figsize=(12, 9.6))
    axr = fig2.add_axes([0.0, 0.0, 1.0, 0.80])
    axr.imshow(mpimg.imread(real))
    axr.set_axis_off()
    axd = fig2.add_axes([0.30, 0.80, 0.40, 0.185])
    axd.imshow(mpimg.imread(diag))
    axd.set_axis_off()
    fig2.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig2)


def root_to_tip(seqs) -> dict:
    """``{tip name: root-to-tip distance}`` in substitutions per site, off the species phylogram.

    What a driven ``substitution`` rate actually produced: every tip sits the same amount of *time*
    from the root, so any spread here is the driver's doing and nothing else."""
    from zombi2.tree import read_newick
    tree, _ = read_newick(seqs.species_phylogram["complete"], assume_extant=True)
    names, out = tree.labels(), {}
    for i in tree.leaves():
        total, node = 0.0, tree.nodes[i]
        while node is not None:
            total += node.end_time - node.birth_time
            node = tree.nodes[node.parent] if node.parent is not None else None
        out[names[i]] = total
    return out


def composite_under_diagram(out: str, diagram_png: str, rows, *, width=12.0, diagram_frac=0.42,
                            pad=0.03, gap=0.30, label=0.36, key=0.40, dpi=182) -> None:
    """The driver·modifier·target diagram on top, then one labelled panel per row.

    ``rows`` is ``[(png, label), ...]``, or ``[(png, label, key), ...]`` to put a colour key under
    the label — see :func:`_draw_key` for the two shapes a key takes. Every panel gets an axes box of
    **its own image's aspect ratio**, so nothing is letterboxed and the diagram sits on the same
    centre line as the panels under it. Saved without ``bbox_inches="tight"``: cropping to content
    pulled the crop in on the row labels, which shifted the panels off centre and clipped the bottom
    one's time axis.

    ``width``, ``gap``, ``label`` and ``key`` are inches; ``pad`` and ``diagram_frac`` are fractions
    of ``width``. The default ``dpi`` is deliberately high: these figures are read by zooming into a
    tree, and at 140 the tip labels went to mush.
    """
    diagram = mpimg.imread(diagram_png)
    imgs = [mpimg.imread(row[0]) for row in rows]
    keys = [row[2] if len(row) > 2 else None for row in rows]
    body = width * (1 - 2 * pad)
    heights = [body * im.shape[0] / im.shape[1] for im in imgs]
    tops = [label + (key if k is not None else 0.0) for k in keys]
    h_diagram = diagram_frac * width * diagram.shape[0] / diagram.shape[1]
    height = gap + h_diagram + sum(h + t + gap for h, t in zip(heights, tops))
    fig = plt.figure(figsize=(width, height))

    def box(x, y, w, h):
        return [x / width, y / height, w / width, h / height]

    y = height - h_diagram
    fig.add_axes(box(width * (1 - diagram_frac) / 2, y, diagram_frac * width,
                     h_diagram)).imshow(diagram)
    for im, h, top, k, row in zip(imgs, heights, tops, keys, rows):
        y -= top + h
        fig.add_axes(box(pad * width, y, body, h)).imshow(im)
        fig.text(pad, (y + h + top - label + 0.08) / height, row[1], fontsize=15, ha="left",
                 va="bottom")
        if k is not None:
            _draw_key(fig.add_axes(box(pad * width, y + h, body, key)), k)
        y -= gap
    for ax in fig.axes:
        ax.set_axis_off()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def composite_two_trees_panel(tree_x_png: str, tree_y_png: str, draw_panel, out: str, *,
                              x_label: str = "trait x", y_label: str = "trait y") -> None:
    """Two trees stacked on the left (one per character), a custom panel spanning both rows on the
    right — the two-trees layout of ``composite_two_trees_scatter`` but with the right panel drawn by
    ``draw_panel(ax)`` (e.g. the model's Markov chain) instead of a scatter."""
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1.25], hspace=0.10, wspace=0.08)
    for row, png, name in ((0, tree_x_png, x_label), (1, tree_y_png, y_label)):
        ax = fig.add_subplot(gs[row, 0])
        ax.imshow(mpimg.imread(png))
        ax.set_axis_off()
        ax.set_title(name, fontsize=15, loc="left")
    axp = fig.add_subplot(gs[:, 1])
    draw_panel(axp)
    axp.set_axis_off()
    fig.savefig(out, dpi=125, bbox_inches="tight")
    plt.close(fig)


def composite_two_trees_scatter(tree_x_png: str, tree_y_png: str, xs, ys, out: str) -> None:
    """Two trees stacked on the left (coloured by trait x, then y), a tip x-vs-y scatter on the right."""
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1.15], hspace=0.10, wspace=0.10)
    for row, png, name in ((0, tree_x_png, "trait x"), (1, tree_y_png, "trait y")):
        ax = fig.add_subplot(gs[row, 0])
        ax.imshow(mpimg.imread(png))
        ax.set_axis_off()
        ax.set_title(name, fontsize=14, loc="left")
    axs = fig.add_subplot(gs[:, 1])
    axs.scatter(xs, ys, c=xs, cmap="viridis", s=26, edgecolors="white", linewidths=0.4)
    axs.set_xlabel("trait x", fontsize=13)
    axs.set_ylabel("trait y", fontsize=13)
    for spine in ("top", "right"):
        axs.spines[spine].set_visible(False)
    fig.savefig(out, dpi=125, bbox_inches="tight")
    plt.close(fig)


def lineages_through_time(complete_tree):
    """(times, counts) step data — the number of lineages alive over time (a diversity skyline).

    Standing diversity starts at 1 (the stem), rises by 1 at each speciation and falls by 1 at each
    extinction; **extant** tips survive to the present and get no event (so the curve does not crash to
    zero at the right edge), and the curve is held flat out to the present."""
    present = max(n.end_time for n in complete_tree.nodes.values())
    events = []
    for node in complete_tree.nodes.values():
        if node.children:                            # speciation: +1 at the split
            events.append((node.end_time, +1))
        elif node.end_time < present - 1e-9:         # extinct tip: -1 at its death
            events.append((node.end_time, -1))
        # extant tips are alive at the present -> no event
    events.sort()
    times, counts, running = [0.0], [1], 1           # the stem lineage
    for t, delta in events:
        running += delta
        times.append(t)
        counts.append(running)
    times.append(present)                            # hold flat to the present
    counts.append(running)
    return times, counts


def viridis_norm(values):
    return cm.viridis, colors.Normalize(min(values), max(values))


# --- ZOMBI2 runs on disk (for the genome / sequence levels, read by Phylustrator) -----------

_DATA = os.path.join(os.path.dirname(__file__), "figures", "_data")


def _zombi(*args) -> None:
    """Run the ZOMBI2 CLI for one pipeline step (located on PATH, else via the module)."""
    exe = shutil.which("zombi2")
    cmd = ([exe] if exe else
           [sys.executable, "-c", "import sys;from zombi2.cli.main import main;sys.exit(main())"])
    subprocess.run(cmd + [str(a) for a in args], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


#: The file a cached run is stamped with. Named for what it holds rather than for the version it
#: used to hold; a directory carrying the old ``.zombi2-version`` simply reads as unstamped and is
#: rebuilt once, which is right — it was built under a rule that ignored the run's own parameters.
_STAMP = ".cache-key"


def _cache_key(depth: int = 2) -> str:
    """What a cached run depends on: the zombi2 that would build it, and **the helper asking for it**.

    The second half is the fix for a trap. The key used to be the version alone, so editing a run —
    a seed, a rate, ``--n-extant`` — left the cache valid and the figure regenerated from the old
    data, silently, with the new parameters sitting in the source unused. Changing the karyotype
    figure's seed from 3 to 71 did exactly that, and the only way out was deleting the directory by
    hand. The helper's source *is* the run's definition, so hashing it invalidates on any edit that
    could change what the run produces. It over-invalidates on an edit that could not — a reworded
    docstring costs one rebuild — which is the side to err on.

    ``depth`` is how far up the stack the helper sits: 2 when called from `_stale` or `_stamp`, which
    is where it is called from, so no ``*_run()`` has to remember to pass anything.
    """
    from zombi2 import __version__

    try:
        source = inspect.getsource(inspect.stack()[depth].frame)
    except (OSError, IndexError):                  # no source (a REPL): fall back to the version
        source = ""
    return f"{__version__} {hashlib.sha256(source.encode()).hexdigest()[:16]}"


def _stale(run: str) -> bool:
    """Is this cached run out of date with what would build it now — and if so, **clear it**?

    The caches under ``figures/_data/`` used to be guarded on "does the file exist", which never
    expires. Four of them survived the one-row-per-event redesign holding the *old* columns
    (``lineage`` · ``copy`` · ``parent`` · ``recipient`` · ``donor``), so the events figure read a
    format zombi2 had not written in months and the rebuild tracebacked. Stamping the cache with
    what produced it is what makes it a cache rather than a fossil — see `_cache_key` for what
    "what produced it" has to mean.

    Clearing it is the other half, and it was missing: the rebuild writes with the CLI, and a run
    directory that already holds a downstream level refuses a re-run of the level above it — "would
    leave them stale" — which is right, and meant every stale cache failed to rebuild and the
    gallery could not be regenerated at all. An expired cache is exactly the thing to throw away, so
    it is thrown away here rather than at each of the five call sites, none of which would have
    remembered.
    """
    try:
        with open(os.path.join(run, _STAMP)) as fh:
            if fh.read().strip() == _cache_key():
                return False
    except OSError:
        pass
    shutil.rmtree(run, ignore_errors=True)
    return True


def _stamp(run: str) -> str:
    with open(os.path.join(run, _STAMP), "w") as fh:
        fh.write(_cache_key())
    return run


def ordered_run() -> str:
    """A cached species + ordered-genomes + sequences run (25 extant genomes, D/T/L + inversions)."""
    run = os.path.join(_DATA, "ordered")
    if _stale(run):
        _zombi("species", run, "--birth", 1.0, "--death", 0.25, "--n-extant", 25, "--seed", 7)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 45,
               "--duplication", 0.22, "--loss", 0.18, "--transfer", 0.1, "--inversion", 0.6,
               "--seed", 19)
        _zombi("sequences", run, "--model", "jc69", "--length", 60, "--seed", 7)
    return _stamp(run)


def karyotype_run() -> str:
    """A cached run for the karyotype figure: seven survivors from three circular chromosomes, with
    enough fission and fusion that the tips end with different karyotypes.

    ``--initial-families`` is high on purpose. A karyotype ring is only legible when a chromosome
    holds tens of genes: at a dozen families spread over several chromosomes each gene takes a third
    of its ring, and the picture is a handful of enormous wedges rather than a genome.

    The seed is chosen too, from a scan of 160. It has to give a spread of chromosome numbers — a run
    where every tip keeps three says nothing about fission and fusion — with no chromosome so short
    it draws as a few dashes, and no two marks close enough on one branch to overlap. This one ends
    with one, two and three chromosomes across the tips, its smallest holds thirty-one genes, and
    each of its four marks accounts for the tips below it."""
    run = os.path.join(_DATA, "karyotype")
    if _stale(run):
        _zombi("species", run, "--birth", 1.0, "--death", 0.0, "--n-extant", 7, "--seed", 71)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 150,
               "--chromosomes", 3, "--topology", "circular",
               "--duplication", 0.05, "--loss", 0.05,
               "--fission", 0.35, "--fusion", 0.35, "--seed", 71)
    return _stamp(run)


def synteny_tree_run() -> str:
    """A cached 30-tip run for the whole-clade synteny figure.

    Rearranged **gently** on purpose. A synteny picture says something only when the blocks are
    mostly collinear and a few have moved: shuffled hard, every row is a fresh permutation and the
    ribbons cross into noise, which reads as decoration rather than as data. At `inversion=0.10` one
    clade comes out strikingly collinear while the other is visibly rearranged, so the panel shows
    structure varying *across the tree* rather than a uniform amount of shuffling."""
    run = os.path.join(_DATA, "synteny_tree")
    if _stale(run):
        _zombi("species", run, "--birth", 1.0, "--death", 0.45, "--n-extant", 30, "--seed", 56)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 14,
               "--duplication", 0.015, "--loss", 0.015,
               "--inversion", 0.10, "--inversion-extent", 3, "--seed", 5)
    return _stamp(run)


def initial_gene_order(run: str) -> list:
    """The gene order the run started with — the ancestral arrangement a synteny figure colours
    against, read from ``initial_genome.tsv`` by column name."""
    with open(os.path.join(run, "genomes", "initial_genome.tsv"), encoding="utf-8") as f:
        ix = {name: k for k, name in enumerate(f.readline().rstrip("\n").split("\t"))}
        return [line.rstrip("\n").split("\t")[ix["family"]] for line in f if line.strip()]


def initial_karyotype(run: str):
    """The genome the run started from, as a Phylustrator genome — the karyotype at the root.

    ``read_genomes`` reads the survivors; this reads ``initial_genome.tsv``, which is the same table
    for time zero. A karyotype figure needs it: without the starting number of chromosomes on the
    page, a tip with three of them could as easily be three fissions as none at all."""
    path = os.path.join(run, "genomes", "initial_genome.tsv")
    by_chrom: dict = {}
    with open(path, encoding="utf-8") as f:
        ix = {name: k for k, name in enumerate(f.readline().rstrip("\n").split("\t"))}
        for line in f:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            by_chrom.setdefault((row[ix["chromosome"]], row[ix["topology"]]), []).append(row)
    chroms = []
    for (cid, topology), rows in by_chrom.items():
        genes = [ph.genomes.Gene(family=r[ix["family"]], strand=int(r[ix["strand"]]), position=i)
                 for i, r in enumerate(rows)]
        chroms.append(ph.genomes.Chromosome(cid, genes, topology=topology))
    return ph.genomes.Genome("initial", chroms)


def events_run() -> str:
    """A cached high-speciation run (40 extant genomes) with plenty of D/T/L for the events figure."""
    run = os.path.join(_DATA, "events")
    if _stale(run):
        _zombi("species", run, "--birth", 1.8, "--death", 0.3, "--n-extant", 40, "--seed", 7)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 45,
               "--duplication", 0.3, "--loss", 0.2, "--transfer", 0.15, "--inversion", 0.5,
               "--seed", 11)
    return _stamp(run)


def aln_run() -> str:
    """A cached 20-species run (species + ordered genomes + JC69 sequences), for the alignment figure."""
    run = os.path.join(_DATA, "aln")
    if _stale(run):
        _zombi("species", run, "--birth", 1.0, "--death", 0.25, "--n-extant", 20, "--seed", 4)
        # no loss so every family survives in every genome — a full one-row-per-tip alignment
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 45,
               "--duplication", 0.04, "--loss", 0.0, "--transfer", 0.0, "--seed", 6)
        # --divergence keeps the alignment from saturating (else every column varies — no signal)
        _zombi("sequences", run, "--model", "jc69", "--length", 60, "--divergence", 0.4, "--seed", 7)
    return _stamp(run)


def protein_run() -> str:
    """A cached 20-species run whose genes hold **proteins**: LG, 60 residues, one family per genome.

    The same shape as `aln_run`, one model apart, so the two cards differ in the alphabet and nothing
    else. Protein models belong to a family or ordered run: a nucleotide genome is measured in base
    pairs and read on either strand, and amino acids have no complement."""
    run = os.path.join(_DATA, "protein")
    if _stale(run):
        _zombi("species", run, "--birth", 1.0, "--death", 0.25, "--n-extant", 20, "--seed", 4)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 45,
               "--duplication", 0.04, "--loss", 0.0, "--transfer", 0.0, "--seed", 6)
        _zombi("sequences", run, "--model", "lg", "--length", 60, "--divergence", 0.7, "--seed", 7)
    return _stamp(run)


def mycoplasma_gff() -> str:
    """Path to a real *Mycoplasma genitalium* G37 GFF3 (NCBI GCF_000027325.1), downloaded + cached."""
    dest = os.path.join(_DATA, "mycoplasma.gff")
    if not os.path.isfile(dest):
        import gzip
        url = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/"
               "GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz")
        gz = dest + ".gz"
        subprocess.run(["curl", "-sSL", "--max-time", "60", url, "-o", gz], check=True)
        with gzip.open(gz, "rb") as f:
            with open(dest, "w") as o:
                o.write(f.read().decode())
        os.remove(gz)
    return dest


def branch_of(token: str) -> str:
    """The branch out of a participant token: ``n3_g467`` -> ``n3``, ``e11_g14`` -> ``e11``.

    The event log names participants, not columns: a gene copy carries the branch it sat on inside
    its own token, which is why there is no ``lineage``/``donor``/``recipient`` column to read."""
    return token.split("_", 1)[0]


def family_events(run: str, tree, kinds=("duplication", "transfer", "loss")) -> tuple:
    """Pick the gene family with the richest history on the tree's lineages, and return
    ``(family_id, events)`` for a branch_events layer. Point events are ``{"kind","node","x"}``;
    transfers are ``{"kind":"transfer","donor","recipient","x"}``.

    Reads the one-row-per-event log: ``time · kind · family · parents · children``, participants
    ``n<species>_g<copy>`` and ``;``-packed. A transfer's children are **donor-side first**, so the
    two ends of the edge come off that one cell."""
    import csv
    from collections import defaultdict

    names = {n.name for n in tree.walk() if n.name}
    byfam = defaultdict(list)
    with open(os.path.join(run, "genomes", "genome_events.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fam = r["family"]
            k = "transfer" if r["kind"].startswith("transfer") else r["kind"]
            if k not in kinds:
                continue
            kids = [t for t in r["children"].split(";") if t]
            born = [t for t in r["parents"].split(";") if t]
            if k == "transfer":
                if len(kids) == 2:
                    donor, recipient = branch_of(kids[0]), branch_of(kids[1])
                    if donor in names and recipient in names:
                        byfam[fam].append({"kind": "transfer", "donor": donor,
                                           "recipient": recipient, "x": float(r["time"])})
            else:
                where = branch_of((kids or born)[0]) if (kids or born) else None
                if where in names:
                    byfam[fam].append({"kind": k, "node": where, "x": float(r["time"])})

    def score(f):
        evs = byfam[f]
        has_transfer = any(e["kind"] == "transfer" for e in evs)   # prefer a family with a transfer
        return (has_transfer, len({e["kind"] for e in evs}), -abs(len(evs) - 16))

    fam = max(byfam, key=score)
    return fam, byfam[fam]


def rearranged_pair(genomes: dict) -> tuple:
    """Two extant genomes that share many families but in the most different order (visible synteny)."""
    import itertools

    def order(g):
        return [str(x.family) for x in g.genes]

    names = sorted(genomes)
    best = None
    for a, b in itertools.combinations(names, 2):
        fa, fb = order(genomes[a]), order(genomes[b])
        shared = set(fa) & set(fb)
        if len(shared) < 8:
            continue
        pos = {f: i for i, f in enumerate([f for f in fb if f in shared])}
        seq = [pos[f] for f in fa if f in shared]
        inv = sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq)) if seq[i] > seq[j])
        if best is None or inv > best[0]:
            best = (inv, a, b)
    return best[1], best[2]


# --- the conditioning diagram (driver · mapping · target), the manual's figure -----------

_INK, _DIM = "#1a1a1a", "#6e6e6e"


#: The conditioning diagram, to one standard. Every one of these figures answers the same three
#: questions in the same three places — what is read, how it is read, what it changes — so two of them
#: side by side can be compared rather than deciphered.
#:
#: DRIVER      the level it comes from, the thing, and what kind of value it offers; its own dynamics
#:             below it, as a state chain, when it has states
#: CONNECTION  the verb, the kind of mapping, and the mapping itself — drawn as a curve when it is one,
#:             because "curve" as a word cannot tell a saturating response from a humped one
#: TARGET      the level, and each parameter it drives with that parameter's kind: a RATE (and its
#:             scope), an EXTENT, or a CHOICE. That line is the one thing the old diagram never said.
#:
#: No colour except a driver's own states, which take the palette the tree beside them uses. Every
#: measurement is in one coordinate system: mixing figure fractions with data units is what made an
#: earlier version drift whenever a height changed.
_BOX_H, _PAD, _MIN_W, _ARROW = 94, 26, 150, 250
_TOP, _HEAD = 56, 30


def _text_w(fig, ax, s, size, **kw):
    """The drawn width of a string in data units — measured, never estimated from character counts.

    Falls back to an estimate when there is no real renderer to ask. ``tests/test_gallery_api.py``
    runs every example with matplotlib mocked out, and a mock cannot be measured or compared; without
    the fallback the canary reports a zombi2 API break that is really a drawing call."""
    try:
        t = ax.text(0, -9999, s, fontsize=size, **kw)
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        width = float(abs(bb.transformed(ax.transData.inverted()).x1
                          - bb.transformed(ax.transData.inverted()).x0))
        t.remove()
        return width
    except (TypeError, AttributeError, ValueError):
        return len(s) * size * 0.6


def _chain(ax, states, arcs, colours, *, cx, y, span):
    """A driver's own dynamics: one circle per state, one arc per transition that exists.

    ``arcs`` is ``[(forward, back), …]`` and either may be ``None`` — a gene family is lost and never
    regained, and that chain has to draw one arrow, not two with a blank label."""
    from matplotlib.patches import Circle, FancyArrowPatch
    step = min(132, max(96, span - 20))
    xs = [cx - step*(len(states)-1)/2 + i*step for i in range(len(states))]
    for st, x, col in zip(states, xs, colours):
        ax.add_patch(Circle((x, y), 15, facecolor=col, edgecolor=_INK, lw=1.2))
        ax.text(x, y+58, st, ha="center", va="center", color=_DIM, fontsize=9)
    for i in range(len(states)-1):
        a, b = xs[i], xs[i+1]
        forward, back = arcs[i]
        if forward is not None:
            ax.add_patch(FancyArrowPatch((a+17, y-4), (b-17, y-4), arrowstyle="-|>", mutation_scale=8,
                                         lw=1.0, color=_DIM, connectionstyle="arc3,rad=-0.42"))
            ax.text((a+b)/2, y-37, forward, ha="center", va="center", color=_DIM, fontsize=8.5)
        if back is not None:
            ax.add_patch(FancyArrowPatch((b-17, y+4), (a+17, y+4), arrowstyle="-|>", mutation_scale=8,
                                         lw=1.0, color=_DIM, connectionstyle="arc3,rad=-0.42"))
            ax.text((a+b)/2, y+36, back, ha="center", va="center", color=_DIM, fontsize=8.5)


def conditioning_png(path, *, driver, connection, target_level, targets,
                     chain=None, curve=None, target_chain=None, returns=False):
    """Draw one conditioning diagram to the standard above, and return the path.

    ``driver``      ``(level, name, kind)`` — "traits", "habitat", "two states"
    ``connection``  ``(verb, mapping kind)`` — "scaled_by", "table"
    ``targets``     ``[(name, kind, mapping)]`` — one entry per parameter this driver reaches.
                    ``kind`` is what it IS: "rate · per copy", "extent · in genes", "choice · a weight".
                    ``mapping`` is drawn under the arrow, because the mapping is part of the
                    *connection*, not of the target.
    ``chain``       ``(states, [(forward, back)], colours)`` for a driver with states
    ``curve``       ``(fn, x label, (lo, hi))`` for a mapping that is a curve
    ``target_chain``the same, when the thing being driven is itself a trait with states
    ``returns``     draw a second arrow, running back from the target to the driver. That one arrow
                    is the whole difference between a conditioned run and a **joint** one: the driver
                    is not grown first and read, it is grown *by* what it drives. Everything else is
                    deliberately the same picture, because everything else is the same.
    """
    from matplotlib.patches import FancyArrowPatch, Rectangle

    n = len(targets)
    verb, mapping_kind = connection

    probe = plt.figure(figsize=(10, 3), dpi=180)          # measure before committing to a canvas
    pax = probe.add_axes([0, 0, 1, 1]); pax.set_xlim(0, 1200); pax.set_ylim(400, 0)
    dw = max(_MIN_W, _PAD*2 + max(_text_w(probe, pax, driver[1], 17),
                                  _text_w(probe, pax, driver[2], 10.5, style="italic"),
                                  _text_w(probe, pax, driver[0], 9, style="italic")))
    if n == 1:
        tw = max(_MIN_W, _PAD*2 + max(_text_w(probe, pax, targets[0][0], 17),
                                      _text_w(probe, pax, targets[0][1], 10.5, style="italic"),
                                      _text_w(probe, pax, target_level, 9, style="italic")))
    else:
        tw = max(_MIN_W, _PAD*2 + max(max(_text_w(probe, pax, a, 15.5),
                                          _text_w(probe, pax, b, 9.5, style="italic"))
                                      for a, b, _ in targets))
    plt.close(probe)

    th = _BOX_H if n == 1 else 30 + 50*n + 14
    mid = _TOP + max(_BOX_H, th)/2
    bottom = mid + max(_BOX_H, th)/2
    W = 30 + dw + _ARROW + tw + 30
    H = max(bottom, (mid + 24 + 74 + 26) if curve else 0) + (118 if (chain or target_chain) else 18)

    fig = plt.figure(figsize=(W/100, H/100), dpi=180)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_axis_off()
    dx, tx = 30, 30 + dw + _ARROW
    dcx, tcx = dx + dw/2, tx + tw/2
    a0, a1 = dx + dw + 14, tx - 14
    ccx = (a0 + a1)/2
    if not returns:
        for cx, label in ((dcx, "DRIVER"), (ccx, "CONNECTION"), (tcx, "TARGET")):
            ax.text(cx, _HEAD, label, ha="center", va="center", color=_DIM, fontsize=11.5,
                    fontweight="bold")
    else:
        # neither box is only a driver or only a target here — each is the other's — so the roles
        # are not labelled, and the pair of arrows says it instead
        ax.text(ccx, _HEAD, "EACH DRIVES THE OTHER", ha="center", va="center", color=_DIM,
                fontsize=11.5, fontweight="bold")

    ax.add_patch(Rectangle((dx, mid-_BOX_H/2), dw, _BOX_H, facecolor="#f2f2f0", edgecolor=_INK, lw=1.5))
    ax.text(dcx, mid-27, driver[0], ha="center", va="center", color=_DIM, fontsize=9, style="italic")
    ax.text(dcx, mid+1, driver[1], ha="center", va="center", color=_INK, fontsize=17)
    ax.text(dcx, mid+29, driver[2], ha="center", va="center", color=_DIM, fontsize=10.5, style="italic")

    out_y = mid - 13 if returns else mid
    ax.add_patch(FancyArrowPatch((a0, out_y), (a1, out_y), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.6, color=_INK))
    if returns:
        ax.add_patch(FancyArrowPatch((a1, mid + 13), (a0, mid + 13), arrowstyle="-|>",
                                     mutation_scale=14, lw=1.6, color=_INK))
    ax.text(ccx, mid-30, verb, ha="center", va="center", color=_INK, fontsize=13.5, family="monospace")
    if not returns:
        ax.text(ccx, mid-12, mapping_kind, ha="center", va="center", color=_DIM, fontsize=10.5,
                style="italic")

    if n == 1:
        name, kind, mapping = targets[0]
        ax.add_patch(Rectangle((tx, mid-_BOX_H/2), tw, _BOX_H, facecolor="#f2f2f0",
                               edgecolor=_INK, lw=1.5))
        ax.text(tcx, mid-27, target_level, ha="center", va="center", color=_DIM, fontsize=9,
                style="italic")
        ax.text(tcx, mid+1, name, ha="center", va="center", color=_INK, fontsize=17)
        ax.text(tcx, mid+29, kind, ha="center", va="center", color=_DIM, fontsize=10.5, style="italic")
        if mapping:
            ax.text(ccx, mid + (30 if returns else 16), mapping, ha="center", va="top", color=_DIM,
                    fontsize=10)
    else:
        ax.add_patch(Rectangle((tx, mid-th/2), tw, th, facecolor="#f2f2f0", edgecolor=_INK, lw=1.5))
        ax.text(tcx, mid-th/2+17, target_level, ha="center", va="center", color=_DIM, fontsize=9,
                style="italic")
        for i, (name, kind, _) in enumerate(targets):
            y = mid - th/2 + 30 + 50*i
            ax.text(tcx, y+16, name, ha="center", va="center", color=_INK, fontsize=15.5)
            ax.text(tcx, y+35, kind, ha="center", va="center", color=_DIM, fontsize=9.5, style="italic")
        for i, (name, _, mapping) in enumerate(targets):
            if mapping:
                ax.text(ccx, mid+14+20*i, f"{name}:  {mapping}", ha="center", va="top",
                        color=_DIM, fontsize=9.5)

    if curve is not None:
        fn, xlabel, (lo_x, hi_x) = curve
        bx, by, bw, bh = ccx-70, mid+24, 140, 74
        ax.add_patch(Rectangle((bx, by), bw, bh, facecolor="white", edgecolor="#bcbcbc", lw=1.0))
        xs = [lo_x + (hi_x-lo_x)*i/199 for i in range(200)]
        vals = [fn(v) for v in xs]
        lo, hi = min(vals), max(vals)
        span = (hi - lo) or 1.0
        ax.plot([bx + (v-lo_x)/(hi_x-lo_x)*bw for v in xs],
                [by + bh - (v-lo)/span*(bh-10) - 5 for v in vals], color=_INK, lw=1.7)
        if lo <= 1.0 <= hi:                     # the line where the driver changes nothing
            y1 = by + bh - (1.0-lo)/span*(bh-10) - 5
            ax.plot([bx, bx+bw], [y1, y1], color=_DIM, lw=0.8, ls=":")
        ax.text(bx+bw/2, by+bh+14, xlabel, ha="center", va="center", color=_DIM, fontsize=9)
        ax.text(bx-10, by+bh/2, "factor", ha="center", va="center", color=_DIM, fontsize=9, rotation=90)

    if chain is not None:
        _chain(ax, *chain, cx=dcx, y=bottom+52, span=dw)
    if target_chain is not None:
        _chain(ax, *target_chain, cx=tcx, y=bottom+52, span=tw)

    fig.savefig(path, dpi=180, transparent=True, bbox_inches="tight")
    plt.close(fig)
    return path
