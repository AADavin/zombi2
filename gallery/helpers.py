"""Shared rendering utilities for the ZOMBI2 examples-gallery preview.

Each example is an :class:`Example` — id, title, caption, tag, and a ``render(out_png)`` function that
simulates with ZOMBI2 and plots with Phylustrator (optionally compositing a matplotlib companion
panel). Level modules (``traits.py``, ``species.py``, …) each expose an ``EXAMPLES`` list; ``build.py``
renders them all and regenerates ``gallery.html``.
"""

from __future__ import annotations

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


def dashed_extinct(tree, extinct_leaf_names: set) -> set:
    """The set of node names whose whole subtree is extinct (their branches draw dashed)."""
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
    """Map the tree image's pixels to time so px(t) == t (see the derivation in the gallery README)."""
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


def composite_beside(tree_png: str, out: str, panel, figsize=(12, 6), ratios=(3, 1.15)) -> None:
    """Tree on the left, a standalone matplotlib panel (its own axes) on the right."""
    img = mpimg.imread(tree_png)
    fig, (axt, axp) = plt.subplots(1, 2, figsize=figsize,
                                   gridspec_kw={"width_ratios": list(ratios), "wspace": 0.15})
    axt.imshow(img)
    axt.set_axis_off()
    panel(axp)
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


def draw_markov(ax, states, palette, switch, rates, *, symbol="λ") -> None:
    """Draw a continuous-time Markov chain: a coloured circle per state (labelled with its driving
    rate), **bidirectional** transition arrows between states, and the switch rate."""
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
    ax.text(0.5, 0.005, f"switch  q = {switch}", ha="center", va="bottom", fontsize=15, color="#666")
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

    pos = {labels[0]: (0.27, 0.72), labels[1]: (0.73, 0.72),
           labels[2]: (0.27, 0.28), labels[3]: (0.73, 0.28)}
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


def ordered_run() -> str:
    """A cached species + ordered-genomes + sequences run (25 extant genomes, D/T/L + inversions)."""
    run = os.path.join(_DATA, "ordered")
    if not os.path.isdir(os.path.join(run, "sequences", "alignments")):
        _zombi("species", run, "--birth", 1.0, "--death", 0.25, "--n-extant", 25, "--seed", 7)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 45,
               "--duplication", 0.22, "--loss", 0.18, "--transfer", 0.1, "--inversion", 0.6,
               "--seed", 19)
        _zombi("sequences", run, "--model", "jc69", "--length", 60, "--seed", 7)
    return run


def events_run() -> str:
    """A cached high-speciation run (40 extant genomes) with plenty of D/T/L for the events figure."""
    run = os.path.join(_DATA, "events")
    if not os.path.isfile(os.path.join(run, "genomes", "genome_events.tsv")):
        _zombi("species", run, "--birth", 1.8, "--death", 0.3, "--n-extant", 40, "--seed", 7)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 45,
               "--duplication", 0.3, "--loss", 0.2, "--transfer", 0.15, "--inversion", 0.5,
               "--seed", 11)
    return run


def phylo_run() -> str:
    """A cached run whose sequences evolve under an **uncorrelated relaxed clock** (ByLineage), so the
    clock tree (branch lengths in substitutions/site) is non-ultrametric. 35 species, for the phylogram."""
    run = os.path.join(_DATA, "phylo")
    if not os.path.isfile(os.path.join(run, "sequences", "clock_species_tree_extant.nwk")):
        _zombi("species", run, "--birth", 1.4, "--death", 0.2, "--n-extant", 35, "--seed", 7)
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 40,
               "--duplication", 0.15, "--loss", 0.12, "--seed", 9)
        _zombi("sequences", run, "--model", "hky85", "--kappa", 2.0, "--length", 500,
               "--substitution", "1.0 * ByLineage(spread=0.6)", "--seed", 7)
    return run


def aln_run() -> str:
    """A cached 20-species run (species + ordered genomes + JC69 sequences), for the alignment figure."""
    run = os.path.join(_DATA, "aln")
    if not os.path.isdir(os.path.join(run, "sequences", "alignments")):
        _zombi("species", run, "--birth", 1.0, "--death", 0.25, "--n-extant", 20, "--seed", 4)
        # no loss so every family survives in every genome — a full one-row-per-tip alignment
        _zombi("genomes", run, "--resolution", "ordered", "--initial-families", 45,
               "--duplication", 0.04, "--loss", 0.0, "--transfer", 0.0, "--seed", 6)
        # --divergence keeps the alignment from saturating (else every column varies — no signal)
        _zombi("sequences", run, "--model", "jc69", "--length", 60, "--divergence", 0.4, "--seed", 7)
    return run


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


def family_events(run: str, tree, kinds=("duplication", "transfer", "loss")) -> tuple:
    """Pick the gene family with the richest history on the tree's lineages, and return
    ``(family_id, events)`` for a branch_events layer. Point events are ``{"kind","node","x"}``;
    transfers are ``{"kind":"transfer","donor","recipient","x"}`` (from the recipient row)."""
    import csv
    from collections import defaultdict

    names = {n.name for n in tree.walk() if n.name}
    byfam = defaultdict(list)
    with open(os.path.join(run, "genomes", "genome_events.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k, fam = r["kind"], r["family"]
            if k not in kinds:
                continue
            if k == "transfer":
                if r["recipient"] and r["donor"] in names and r["recipient"] in names:
                    byfam[fam].append({"kind": "transfer", "donor": r["donor"],
                                       "recipient": r["recipient"], "x": float(r["time"])})
            elif r["lineage"] in names:
                byfam[fam].append({"kind": k, "node": r["lineage"], "x": float(r["time"])})

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


# --- the conditioning diagram (driver · modifier · target), the manual's figure -----------

def draw_conditioning(ax, *, driver, states, switch, mapping, target, target_base=None,
                      target_sub=None, symbol="×", target_run="genome", state_colors=None):
    """Reproduce the manual's driver→modifier→target diagram. ``switch`` is a {"a->b": rate} dict (an
    arrow is drawn for each positive rate, so an irreversible trait shows one arrow); ``mapping`` is the
    per-state multiplier; ``state_colors`` tints the state nodes to match the tree's palette. The state
    names sit *outside* (below) their circles. ``target_base`` is the rate's base value (``None`` for a
    target that is not a rate, e.g. a recipient-choice slot); ``target_sub`` overrides the italic
    caption under TARGET."""
    from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

    ax.set_xlim(0, 660)
    ax.set_ylim(250, 0)                       # y grows downward, like the SVG
    ax.set_aspect("auto")
    ax.set_axis_off()
    ink, dim = "#1a1a1a", "#6e6e6e"

    for x, t in ((120, "DRIVER"), (345, "MODIFIER"), (566, "TARGET")):
        ax.text(x, 30, t, ha="center", va="center", color=dim, fontsize=12.5, fontweight="bold")

    ax.add_patch(Rectangle((45, 96), 150, 60, fill=True, facecolor="#f2f2f0", edgecolor=ink,
                           lw=1.6, joinstyle="round"))
    ax.text(120, 126, driver, ha="center", va="center", color=ink, fontsize=16)

    n = len(states)
    r_st = 15
    xs = [72 + (168 - 72) * (k / (n - 1) if n > 1 else 0.5) for k in range(n)]
    pos = dict(zip(states, xs))
    for k, s in enumerate(states):
        x = xs[k]
        col = (state_colors or {}).get(s, "#c9c9c9")
        ax.add_patch(Circle((x, 202), r_st, facecolor=col, edgecolor=ink, lw=1.2))
        ax.text(x, 202 + r_st + 11, s, ha="center", va="top", color=ink, fontsize=10.5)
    for key, rate in switch.items():
        if rate <= 0:
            continue
        a, b = key.split("->")
        xa, xb = pos[a], pos[b]
        inner_l, inner_r = min(xa, xb) + r_st, max(xa, xb) - r_st
        if xa < xb:
            ax.add_patch(FancyArrowPatch((inner_l, 196), (inner_r, 196),
                                         connectionstyle="arc3,rad=-0.6", arrowstyle="-|>",
                                         mutation_scale=9, lw=1.1, color=ink))
        else:
            ax.add_patch(FancyArrowPatch((inner_r, 208), (inner_l, 208),
                                         connectionstyle="arc3,rad=-0.6", arrowstyle="-|>",
                                         mutation_scale=9, lw=1.1, color=ink))

    ax.add_patch(FancyArrowPatch((203, 126), (486, 126), arrowstyle="-|>", mutation_scale=15,
                                 lw=1.7, color=ink))
    ax.text(345, 112, "DrivenBy", ha="center", va="center", color=ink, fontsize=14.5, style="italic")
    for i, s in enumerate(states):
        ax.text(345, 152 + i * 19, f"{s} {symbol} {mapping.get(s, 1)}", ha="center", va="center",
                color=dim, fontsize=13.5)

    ax.add_patch(Rectangle((496, 96), 140, 60, fill=True, facecolor="#f2f2f0", edgecolor=ink,
                           lw=1.6, joinstyle="round"))
    if target_base is None:
        ax.text(566, 126, target, ha="center", va="center", color=ink, fontsize=15)
    else:
        ax.text(566, 120, target, ha="center", va="center", color=ink, fontsize=15)
        ax.text(566, 142, f"base {target_base}", ha="center", va="center", color=dim, fontsize=13)


def conditioning_png(path, **kw):
    """Render :func:`draw_conditioning` to its own (transparent) PNG, so it can be placed small and
    undistorted above a realization."""
    fig, ax = plt.subplots(figsize=(9.5, 3.5))
    draw_conditioning(ax, **kw)
    fig.savefig(path, dpi=180, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return path
