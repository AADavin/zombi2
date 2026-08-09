"""Sequence-level examples: the tree the sequences evolved down, and an alignment beside it.

The tree is drawn with ``ph.trees``; the alignment grid with ``ph.genomes`` — both from Phylustrator.
"""

from __future__ import annotations

import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import helpers as h
from zombi2.params import Drift, LogNormal, PerSite
from helpers import Example

import phylustrator as ph
from zombi2.params.distributions import Gamma


# --- the relaxed clocks, one card each: one species tree, one shared colour scale ------------

_CLOCK_SEED, _CLOCK_TIPS, _CLOCK_DIV = 11, 70, 1.0

#: name -> the clock's shape. `divergence` solves for the base, so all four are calibrated to the
#: same root-to-tip divergence and only their *pattern* of rate variation differs — which is the
#: whole comparison. `bins` is the discrete-bin (rate-category) form of the autocorrelated clock.
_CLOCKS = {
    "ucln":  PerSite().varying_among('lineages', LogNormal(0.0, 0.55)),
    "ugam":  PerSite().varying_among('lineages', Gamma(shape=3.31, scale=0.302)),
    "auto":  PerSite().varying_among('lineages', Drift(LogNormal(0.0, 0.4))),
    "bins":  PerSite().varying_among('lineages', Drift(LogNormal(0.0, 0.45), bins=6)),
}

_CLOCK_CACHE: dict = {}


def _clock_panels() -> dict:
    """Simulate all four clocks down ONE species tree and return
    ``{name: (tree, {branch: log10 rate})}`` plus the shared colour ``limits``.

    Built together and cached because the four cards must share a colour scale: coloured
    independently, each would normalise to its own min and max and the same green would mean a
    different rate on every card. The range is clipped to the 2nd-98th percentile so one fast branch
    cannot flatten the rest, and it is logarithmic because these clocks are lognormal."""
    import math
    from zombi2.genomes import simulate_genomes_family
    from zombi2.sequences import hky85, simulate_sequences
    from zombi2.species import simulate_species_tree

    if _CLOCK_CACHE:
        return _CLOCK_CACHE
    sp = simulate_species_tree(birth=1.0, n_extant=_CLOCK_TIPS, seed=_CLOCK_SEED)
    ct = sp.complete_tree                       # pure birth, so clock and time trees share node ids
    g = simulate_genomes_family(ct, initial_families=1, duplication=0.0, loss=0.0, seed=_CLOCK_SEED)
    panels, every = {}, []
    for name, shape in _CLOCKS.items():
        seqs = simulate_sequences(g, model=hky85(kappa=2.0), length=600,
                                  substitution=shape, divergence=_CLOCK_DIV, seed=_CLOCK_SEED)
        tree = ph.trees.loads(seqs.species_phylogram["extant"])
        logr = {}
        for n in tree.walk():
            if n.name and n.name.startswith("n") and n.length:
                nd = ct.nodes[int(n.name[1:])]
                dt = nd.end_time - nd.birth_time
                if dt > 0:
                    logr[n.name] = math.log10(n.length / dt)
        panels[name] = (tree, logr)
        every += list(logr.values())
    every.sort()
    lo, hi = every[int(0.02 * len(every))], every[int(0.98 * len(every)) - 1]
    _CLOCK_CACHE.update(panels=panels, limits=(lo, hi))
    return _CLOCK_CACHE


def _clock_figure(name, out):
    """One clock's phylogram: branches coloured by lineage rate, on the scale all four share.

    No time axis — the card is about the *pattern* of rate variation, and the four trees are drawn on
    one branch-length scale, so their widths stay comparable without one."""
    got = _clock_panels()
    tree, logr = got["panels"][name]
    style = ph.Style(width=1250, height=1150, margin=70, branch_width=1.9)
    (ph.trees.plot(tree, style=style)
     + ph.trees.color_branches(logr, cmap="viridis", limits=got["limits"])
     + ph.trees.colorbar("lineage rate  (subs/site per unit time)", loc="bottom-left",
                         width=210, height=14, size=20,
                         labels=tuple(f"{10 ** v:.2f}" for v in got["limits"]))).save(out)


def clock_ucln(out):
    _clock_figure("ucln", out)


def clock_ugam(out):
    _clock_figure("ugam", out)


def clock_autocorrelated(out):
    _clock_figure("auto", out)


def clock_discrete_bin(out):
    _clock_figure("bins", out)


def _best_single_copy(M) -> str:
    """The family that is single-copy (exactly one gene) in the most genomes — the cleanest
    one-row-per-tip alignment."""
    best, best_n = M.cols[0], -1
    for j, c in enumerate(M.cols):
        n = sum(1 for row in M.values if row[j] == 1)
        if n > best_n:
            best, best_n = c, n
    return best


def alignment_beside_tree(out):
    run = h.aln_run()                                    # 20 species
    tree = ph.trees.read(run + "/species/species_extant.nwk")
    fam = _best_single_copy(ph.zombi.read_profiles(run))
    aln = ph.zombi.read_alignment(run, fam)
    fig = ph.trees.plot(tree, style=h.style())                            # no leaf labels
    ph.beside(fig, ph.genomes.alignment(aln, letters=False, legend=True),   # no letters; nucleotide key
              width=1150, tree_fraction=0.30, footer=70).save(out)  # no title


# --- a small tree with numbered nodes, beside the ancestral sequence at each -----------------

_NT = {"A": "#4E9F50", "C": "#3B7DD8", "G": "#F2A93B", "T": "#D75455", "-": "#dddddd", "N": "#bbbbbb"}


def _seq_row(ax, y, seq, cell_w, cell_h):
    for j, ch in enumerate(seq):
        ax.add_patch(mpatches.Rectangle((j * cell_w, y), cell_w, cell_h,
                                        facecolor=_NT.get(ch, "#cccccc"), edgecolor="white", lw=0.4))


def _root_pixel(tree, style, margin):
    """The stem's left-end pixel (the initial genome): x = margin, y = the root's pixel-y. ph places
    each internal node at the mean of its children, so averaging tip pixel-ys up the tree gives it."""
    fig = ph.trees.plot(tree, style=style)
    tip_y = {t.name: t.y for t in fig.geometry().tips}
    kids = {n.name: [c.name for c in n.children] for n in tree.walk()}

    def yof(name):
        return tip_y[name] if name in tip_y else sum(yof(c) for c in kids[name]) / len(kids[name])

    return margin, yof(next(n for n in tree.walk() if n.is_root).name)


def numbered_ancestral(out):
    from zombi2.genomes import simulate_genomes_family
    from zombi2.sequences import jc69, simulate_sequences
    from zombi2.species import simulate_species_tree

    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=20, seed=12)   # a longer stem
    ct = sp.complete_tree
    g = simulate_genomes_family(ct, initial_families=1, duplication=0.0, loss=0.0, seed=5)
    seqs = simulate_sequences(g, model=jc69(), length=24, divergence=0.35, seed=7)  # keep it unsaturated
    fam = next(iter(seqs.ancestral))

    internal = sorted((i for i in ct.nodes if ct.nodes[i].children),
                      key=lambda i: (ct.nodes[i].birth_time, i))       # root (crown) first
    numbering = {i: k + 1 for k, i in enumerate(internal)}             # node id -> number
    tree = ph.trees.loads(ct.to_newick())
    for n in tree.walk():
        if n.name and n.name.startswith("n"):
            n.name = "" if n.is_leaf else f"{numbering[int(n.name[1:])]}"

    margin = 92
    style = ph.Style(width=900, height=1500, margin=margin, branch_width=2.6)
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(tree, style=style)
     + ph.trees.node_labels(size=28, color="#333333", offset=15)
     + ph.trees.time_axis("time", tick_size=18, label_size=22)).save(tree_png)
    mx, my = _root_pixel(ph.trees.loads(ct.to_newick()), style, margin)

    rows = [("0", "initial genome", seqs.founding[fam], "#7a5cc0")]
    for i in internal:
        rows.append((str(numbering[i]), "crown" if i == ct.root else "",
                     seqs.ancestral[fam][f"n{i}_g{i}"], "#111111"))

    fig = plt.figure(figsize=(12.5, 11.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.95], wspace=0.03)
    axt = fig.add_subplot(gs[0, 0])
    img = mpimg.imread(tree_png)
    axt.imshow(img)
    sx, sy = img.shape[1] / style.width, img.shape[0] / style.height
    mx, my = mx * sx, my * sy
    axt.text(mx, my, "0", ha="center", va="center", color="white", fontsize=17, fontweight="bold",
             zorder=6, bbox=dict(boxstyle="circle,pad=0.5", fc="#7a5cc0", ec="white", lw=2))
    axt.text(mx, my - 104, "initial\ngenome", ha="center", va="bottom", fontsize=13,
             color="#7a5cc0", fontweight="bold")
    axt.set_xlim(0, img.shape[1])
    axt.set_ylim(img.shape[0], 0)
    axt.set_axis_off()

    axr = fig.add_subplot(gs[0, 1])
    L = len(rows[0][2])
    cell_w, cell_h, gap = 1.0, 0.72, 0.34
    for r, (num, name, seq, col) in enumerate(rows):
        y = (len(rows) - 1 - r) * (cell_h + gap)
        cy = y + cell_h / 2
        _seq_row(axr, y, seq, cell_w, cell_h)
        axr.text(-1.3, cy, num, ha="center", va="center", color="white", fontsize=11,
                 fontweight="bold", zorder=3, bbox=dict(boxstyle="circle,pad=0.32", fc=col, ec="none"))
        if name:
            axr.text(-2.9, cy, name, ha="right", va="center", fontsize=12, color=col, fontweight="bold")
    axr.set_xlim(-11, L * cell_w + 0.5)
    axr.set_ylim(-1.9, len(rows) * (cell_h + gap))
    for k, b in enumerate("ACGT"):                                     # nucleotide key (no label)
        axr.add_patch(mpatches.Rectangle((k * 2.0, -1.75), 0.7, 0.5, facecolor=_NT[b],
                                         edgecolor="white", lw=0.4))
        axr.text(k * 2.0 + 0.95, -1.5, b, ha="left", va="center", fontsize=11, color="#444")
    axr.set_axis_off()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


_C_ANCESTRAL = '''\
### simulate  —  a small tree (20 tips), one gene family, JC69 sequences
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import simulate_sequences, jc69

sp = simulate_species_tree(birth=1.0, n_extant=20, seed=12)
ct = sp.complete_tree
g = simulate_genomes_family(ct, initial_families=1, seed=5)
seqs = simulate_sequences(g, model=jc69(), length=24, divergence=0.35, seed=7)   # unsaturated
fam = next(iter(seqs.ancestral))
seqs.founding[fam]              # the INITIAL genome (node 0), at the start of the stem
seqs.ancestral[fam]["n0_g0"]    # the CROWN (node 1); n1_g1, n2_g2, ... the other internal nodes

### plot  —  the tree with internal nodes numbered, beside the ancestral sequence at each
import phylustrator as ph

internal = sorted((i for i in ct.nodes if ct.nodes[i].children), key=lambda i: ct.nodes[i].birth_time)
number = {i: k + 1 for k, i in enumerate(internal)}     # root (crown) = 1, then by depth
tree = ph.trees.loads(ct.to_newick())
for n in tree.walk():                                   # relabel internal nodes with their number
    if n.name and n.name.startswith("n"):
        n.name = "" if n.is_leaf else str(number[int(n.name[1:])])
(ph.trees.plot(tree) + ph.trees.node_labels() + ph.trees.time_axis("time")).save("tree.png")
# beside it, a matplotlib colour grid: one free-floating row per numbered node — the founding sequence
# (node 0) and each seqs.ancestral[fam][...] — NOT tip-aligned (they belong to internal nodes)'''

_C_ALN = '''\
### simulate  —  20 species, JC69 sequences (no loss, so a full alignment)
zombi2 species   run --birth 1.0 --death 0.25 --n-extant 20 --seed 4
zombi2 genomes   run --resolution ordered --initial-families 45 \\
                     --duplication 0.04 --loss 0.0 --transfer 0.0 --seed 6
zombi2 sequences run --model jc69 --length 60 --divergence 0.4 --seed 7   # --divergence: don't saturate

### plot  —  Phylustrator (ph.trees for the tree, ph.genomes for the alignment)
import phylustrator as ph

tree = ph.trees.read("run/species/species_extant.nwk")
aln = ph.zombi.read_alignment("run", family=2)     # single-copy in every genome (one row per tip)
fig = ph.trees.plot(tree)                          # no leaf labels
ph.beside(fig, ph.genomes.alignment(aln, letters=False),   # colour blocks + a nucleotide key
          footer=70).save("alignment.png")'''


_C_CLOCKS = '''\
### simulate  —  one species tree, four clocks, all calibrated to the same divergence
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import simulate_sequences, hky85
from zombi2.params import Drift, Gamma, LogNormal, PerSite

sp = simulate_species_tree(birth=1.0, n_extant=70, seed=11)     # pure birth: no extinction
ct = sp.complete_tree
g = simulate_genomes_family(ct, initial_families=1, seed=11)

clocks = {"uncorrelated lognormal": PerSite().varying_among('lineages', LogNormal(0.0, 0.55)),
          "uncorrelated gamma":     PerSite().varying_among('lineages', Gamma(shape=3.31, scale=0.302)),
          "autocorrelated":         PerSite().varying_among('lineages', Drift(LogNormal(0.0, 0.4))),
          "discrete-bin":           PerSite().varying_among('lineages', Drift(LogNormal(0.0, 0.45), bins=6))}

# the clock's SHAPE alone, with divergence solving for the base — so the four differ only in
# their pattern of rate variation, not in how far the sequences ran
seqs = simulate_sequences(g, model=hky85(kappa=2), length=600,
                          substitution=clocks["discrete-bin"], divergence=1.0, seed=11)

### plot  —  the clock tree, branches coloured by lineage rate (clock length / elapsed time)
import math
import phylustrator as ph

tree = ph.trees.loads(seqs.species_phylogram["extant"])
rate = {n.name: math.log10(n.length / (ct.nodes[int(n.name[1:])].end_time
                                       - ct.nodes[int(n.name[1:])].birth_time))
        for n in tree.walk() if n.name and n.name.startswith("n") and n.length}
# limits= fixes the colour range, so several of these share one scale and a colour means
# the same rate on each; without it every figure normalises to its own min and max
(ph.trees.plot(tree, style=ph.Style(width=1250, height=1150, branch_width=1.9))
 + ph.trees.color_branches(rate, cmap="viridis", limits=(-1.31, -0.33))
 + ph.trees.colorbar("lineage rate  (subs/site per unit time)", loc="bottom-left",
                     labels=("0.05", "0.47"))).save("clock.png")'''


EXAMPLES = [
    Example("clock_ucln", "Uncorrelated lognormal clock",
            "Every lineage draws its own rate, with no memory of its parent, so the colour is "
            "salt-and-pepper. <code>substitution&nbsp;=&nbsp;PerSite().varying_among('lineages',&nbsp;LogNormal(0.0, 0.55))</code>.",
            "phylustrator · clocks", clock_ucln, code=_C_CLOCKS),
    Example("clock_ugam", "Uncorrelated gamma clock",
            "The same independent draw with a gamma instead of a lognormal. "
            "<code>varying_among('lineages',&nbsp;Gamma(shape=3.31,&nbsp;scale=0.302))</code>.",
            "phylustrator · clocks", clock_ugam, code=_C_CLOCKS),
    Example("clock_autocorrelated", "Autocorrelated clock",
            "A daughter starts at its parent's rate and is nudged, so the colour moves in <b>clades</b> "
            "rather than branch to branch. <code>substitution&nbsp;=&nbsp;PerSite().varying_among('lineages',&nbsp;Drift(LogNormal(0.0, 0.4)))</code>.",
            "phylustrator · clocks", clock_autocorrelated, code=_C_CLOCKS),
    Example("clock_discrete_bin", "Discrete-bin clock",
            "The same inherited drift in <b>steps</b>: the rate takes one of a few values and a daughter "
            "moves to a neighbouring one. <code>varying_among('lineages',&nbsp;Drift(LogNormal(0.0, 0.45),&nbsp;bins=6))</code>.",
            "phylustrator · clocks", clock_discrete_bin, code=_C_CLOCKS),
    Example("seq_ancestral", "Ancestral sequences at the nodes",
            "A small tree with its internal nodes numbered, and beside it the sequence at each. The rows "
            "are one per node, <i>not</i> aligned to the tips.",
            "phylustrator · ancestral", numbered_ancestral, code=_C_ANCESTRAL),
    Example("seq_alignment", "Alignment beside the tree",
            "A single-copy family across 20 species, residues coloured (with a nucleotide key), each "
            "row locked to its tip. <code>beside(tree,&nbsp;alignment(aln))</code>.",
            "phylustrator", alignment_beside_tree, code=_C_ALN),
]
