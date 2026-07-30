"""Sequence-level examples: the tree the sequences evolved down, and an alignment beside it.

The tree is drawn with ``ph.trees``; the alignment grid with ``ph.genomes`` — both from Phylustrator.
"""

from __future__ import annotations

import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import helpers as h
from helpers import Example

import phylustrator as ph


def species_phylogram(out):
    run = h.phylo_run()
    # the clock tree: branch lengths in substitutions/site (non-ultrametric under a relaxed clock)
    tree = ph.trees.read(run + "/sequences/clock_species_tree_extant.nwk")
    style = ph.Style(width=1250, height=1050, margin=82, branch_width=1.6)
    (ph.trees.plot(tree, style=style)
     + ph.trees.tip_labels()
     + ph.trees.note("uncorrelated relaxed clock (ByLineage)", loc="top-left", size=22)
     + ph.trees.time_axis("substitutions / site", tick_size=20, label_size=26)).save(out)


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


# --- an autocorrelated-clock phylogram, branches coloured by lineage rate ------------------

def _clock_rates(ct, clock_tree) -> dict:
    """Per-branch lineage clock rate = (branch length in subs/site) / (branch duration in time). With
    no extinction the clock tree and the time tree share node ids, so this is exact per branch."""
    rates = {}
    for n in clock_tree.walk():
        if not n.name or not n.name.startswith("n"):
            continue
        nid = int(n.name[1:])
        dt = ct.nodes[nid].end_time - ct.nodes[nid].birth_time
        if dt > 0:
            rates[n.name] = n.length / dt
    return rates


def autocorrelated_phylogram(out):
    from zombi2.genomes import simulate_genomes_family
    from zombi2.rates import modifiers as mod
    from zombi2.sequences import hky85, simulate_sequences
    from zombi2.species import simulate_species_tree

    # a rich Yule tree (pure birth) so the clock tree and time tree share node ids exactly
    sp = simulate_species_tree(birth=1.0, n_extant=120, seed=7)
    ct = sp.complete_tree
    g = simulate_genomes_family(ct, initial_families=1, duplication=0.0, loss=0.0, seed=9)
    seqs = simulate_sequences(g, model=hky85(kappa=2.0), length=600,
                              substitution=1.0 * mod.FromParent(spread=0.6), seed=7)
    tree = ph.trees.loads(seqs.species_phylogram["extant"])
    rates = _clock_rates(ct, tree)
    style = ph.Style(width=1500, height=784, margin=80, branch_width=1.5)
    (ph.trees.plot(tree, style=style)
     + ph.trees.color_branches(rates, cmap="viridis")
     + ph.trees.colorbar("subs/site per unit time", loc="bottom-left", size=18)
     + ph.trees.note("autocorrelated relaxed clock", loc="top-right", size=22)
     + ph.trees.time_axis("substitutions / site", tick_size=20, label_size=26)).save(out)


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
    seqs = simulate_sequences(g, model=jc69(), length=24, seed=7)
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


_C_AUTOCORR = '''\
### simulate  —  a Yule tree, sequences under the AUTOCORRELATED clock (FromParent)
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import simulate_sequences, hky85
from zombi2.rates import modifiers as mod

sp = simulate_species_tree(birth=1.0, n_extant=120, seed=7)          # pure birth: no extinction
ct = sp.complete_tree
g = simulate_genomes_family(ct, initial_families=1, seed=9)
seqs = simulate_sequences(g, model=hky85(kappa=2), length=600,
                          substitution=1.0 * mod.FromParent(spread=0.6), seed=7)  # rate drifts parent->child

### plot  —  the clock tree, branches coloured by lineage rate (= clock length / time)
import phylustrator as ph

tree = ph.trees.loads(seqs.species_phylogram["extant"])
rate = {f"n{i}": ct... }        # subs/site per unit time, per branch
(ph.trees.plot(tree)
 + ph.trees.color_branches(rate, cmap="viridis")
 + ph.trees.colorbar("subs/site per unit time", loc="bottom-left")
 + ph.trees.time_axis("substitutions / site")).save("phylogram.png")'''

_C_ANCESTRAL = '''\
### simulate  —  a small tree (20 tips), one gene family, JC69 sequences
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences import simulate_sequences, jc69

sp = simulate_species_tree(birth=1.0, n_extant=20, seed=12)
g = simulate_genomes_family(sp.complete_tree, initial_families=1, seed=5)
seqs = simulate_sequences(g, model=jc69(), length=24, seed=7)

### the reconstructed sequence at every internal node, keyed to the phylogram
seqs.founding[fam]              # the INITIAL genome (node 0), at the start of the stem
seqs.ancestral[fam]["n0_g0"]    # the CROWN (node 1); n1_g1, n2_g2, ... the other internal nodes
# plotted: the tree with numbered internal nodes, beside one coloured row per node
# — the rows are NOT tip-aligned, they belong to internal nodes'''

_C_PHYLO = '''\
### simulate  —  a relaxed clock, so branch lengths are substitutions (non-ultrametric)
zombi2 species   run --birth 1.4 --death 0.2 --n-extant 35 --seed 7
zombi2 genomes   run --resolution ordered --initial-families 40 --duplication 0.15 --loss 0.12 --seed 9
zombi2 sequences run --model hky85 --kappa 2 --length 500 \\
                     --substitution "1.0 * ByLineage(spread=0.6)" --seed 7

### plot  —  the clock tree (branch lengths in substitutions/site)
import phylustrator as ph

tree = ph.trees.read("run/sequences/clock_species_tree_extant.nwk")
(ph.trees.plot(tree)
 + ph.trees.tip_labels()
 + ph.trees.note("uncorrelated relaxed clock (ByLineage)", loc="top-left", size=22)
 + ph.trees.time_axis("substitutions / site", tick_size=20, label_size=26)).save("phylogram.png")'''

_C_ALN = '''\
### simulate  —  20 species, JC69 sequences (no loss, so a full alignment)
zombi2 species   run --birth 1.0 --death 0.25 --n-extant 20 --seed 4
zombi2 genomes   run --resolution ordered --initial-families 45 \\
                     --duplication 0.04 --loss 0.0 --transfer 0.0 --seed 6
zombi2 sequences run --model jc69 --length 60 --seed 7

### plot  —  Phylustrator (ph.trees for the tree, ph.genomes for the alignment)
import phylustrator as ph

tree = ph.trees.read("run/species/species_extant.nwk")
aln = ph.zombi.read_alignment("run", family=0)     # single-copy in every genome (one row per tip)
fig = ph.trees.plot(tree)                          # no leaf labels
ph.beside(fig, ph.genomes.alignment(aln, letters=False),   # colour blocks + a nucleotide key
          footer=70).save("alignment.png")'''


EXAMPLES = [
    Example("seq_phylogram", "Sequence phylogram",
            "The clock tree the sequences evolve down — branch lengths are substitutions/site under an "
            "uncorrelated relaxed clock, so the tips are <i>not</i> level.",
            "phylustrator · phylogram", species_phylogram, code=_C_PHYLO),
    Example("seq_phylogram_autocorr", "Autocorrelated-clock phylogram",
            "The other clock we ship: under the <b>autocorrelated</b> clock the rate drifts parent→child, "
            "so related lineages share a rate — branches coloured by lineage rate move in blocks, not "
            "salt-and-pepper. <code>substitution&nbsp;=&nbsp;FromParent(spread)</code>.",
            "phylustrator · phylogram", autocorrelated_phylogram, code=_C_AUTOCORR),
    Example("seq_ancestral", "Ancestral sequences at the nodes",
            "A small tree with its internal nodes numbered (0&nbsp;=&nbsp;initial genome, 1&nbsp;=&nbsp;crown, "
            "…); beside it the reconstructed sequence at each — one free-floating row per node, "
            "<i>not</i> aligned to the tips. <code>seqs.ancestral</code>.",
            "phylustrator · ancestral", numbered_ancestral, code=_C_ANCESTRAL),
    Example("seq_alignment", "Alignment beside the tree",
            "A single-copy family across 20 species, residues coloured (with a nucleotide key), each "
            "row locked to its tip. <code>beside(tree,&nbsp;alignment(aln))</code>.",
            "phylustrator", alignment_beside_tree, code=_C_ALN),
]
