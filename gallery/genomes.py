"""Genome-level examples: gene arrangements, synteny, and gene-family history on the species tree.

Trees are drawn with ``ph.trees``; genomes, synteny and the copy-number heatmap with ``ph.genomes``
— both from Phylustrator. The ``code`` on each example is copy-paste reproducible: the exact ZOMBI2 CLI run (seeds and all),
then the plotting, split by a ``###`` divider.
"""

from __future__ import annotations

import matplotlib
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2.params import LogNormal, PerCopy, Random

# a square style for the ring figures — the classic thin genes on a solid backbone (Adrián's preference
# for the dense real-genome-style rings; the chunky "arrow" style is reserved for the sparse inversion figure)
_RING = ph.Style(width=560, height=560, margin=42, gene_stroke_width=0.8, gene_style="wedge")


def _fullest(genomes: dict):
    return max(genomes.values(), key=lambda g: len(g.genes))


def circular_ordered(out):
    G = ph.zombi.read_genomes(h.ordered_run() + "/genomes")
    (ph.genomes.plot(_fullest(G), layout="circular", style=_RING) + ph.genomes.genes(by="family")).save(out)


def synteny_pair(out):
    G = ph.zombi.read_genomes(h.ordered_run() + "/genomes")
    a, b = h.rearranged_pair(G)
    (ph.genomes.stack([G[a], G[b]]) + ph.genomes.genes(by="family") + ph.genomes.synteny(opacity=0.42)).save(out)


def synteny_on_the_tree(out):
    """Every tip's gene order beside the tree, homologues ribboned between neighbouring genomes.

    The pair figure above shows two genomes; this shows the whole clade at once, which is where a
    synteny picture starts saying something about *history* — a collinear block that survives down
    one clade and breaks in another.

    Genes are coloured by their position in the **ancestral** order, so a genome still in that order
    reads as a clean gradient and every rearrangement is a break in it. Categorical colours would
    make the reader match hues instead of reading the shape."""
    run = h.synteny_tree_run()
    tree = ph.trees.read(run + "/species/species_extant.nwk")
    genomes = ph.zombi.read_genomes(run)
    style = ph.Style(width=1240, height=840, margin=72, branch_width=2.4)
    panel = ph.genomes.tracks(list(genomes.values()), reference=h.initial_gene_order(run),
                              cmap="viridis", opacity=0.38)
    ph.beside(ph.trees.plot(tree, style=style), panel,
              width=1240, height=840, tree_fraction=0.27, gap=22, pad=34).save(out)


def tree_with_events(out):
    run = h.events_run()                                          # 40 species, high speciation
    tree = ph.trees.read(run + "/species/species_extant.nwk")
    fam, events = h.family_events(run, tree)
    style = ph.Style(width=1500, height=1150, margin=95, branch_width=1.5)
    (ph.trees.plot(tree, style=style)
     + ph.trees.branch_events(events, legend_title="events",
                        legend_loc="top-left", legend_size=23, size=8)
     + ph.trees.tip_labels() + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def tree_with_profiles(out):
    run = h.ordered_run()                                        # 25 species
    tree = ph.trees.read(run + "/species/species_extant.nwk")
    M = ph.zombi.read_profiles(run)
    var = [j for j, c in enumerate(M.cols) if len({row[j] for row in M.values}) > 1]
    Mv = ph.genomes.Matrix(rows=M.rows, cols=[M.cols[j] for j in var],
                   values=[[row[j] for j in var] for row in M.values])
    fig = ph.trees.plot(tree, style=h.style())                          # no tip labels
    ph.beside(fig, ph.genomes.heatmap(Mv), width=1320, tree_fraction=0.30,
              footer=64).save(out)                                # no title; footer holds the colour key


def circular_nucleotide(out):
    genome = next(iter(ph.genomes.read_gff(h.mycoplasma_gff()).values()))   # real M. genitalium, 546 genes
    style = ph.Style(width=660, height=660, margin=50, gene_stroke_width=0.3,
                     gene_style="wedge")            # thin — 546 chunky arrows would overlap

    (ph.genomes.plot(genome, layout="circular", coordinates="nucleotide", style=style)
     + ph.genomes.genes(by="strand", palette={"1": "#3a7ca5", "-1": "#c1443c"})
     + ph.genomes.position_axis()).save(out)                                # inner bp coordinate ring


# --- one inversion, before -> after (Phylustrator's chunky circular arrows) ------------------

def _family_palette(n):
    cmap = matplotlib.colormaps["viridis"]
    return {str(i): matplotlib.colors.to_hex(cmap(i / max(n - 1, 1))) for i in range(n)}


def _to_ph_genome(zchrom, name):
    genes = [ph.genomes.Gene(family=str(g.family), strand=g.strand, position=i)
             for i, g in enumerate(zchrom.genes)]
    return ph.genomes.Genome(name, [ph.genomes.Chromosome(str(zchrom.id), genes,
                                                          topology=zchrom.topology)])


def _apply_inversion(zchrom, start, length):
    """A copy of the chromosome with genes [start, start+length) reversed and their strands flipped —
    exactly ZOMBI2's ordered inversion operator (non-wrapping blocks only)."""
    import copy

    class _G:
        def __init__(self, family, strand):
            self.family, self.strand = family, strand
    genes = [_G(g.family, g.strand) for g in zchrom.genes]
    block = [_G(g.family, -g.strand) for g in reversed(genes[start:start + length])]
    out = copy.copy(zchrom)
    out.genes = genes[:start] + block + genes[start + length:]
    return out


def _first_clean_inversion(n):
    """Grow a small ordered genome (inversions only) and return (before_chrom, first_inversion) where
    the first recorded inversion on the stem lineage is a visible, non-wrapping 3–5 gene block."""
    from zombi2.genomes import simulate_genomes_ordered
    from zombi2.species import simulate_species_tree

    ct = simulate_species_tree(birth=1.0, death=0.0, n_extant=2, seed=1).complete_tree
    for seed in range(1, 400):
        r = simulate_genomes_ordered(ct, initial_families=n, inversion=0.5, inversion_extent=4,
                                     topology="circular", seed=seed)
        invs = [rr for rr in r.rearrangements if rr.lineage == 0]
        if invs and 3 <= invs[0].length <= 5 and invs[0].start + invs[0].length <= n:
            return r.initial_genome[0], invs[0]
    raise RuntimeError("no clean first inversion found")


def _ring_png(genome, path, palette, *, highlight=None, size=620):
    # slimmer bodies so the flared head reads: each gene is an arrow, not a chunky pentagon
    style = ph.Style(width=size, height=size, margin=int(size * 0.14), gene_stroke_width=1.0,
                     ring_gene_frac=0.18)
    fig = ph.genomes.plot(genome, layout="circular", style=style)
    if highlight is not None:
        fig = fig + ph.genomes.highlight(genome, start=highlight[0], end=highlight[1])  # behind the genes
    (fig + ph.genomes.genes(by="family", palette=palette)).save(path)
    return path


def genome_inversion(out):
    n = 12
    pal = _family_palette(n)
    zchrom, inv = _first_clean_inversion(n)
    start, length = inv.start, inv.length
    before = _to_ph_genome(zchrom, "before")
    after = _to_ph_genome(_apply_inversion(zchrom, start, length), "after")
    seg = (start, start + length - 1)
    p_before = _ring_png(before, out.replace(".png", "_before.png"), pal, highlight=seg)
    p_after = _ring_png(after, out.replace(".png", "_after.png"), pal, highlight=seg)

    fig, ax = plt.subplots(figsize=(12.5, 6.4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 50)
    ax.set_aspect("equal")
    ax.set_axis_off()

    def place(png, cx, cy, w):
        img = mpimg.imread(png)
        hgt = w * img.shape[0] / img.shape[1]
        ax.imshow(img, extent=[cx - w / 2, cx + w / 2, cy - hgt / 2, cy + hgt / 2], zorder=2)

    place(p_before, 23, 25, 44)
    place(p_after, 77, 25, 44)
    ax.add_patch(mpatches.FancyArrowPatch((46, 25), (54, 25), arrowstyle="-|>", mutation_scale=24,
                                          lw=3.0, color="#333", zorder=1))
    ax.text(50, 34, "Inversion", ha="center", fontsize=17, color="#111", fontweight="bold")
    ax.text(23, 2.0, "before", ha="center", fontsize=14, color="#333")
    ax.text(77, 2.0, "after", ha="center", fontsize=14, color="#333")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --- an HGT "highway": transfers steered between two clades, counted in a barplot -------------

_CLA, _CLB, _REST = "#2C6E9E", "#D1642F", "#c9c9c9"


def _extant_tips_under(ct, i, present):
    out, st = set(), [i]
    while st:
        j = st.pop()
        ch = ct.nodes[j].children
        if ch:
            st.extend(ch)
        elif ct.nodes[j].end_time >= present - 1e-9:
            out.add(j)
    return out


def _pick_two_clades(ct, lo=6, hi=9):
    present = max(n.end_time for n in ct.nodes.values())
    internal = [i for i in ct.nodes if ct.nodes[i].children]
    sized = [(i, _extant_tips_under(ct, i, present)) for i in internal]
    cand = [(i, t) for i, t in sized if lo <= len(t) <= hi]
    a_id, a_tips = cand[0]
    b_id, _ = next((i, t) for i, t in cand if not (t & a_tips))
    return a_id, b_id


def transfer_highway(out):
    from collections import Counter

    from zombi2.genomes import Clades, simulate_genomes_family
    from zombi2.genomes._transfer import resolve_groups
    from zombi2.params.mapping import Between
    from zombi2.species import simulate_species_tree

    sp = simulate_species_tree(birth=1.0, n_extant=30, seed=11)
    ct = sp.complete_tree
    A, B = _pick_two_clades(ct)
    highway = Clades({"A": A, "B": B}, Between({("A", "B"): 10.0, ("B", "A"): 10.0}, default=1.0))
    g = simulate_genomes_family(ct, initial_families=50, transfer=0.4, duplication=0.1, loss=0.1,
                                transfer_to=highway, seed=7)

    grp = resolve_groups(ct, {"A": A, "B": B})
    trans = [e for e in g.edges if e.kind == "transfer" and e.recipient is not None]
    pair = Counter((grp[e.donor], grp[e.recipient]) for e in trans)
    counts = {"within A": pair[("A", "A")], "A → B": pair[("A", "B")],
              "B → A": pair[("B", "A")], "within B": pair[("B", "B")]}

    present = max(n.end_time for n in ct.nodes.values())
    tree = ph.trees.loads(ct.to_newick())
    labels = {f"n{i}": grp[i] for i in ct.nodes}
    style = ph.Style(width=1150, height=720, margin=72, branch_width=2.6)
    tree_png = out.replace(".png", "_tree.png")
    fig_t = (ph.trees.plot(tree, style=style)
             + ph.trees.color_branches(labels, palette={"A": _CLA, "B": _CLB, "rest": _REST})
             + ph.trees.highlight_clade(f"n{A}", color=_CLA, opacity=0.13)
             + ph.trees.highlight_clade(f"n{B}", color=_CLB, opacity=0.13)
             + ph.trees.time_axis("time", tick_size=20, label_size=26))
    fig_t.save(tree_png)
    tip_y = {t.name: t.y for t in fig_t.geometry().tips}
    tip_x = fig_t.geometry().tip_x

    def clade_y(node):
        tips = _extant_tips_under(ct, node, present)
        return sum(tip_y[f"n{t}"] for t in tips) / len(tips)

    img = mpimg.imread(tree_png)
    sx, sy = img.shape[1] / style.width, img.shape[0] / style.height
    fig = plt.figure(figsize=(13.5, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[3, 1.55], wspace=0.14)
    axt = fig.add_subplot(gs[0, 0])
    axt.imshow(img)
    for node, col, name in ((A, _CLA, "clade A"), (B, _CLB, "clade B")):
        axt.text((tip_x + 26) * sx, clade_y(node) * sy, name, ha="left", va="center",
                 fontsize=12, color=col, fontweight="bold")
    axt.set_xlim(0, img.shape[1] * 1.14)
    axt.set_ylim(img.shape[0], 0)
    axt.set_axis_off()

    axp = fig.add_subplot(gs[0, 1])
    keys = list(counts)
    vals = [counts[k] for k in keys]
    axp.bar(range(4), vals, color=[_CLA, "#7B4FA0", "#7B4FA0", _CLB], edgecolor="white", width=0.78)
    axp.set_xticks(range(4))
    axp.set_xticklabels(keys, fontsize=13)
    axp.set_ylabel("gene transfers", fontsize=15)
    top = max(vals)
    for i, v in enumerate(vals):
        axp.text(i, v + top * 0.015, str(v), ha="center", va="bottom", fontsize=12)
    axp.plot([1, 1, 2, 2], [top * 1.08, top * 1.13, top * 1.13, top * 1.08], color="#7B4FA0", lw=1.6)
    axp.text(1.5, top * 1.16, "highway  (A ↔ B)", ha="center", va="bottom", fontsize=13,
             color="#7B4FA0", fontweight="bold")
    axp.set_ylim(0, top * 1.28)
    axp.tick_params(labelsize=12)
    for spine in ("top", "right"):
        axp.spines[spine].set_visible(False)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


# --- copy-paste-reproducible snippets shown on the detail view --------------
# The ordered run these four share:
_SIM_ORDERED = '''\
### simulate  —  ZOMBI2 CLI (deterministic given the seeds)
zombi2 species run --birth 1.0 --death 0.25 --n-extant 25 --seed 7
zombi2 genomes run --resolution ordered --initial-families 45 \\
                   --duplication 0.22 --loss 0.18 --transfer 0.1 --inversion 0.6 --seed 19'''

_C_CIRC_ORD = _SIM_ORDERED + '''

### plot  —  Phylustrator
import phylustrator as ph

g = ph.zombi.read_genomes("run/genomes")["n50"]            # the genome with the most genes
(ph.genomes.plot(g, layout="circular")
 + ph.genomes.genes(by="family")).save("ring.png")'''

_C_SYNTENY = _SIM_ORDERED + '''

### plot  —  Phylustrator
import phylustrator as ph

G = ph.zombi.read_genomes("run/genomes")
(ph.genomes.stack([G["n50"], G["n55"]])              # one genome per row
 + ph.genomes.genes(by="family")
 + ph.genomes.synteny()).save("synteny.png")         # ribbons link same-family genes'''

_C_SYNTENY_TREE = '''\
### simulate  —  30 species, 14 families on one chromosome, gently rearranged
zombi2 species run --birth 1.0 --death 0.45 --n-extant 30 --seed 56
zombi2 genomes run --resolution ordered --initial-families 14 \\
                   --duplication 0.015 --loss 0.015 \\
                   --inversion 0.10 --inversion-extent 3 --seed 5

### plot  —  Phylustrator
import csv
import phylustrator as ph

tree = ph.trees.read("run/species/species_extant.nwk")
G = ph.zombi.read_genomes("run")
with open("run/genomes/initial_genome.tsv") as fh:            # the ancestral arrangement
    ancestral = [r["family"] for r in csv.DictReader(fh, delimiter="\\t")]

style = ph.Style(width=1240, height=840, margin=72, branch_width=2.4)
panel = ph.genomes.tracks(list(G.values()),      # one gene track per genome
                          reference=ancestral,   # colour by ANCESTRAL position…
                          cmap="viridis")        # …so a rearrangement is a break in the gradient
ph.beside(ph.trees.plot(tree, style=style), panel,
          width=1240, height=840, tree_fraction=0.27).save("synteny_tree.png")'''

_C_EVENTS = '''\
### simulate  —  high speciation, 40 species, with D/T/L (a different tree)
zombi2 species run --birth 1.8 --death 0.3 --n-extant 40 --seed 7
zombi2 genomes run --resolution ordered --initial-families 45 \\
                   --duplication 0.3 --loss 0.2 --transfer 0.15 --inversion 0.5 --seed 11

### plot  —  Phylustrator (duplication=square, loss=cross, transfer=arrow)
import csv
import phylustrator as ph

tree = ph.trees.read("run/species/species_extant.nwk")
events = []
branch = lambda tok: tok.split("_", 1)[0]      # n3_g467 -> n3: the copy carries its own branch
for r in csv.DictReader(open("run/genomes/genome_events.tsv"), delimiter="\\t"):
    if r["family"] != "1":
        continue
    kids = [t for t in r["children"].split(";") if t]
    if r["kind"].startswith("transfer"):               # children are donor-side first
        events.append({"kind": "transfer", "donor": branch(kids[0]),
                       "recipient": branch(kids[1]), "x": float(r["time"])})
    elif r["kind"] == "duplication":
        events.append({"kind": "duplication", "node": branch(kids[0]), "x": float(r["time"])})
    elif r["kind"] == "loss":
        events.append({"kind": "loss", "node": branch(r["parents"]), "x": float(r["time"])})
(ph.trees.plot(tree, style=ph.Style(width=1500, height=1150, branch_width=1.5))
 + ph.trees.branch_events(events, legend_title="events", legend_loc="top-left", legend_size=23, size=8)
 + ph.trees.tip_labels() + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("events.png")'''

_C_PROFILES = _SIM_ORDERED + '''

### plot  —  Phylustrator (ph.trees for the tree, ph.genomes for the heatmap)
import phylustrator as ph

tree = ph.trees.read("run/species/species_extant.nwk")
p = ph.zombi.read_profiles("run")                  # genomes x families copy number
cols = [j for j, c in enumerate(p.cols)      # keep the families that vary
        if len({row[j] for row in p.values}) > 1]
p = ph.genomes.Matrix(p.rows, [p.cols[j] for j in cols],
              [[row[j] for j in cols] for row in p.values])
fig = ph.trees.plot(tree)                          # heatmap rows lock to the tips (no labels needed)
ph.beside(fig, ph.genomes.heatmap(p), footer=64).save("profiles.png")'''

_C_CIRC_NUC = '''\
### a real genome  —  Mycoplasma genitalium G37 (NCBI GCF_000027325.1)
# fetch the GFF3, then plot it directly — no simulation
curl -sL <ncbi>/GCF_000027325.1_..._genomic.gff.gz | gunzip > mycoplasma.gff

### plot  —  Phylustrator
import phylustrator as ph

G = ph.genomes.read_gff("mycoplasma.gff")            # real genes at their real base positions
genome = next(iter(G.values()))              # one circular chromosome, 546 genes
(ph.genomes.plot(genome, layout="circular", coordinates="nucleotide")
 + ph.genomes.genes(by="strand", palette={"1": "#3a7ca5", "-1": "#c1443c"})).save("ring.png")'''


_C_INVERSION = '''\
### simulate  —  a small ordered genome, inversions only; take one recorded inversion
import copy, dataclasses
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_ordered

ct = simulate_species_tree(birth=1.0, n_extant=2, seed=1).complete_tree
r = simulate_genomes_ordered(ct, initial_families=12, inversion=0.5, inversion_extent=4, seed=3)
chrom = r.initial_genome[0]             # the starting ring
inv = r.rearrangements[0]               # its first inversion: a (start, length) block
s, L = inv.start, inv.length
# "after" = that block reversed with strands flipped — ZOMBI2's inversion operator
flipped = copy.copy(chrom)
flipped.genes = (chrom.genes[:s]
                 + [dataclasses.replace(g, strand=-g.strand) for g in reversed(chrom.genes[s:s + L])]
                 + chrom.genes[s + L:])

### plot  —  Phylustrator circular; highlight() before genes() sits the band behind them
import phylustrator as ph

def as_genome(c, name):                 # a ZOMBI2 chromosome -> a Phylustrator genome
    genes = [ph.genomes.Gene(family=str(g.family), strand=g.strand, position=i)
             for i, g in enumerate(c.genes)]
    return ph.genomes.Genome(name, [ph.genomes.Chromosome(str(c.id), genes, topology=c.topology)])

for genome in (as_genome(chrom, "before"), as_genome(flipped, "after")):
    (ph.genomes.plot(genome, layout="circular")
     + ph.genomes.highlight(genome, start=s, end=s + L - 1)
     + ph.genomes.genes(by="family")).save(genome.name + ".png")'''

_C_HIGHWAY = '''\
### simulate  —  transfers steered to run BETWEEN two clades (a HGT highway)
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family, Clades
from zombi2.params.mapping import Between

ct = simulate_species_tree(birth=1.0, n_extant=30, seed=11).complete_tree
# clades are a FACT OF THE TREE; Between weights WHO RECEIVES by (donor clade, recipient clade)
g = simulate_genomes_family(ct, initial_families=50, transfer=0.4, duplication=0.1, loss=0.1,
        transfer_to=Clades({"A": nodeA, "B": nodeB},
                           Between({("A","B"): 10, ("B","A"): 10}, default=1)), seed=7)

### plot  —  tree coloured by clade (A, B, rest); a barplot of transfer counts by clade pair
import phylustrator as ph
from zombi2.genomes._transfer import resolve_groups
from collections import Counter

grp = resolve_groups(ct, {"A": nodeA, "B": nodeB})              # each lineage's clade (A / B / rest)
labels = {f"n{i}": grp[i] for i in ct.nodes}
(ph.trees.plot(ph.trees.loads(ct.to_newick()))
 + ph.trees.color_branches(labels, palette={"A": "#2C6E9E", "B": "#D1642F", "rest": "#c9c9c9"})
 + ph.trees.highlight_clade(f"n{nodeA}") + ph.trees.highlight_clade(f"n{nodeB}")
 + ph.trees.time_axis("time")).save("tree.png")
# beside it, a matplotlib barplot of transfer counts by (donor, recipient) clade — A<->B towers:
transfers = [e for e in g.edges if e.kind == "transfer" and e.recipient is not None]
counts = Counter((grp[e.donor], grp[e.recipient]) for e in transfers)'''


# --- how much families differ from one another ----------------------------------------------------

_PANGENOME_SPREAD = 1.4
_ABSENT, _PRESENT = "#F1EFE9", "#26565B"


def _pangenome_runs():
    """Two genome runs on one species tree: every family alike, then each rate varying by family.

    Same tree, same seed, and the same mean duplication / transfer / loss rate in both — only how
    much families differ from one another changes. The Python API rather than the CLI, because the
    per-family draw is an object and the point of the figure is the one line that differs."""
    from zombi2 import genomes as zg
    from zombi2 import species as zs

    tree = zs.simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=11)
    base = dict(duplication=0.06, transfer=0.10, loss=0.30, origination=0.30,
                initial_families=200, max_family_size=6, seed=7)
    spread = Random('families', LogNormal(0.0, _PANGENOME_SPREAD))
    varied = dict(base, duplication=PerCopy(0.06).varying_among(spread),
                  transfer=PerCopy(0.10).varying_among(spread),
                  loss=PerCopy(0.30).varying_among(spread))
    return tree, [("Every family alike", "duplication=0.06\ntransfer=0.10\nloss=0.30",
                   zg.simulate_genomes_family(tree, **base)),
                  ("Each rate varies by family",
                   f"duplication=PerCopy(0.06).varying_among('families', LogNormal(0.0, {_PANGENOME_SPREAD}))\n"
                   f"transfer=PerCopy(0.10).varying_among('families', LogNormal(0.0, {_PANGENOME_SPREAD}))\n"
                   f"loss=PerCopy(0.30).varying_among('families', LogNormal(0.0, {_PANGENOME_SPREAD}))",
                   zg.simulate_genomes_family(tree, **varied))]


def pangenome_by_family(out):
    """Whether a pangenome has a core at all, decided by one parameter.

    The top row is the profile matrix with families sorted by how many genomes carry them; the
    bottom row is the gene-frequency spectrum, which is the quantity a pangenome paper plots. With
    every family alike the spectrum is a hump in the middle and **no** family is in every genome.
    Let families differ and it goes bimodal — a spike of universal families over a flat cloud, the
    U-shape real pangenomes show.

    The genome sizes differ too (90 against 120 genes), and that is not a flaw in the comparison: a
    per-copy rate compounds, so spreading it around a fixed mean raises the expected copy number."""
    import numpy as np

    tree, panels = _pangenome_runs()
    n = len(tree.complete_tree.extant_leaves())
    grid_style = ph.Style(width=760, height=1000, margin=0, background=None)
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2),
                             gridspec_kw=dict(height_ratios=[2.6, 1.15], hspace=0.40, wspace=0.20))

    for col, (title, params, result) in enumerate(panels):
        m = result.profiles.matrix
        prev = (m > 0).sum(axis=1)
        order = np.argsort(-prev)                       # commonest family at the top
        M = ph.genomes.Matrix(rows=[str(result.profiles.families[i]) for i in order],
                              cols=[str(c) for c in result.profiles.species],
                              values=[[int(v > 0) for v in m[i]] for i in order])
        png = out.replace(".png", f"_grid{col}.png")
        ph.genomes.grid(M, palette={0: _ABSENT, 1: _PRESENT}, borders=False,
                        style=grid_style).save(png)

        ax = axes[0, col]
        ax.imshow(mpimg.imread(png), aspect="auto", interpolation="antialiased")
        ax.set_title(title, fontsize=12.5, color="#16191C", pad=44, fontweight="semibold")
        ax.text(0.5, 1.012, params, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=7.4, color="#6C6F6A", family="monospace", linespacing=1.5)
        ax.set_xlabel("genomes", fontsize=9.5, color="#6C6F6A", labelpad=4)
        if col == 0:
            ax.set_ylabel("gene families\n(sorted by how many genomes carry them)",
                          fontsize=9.5, color="#6C6F6A")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#D6D5CC")
        core = int((prev >= 0.95 * n).sum())
        ax.text(0.97, 0.035, f"{core} core  \u00b7  {m.sum() / n:.0f} genes per genome",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5,
                color="#16191C" if core else "#9E3C29", fontweight="semibold",
                bbox=dict(facecolor="white", edgecolor="#D6D5CC", boxstyle="round,pad=0.34",
                          linewidth=0.7))

        ax = axes[1, col]
        hist, _ = np.histogram(prev, bins=np.arange(1, n + 2))
        ax.bar(np.arange(1, n + 1), hist / hist.sum(), width=0.86, color=_PRESENT, linewidth=0)
        ax.set_xlabel("in how many genomes", fontsize=9.5, color="#6C6F6A", labelpad=3)
        if col == 0:
            ax.set_ylabel("share of families", fontsize=9.5, color="#6C6F6A")
        ax.set_xlim(0.3, n + 0.7); ax.set_ylim(0, 0.32)
        ax.set_yticks([0, 0.15, 0.30])
        ax.tick_params(labelsize=8.5, colors="#6C6F6A", length=3)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#D6D5CC")

    fig.subplots_adjust(top=0.885, bottom=0.095, left=0.105, right=0.975)
    fig.savefig(out, dpi=170, facecolor="white")
    plt.close(fig)


_C_PANGENOME = '''### simulate  —  two runs, one species tree, the same mean rates
from zombi2 import genomes, species
from zombi2.params import LogNormal, PerCopy, Random

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=11)
base = dict(duplication=0.06, transfer=0.10, loss=0.30, origination=0.30,
            initial_families=200, max_family_size=6, seed=7)

alike  = genomes.simulate_genomes_family(tree, **base)

spread = Random('families', LogNormal(0.0, 1.4))            # one draw per family, mean unchanged
varied = genomes.simulate_genomes_family(
    tree, **dict(base, duplication=PerCopy(0.06).varying_among(spread),
                 transfer=PerCopy(0.10).varying_among(spread),
                 loss=PerCopy(0.30).varying_among(spread)))

### plot  —  the profile matrix (ph.genomes.grid) over its frequency spectrum
import numpy as np, phylustrator as ph

p = varied.profiles
prev = (p.matrix > 0).sum(axis=1)                  # genomes carrying each family
order = np.argsort(-prev)                          # commonest at the top
M = ph.genomes.Matrix(rows=[str(p.families[i]) for i in order],
                      cols=[str(c) for c in p.species],
                      values=[[int(v > 0) for v in p.matrix[i]] for i in order])
ph.genomes.grid(M, palette={0: "#F1EFE9", 1: "#26565B"}, borders=False).save("pangenome.png")

np.histogram(prev, bins=np.arange(1, len(p.species) + 2))   # the frequency spectrum'''


# --- a transition partway through the run ---------------------------------------------------------
#
# The clade is named off the TREE, and its factor is written as a schedule, so the change is scoped
# to a group of lineages AND to a time. Chaining two verbs cannot say this: scaled_by(clade, ...)
# .changing_at(...) multiplies two factors that each apply to every lineage, so the window would fall
# on the whole tree instead of on the clade.

_SW_SEED, _SW_TIPS, _SW_T0 = 2, 34, 2.0      # _SW_T0: when the clade's loss rate changes
_SW_FACTOR = 20                              # what it changes to
_SW_AFTER, _SW_BEFORE = "#b2182b", "#4393c3"


def _switch_run():
    """One tree, one genome run: a clade whose loss rate jumps twenty-fold at ``_SW_T0`` and which
    stops receiving transfers at the same moment. Nothing else in the run changes, and nothing about
    those lineages differs before that time."""
    from zombi2.genomes import simulate_genomes_family
    from zombi2.params import Clade, PerCopy, PerLineage, Recipients
    from zombi2.species import simulate_species_tree

    sp = simulate_species_tree(birth=1.0, n_extant=_SW_TIPS, seed=_SW_SEED)
    ct = sp.complete_tree
    clade = Clade({"selected": ["n13", "n33"]})
    on = {0: 1.0, _SW_T0: 20.0}                       # the schedule: one factor, from t on
    off = {0: 1.0, _SW_T0: 0.0}
    g = simulate_genomes_family(
        ct, initial_families=220, duplication=PerCopy(0.02), origination=PerLineage(1.0),
        transfer=PerCopy(0.05),
        loss=PerCopy(0.02).scaled_by(clade, {"selected": on, "rest": 1.0}),
        transfer_to=Recipients().weighted_by(clade, {"selected": off, "rest": 1.0}),
        seed=_SW_SEED)
    return ct, clade, g


def clade_transition(out):
    """One clade's loss rate changes, and only from a given time.

    The shading names the clade the run selected; the colour is each lineage's loss rate through
    time, so it changes *along* the branches the threshold crosses; the bars are what each genome is
    left with."""
    import numpy as np

    ct, clade, g = _switch_run()
    inside = set(clade.resolve(ct)["selected"])
    labels = ct.labels()
    # the clade's own root — the shading marks which lineages the run selected, before any of them
    # behaves differently. Its OWN branch is inside the clade, which is what `Clade.resolve` shows.
    clade_root = labels[min(inside, key=lambda i: ct.nodes[i].birth_time)]

    # The colour is the lineage's STATE THROUGH TIME, not which group it is in. A selected lineage is
    # an ordinary one until _SW_T0 — same loss rate, same transfers — so painting it red
    # from the clade's origin would show membership where the figure is about a change. Every branch
    # is a mosaic: the one the switch falls inside is blue up to it and red after.
    history = {}
    for i, nd in ct.nodes.items():
        a, b = nd.birth_time, nd.end_time
        if i not in inside or b <= _SW_T0:
            history[labels[i]] = [("base rate", b - a)]
        elif a >= _SW_T0:
            history[labels[i]] = [("loss x20", b - a)]
        else:                                        # the switch falls inside this branch
            history[labels[i]] = [("base rate", _SW_T0 - a),
                                  ("loss x20", b - _SW_T0)]
    palette = {"loss x20": _SW_AFTER, "base rate": _SW_BEFORE}

    fig = (ph.trees.plot(ph.trees.loads(ct.to_newick()),
                         style=ph.Style(width=1150, height=980, margin=70, branch_width=3.0))
           + ph.trees.color_history(history, palette=palette)
           + ph.trees.highlight_clade(clade_root, color="#c8b9a0", opacity=0.22)
           + ph.trees.time_marker(_SW_T0, label=f"loss x{_SW_FACTOR} from t = {_SW_T0:g}")
           + ph.trees.time_axis("time")
           + ph.trees.legend("loss rate"))
    tmp = out.replace(".png", "_tree.png")
    fig.save(tmp)

    # bars in the tree's own tip order, so row k of the panel is tip k of the figure
    order = [t.name for t in fig.geometry().tips]
    sizes = [len(g.genomes[name]) for name in order]
    colours = [palette[history[name][-1][0]] for name in order]

    geo = fig.geometry()
    tips = geo.tips
    span = (tips[-1].y - tips[0].y) / max(len(tips) - 1, 1)      # one tip's share of the height

    def panel(ax):
        ax.barh([t.y for t in tips], sizes, color=colours, height=span * 0.72, linewidth=0)
        ax.set_xlabel("genes per genome")
        ax.axvline(np.mean([n for n, c in zip(sizes, colours) if c == _SW_BEFORE]),
                   color="0.45", lw=1.0, ls=":")

    body = out.replace(".png", "_body.png")
    h.composite_beside(tmp, body, panel, figsize=(12.5, 7.4), ratios=(3, 1.0),
                       geometry=geo, wspace=0.02)
    diag = h.conditioning_png(
        out.replace(".png", "_diag.png"),
        # a clade is read off the tree rather than grown, and a schedule is its connection: the
        # switch happens to the whole clade at a moment the run fixes, not at a rate per lineage
        driver=("the tree", "one clade", "a named group"),
        connection=("changing_at", "a schedule"),
        target_level="genomes",
        targets=[("loss", "rate · per copy", f"× {_SW_FACTOR} inside the clade")],
        chain=(("base rate", f"× {_SW_FACTOR}"),
               [(f"at t = {_SW_T0:g}", None)],
               (_SW_BEFORE, _SW_AFTER)))
    h.composite_under_diagram(out, diag, [(body, "")], diagram_frac=0.46)


_C_TRANSITION = '''\
from zombi2.params import Clade, PerCopy, Recipients

clade = Clade({"selected": ["n13", "n33"]})

# one factor, scoped to the clade AND to a time: x20 loss from t=2, and no transfer in from then
my_genomes = simulate_genomes_family(
    tree, initial_families=220, duplication=PerCopy(0.02), transfer=PerCopy(0.05),
    loss=PerCopy(0.02).scaled_by(clade, {"endo": {0: 1.0, 2.0: 20.0}, "rest": 1.0}),
    transfer_to=Recipients().weighted_by(clade, {"endo": {0: 1.0, 2.0: 0.0}, "rest": 1.0}),
    seed=2)'''


EXAMPLES = [
    Example("genome_circular_ordered", "Circular genome (ordered)",
            "A genome as a ring — genes evenly spaced by rank, coloured by family, arrows by strand. "
            "<code>plot(g,&nbsp;layout=&quot;circular&quot;)&nbsp;+&nbsp;genes()</code>.",
            "phylustrator · circular", circular_ordered, code=_C_CIRC_ORD),
    Example("genome_synteny", "Synteny between two genomes",
            "Two genomes, one per row; ribbons link same-family genes and cross where the order was "
            "rearranged. <code>stack([a,b])&nbsp;+&nbsp;synteny()</code>.",
            "phylustrator · synteny", synteny_pair, code=_C_SYNTENY),
    Example("genome_synteny_tree", "Synteny across a whole clade",
            "Every tip's gene order beside the tree. Genes are coloured by their <em>ancestral</em> "
            "position, so each rearrangement is a break in the gradient.",
            "phylustrator · synteny", synteny_on_the_tree, code=_C_SYNTENY_TREE),
    Example("genome_tree_events", "Gene-family events on the tree",
            "One family's history on the species tree: duplications (squares), losses (crosses), "
            "transfers (arrows, donor→recipient).",
            "phylustrator · events", tree_with_events, code=_C_EVENTS),
    Example("genome_tree_profiles", "Profile copy-number",
            "A family&nbsp;×&nbsp;genome copy-number heatmap, its rows locked to the tips. "
            "<code>beside(tree,&nbsp;heatmap(profiles))</code>.",
            "phylustrator", tree_with_profiles, code=_C_PROFILES),
    Example("genome_circular_nucleotide", "Real genome (Mycoplasma)",
            "A real bacterium: <i>Mycoplasma genitalium</i>, 546 genes at their true base positions. "
            "The forward/reverse switch marks the replication origin.",
            "phylustrator · real GFF", circular_nucleotide, code=_C_CIRC_NUC),
    Example("genome_inversion", "An inversion, before → after",
            "One inversion on a circular genome: the segment is reversed and its strands flip. The band "
            "marks it in both rings.",
            "phylustrator · circular", genome_inversion, code=_C_INVERSION),
    Example("genome_transfer_highway", "A transfer highway between clades",
            "Transfers steered to run <i>between</i> two clades, by topology rather than by a trait. The "
            "barplot counts them by clade pair, so A↔B towers over within-clade.",
            "clades · transfer_to", transfer_highway, code=_C_HIGHWAY),
    Example("genome_pangenome_by_family", "Core and accessory, from one parameter",
            "Two runs at the same mean rates. Every family alike gives no core at all; letting families "
            "differ gives 28 core families and a U-shaped spectrum.",
            "phylustrator · heterogeneity", pangenome_by_family, code=_C_PANGENOME),
    Example("genome_clade_transition", "One clade's loss rate changes, from a given time",
            "The shading marks the clade the run selected. Its lineages lose genes at the base rate "
            "up to the dashed line and twenty times faster after it, so the colour changes <i>along</i> "
            "the branches the line crosses — nothing about them differs beforehand. The bars are what "
            "each genome is left with: about 140 genes against 270 outside, all of it lost in the last "
            "third of the run. One factor, scoped to a group <b>and</b> to a time — chaining "
            "<code>scaled_by</code> with <code>changing_at</code> cannot say this, because the two "
            "factors would each apply to every lineage. "
            "<code>scaled_by(clade,&nbsp;{'selected':&nbsp;{0:&nbsp;1.0,&nbsp;2.0:&nbsp;20.0},&nbsp;'rest':&nbsp;1.0})</code>.",
            "clades · schedules", clade_transition, code=_C_TRANSITION),
]
