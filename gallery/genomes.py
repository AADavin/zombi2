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
    cmap = matplotlib.colormaps["tab20"]
    return {str(i): matplotlib.colors.to_hex(cmap((i % 20) / 19)) for i in range(n)}


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
    style = ph.Style(width=size, height=size, margin=int(size * 0.14), gene_stroke_width=1.0)
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
    from zombi2.rates.mapping import Between
    from zombi2.species import simulate_species_tree

    sp = simulate_species_tree(birth=1.0, n_extant=30, seed=11)
    ct = sp.complete_tree
    A, B = _pick_two_clades(ct)
    highway = Clades({"A": A, "B": B}, Between({("A", "B"): 10.0, ("B", "A"): 10.0}, default=1.0))
    g = simulate_genomes_family(ct, initial_families=50, transfer=0.4, duplication=0.1, loss=0.1,
                                transfer_to=highway, seed=7)

    grp = resolve_groups(ct, {"A": A, "B": B})
    trans = [e for e in g.events if e.kind == "transfer" and e.recipient is not None]
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

g = ph.zombi.read_genomes("run/genomes")["n55"]            # the genome with the most genes
(ph.genomes.plot(g, layout="circular")
 + ph.genomes.genes(by="family")).save("ring.png")'''

_C_SYNTENY = _SIM_ORDERED + '''

### plot  —  Phylustrator
import phylustrator as ph

G = ph.zombi.read_genomes("run/genomes")
(ph.genomes.stack([G["n37"], G["n54"]])              # one genome per row
 + ph.genomes.genes(by="family")
 + ph.genomes.synteny()).save("synteny.png")         # ribbons link same-family genes'''

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
for r in csv.DictReader(open("run/genomes/genome_events.tsv"), delimiter="\\t"):
    if r["family"] != "1":
        continue
    if r["kind"] == "transfer" and r["recipient"]:      # the recipient row carries donor + recipient
        events.append({"kind": "transfer", "donor": r["donor"],
                       "recipient": r["recipient"], "x": float(r["time"])})
    elif r["kind"] in ("duplication", "loss"):
        events.append({"kind": r["kind"], "node": r["lineage"], "x": float(r["time"])})
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
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_ordered

ct = simulate_species_tree(birth=1.0, n_extant=2, seed=1).complete_tree
r = simulate_genomes_ordered(ct, initial_families=12, inversion=0.5, inversion_extent=4, seed=SEED)
before = r.initial_genome[0]            # the starting ring
inv = r.rearrangements[0]               # its first inversion: a (start, length) block
# "after" = that block reversed with strands flipped — ZOMBI2's inversion operator

### plot  —  Phylustrator circular; highlight() before genes() sits the band behind them
import phylustrator as ph

for genome in (before, after):
    (ph.genomes.plot(genome, layout="circular")
     + ph.genomes.highlight(genome, start=inv.start, end=inv.start + inv.length - 1)
     + ph.genomes.genes(by="family")).save(...)'''

_C_HIGHWAY = '''\
### simulate  —  transfers steered to run BETWEEN two clades (a HGT highway)
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family, Clades
from zombi2.rates.mapping import Between

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
transfers = [e for e in g.events if e.kind == "transfer" and e.recipient is not None]
counts = Counter((grp[e.donor], grp[e.recipient]) for e in transfers)'''


EXAMPLES = [
    Example("genome_circular_ordered", "Circular genome (ordered)",
            "A genome as a ring — genes evenly spaced by rank, coloured by family, arrows by strand. "
            "<code>plot(g,&nbsp;layout=&quot;circular&quot;)&nbsp;+&nbsp;genes()</code>.",
            "phylustrator · circular", circular_ordered, code=_C_CIRC_ORD),
    Example("genome_synteny", "Synteny between two genomes",
            "Two genomes, one per row; ribbons link same-family genes and cross where the order was "
            "rearranged. <code>stack([a,b])&nbsp;+&nbsp;synteny()</code>.",
            "phylustrator · synteny", synteny_pair, code=_C_SYNTENY),
    Example("genome_tree_events", "Gene-family events on the tree",
            "One family's history on the species tree: duplications (squares), losses (crosses) and "
            "transfers (arrows, donor→recipient). <code>plot(tree)&nbsp;+&nbsp;branch_events(…)</code>.",
            "phylustrator · events", tree_with_events, code=_C_EVENTS),
    Example("genome_tree_profiles", "Profile copy-number",
            "A family&nbsp;×&nbsp;genome copy-number heatmap, its rows locked to the tips. "
            "<code>beside(tree,&nbsp;heatmap(profiles))</code>.",
            "phylustrator", tree_with_profiles, code=_C_PROFILES),
    Example("genome_circular_nucleotide", "Real genome (Mycoplasma)",
            "A real bacterium — <i>Mycoplasma genitalium</i>, 546 genes at their true base positions, "
            "coloured by strand; the forward/reverse switch marks the replication origin. "
            "<code>read_gff(…)</code>.",
            "phylustrator · real GFF", circular_nucleotide, code=_C_CIRC_NUC),
    Example("genome_inversion", "An inversion, before → after",
            "One inversion on a circular genome: the affected segment is reversed and its strands flip "
            "(the arrows turn round). The band marks the segment in both rings. "
            "<code>highlight(g,&nbsp;start,&nbsp;end)&nbsp;+&nbsp;genes()</code>.",
            "phylustrator · circular", genome_inversion, code=_C_INVERSION),
    Example("genome_transfer_highway", "A transfer highway between clades",
            "Transfers steered to run <i>between</i> two clades (a <code>Clades&nbsp;+&nbsp;Between</code> "
            "kernel — topology, not a trait). Tree coloured by clade; the barplot counts transfers by "
            "clade pair, so A↔B towers over within-clade.",
            "clades · transfer_to", transfer_highway, code=_C_HIGHWAY),
]
