"""Genome-level examples: gene arrangements, synteny, and gene-family history on the species tree.

Trees are drawn with ``ph.trees``; genomes, synteny and the copy-number heatmap with ``ph.genomes``
— both from Phylustrator. The ``code`` on each example is copy-paste reproducible: the exact ZOMBI2 CLI run (seeds and all),
then the plotting, split by a ``###`` divider.
"""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph

# a square style for the ring figures
_RING = ph.Style(width=560, height=560, margin=42, gene_stroke_width=0.8)


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
    style = ph.Style(width=660, height=660, margin=50, gene_stroke_width=0.3)
    (ph.genomes.plot(genome, layout="circular", coordinates="nucleotide", style=style)
     + ph.genomes.genes(by="strand", palette={"1": "#3a7ca5", "-1": "#c1443c"})
     + ph.genomes.position_axis()).save(out)                                # inner bp coordinate ring


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
]
