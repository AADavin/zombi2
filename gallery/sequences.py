"""Sequence-level examples: the tree the sequences evolved down, and an alignment beside it.

The tree is drawn with ``ph.trees``; the alignment grid with ``ph.genomes`` — both from Phylustrator.
"""

from __future__ import annotations

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
    Example("seq_alignment", "Alignment beside the tree",
            "A single-copy family across 20 species, residues coloured (with a nucleotide key), each "
            "row locked to its tip. <code>beside(tree,&nbsp;alignment(aln))</code>.",
            "phylustrator", alignment_beside_tree, code=_C_ALN),
]
