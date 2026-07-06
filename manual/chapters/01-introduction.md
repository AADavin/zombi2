# Introduction

Evolutionary studies rely on statistical analyses of data observed in the present to make
inferences about what happened in the past. This creates a fundamental problem: we rarely have
access to the real history, so we cannot easily confirm that our methods are recovering it
correctly.

Simulation offers a way out. By generating artificial datasets under an explicit model, we build
worlds in which we know *everything* about how the data came to be — a solid ground truth against
which reconstruction methods can be tested, calibrated, and compared.

**ZOMBI2** is a phylogenetic simulator — a ground-up redesign of ZOMBI [@davin2020zombi] built
around a fast Rust engine, a simple command-line interface, and a composable Python library. It
offers a collection of models for simulating species trees, gene trees, and traits, independently
or jointly, and records the true history behind every dataset it produces. For example, ZOMBI2
can:

- simulate very large species trees under a coalescent point-process model, and decorate them
  *a posteriori* with extinct ("ghost") lineages;
- evolve gene families along any tree by duplication, transfer, loss and origination, recovering
  their gene trees, reconciliations, and copy-number profiles;
- descend to nucleotide resolution — structural rearrangements from which genes emerge as
  conserved "blocks", and real DNA sequences under standard substitution models;
- evolve phenotypic traits (Brownian motion, Ornstein–Uhlenbeck, Mk, threshold, biogeographic
  ranges) along the same tree;
- couple these processes into models of *coevolution* — gene content that shapes traits, traits
  that shape gene content, or gene families that evolve non-independently.

## The levels of evolution

Evolution is a single process, but it leaves its trace at several levels at once: a lineage
diversifies into a species tree; along that tree, genomes gain and lose gene families; the
organisms carry traits that drift and adapt; and inside every gene, molecular sequences
accumulate substitutions. Reality draws no line between these levels. We do — because a
tractable model usually captures only one of them at a time.

ZOMBI2 is organised around that separation. It simulates four levels: **species** (the dated
tree), **genomes** (the gene families along it), **traits** (phenotypes on the tree), and
**sequences** (the nucleotides or amino acids in each gene). There are two ways to put them
together.

![The four levels of evolution and the two ways ZOMBI2 composes them. In the hierarchical pipeline (left), each level is simulated conditioned on the one above; in a joint `coevolve` run (right), species, genomes and traits are mutually coupled. Sequences are a downstream layer in both. You always choose which levels to simulate — you need not run them all.](figures/levels.pdf){width=100%}

The usual way is **hierarchical**: simulate one level, then the next *conditioned on* it — a
species tree, then gene families along its branches, then sequences along the resulting gene
trees. Each level treats the one above as a fixed backbone. This keeps the levels independent
and easy to reason about, and lets a single species tree seed many genome, trait, or sequence
runs, so one backbone underlies a whole benchmark. It is the design ZOMBI has always had:
species, then genomes, then sequences, in that order.

Sometimes the levels are *not* independent — a trait may change the rate at which genes are
gained or lost; gene content may decide which lineages survive and diversify. When influence
runs *between* levels, a one-directional pipeline is no longer enough and the levels must be
simulated **jointly**. ZOMBI2's `coevolve` mode does this, coupling species, genomes, and
traits along directed links you choose. Sequences remain a downstream layer for now, though one
can imagine them feeding back on the tree too.

## For the impatient

A minimal run can be just two commands, for example: build a species tree, then evolve gene
families along it.

```bash
# 1. a species tree: 20 extant species from a birth–death process
zombi2 species --tips 20 --age 5 --birth 1.0 --death 0.3 --seed 1 -o my_tree

# 2. gene families evolving along it (duplication, transfer, loss, origination)
zombi2 genomes -t my_tree/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o my_genomes
```

The first writes `my_tree/species_tree.nwk`; the second writes gene trees and the presence/absence
**profile matrix** (gene families × species) into `my_genomes/`. You can do the same run from
Python:

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_families=20, seed=42)

genomes.profiles.matrix     # gene families × extant species (copy numbers)
genomes.write("out/")       # trees, event tables, transfers, profiles
```
