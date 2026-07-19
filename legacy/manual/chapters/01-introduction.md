```{=latex}
\part{Getting started}
```

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

ZOMBI2 simulates four **levels** — species, genomes, traits and sequences — and composes them
either as a hierarchical pipeline or as a jointly-coupled `coevolve` run. The next chapter,
[A tour of ZOMBI2](#a-tour-of-zombi2), lays out those levels, the two ways to compose them, how
every rate works, and the vocabulary the rest of the book uses. What follows is the shortest
possible taste.

## For the impatient

A minimal run can be just two commands, for example: build a species tree, then evolve gene
families along it.

```bash
# 1. a species tree: 20 extant species from a birth–death process
zombi2 species --tips 20 --age 5 --birth 1.0 --death 0.3 --seed 1 -o my_tree

# 2. gene families evolving along it (duplication, transfer, loss, origination)
zombi2 genomes -t my_tree/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --initial-families 20 --seed 42 -o my_genomes
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
