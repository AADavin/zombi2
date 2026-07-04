# Introduction

**ZOMBI2** is a simulator of genome and species-tree evolution for phylogenetics — a
ground-up redesign of ZOMBI [@davin2020zombi], rebuilt around a fast native engine, a simple
command-line interface, and a composable Python library. Its purpose is to generate
**benchmark datasets** — gene trees, reconciliations, and copy-number profiles — under
explicit, reproducible evolutionary models, so that phylogenetic methods can be tested against
data whose true history is known.

## The two-step workflow

ZOMBI2 simulates evolution in two symmetric steps:

1. build a **species tree** — backward in time as a reconstructed birth–death process
   conditioned on the number of extant species, or forward in time as a complete tree that
   keeps extinct lineages;
2. evolve **gene families** forward along that tree — duplication, transfer, loss and
   origination (DTL), optionally at nucleotide resolution or with gene order.

::: note
The two halves are independent. Any Newick tree can feed the gene-family step, and one species
tree can be reused across many genome simulations — so the same backbone underlies every
dataset in a benchmark.
:::

From those two steps you obtain the reconstructed **complete and pruned gene trees**, their
reconciliations against the species tree, and the presence/absence **profile matrix** — the raw
material for phylogenetic-profiling, reconciliation, and molecular-dating analyses.

![One gene family evolving along a species tree: a duplication ($\square$), a loss ($\circ$) and a transfer (arrow) placed on the branches where the continuous-time (Gillespie) process fires them.](figures/species_tree_events.pdf)

A complete simulation is a handful of lines:

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, seed=42)

genomes.profiles.matrix     # gene families × extant species (copy numbers)
genomes.write("out/")       # trees, event tables, transfers, profiles
```

## Why it matters

Because ZOMBI2 records the *true* history behind every dataset, it is a testing ground for
methods that must infer that history from incomplete data — reconciliation, gene-tree/species-tree
inference, and divergence-time estimation. It also generates data for questions where genome
content and organismal traits are coupled: for example, dating the bacterial tree of life from
the gene-family signature of the Great Oxidation Event [@davin2025goe], where the acquisition of
oxygen-using gene families marks a datable transition. Models of that kind — traits that shape
gene content, and gene content that shapes diversification — are the subject of the coevolution
chapters.

## How this manual is organised

The manual is a **concepts-and-tutorial** companion to the software. Part I gets you installed
and running end to end. Parts II–VI then cover the models in the order you meet them in a
simulation: species trees, gene families, traits, coevolution, and sequence evolution — each
with worked examples and figures. The exhaustive command-line and Python API reference lives in
the online documentation rather than here, so this book stays readable front to back.
