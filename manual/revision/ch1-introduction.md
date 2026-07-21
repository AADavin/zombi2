# Introduction

## Why simulate

Evolutionary biology infers the past from what survives into the present: a gene tree from an alignment, a rate of gene loss from a set of genomes, an ancestral body size from the sizes of living species. The true history is gone, so there is nothing to check the answer against.

Simulation is the way around this. You choose a model, run it forward in time, and get a dataset whose history you already know: which lineages went extinct, which gene was transferred and when, what the sequence at each internal node was. Run a method on that dataset and you can measure how much of the history it recovers. This is how phylogenetic methods are tested, calibrated and compared.

A simulator is only useful for this if it says what it did. ZOMBI2 writes the full event log of every run, in the order the events fired, and records the state at every internal node of the tree, not just at the tips. The ground truth is part of the output.

## What ZOMBI2 is

ZOMBI2 simulates four levels of evolution: **Species**, **Genomes**, **Sequences** and **Traits**. One program covers all four. You can run a single level, run several in sequence, or let one level drive another.

It comes as a Python library and as a command-line tool. Both drive the same engine and take the same parameters, so a run can be written either way.

The engine is fast. A birth–death tree with 100,000 surviving species takes about a second on a laptop, so the size of a simulation is rarely a constraint.

## What it can do

- **Species** — dated species trees from a birth–death process. The rates can be constant, change through time, slow down as the clade fills up, or drift from parent to daughter lineage. Mass-extinction pulses, incomplete sampling and fossil recovery are all available. You get both the complete tree, with the extinct lineages still on it, and the tree of the survivors.

- **Genomes** — gene families evolving along a species tree by duplication, transfer, loss and origination. The event log is a full genealogy, so every family's true gene tree comes out of it, along with the copy-number profile of each species. At the ordered resolution the genes sit at positions on chromosomes, and inversions, transpositions, translocations, fissions and fusions rearrange them.

- **Sequences** — a nucleotide alignment for each gene family, evolved down that family's own gene tree under JC69, K80, HKY85 or GTR. The clock can be strict or relaxed. You get the alignment, the phylogram in substitutions per site, and the sequence at every internal node.

- **Traits** — a phenotype evolving along a tree. A continuous trait diffuses, either freely or pulled toward an optimum; a discrete trait switches between the states you name. The value is recorded at every node, and a discrete trait also records every change along every branch.

- **Coupling** — a rate at one level can read its value from another level instead of being a number you type. When the driver can be simulated first, it is written to a file and handed to the level it drives. When the two are entangled, as when a trait drives speciation and so shapes the tree it is evolving on, both are grown together in a single run from the Python library.

## For the impatient

Two commands: build a species tree, then evolve gene families along it.

```bash
# 1. a species tree: 20 surviving species from a birth-death process
zombi2 species --birth 1 --death 0.3 --n-extant 20 --seed 1 -o out/

# 2. gene families along it, by duplication, transfer, loss and origination
zombi2 genomes -t out/species_complete.nwk \
    --duplication 0.2 --transfer 0.1 --loss 0.25 --origination 0.5 --seed 42 -o out/
```

`out/` now holds the two trees (`species_complete.nwk`, which keeps the extinct lineages, and `species_extant.nwk`, which does not), the event log of each level (`species_events.tsv` and `genome_events.tsv`) and the matrix of gene-family copy numbers per species (`profiles.tsv`).

Here is the same run from Python:

```python
from zombi2 import species, genomes

sp  = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
gen = genomes.simulate_genomes_unordered(sp.complete_tree,
                                         duplication=0.2, transfer=0.1,
                                         loss=0.25, origination=0.5, seed=42)

sp.write("out/")
gen.write("out/")
```

Each CLI flag is the Python argument of the same name, hyphenated, and `.write()` does what `-o` does. In a session you can also read a result straight off, without writing anything: `gen.profiles` is the copy-number matrix, and `gen.gene_trees` gives each gene family's true gene tree.

The next chapter, *A tour of ZOMBI2*, lays out the four levels, the three ways they can relate, how every rate is written, and the vocabulary the rest of the book uses.
