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

Some questions ZOMBI2 is built to answer, each of which is one run and a comparison:

- **How well does a reconciliation method recover the truth?** Evolve gene families under duplication, transfer and loss, and you get every family's true gene tree together with the event behind every node. Reconcile the gene trees against the species tree and score what the method found against what actually happened.

- **Can a transfer be detected when the donor is gone?** Transfers in ZOMBI2 can come from lineages that later go extinct, so a gene arrives in a survivor from a donor that leaves no other trace. The event log names that donor, so you can ask how often a detection method finds a transfer whose source is no longer on the tree.

- **What does genome reduction look like in host-restricted bacteria?** Evolve a lifestyle trait, free-living or host-restricted, then let it drive the loss rate. Lineages that move inside a host shed genes faster, and you can measure how much of the resulting genome-size pattern a method attributes to lifestyle rather than to shared ancestry.

- **Does a dating method survive a clock that is not strict?** Give the substitution rate a relaxed clock, so lineages evolve at different paces, and compare the dates a method infers from the alignment against the true node ages.

- **How accurate is ancestral sequence reconstruction?** A run records the sequence at every internal node, not just the tips, so a reconstruction can be compared residue by residue against the sequence that really sat there.

- **Is a trait correlation real, or an artefact of the tree?** Two traits evolving on the same tree look correlated at the tips whether or not either drives the other, because they share ancestry. Simulate with the correlation switched off, on the same tree, and you have the baseline any comparative method has to beat.

- **Does a trait actually drive diversification?** Let a trait set the speciation rate, so the tree and the trait grow together, and test whether a state-dependent method recovers the effect — or reports one when the trait is inert.

- **What signal survives in gene order?** At the ordered resolution, inversions, transpositions and translocations rearrange genes on chromosomes, and fissions and fusions change the karyotype itself, so synteny and rearrangement methods can be tested against the moves that were actually made.

## Installing it

ZOMBI2 needs Python 3.10 or newer and depends only on NumPy, so there is no compiler and no toolchain to set up.

```bash
pip install zombi2
```

`zombi2 --version` confirms the install, and `zombi2 -h` lists the commands, one per level.

## For the impatient

Two commands: build a species tree, then evolve gene families along it.

```bash
# 1. a species tree: 20 surviving species from a birth-death process
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1

# 2. gene families along it, by duplication, transfer, loss and origination
zombi2 genomes out/ \
    --duplication 0.2 --transfer 0.1 --loss 0.25 --origination 0.5 --seed 42
```

`out/` now holds one directory per level:

```
out/species/    species_complete.nwk   the tree, extinct lineages kept
                species_extant.nwk     the survivors only
                species_events.tsv     every speciation and extinction, with its time
out/genomes/    genome_events.tsv      every duplication, transfer, loss and origination
                profiles.tsv           gene-family copy numbers per species
out/logs/       species.log            what was run, and with which parameters
                genomes.log
```

Each level keeps to its own directory, and the outputs that run to one file per gene family get a directory of their own inside it, so a run of a few hundred families stays legible. Pass `--flat` to any command and it writes everything straight into `out/` instead.

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

Each CLI flag is the Python argument of the same name, hyphenated, and `.write()` does what the run directory does. In a session you can also read a result straight off, without writing anything: `gen.profiles` is the copy-number matrix, and `gen.gene_trees` gives each gene family's true gene tree.

The next chapter, *A tour of ZOMBI2*, lays out the four levels, the three ways they can relate, how every rate is written, and the vocabulary the rest of the book uses.
