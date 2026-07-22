# A tour of ZOMBI2

ZOMBI2 simulates evolution at four levels: Species Tree, Genomes, Sequences and Traits. This chapter introduces the four levels, the three ways they can relate, and the single shape every rate takes. It is the vocabulary the rest of the book uses.

## The four levels of ZOMBI2

ZOMBI2 simulates evolution at four different levels:

- **Species** — the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time. Every ZOMBI2 workflow has some species tree, so this is the first thing to run.
- **Genomes** — the genes that exist in each lineage of the tree. Genomes can be simulated at different levels of resolution: a full genome represented with nucleotides, or just a set of gene families. Genomes are always simulated within a species tree; you always need one to obtain them.
- **Sequences** — the nucleotides or amino acids inside each gene. Sequences *always* evolve on gene trees, so you need to simulate genomes before using this level.
- **Traits** — phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum. Traits can influence the evolution of the other levels in different ways that we cover in this chapter.

The most general way to connect the four levels is shown in Figure 1.

![The four levels of ZOMBI2. Species, genomes and sequences form a chain of ancestry: a genome lives on the species tree, a sequence inside a gene. Traits branch to the side, because a trait can ride any species tree.](figures/fig-2-1-four-levels_print.png){width=40%}

You can run a simulation in which every level is simulated. We write it in this notation:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species}) \cdot P(\text{Sequences} \mid \text{Genomes}) \cdot P(\text{Traits} \mid \text{Species})$$

But you can always choose which levels to run; you need not run them all. For example, you do not need to simulate sequences if you are only interested in gene trees, which the genome level already produces:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species})$$

And you do not have to simulate genomes if all you want is a species tree with some traits on it:

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species})$$

In ZOMBI2 everything depends on a species tree, and in most cases you begin a workflow by simulating the tree alone. There are a few exceptions, which we cover a bit later.

## Time

ZOMBI2 is a forward simulator: evolution runs from an ancestral state at time 0 to the present.

Time is imposed by the species tree, and every rate is measured against that time scale. If your tree runs from 0 at the root to 1 at the tips, your simulation lasts one unit of time. Time 0 is the origin of the founding lineage, and every time you give ZOMBI2 — the moment of a mass extinction, the breakpoints of a rate that changes through time — is measured on that scale.

The founding lineage lives for a while before it first splits, so it has a duration of its own.

## Rates

In ZOMBI2 everything is driven by events that fire over time. The kind of event depends on the level being simulated. Some of the basic events are:

- **Species** — speciations and extinctions
- **Genomes** — duplications, transfers, losses, originations, inversions, transpositions.
- **Sequences** — mutations
- **Traits** — phenotypic changes

The frequency at which an event fires depends on its **effective rate**:

$$\text{effective rate} = \text{scope}(\text{base}) \times \text{modifiers}$$

The **base** is the speed of a single event (how fast), in units of inverse time. The **scope** wraps that base to say how many independent chances the event has: per lineage, per copy, or per site. The **modifiers** are dimensionless context multipliers that make a rate faster or slower depending on some factor — the lineage where the event happens, the gene family affected, or the total diversity present in the simulation.

Most of the time you do not need to touch either the default scope or the modifiers — but you can, and it is this flexibility that lets ZOMBI2 reach a wide range of scenarios. With them you could simulate, for example:

- a burst of change concentrated early and then tapering off — an early radiation, or a trait that diversifies fast at first and settles — by giving the rate a schedule that starts high and drops later;
- a molecular clock that speeds up and slows down along the tree, with closely related lineages ticking at similar rates, by letting each lineage inherit its rate from its parent and drift a little at every split;
- a radiation that starts fast and then eases off as the clade fills up toward a carrying capacity, by having the speciation rate read the total diversity present at each moment.

The full rate reference — how these units work in detail, the default scope at each level, and the complete catalogue of modifiers — is **Appendix A**, which also covers the Gillespie algorithm that turns these rates into the events of a simulation.

## Going beyond the basic simulation: conditioning and joining levels

Some evolutionary scenarios cannot be simulated by the paradigm we have described so far. For example:

- A trait evolves along a tree and itself controls how fast that tree speciates.
- Gene content decides survival: lineages that acquire a key gene diversify faster than the rest.
- Two levels feed back on each other: a trait raises a gene family's loss rate, while carrying that family in turn pulls on the trait, so the two co-evolve.

ZOMBI2 adds two new connections between levels to reach scenarios like these. One is **conditioning**, the other is **joining**. Because they need some care, each has its own chapter; here we give only a brief overview.

When we **condition** one level on another, a parameter of the second level stops being a fixed number you set and instead reads the state of the first. In the run below, the gene-loss rate is no longer constant: it depends on the trait's value on each branch. For example, you can use this to simulate that aquatic lineages lose their olfactory genes faster. Because the trait can be simulated first and then held fixed, conditioning is still two ordinary runs in order: you simulate the driver level, write it out, and feed it to the second, exactly as you already feed a species tree to a genome run.

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species}) \cdot P(\text{Genomes} \mid \text{Species}, \text{Traits})$$

When we **join** two levels, neither can be simulated first, because each depends on the other as it unfolds. If a trait speeds up speciation, then faster-speciating lineages leave more descendants, so the tree's shape depends on the trait, but the trait is evolving along that very tree at the same time. No order works, so the two levels are grown together in a single run, step by step. When the coupling reaches back into the species tree, the tree itself becomes an output rather than an input.

$$P(\text{Species}, \text{Traits})$$

We **join** whenever a coupling would form a loop: when one level shapes another and is shaped back, directly or through the tree. If the influence runs only one way, we condition; if it runs in a loop, we join. In the directional case we can still name the variable the driver sets on its target: for a trait driving speciation, that variable is the speciation rate.

## Using ZOMBI2 in Python

Each level is a function in its own subpackage, and they compose by feeding one level's result into the next. A run returns a *result object*; you read it directly in the session or write it to disk with `.write()`. A whole workflow is a short script:

```python
from zombi2 import species, genomes, sequences, traits
from zombi2.sequences import substitution_models as sm

sp   = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
gen  = genomes.simulate_genomes_unordered(sp, duplication=0.2, loss=0.25, origination=0.5, seed=42)
seqs = sequences.simulate_sequences(gen, model=sm.hky85(kappa=2.0), length=300, seed=7)
bm   = traits.simulate_continuous(sp, rate=1.0, seed=1)
```

Each call takes the object it depends on — genomes and traits read the species result, sequences reads the genomes result — so the script reads top to bottom in exactly the `P(·)` order from the start of this chapter. Every level's function, its arguments, and its result object are covered in that level's own chapter.

## Using ZOMBI2 from the CLI

The same simulations run from the command line. Each level is a subcommand of `zombi2`, and its flags are the long-form names of the Python arguments:

```bash
# a dated species tree (20 extant tips)
zombi2 species --birth 1 --death 0.3 --n-extant 20 --seed 1 -o out/

# gene families along it
zombi2 genomes -t out/species_complete.nwk --duplication 0.2 --loss 0.25 --origination 0.5 --seed 42 -o out/
```

A rate flag takes a rate **written exactly as you would write it in Python** — a bare number, or a scope wrapper and modifiers composed with `*`, quoted so the shell keeps it in one piece:

```bash
# speciation drops to a third of its rate at time 3 (a skyline)
zombi2 species --birth "1.0 * OnTime({0: 1.0, 3: 0.3})" --death 0.3 --total-time 5 --seed 1 -o out/
```

`-o` sets the output directory and `-t` feeds one level's tree into the next, so a pipeline is a sequence of commands sharing a directory; a `--params` TOML file can hold the settings for a whole pipeline at once. On the clean core the CLI covers all four levels — **species**, **genomes**, **sequences**, **traits**; the coupled models are run from Python until their commands land.

## Output in ZOMBI2

Every run can be written to disk with `result.write("out/", outputs=[...])`; with no `outputs` argument it writes that level's **default** set of files. The formats are uniform across levels — **trees** in Newick, **tables and event logs** in TSV, **sequences** in FASTA — and branch lengths are in units of time everywhere except the sequence phylograms, which are in substitutions per site. At every level the **event log** (`*_events.tsv`) is the true, ordered history the run actually followed: the source of truth from which the summaries are derived. Appendix B lists every file, level by level.
