# A tour of ZOMBI2

This chapter is the whole tool in one pass. The most important concepts are explained here, but you might not need all of them: to simulate something simple, like a species tree, you can go straight to the dedicated chapter.

## The four levels of ZOMBI2

ZOMBI2 simulates four levels of evolution.

- **Species**, the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time.
- **Genomes**, the genes that exist in each lineage. Genomes can be simulated at three **resolutions**: gene families alone, genes placed on chromosomes, or a full nucleotide genome. A genome is always simulated within a species tree.
- **Sequences**, the nucleotides or amino acids inside each gene. Sequences evolve on gene trees (which are generated during the genome simulation), so genomes must be simulated first.
- **Traits**, phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum.

![The four levels of ZOMBI2. Everything starts from the species tree, which is the general backbone of a simulation run. Then, in that species tree, you can simulate genomes, or traits. If you simulate genomes, you can also simulate their sequences.](figures/fig-2-1-four-levels_print.png){width=45%}

A run in which every level is simulated:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species}) \cdot P(\text{Sequences} \mid \text{Genomes}) \cdot P(\text{Traits} \mid \text{Species})$$

You need not run them all. Skip sequences if you only want gene trees, which the genome level already produces; skip genomes if you want a species tree with traits on it. Everything depends on a species tree, so a workflow almost always begins by simulating one alone. The exception is a **joint** model, below, where the tree is an output rather than an input.

## Time

ZOMBI2 is a forward simulator: evolution runs from an ancestral state at time 0 to the present.

Time is imposed by the species tree, and every rate is measured against that scale. If your tree runs from 0 at the origin to 1 at the tips, the simulation lasts one unit of time. Time 0 is the origin of the founding lineage, and every time you give ZOMBI2 is measured on it: the moment of a mass extinction, the breakpoints of a rate that changes through time.

The founding lineage lives for a while before it first splits. That stretch is the tree's **stem**, and the first split is its **crown**.

![The stem. A run starts from one lineage at time 0, the origin, which lives for a while before it first splits. Everything up to that split is the stem, and it is ordinary simulated time: genes are gained and lost along it, traits drift along it. Every tree ZOMBI2 writes therefore gives its root a branch length.](figures/stem.pdf){width=70%}

Some tree viewers do not draw the stem by default, so be careful. We recommend plotting trees with Phylustrator, a companion package to ZOMBI2 maintained by the same author.

## Rates

In ZOMBI2 different events happen over time with different probability: speciations and extinctions at the species level; duplications, transfers, losses, originations and rearrangements at the genome level; substitutions in a sequence; changes in a trait.

How often an event occurs is controlled by its **effective rate**:

$$\text{effective rate} = \text{scope}(\text{base}) \times \text{modifiers}$$

The **base** is the speed of a single event, in units of inverse time. The **scope** says how we measure the event: per lineage, per copy, per site. The **modifiers** are the dimensionless multipliers a rate picks up from its context. Appendix A is the full reference: the units, each level's default scope, every modifier and which levels accept it, and the Gillespie algorithm that turns rates into events.

## Conditioning

ZOMBI2 includes options for dependencies between different parts of the simulation, so that complex scenarios can be simulated. In some cases we want to run a model that depends on a different model — genome evolution that depends on some trait, for instance. If both models can strictly be run in sequence, we talk about **conditioning**.

This is best explained with an example. Imagine we are simulating the evolution of mammals and their olfactory genes. A habitat trait switches between aquatic and terrestrial along the tree, and aquatic lineages lose those genes four times faster.

![Conditioning: a habitat trait is grown first and held fixed, and a gene loss rate reads it. The driver, the connection and the target are the three parts every conditioned run has, and Chapter 9 takes them one at a time.](figures/conditioning_print.png){width=95%}

We can write this run like this:

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species}) \cdot P(\text{Genomes} \mid \text{Species}, \text{Traits})$$

When we condition there are three things to pay attention to:

- the **driver**: the variable controlling the run.
- the **target**: the variable that the driver modifies.
- the **connection**: how the driver controls the target.

Not everything can be connected in ZOMBI2, but it is flexible enough to allow very specific rules — a trait can set how often a lineage loses genes, and a gene family's presence can set how fast a trait changes. There is a full chapter devoted to conditioning, and a user should read it to study the different cases.

## Joining

**Joining** simulates two levels **at the same time**, for the scenarios where each one shapes the other:

- A trait controls how fast a lineage speciates — body size, or a habitat.
- Gene content decides survival: lineages that acquire a key gene diversify faster.

Take the first. Lineages in the fast state split more often, so the trait decides the shape of the tree; and the trait evolves along that very tree, so the tree decides where the trait can go. Neither can be grown first and handed over, so one run grows both, and the tree comes out as a result rather than going in as an input.

Because neither is finished first, neither can be written out and handed over. There is no conditional probability to write, only a joint one:

$$P(\text{Species}, \text{Traits})$$

The test is one question: can the driver be grown first, on its own, and handed over? If it can, condition. If it cannot, join. Chapter 9 works through conditioning and Chapter 10 through joining.
