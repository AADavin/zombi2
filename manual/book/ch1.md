# Introduction

## Why simulate

Evolutionary biology infers the past from what survives into the present: a gene tree from an alignment, a rate of gene loss from a set of genomes, an ancestral body size from the sizes of living species. The true history is gone, so there is nothing to check the answer against.

Simulation is the way around this. You create a dataset whose history you already know: which lineages went extinct, which gene was transferred and when, what the sequence at each internal node was. Run a method on that dataset and you can measure how much of the history it recovers. This is how phylogenetic methods are tested, calibrated and compared.

## What ZOMBI2 is

ZOMBI2 simulates four levels of evolution: **Species**, **Genomes**, **Sequences** and **Traits**. You can run one level, several in sequence, or let one drive another. It is a Python library and a command-line tool over the same engine, taking the same parameters, so a run can be written either way.

The next sections are the whole tool in one pass — the concepts every level shares. You might not need all of them: to simulate something simple, like a species tree, you can go straight to its chapter.

## The four levels of ZOMBI2

- **Species**, the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time.
- **Genomes**, the genes that exist in each lineage. Genomes can be simulated at three **resolutions** — **family**, the gene families alone; **ordered**, the genes placed on chromosomes; **nucleotide**, a genome as DNA coordinates (Chapters 3 to 5). A genome is always simulated within a species tree — even in a joint run, below, where the two grow together.
- **Sequences**, the nucleotides or amino acids inside each gene: the letters themselves, which even a nucleotide-resolution genome does not carry. Sequences evolve on gene trees (which are generated during the genome simulation), so genomes must be simulated first.
- **Traits**, phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum.

![The four levels of ZOMBI2. Everything starts from the species tree, which is the general backbone of a simulation run. Then, in that species tree, you can simulate genomes, or traits. If you simulate genomes, you can also simulate their sequences.](figures/fig-2-1-four-levels_print.png){width=45%}

A run in which every level is simulated:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species}) \cdot P(\text{Sequences} \mid \text{Genomes}) \cdot P(\text{Traits} \mid \text{Species})$$

You need not run them all. Skip sequences if you only want gene trees, which the genome level already produces; skip genomes if you want a species tree with traits on it. Everything depends on a species tree, so a workflow almost always begins by simulating one alone. The exception is a **joint** model, below, where the tree is an output rather than an input.

## Time

ZOMBI2 is a forward simulator: evolution runs from an ancestral state at time 0 to the present.

Time is imposed by the species tree, and every rate is measured against that scale. How much time a run covers is set by its stop condition: `total_time` fixes the height directly, and `n_extant` stops at a tip count, so the height comes out of the run. Either way, if the tree runs from 0 at the origin to 1 at the tips, the simulation lasted one unit of time. Time 0 is the origin of the founding lineage, and every time you give ZOMBI2 is measured on it: the moment of a mass extinction, the breakpoints of a rate that changes through time.

The founding lineage lives for a while before it first splits. That stretch is the tree's **stem**, and the first split is its **crown**. The stem is ordinary simulated time, so every tree ZOMBI2 writes gives its root a branch length.

![The stem. A run starts from one lineage at time 0, the origin, which lives for a while before it first splits. Everything up to that split is the stem, and it is ordinary simulated time: genes are gained and lost along it, traits drift along it. Every tree ZOMBI2 writes therefore gives its root a branch length.](figures/stem.pdf){width=70%}

Some tree viewers do not draw the stem by default, so a plotted tree can look shorter than the height the run reports. We recommend plotting trees with Phylustrator (`pip install phylustrator`), a companion package to ZOMBI2 maintained by the same author.

## Rates

In ZOMBI2 different events happen over time with different probability: speciations and extinctions at the species level; duplications, transfers, losses, originations — the birth of a new gene family — and rearrangements at the genome level; substitutions in a sequence; changes in a trait.

How often an event occurs is controlled by its **effective rate**:

$$\text{effective rate} = \text{scope}(\text{base}) \times \text{modifiers}$$

The **base** is the expected number of events per unit time, for one unit of whatever the scope counts. The **scope** names that unit — a lineage, a gene copy, a site — and so decides how many independent chances run at once: a loss rate of 0.25 counted per copy, in a genome of forty copies, is a total loss rate of 40 × 0.25 = 10 per unit time; counted per lineage, it stays 0.25 however large the genome grows. The **modifiers** are the dimensionless multipliers a rate picks up from its context; they are written as **verbs** chained onto the rate — `.changing_at(…)`, `.scaled_by(…)` — and that is what this book calls them. Appendix A is the full reference: the units, each level's default scope, every modifier and which levels accept it, and what the engine does with a rate once it has one.

## Conditioning

ZOMBI2 includes options for dependencies between different parts of the simulation, so that complex scenarios can be simulated. In some cases we want to run a model that depends on a different model — genome evolution that depends on some trait, for instance. If both models can strictly be run in sequence, we talk about **conditioning**.

This is best explained with an example. Imagine we are simulating the evolution of mammals and their olfactory genes. A habitat trait switches between aquatic and terrestrial along the tree, and aquatic lineages lose those genes four times faster.

We can write this run like this:

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species}) \cdot P(\text{Genomes} \mid \text{Species}, \text{Traits})$$

When we condition there are three things to pay attention to:

- the **driver**: the variable controlling the run.
- the **target**: the variable that the driver modifies.
- the **connection**: how the driver controls the target.

![Conditioning: a habitat trait is grown first and held fixed, and a gene loss rate reads it. The driver, the connection and the target are the three parts every conditioned run has, and Chapter 8 takes them one at a time.](figures/conditioning_print.png){width=95%}

Not everything can be connected in ZOMBI2, but it is flexible enough to allow very specific rules — a trait can set how often a lineage loses genes, and a gene family's presence can set how fast a trait changes. There is a full chapter devoted to conditioning, and a user should read it to study the different cases.

## Joining

**Joining** simulates two levels **at the same time**, for the scenarios where each one shapes the other:

- A trait controls how fast a lineage speciates — body size, or a habitat.
- Gene content decides diversification: lineages that acquire a key gene split more often.

Take the first. Lineages in the fast state split more often, so the trait decides the shape of the tree; and the trait evolves along that very tree, so the tree decides where the trait can go. Neither can be grown first and handed over, so one run grows both, and the tree comes out as a result rather than going in as an input.

Because neither is finished first, neither can be written out and handed over. There is no conditional probability to write, only a joint one:

$$P(\text{Species}, \text{Traits})$$

The test is one question: can the driver be grown first, on its own, and handed over? If it can, condition. If it cannot, join. Chapter 8 works through conditioning and Chapter 9 through joining.

## What it can do

Some questions ZOMBI2 is built to answer:

- **How well does a reconciliation method recover the truth?** Evolve gene families under duplication, transfer and loss, and you get every family's true gene tree with the event behind every node. Reconcile against the species tree and score what the method found against what happened.

- **Can a transfer be detected when the donor is gone?** Transfers can come from lineages that later go extinct, so a gene arrives in a survivor from a donor that leaves no other trace. The event log names that donor, so you can ask how often a detection method finds a transfer whose donor is no longer on the tree.

- **What does genome reduction look like in host-restricted bacteria?** Evolve a lifestyle trait, free-living or host-restricted, and let it drive the loss rate. Lineages that move inside a host shed genes faster, so you can measure how much of the genome-size pattern a method attributes to lifestyle rather than to shared ancestry.

- **Can Bayesian methods of inference (relaxed molecular clocks) recover the true ages of a tree?** Give the substitution rate a relaxed clock, so lineages evolve at different paces, and compare the dates a method infers from the alignment against the true node ages.

- **How accurate is ancestral sequence reconstruction?** A run records the sequence at every internal node, not just the tips, so a reconstruction can be compared residue by residue against the one that really sat there.

- **Is a trait correlation real, or an artefact of the tree?** Two traits evolving on the same tree look correlated at the tips whether or not either drives the other, because they share ancestry. Simulate with the correlation switched off, on the same tree, and you have the baseline any comparative method has to beat.

- **Does a trait actually drive diversification?** Let a trait set the speciation rate, so the tree and the trait grow together, and test whether a state-dependent method recovers the effect.

- **What signal survives in gene order?** With genes placed on chromosomes — the **ordered** genome resolution (Chapter 4) — inversions, transpositions and translocations rearrange them, and fissions and fusions change the karyotype itself, so synteny and rearrangement methods can be tested against the moves that were actually made.

## Installing it

ZOMBI2 needs Python 3.10 or newer and depends on NumPy and tqdm, plus the `tomli` backport on Python 3.10:

```bash
pip install zombi2
```

`zombi2 --version` confirms the install, and `zombi2 -h` lists the commands: one per level (Chapters 2 to 7 cover them in order), plus `joint`, which grows the species tree and a level that drives it in one run (Chapter 9), and `tools`, analyses that read a finished run (Appendix C).

ZOMBI2 is pure Python over NumPy, with no compiled part to build, and the test suite runs on **Linux,
macOS and Windows** on every change, so the same command does the same thing on all three. The one
thing that differs on Windows — how backslash paths are read inside a rate expression or a
`--params` file — is in Appendix A, under "Paths on Windows".

## The examples gallery

Every figure in this book is drawn from a run, and the code behind a great many of them lives in the
[examples gallery](https://aadavin.github.io/zombi2/gallery.html) — a page of worked examples, each
with the exact commands that made it. The gallery numbers them by section: `Sp` species, `Ge` genomes,
`Sq` sequences and `Tr` traits for the four levels, then `Co` conditioning and `Jo` joining for the
two ways of coupling them (Chapters 8 and 9). This book cites them the way it cites a
figure, so a paragraph that describes what an inversion does to a chromosome ends **(Ge7)**.

## For the impatient

Two commands: build a species tree, then evolve gene families along it.

```bash
# 1. a species tree — by default 20 extant tips at birth 1.0, every tip sampled
zombi2 species out/

# 2. gene families along it, by duplication, transfer, loss and origination
#    (defaults 0.2, 0.1, 0.25 and 0.5)
zombi2 genomes out/
```

`out/` now holds a report of the whole run, and one directory per level:

```
out/run.zombi2                         every file the run wrote, and how to reproduce it
out/species/    species_complete.nwk   the tree, extinct lineages kept
                species_extant.nwk     the survivors only
                species_events.tsv     every speciation and extinction, one row each
                species_fates.tsv      which tips survived, died out, or went unsampled
                species_summary.json   what the run produced, in numbers
                species.log            what was run, and with which parameters
out/genomes/    genome_events.tsv      every duplication, transfer, loss and origination
                genomes.tsv            the genes each species ends up with, ancestors included
                profiles.tsv           gene-family copy numbers, extant species only
                initial_genome.tsv     the genome the run started from
                gene_trees/            one true gene tree per family
                genome_summary.json    what the run produced, in numbers
                genomes.log
```

Open `run.zombi2` first — plain text, one page: it names every file the run wrote and carries the commands that rebuild it. Below it each level keeps to its own directory, run log included, and an output with one file per gene family gets a directory of its own inside the level's, so a few hundred families stay legible. `--flat` writes every file straight into `out/` instead, and writes no run report; `--quiet` turns off the progress bar.

Here is the same model from Python, with seeds pinned so it is repeatable. The rates are spelled out because they are exactly the command line's defaults, and only the command line has defaults: in Python every genome rate starts at 0, so a call that names none of them still runs, and gives a history with no duplication, transfer or loss — just the families it started with, copied down the tree.

```python
from zombi2 import species, genomes

sp  = species.simulate_species_tree(birth=1.0, n_extant=20, seed=1)
gen = genomes.simulate_genomes_family(sp, duplication=0.2, transfer=0.1,
                                      loss=0.25, origination=0.5, seed=42)

sp.write("out/species/")
gen.write("out/genomes/")
```
