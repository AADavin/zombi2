# Introduction

## Why simulate

Evolutionary biology infers the past from what survives into the present: a gene tree from an alignment, a rate of gene loss from a set of genomes, an ancestral body size from the sizes of living species. The true history is gone, so there is no way to know whether these inferences are correct.

Simulation is the way around this. You create a dataset whose history you already know: which lineages went extinct, which gene was transferred and when, what the sequence at each internal node was. Run a method on that dataset and you can measure how much of the history it recovers. This is how phylogenetic methods are tested, calibrated and compared.

## What ZOMBI2 is

ZOMBI2 simulates four levels of evolution: **Species**, **Genomes**, **Sequences** and **Traits**. It is a Python library and a command-line tool over the same engine, taking the same parameters, so a run can be written either way. The next sections are the whole tool in one pass, the concepts every level shares. You might not need all of them: to simulate something simple, like a species tree, you can go straight to its chapter.

## The four levels of ZOMBI2

- **Species**, the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time.
- **Genomes**, the genes that exist in each lineage. Genomes can be simulated at three **resolutions**: **family**, the gene families alone; **ordered**, the genes placed on chromosomes; **nucleotide**, a genome as DNA coordinates (Chapters 3 to 5).
- **Sequences**, the nucleotides or amino acids of each gene. Sequences evolve on gene trees (which are generated during the genome simulation), so genomes must be simulated first.
- **Traits**, phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum. They can be continuous or discrete.

![The four levels of ZOMBI2. Everything starts from the species tree, which is the general backbone of a simulation run. Then, in that species tree, you can simulate genomes, or traits. If you simulate genomes, you can also simulate their sequences.](figures/fig-2-1-four-levels_print.png){width=45%}

A run in which every level is simulated:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species}) \cdot P(\text{Sequences} \mid \text{Genomes}) \cdot P(\text{Traits} \mid \text{Species})$$

You can run only the levels you need. Skip sequences if you only want gene trees, which the genome level already produces; skip genomes if you want a species tree with traits on it. Everything depends on a species tree, so a workflow almost always begins by simulating one alone. The exception is a **joint** model, in which traits and species, or genomes and species, are simulated simultaneously.

## Time

ZOMBI2 is a forward simulator: evolution runs from an ancestral state at time 0 to the present.

Time is imposed by the species tree, and every rate is measured against that scale. How much time a run covers is set by its stop condition: `total_time` fixes the height directly, and `n_extant` stops at a tip count, so the height comes out of the run. Either way, if the tree runs from 0 at the origin to 1 at the tips, the simulation lasted one unit of time. Time 0 is the origin of the founding lineage, and every time you give ZOMBI2 is measured on it: the moment of a mass extinction, the breakpoints of a rate that changes through time.

The founding lineage lives for a while before it first splits. That stretch is the tree's **stem**, and the first split is its **crown**. The stem is ordinary simulated time, so every tree ZOMBI2 writes gives its root a branch length.

![The stem. A run starts from one lineage at time 0, the origin, which lives for a while before it first splits. Everything up to that split is the stem, and it is ordinary simulated time: genes are gained and lost along it, traits drift along it. Every tree ZOMBI2 writes therefore gives its root a branch length.](figures/stem.pdf){width=70%}

Some tree viewers do not draw the stem by default, so a plotted tree can look shorter than the height the run reports. We recommend plotting trees with Phylustrator (`pip install phylustrator`), a companion package to ZOMBI2 maintained by me.

## Rates

In ZOMBI2 different events happen over time with different probability: speciations and extinctions at the species level; duplications, transfers, losses, originations (the birth of a new gene family) and rearrangements at the genome level; substitutions in a sequence; changes in a trait.

How often an event occurs is controlled by its **effective rate**:

$$\text{effective rate} = \text{scope}(\text{base}) \times \text{modifiers}$$

The **base** is the expected number of events per unit time, for one unit of whatever the scope counts. The **scope** names that unit (a lineage, a gene copy, a site) and so decides how many independent chances run at once: a loss rate of 0.25 counted per copy, in a genome of forty copies, is a total loss rate of 40 × 0.25 = 10 per unit time; counted per lineage, it stays 0.25 however large the genome grows. The **modifiers** are dimensionless multipliers. There are different types, check Appendix A for the full reference: the units, each level's default scope, every modifier and which levels accept it, and what the engine does with a rate once it has one.

## Dependent runs

ZOMBI2 can also simulate dependencies between the levels, so that complex scenarios can be simulated. For example, we can simulate a scenario like the evolution of mammals and their olfactory genes. A habitat trait switches between aquatic and terrestrial along the tree, and aquatic lineages lose those genes four times faster.

Every dependency is written as a **connection**, with three parts:

- the **driver**: the value the connection takes as input, the habitat of each lineage.
- the **target**: the parameter that depends on it, the loss rate.
- the **link**: what joins them. It takes the driver, applies some transformation and returns the value for the target.

![A connection: a habitat trait is simulated first and held fixed, and a gene loss rate depends on it. The driver, the link and the target are the three parts every connection has. More details can be found in Chapter 8.](figures/conditioning_print.png){width=95%}

In this example the habitat can be simulated first, and the genome run then depends on its finished history. The run is **conditioned**: two ordinary runs, in order.

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species}) \cdot P(\text{Genomes} \mid \text{Species}, \text{Traits})$$

Sometimes neither level can be simulated first. A trait controls how fast a lineage speciates, and the trait evolves along the very tree those speciations build; or a key gene decides diversification, and the tree decides which genomes exist. Then one run simulates both levels at the same time, the run is **joint**, and the tree comes out as a result rather than going in as an input. There is no conditional probability to write, only a joint one:

$$P(\text{Species}, \text{Traits})$$

The test is one question: can the driver be simulated first, on its own, and handed over? If it can, the run is conditioned. If it cannot, it is joint. Chapter 8 works through both, and Appendix C is the full map of what can drive what.

## What it can do

Some questions that can be explored with ZOMBI2:

- **How well does a reconciliation method recover the truth?** Evolve gene families under duplication, transfer and loss, and you get every family's true gene tree with the event behind every node. Reconcile against the species tree and score what the method found against what happened.

- **Can a transfer be detected when the donor is extinct?** Transfers can come from lineages that later go extinct, so a gene arrives in a survivor from a donor that leaves no other trace. The event log registers the donor, so you can ask how often a detection method finds a transfer whose donor is no longer on the tree.

- **What does genome reduction look like in parasitic bacteria?** Evolve a lifestyle trait, free-living or parasitic, and make the loss rate dependent on this trait. Lineages that become parasitic shed genes faster, so you can measure how much of the genome-size reconstruction a method attributes to lifestyle rather than to shared ancestry.

- **Can Bayesian methods of inference (relaxed molecular clocks) recover the true ages of a tree?** Give the substitution rate a relaxed clock, so lineages evolve at different paces, and compare the dates a method infers from the alignment against the true node ages.

- **How accurate is ancestral sequence reconstruction?** A run records the sequence at every internal node, not just the tips, so a reconstruction can be compared residue by residue against the real one.

- **Is a trait correlation real, or an artefact of the tree?** Two traits evolving on the same tree look correlated at the tips whether or not there is a connection between the two, because they share ancestry. Simulate without the connection, on the same tree, and you have the baseline any comparative method has to beat.

- **Does a trait actually drive diversification?** Let a trait set the speciation rate, so the tree and the trait grow together, and test whether a state-dependent method recovers the effect.

- **What signal survives in gene order?** With genes placed on chromosomes at the **ordered** genome resolution (Chapter 4), inversions, transpositions and translocations rearrange them, and fissions and fusions change the karyotype itself, so synteny and rearrangement inference methods can be tested.

## Installing it

ZOMBI2 needs Python 3.10 or newer and depends on NumPy and tqdm, plus the `tomli` backport on Python 3.10:

```bash
pip install zombi2
```

`zombi2 --version` confirms the install, and `zombi2 -h` lists the commands.

ZOMBI2 is pure Python over NumPy, with no compiled part to build, and the test suite runs on **Linux,
macOS and Windows** on every change, so the same command does the same thing on all three. Two things differ on
Windows and are worth knowing before they bite:

- **Paths in a rate expression.** A driver path goes inside the rate. On the command line ZOMBI2
  reads the backslashes as written, so `PerCopy(0.25).scaled_by('C:\Users\me\trait_events.tsv',
  {...})` works as pasted. In a Python script the same line is Python source, and Python reads `\U`
  as an escape, so it is a `SyntaxError`. Worse, `C:\temp` is read quietly as `C:` followed by a
  tab. Write the path as a raw string there: `r'C:\Users\me\trait_events.tsv'`. Forward slashes work
  in all three places, the command line, a Python script and a `--params` file, and Windows
  accepts them.
- **Paths in a `--params` file**, the TOML file that can carry a run's flags (`--params run.toml`).
  TOML, not ZOMBI2, reads that file, and TOML's ordinary `"…"` string
  treats a backslash as an escape, so `C:\Users` fails there with a message about a hex value. Put a
  value containing a path in a TOML **literal** string instead, `'''…'''`, which is taken exactly as
  written:

  ```toml
  transfer-to = '''Recipients().weighted_by('C:\Users\me\trait_events.tsv', {'competent': 2.0})'''
  ```


## The examples gallery

ZOMBI2 includes an [examples gallery](https://aadavin.github.io/zombi2/gallery.html) on how to
simulate many different scenarios: a page of worked examples, each with the exact commands that made
it. The gallery numbers them by section: `Sp` species, `Ge` genomes, `Sq` sequences and `Tr` traits
for the four levels, then `Co` conditioning and `Jo` joining, the two kinds of dependent run
(Chapter 8). This book cites them the way it cites a figure. For example, a paragraph that describes
what an inversion does to a chromosome has a developed example in the gallery **(Ge9)**.

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

Open `run.zombi2` first. It is plain text, one page: it names every file the run wrote. `--flat` writes every file straight into `out/` instead, and writes no run report; `--quiet` turns off the progress bar.

Here is the same model in Python. Every genome rate starts at 0, so a call that names none of them still runs, and gives a history with no duplication, transfer or loss, just the families it started with, copied down the tree.

```python
from zombi2 import species, genomes

sp  = species.simulate_species_tree(birth=1.0, n_extant=20, seed=1)
gen = genomes.simulate_genomes_family(sp, duplication=0.2, transfer=0.1,
                                      loss=0.25, origination=0.5, seed=42)

sp.write("out/species/")
gen.write("out/genomes/")
```
