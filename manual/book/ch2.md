# A tour of ZOMBI2

This chapter introduces the four levels ZOMBI2 simulates, the three ways they can relate, and the single shape every rate takes.

## The four levels of ZOMBI2

- **Species** — the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time. Every workflow has a species tree, so this is the first thing to run.
- **Genomes** — the genes that exist in each lineage. Genomes can be simulated at three **resolutions**: gene families alone, genes placed on chromosomes, or a full nucleotide genome. A genome is always simulated within a species tree.
- **Sequences** — the nucleotides or amino acids inside each gene. Sequences evolve on gene trees, so genomes must be simulated first.
- **Traits** — phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum.

![The four levels of ZOMBI2. Species, genomes and sequences form a chain of ancestry: a genome lives on the species tree, a sequence inside a gene. Traits branch to the side, because a trait can ride any species tree.](figures/fig-2-1-four-levels_print.png){width=40%}

A run in which every level is simulated:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species}) \cdot P(\text{Sequences} \mid \text{Genomes}) \cdot P(\text{Traits} \mid \text{Species})$$

You need not run them all. Skip sequences if you only want gene trees, which the genome level already produces:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species})$$

Skip genomes if you want a species tree with traits on it:

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species})$$

Everything depends on a species tree, so a workflow almost always begins by simulating one alone. The exception is a **joint** model, below, where the tree is an output rather than an input.

## Time

ZOMBI2 is a forward simulator: evolution runs from an ancestral state at time 0 to the present.

Time is imposed by the species tree, and every rate is measured against that scale. If your tree runs from 0 at the root to 1 at the tips, the simulation lasts one unit of time. Time 0 is the origin of the founding lineage, and every time you give ZOMBI2 — the moment of a mass extinction, the breakpoints of a rate that changes through time — is measured on it.

The founding lineage lives for a while before it first splits. That stretch is the tree's **stem**, and the first split is its **crown**.

![The stem. A run starts from one lineage at time 0 — the origin — which lives for a while before it first splits. Everything up to that split is the stem, and it is ordinary simulated time: genes are gained and lost along it, traits drift along it. Every tree ZOMBI2 writes therefore gives its root a branch length.](figures/stem.pdf){width=70%}

## Rates

Everything is driven by events that fire over time:

- **Species** — speciations and extinctions
- **Genomes** — duplications, transfers, losses, originations, inversions, transpositions
- **Sequences** — substitutions
- **Traits** — phenotypic changes

How often an event fires is its **effective rate**:

$$\text{effective rate} = \text{scope}(\text{base}) \times \text{modifiers}$$

The **base** is the speed of a single event (how fast), in units of inverse time. The **scope** wraps it to say how many independent chances the event has: per lineage, per copy, per site. The **modifiers** are dimensionless multipliers that make a rate faster or slower depending on context — the lineage, the gene family, the total diversity present.

You rarely need to touch the default scope or the modifiers. With them you can simulate:

- an early burst that tapers off, by giving the rate a schedule that starts high and drops;
- a molecular clock that speeds and slows along the tree, with close relatives ticking alike, by letting each lineage inherit its rate from its parent and drift at every split;
- a radiation that eases off as the clade fills up, by having the speciation rate read the diversity present at each moment.

**Appendix A** is the full rate reference — the units, each level's default scope, the catalogue of modifiers — and the Gillespie algorithm that turns rates into events.

## Going beyond: conditioning and joining levels

Some scenarios need more than the levels running in sequence:

- A trait evolves along a tree and controls how fast that tree speciates.
- Gene content decides survival: lineages that acquire a key gene diversify faster.
- Two levels feed back: a trait raises a gene family's loss rate, and carrying that family pulls on the trait.

ZOMBI2 adds two connections for these — **conditioning** and **joining**. Each has its own chapter; this is the shape of them.

**Conditioning** makes a parameter of one level read the state of another instead of being a fixed number. Aquatic lineages lose their olfactory genes faster: the gene-loss rate depends on the trait's value on each branch. Because the trait can be simulated first and then held fixed, this is still two ordinary runs in order — simulate the driver, write it out, feed it to the second, exactly as you already feed a species tree to a genome run.

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species}) \cdot P(\text{Genomes} \mid \text{Species}, \text{Traits})$$

**Joining** is for when neither level can go first. If a trait speeds up speciation, faster-speciating lineages leave more descendants, so the tree's shape depends on the trait — while the trait is evolving along that very tree. No order works, so both are grown together in one run, and the tree becomes an output rather than an input.

$$P(\text{Species}, \text{Traits})$$

The test is whether the influence runs one way or in a loop. One way: condition. In a loop — one level shaping another and being shaped back, directly or through the tree — join.

## Using ZOMBI2 in Python

Each level is a function in its own subpackage, and they compose by feeding one level's result into the next. A run returns a *result object*; read it in the session or write it to disk with `.write()`. A whole workflow is a short script:

```python
from zombi2 import species, genomes, sequences, traits
from zombi2.sequences import substitution_models as sm

sp   = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
gen  = genomes.simulate_genomes_family(sp, duplication=0.2, loss=0.25, origination=0.5, seed=42)
seqs = sequences.simulate_sequences(gen, model=sm.hky85(kappa=2.0), length=300, seed=7)
bm   = traits.simulate_continuous(sp, rate=1.0, seed=1)
```

Each call takes the object it depends on — genomes and traits read the species result, sequences reads the genomes result — so the script reads top to bottom in the `P(·)` order above.

## Using ZOMBI2 from the CLI

The same simulations run from the command line. Each level is a subcommand of `zombi2`, and its flags are the long-form names of the Python arguments:

```bash
# a dated species tree (20 extant tips)
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1

# gene families along it
zombi2 genomes out/ --duplication 0.2 --loss 0.25 --origination 0.5 --seed 42
```

A rate flag takes a rate **written exactly as you would write it in Python** — a bare number, or a scope wrapper and modifiers composed with `*`, quoted so the shell keeps it in one piece:

```bash
# speciation drops to a third of its rate at time 3 (a skyline)
zombi2 species out/ --birth "1.0 * OnTime({0: 1.0, 3: 0.3})" --death 0.3 --total-time 5 --seed 1
```

Every command takes one positional argument, the **run directory**. It is both where that command writes and where it reads the level before it, so a pipeline is the same directory named once per command and nothing is passed by hand. `--from` overrides the reading half, for a tree from elsewhere or a run you would rather not write into; a `--params` TOML file can hold a whole pipeline's settings.

Because the levels share one directory, a command refuses to re-run a level in place when a later level was built from it — that would leave the later output out of step. `--force` re-runs anyway and removes the now-stale downstream. The CLI covers all four levels; the coupled models are run from Python until their commands land.

## Output in ZOMBI2

Every run can be written with `result.write("out/", outputs=[...])`; with no `outputs` it writes that level's **default** set. The formats are uniform — **trees** in Newick, **tables and event logs** in TSV, **sequences** in FASTA. Branch lengths are in time everywhere except the sequence phylograms, which are in substitutions per site.

At every level the **event log** (`*_events.tsv`) is the true, ordered history the run followed: the source of truth the summaries are derived from. Appendix B lists every file, level by level.
