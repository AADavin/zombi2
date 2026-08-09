# A tour of ZOMBI2

This chapter is the whole tool in one pass: the levels, how time and rates work, how one level drives another, and how a run is written in Python and on the command line.

## The four levels of ZOMBI2

ZOMBI2 simulates four levels of evolution.

- **Species**, the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time. Every workflow has a species tree, so this is the first thing to run.
- **Genomes**, the genes that exist in each lineage. Genomes can be simulated at three **resolutions**: gene families alone, genes placed on chromosomes, or a full nucleotide genome. A genome is always simulated within a species tree.
- **Sequences**, the nucleotides or amino acids inside each gene. Sequences evolve on gene trees, so genomes must be simulated first.
- **Traits**, phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum.

![The four levels of ZOMBI2. Everything starts from the species tree, which forks: a genome evolves along it, and so does a trait. Sequences continue below genomes, because a sequence lives inside a gene.](figures/fig-2-1-four-levels_print.png){width=45%}

A run in which every level is simulated:

$$P(\text{Species}) \cdot P(\text{Genomes} \mid \text{Species}) \cdot P(\text{Sequences} \mid \text{Genomes}) \cdot P(\text{Traits} \mid \text{Species})$$

You need not run them all. Skip sequences if you only want gene trees, which the genome level already produces; skip genomes if you want a species tree with traits on it. Everything depends on a species tree, so a workflow almost always begins by simulating one alone. The exception is a **joint** model, below, where the tree is an output rather than an input.

## Time

ZOMBI2 is a forward simulator: evolution runs from an ancestral state at time 0 to the present.

Time is imposed by the species tree, and every rate is measured against that scale. If your tree runs from 0 at the origin to 1 at the tips, the simulation lasts one unit of time. Time 0 is the origin of the founding lineage, and every time you give ZOMBI2 is measured on it: the moment of a mass extinction, the breakpoints of a rate that changes through time.

The founding lineage lives for a while before it first splits. That stretch is the tree's **stem**, and the first split is its **crown**.

![The stem. A run starts from one lineage at time 0, the origin, which lives for a while before it first splits. Everything up to that split is the stem, and it is ordinary simulated time: genes are gained and lost along it, traits drift along it. Every tree ZOMBI2 writes therefore gives its root a branch length.](figures/stem.pdf){width=70%}

Some tree viewers do not draw the stem by default. It is there, and evolution takes place in it.

## Rates

Everything is driven by events that fire over time: speciations and extinctions at the species level; duplications, transfers, losses, originations and rearrangements at the genome level; substitutions in a sequence; changes in a trait.

A rate is written from its **scope**, which answers *per what?*:

```python
from zombi2.params.scope import PerCopy, PerLineage

loss = PerCopy(0.25)      # each gene copy is lost at 0.25 — a big genome loses often
loss = PerLineage(0.25)   # the lineage loses at 0.25, whatever its genome holds
```

Those are two different models, not two spellings of one. Multiply the first by a genome of a thousand genes and it fires a thousand times as often as with one; the second never notices the genome's size. A bare number takes the event's default scope, which is why most runs write no scope at all — but where an event offers a choice, write it. Appendix A lists which rates have one.

How often an event fires is its **effective rate**:

$$\text{effective rate} = \text{scope}(\text{base}) \times \text{modifiers}$$

The **base** is the speed of a single event, in units of inverse time. The **scope** wraps it to say how many independent chances the event has right now: per lineage, per copy, per site. The **modifiers** are the dimensionless multipliers a rate picks up from its context.

You do not write a modifier. You chain a **verb** onto the rate, and the verb says what the number does to it. There are three: `scaled_by` multiplies the base, `set_by` replaces it in the rate's own units so nothing is written in front, and `weighted_by` compares the candidates of a **choice** — `transfer_to` is the only one. Only `scaled_by` multiplies, so only it is in the equation above.

The verb's first argument is the **driver**, the thing the rate reads. Each kind of name has its own module — scopes in `params.scope`, drivers in `params.driver`, laws in `params.law`, distributions in `params.distributions` — and all of them are also importable straight from `zombi2.params` once you know which is which. One row here per thing you can depend on, not per verb:

| The rate depends on | Written |
|---|---|
| the clock | `PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})` |
| how many lineages are standing | `PerLineage(1.0).scaled_by(TotalDiversity(cap=100))` |
| chance, once per lineage | `PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3))` |
| chance, once per gene family | `PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))` |
| the parent's value, drifting | `PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.3)))` |
| where in the tree a lineage sits | `PerCopy(0.25).scaled_by(Clade({'cave': ['n12']}), {'cave': 3.0})` |
| a value another level evolved | `PerCopy(0.25).scaled_by(habitat, {'aquatic': 4.0})` |

Two of those drivers are written so often that each has a verb of its own — `changing_at` for the clock and `varying_among` for chance — and `scaled_by` refuses both by name, so there is one spelling for each rather than two.

Verbs chain, and their factors multiply:

```python
from zombi2.params.distributions import LogNormal

# loss triples after time 2, and varies from family to family on top of that
loss = PerCopy(0.25).changing_at({0: 1.0, 2: 3.0}).varying_among('families', LogNormal(0.0, 0.5))
```

Much of what the literature names as a model is one verb on one rate:

| What it is usually called | What it is here |
|---|---|
| skyline / episodic birth–death | `changing_at` on `birth` |
| diversity-dependent diversification | `scaled_by(TotalDiversity(...))` on `birth` |
| uncorrelated ("relaxed") molecular clock | `varying_among('lineages', ...)` on `substitution` |
| autocorrelated clock | `varying_among('lineages', Drift(...))` on `substitution` |
| rate heterogeneity across gene families | `varying_among('families', ...)` on a genome rate |
| state-dependent diversification (BiSSE and kin) | `scaled_by` on `birth`, in a joint run |

None of those is a separate code path with its own function and its own parameters. They are the same grammar pointed at different rates, which is why they combine: a relaxed clock *and* an early burst is one rate with two verbs, not two models.

Not every driver is available at every level. Some combinations mean nothing — a species tree has no gene families for `varying_among('families', ...)` to vary over — and others are simply not implemented yet. Either way the level refuses the rate and says what it does take, rather than ignoring it. `zombi2 <command> -h` lists them, and Appendix A is the full reference: the units, each level's default scope, what each level accepts, how to write a driver of your own, and the Gillespie algorithm that turns rates into events.

## Conditioning

Running the levels in sequence already makes each one depend on the tree. Sometimes you want more: a rate that reads a value which varies from lineage to lineage, rather than a fixed number. That is **one thing driving another**, and writing it is three questions.

**What depends?** A rate, or — at the genome level — an **extent** (how much an event takes) or the **transfer recipient** (which lineage a transfer goes to).

**On what?** The driver: a trait's state, a gene family's presence, a sequence's GC content.

**How?** The mapping, which turns each value of the driver into a number: a table over named states, a curve over a number.

The verb is not a fourth thing to learn. It follows from the first answer: `scaled_by` on a rate or an extent, `weighted_by` on a recipient, and `set_by` when you mean to replace the base rather than multiply it.

Take olfactory genes. A habitat trait switches between aquatic and terrestrial along the tree, and aquatic lineages lose those genes four times faster. Gene loss is what depends, the habitat is what it depends on, and the mapping turns one into the other: on an aquatic branch the loss rate is `0.25 × 4`, elsewhere `0.25 × 1`.

![Conditioning. The **driver** is a value already simulated — a habitat trait here, with the two states each lineage switches between shown below it. What depends is a rate in the run that comes next, gene loss here, and the arrow between them carries the **mapping**: one multiplier per state of the driver, so a branch's habitat sets that branch's loss rate. The driver is finished and written to a file before the second run starts, which is what makes this two ordinary commands.](figures/conditioning_print.png){width=95%}

What makes this **conditioning** is that the driver can be finished before the target starts. The habitat is unaffected by how many genes a lineage has, so the trait run completes and writes its event log, and the genome run reads that log and looks up a multiplier per branch. Nothing about it needs a special engine, and the run factorises:

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species}) \cdot P(\text{Genomes} \mid \text{Species}, \text{Traits})$$

Driver and target need not be at different levels. One gene family can drive another's rates, one character can drive another's, one gene's sequence can drive another's, because a level that holds many separable objects can be run twice: the first run is finished, and the second reads it. Chapter 9 works through the whole of it, within a level and across.

## Joining

**Joining** is for when the driver cannot be finished first. Two scenarios need it:

- A trait evolves along a tree and controls how fast that tree speciates.
- Gene content decides survival: lineages that acquire a key gene diversify faster.

If a trait speeds up speciation, faster-speciating lineages leave more descendants, so the tree's shape depends on the trait while the trait is evolving along that very tree. No order works, so both are grown together in one run, and the tree becomes an output rather than an input.

![Joining, with the same parts and one thing that changes everything. Body size drives speciation through the same `scaled_by`, and the speciation rate makes the tree — but that tree is the one the trait is evolving along, so the loop closes. Note where: not from the rate back to the trait, since a speciation rate does not change a body size. It decides which lineages split, and every split hands the parent's state to both daughters. The tree is the third thing in the cycle, and it is the whole difference between the two figures — a fixed input there, an output here.](figures/joining_print.png){width=95%}

Because neither is finished first, neither can be written out and handed over. There is no conditional probability to write, only a joint one:

$$P(\text{Species}, \text{Traits})$$

The test is one question: can the driver be grown first, on its own, and handed over? If it can, condition. If it cannot, join. Chapter 9 works through conditioning and Chapter 10 through joining.

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

Each call takes the object it depends on: genomes and traits read the species result, sequences reads the genomes result. So the script reads top to bottom in the `P(·)` order above.

## Using ZOMBI2 from the CLI

The same simulations run from the command line. Each level is a subcommand of `zombi2`, and its flags are the long-form names of the Python arguments:

```bash
# a dated species tree (20 extant tips)
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1

# gene families along it
zombi2 genomes out/ --duplication 0.2 --loss 0.25 --origination 0.5 --seed 42
```

A rate flag takes a rate **written exactly as you would write it in Python**: a bare number, or a scope with verbs chained onto it, quoted so the shell keeps it in one piece.

```bash
# speciation drops to a third of its rate at time 3 (a skyline). --n-extant is conditioned on
# survival, so it runs for any seed; a --total-time run can go extinct and stop with an error.
zombi2 species skyline/ --birth "PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})" \
    --death 0.3 --n-extant 30 --seed 1
```

Every command takes one positional argument, the **run directory** — both where it writes and where it reads the level before it, so a pipeline names one directory per command and nothing is passed by hand. The CLI covers all four levels, and `zombi2 joint` grows a tree and its driver in one pass. Conditioning has no command of its own: a driven rate is written on the level command that reads it, like any other rate. Appendix C is the command reference.

## Output in ZOMBI2

Every run can be written with `result.write("out/", outputs=[...])`; with no `outputs` it writes that level's **default** set. Trees are Newick, tables and event logs are TSV, sequences are FASTA. Branch lengths are in time everywhere except the sequence phylograms, which are in substitutions per site. Appendix B lists every file, level by level.

Each level draws its own stream of random numbers, started by that level's `seed`, so the same seed, parameters and version give the same run event for event. A run given no seed draws one and records it, so it can still be repeated.
