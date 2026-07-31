# A tour of ZOMBI2

This chapter introduces the most important concepts to understand how ZOMBI2 works and what it can do.

## The four levels of ZOMBI2

ZOMBI2 can simulate different levels of evolution, namely:

- **Species** — the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time. Every workflow has a species tree, so this is the first thing to run.
- **Genomes** — the genes that exist in each lineage. Genomes can be simulated at three **resolutions**: gene families alone, genes placed on chromosomes, or a full nucleotide genome. A genome is always simulated within a species tree.
- **Sequences** — the nucleotides or amino acids inside each gene. Sequences evolve on gene trees, so genomes must be simulated first.
- **Traits** — phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum.

![The four levels of ZOMBI2. Everything starts from the species tree, which forks: a genome evolves along it, and so does a trait. Sequences continue below genomes, because a sequence lives inside a gene.](figures/fig-2-1-four-levels_print.png){width=45%}

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

The founding lineage lives for a while before it first splits. That stretch is the tree's **stem**, and the first split is its **crown** — the figure below shows both.

![The stem. A run starts from one lineage at time 0 — the origin — which lives for a while before it first splits. Everything up to that split is the stem, and it is ordinary simulated time: genes are gained and lost along it, traits drift along it. Every tree ZOMBI2 writes therefore gives its root a branch length.](figures/stem.pdf){width=70%}

Be careful because some tree viewers do not show the stem by default: trees simulated with ZOMBI2 do include a stem in which evolution takes place. 

## Rates

Everything is driven by events that fire over time:

- **Species** — speciations and extinctions
- **Genomes** — duplications, transfers, losses, originations, inversions, transpositions
- **Sequences** — substitutions
- **Traits** — phenotypic changes

How often an event fires is its **effective rate**:

$$\text{effective rate} = \text{scope}(\text{base}) \times \text{modifiers}$$

The **base** is the speed of a single event (how fast), in units of inverse time. The **scope** wraps it to say how many independent chances the event has: per lineage, per copy, per site. The **modifiers** are dimensionless multipliers that make a rate faster or slower depending on context — the lineage (which branch of the species tree), the gene family, the total diversity present.

A bare number is a valid rate, and the default scope is the right one for each event, so most runs write neither. The modifiers are where the range of the tool lives, and they get a section of their own below.

**Appendix A of the PDF manual** is the full rate reference — the units, each level's default scope, the catalogue of modifiers — and the Gillespie algorithm that turns rates into events.

## Modifiers

The base of a rate says how fast. A modifier says **what it depends on**:

| Modifier | The rate depends on | Written |
|---|---|---|
| `OnTime` | the clock — a schedule of intervals | `OnTime({0: 1.0, 3: 0.3})` |
| `OnTotalDiversity` | how many lineages are standing right now | `OnTotalDiversity(cap=100)` |
| `FromParent` | the parent's value, drifting at each split | `FromParent(spread=0.3)` |
| `ByLineage` | the lineage, drawn independently | `ByLineage(spread=0.3)` |
| `ByFamily` | the gene family, drawn independently | `ByFamily(spread=0.5)` |
| `DrivenBy` | **a driver** — a trait's state, a gene's presence | `DrivenBy('trait', {'hot': 4.0})` |

Each is a dimensionless multiplier, so they multiply, and a rate can carry several:

```python
from zombi2.rates import modifiers as mod

# loss triples after time 2, and varies from family to family on top of that
loss = 0.25 * mod.OnTime({0: 1.0, 2: 3.0}) * mod.ByFamily(spread=0.5)
```

Much of what the literature names as a model is one modifier on one rate:

| What it is usually called | What it is here |
|---|---|
| skyline / episodic birth–death | `OnTime` on `birth` |
| diversity-dependent diversification | `OnTotalDiversity` on `birth` |
| uncorrelated ("relaxed") molecular clock | `ByLineage` on `substitution` |
| autocorrelated clock | `FromParent` on `substitution` |
| rate heterogeneity across gene families | `ByFamily` on a genome rate |
| state-dependent diversification (BiSSE and kin) | `DrivenBy` on `birth` |

None of those is a separate code path with its own function and its own parameters. They are the same grammar pointed at different rates, which is why they combine: a relaxed clock *and* an early burst is one rate with two modifiers, not two models.

Not every modifier is available at every level. Some combinations mean nothing — a species tree has no gene families for `ByFamily` to vary over — and others are simply not implemented yet. Either way the level refuses the rate and says which modifiers it does take, rather than ignoring it. `zombi2 <command> -h` lists them for that level, and Appendix A of the PDF manual gives the full table.

## Conditioning

Some scenarios need more than the levels running in sequence. A rate can read something that varies from lineage to lineage instead of being a fixed number: a habitat trait makes gene loss faster on some branches than others, a gene's presence makes a lineage speciate faster. That is **one thing driving another**, and it has three parts:

- the **driver** — the value that is read: a trait's state, a gene family's presence.
- the **target** — the rate that reads it.
- the **modifier** — `DrivenBy`, which says what each value of the driver does to that rate.

Take olfactory genes. A habitat trait switches between aquatic and terrestrial along the tree, and aquatic lineages lose those genes four times faster. The habitat is the driver, gene loss is the target, and the modifier is what turns one into the other: on a branch that is aquatic the loss rate is `0.25 × 4`, and elsewhere it is `0.25 × 1`.

![Conditioning has three parts. The **driver** is a value that was already simulated — here a habitat trait, whose two states each lineage switches between are shown below it. The **target** is a rate in a later run, here gene loss. The **modifier**, `DrivenBy`, carries a multiplier for each state of the driver, so a branch's habitat sets that branch's loss rate. The driver is finished before the target's run begins, which is what makes this two ordinary runs rather than one joint one.](figures/conditioning_print.png){width=95%}

What decides how the run is organised is whether the driver can be finished before the target starts. Here it can, and that is **conditioning**: the trait run finishes and writes its event log, and the genome run reads that log and looks up a multiplier per branch. Nothing about it needs a special engine — it is the same handoff as feeding a species tree to a genome run, one level further along, which is why the whole thing factorises:

$$P(\text{Species}) \cdot P(\text{Traits} \mid \text{Species}) \cdot P(\text{Genomes} \mid \text{Species}, \text{Traits})$$

## Joining

**Joining** is for when the driver cannot be finished first. Two scenarios need it:

- A trait evolves along a tree and controls how fast that tree speciates.
- Gene content decides survival: lineages that acquire a key gene diversify faster.

If a trait speeds up speciation, faster-speciating lineages leave more descendants, so the tree's shape depends on the trait — while the trait is evolving along that very tree. No order works, so both are grown together in one run, and the tree becomes an output rather than an input.

![Joining, with the same three parts and a fourth thing that changes everything. Body size drives speciation through the same `DrivenBy` modifier, and the speciation rate creates the tree — but that tree is the one the trait is evolving along, so the loop closes. Note where it closes: not from the rate back to the trait, since a speciation rate does not change a body size. It decides which lineages split, and every split hands the parent's trait state to both daughters. The tree is the third thing in the cycle, and it is the difference between the two figures: in conditioning it is a fixed input, given before either level runs, and here it is an output.](figures/joining_print.png){width=95%}

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
# speciation drops to a third of its rate at time 3 (a skyline). --n-extant is conditioned on
# survival, so it runs for any seed; a --total-time run can go extinct and stop with an error.
zombi2 species skyline/ --birth "1.0 * OnTime({0: 1.0, 3: 0.3})" --death 0.3 --n-extant 30 --seed 1
```

Every command takes one positional argument, the **run directory**. It is both where that command writes and where it reads the level before it, so a pipeline is the same directory named once per command and nothing is passed by hand. `--from` overrides the reading half, for a tree from elsewhere or a run you would rather not write into; a `--params` TOML file can hold a whole pipeline's settings.

Because the levels share one directory, a command refuses to re-run a level in place when a later level was built from it — that would leave the later output out of step. `--force` re-runs anyway and removes the now-stale downstream. The CLI covers all four levels, and `zombi2 joint` grows a tree and its driver in one pass. Conditioning has no command of its own: a driven rate is written on the level command that reads it, like any other rate.

## Output in ZOMBI2

Every run can be written with `result.write("out/", outputs=[...])`; with no `outputs` it writes that level's **default** set. The formats are uniform:

| Output | Format |
|---|---|
| trees | Newick |
| tables and event logs | TSV |
| sequences | FASTA |

Branch lengths are in time everywhere except the sequence phylograms, which are in substitutions per site. At every level the **event log** (`*_events.tsv`) is the true, ordered history the run followed, and Appendix B lists every file, level by level.

## Reproducing a run

One pseudo-random stream drives a whole run, and `seed` starts it: the same seed with the same parameters gives the same run, event for event. A run given no seed draws one and records it in the log, so an unseeded run can still be repeated. Three other things have to match:

- **The version.** A seed names a position in a stream of draws, not a history, so anything that changes how many draws a run makes, or in what order, moves everything after it. Every run log records the version it was made with, and `pip install zombi2==<that version>` gets the old run back.
- **The inputs, by content.** A run that reads a species tree, a `--tip-fates` table or a `DrivenBy` driver file is reproducible only with those files; the log records each by SHA-256 as well as by name.
- **Serial or parallel.** `parallel=` is a separate engine, so a parallel run and a serial run of the same seed differ; both are valid draws of the same process. A parallel run is identical for any worker count (Chapter 4).
