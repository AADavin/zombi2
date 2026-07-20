# A tour of ZOMBI2

ZOMBI2 simulates evolution at four levels, and lets you either run them one after another or grow them together. This chapter introduces the four levels, the three ways they can relate, and the single shape every rate takes. It is the vocabulary the rest of the book uses.

## The four levels of ZOMBI2

ZOMBI2 simulates evolution at four different levels:

- **Species** — the tree of lineages: a strictly bifurcating rooted tree, with branches measured in time. Every ZOMBI2 workflow has some species tree, so this is the first thing to run.
- **Genomes** — the genes that exist in each lineage of the tree. Genomes can be simulated at different levels of resolution: a full genome represented with nucleotides, or just a set of gene families. Genomes are always simulated within a species tree; you always need one to obtain them.
- **Sequences** — the nucleotides or amino acids inside each gene. Sequences *always* evolve on gene trees, so you need to simulate genomes before using this level.
- **Traits** — phenotypes evolving along a tree: body size, optimal growth temperature, the presence or absence of a flagellum. Traits can influence the evolution of the other levels in different ways that we cover in this chapter.

The most general way to connect the four levels is shown in Figure 1.

![The four levels of ZOMBI2. Species, genomes and sequences form a chain of ancestry: a genome lives on the species tree, a sequence inside a gene. Traits branch to the side, because a trait can ride any species tree.](figures/fig-2-1-four-levels_print.png){width=40%}

You can run a simulation in which every level is simulated. We write it in this notation:

\begin{center}
P(Species) · P(Genomes | Species) · P(Sequences | Genomes) · P(Traits | Species)
\end{center}

But you can always choose which levels to run; you need not run them all. For example, you do not need to simulate sequences if you are only interested in gene trees, which the genome level already produces:

\begin{center}
P(Species) · P(Genomes | Species)
\end{center}

And you do not have to simulate genomes if all you want is a species tree with some traits on it:

\begin{center}
P(Species) · P(Traits | Species)
\end{center}

In ZOMBI2 everything depends on a species tree, and in most cases you begin a workflow by simulating the tree alone. There are a few exceptions, which we cover a bit later.

## Rates

In ZOMBI2 everything is driven by events that fire over time. The kind of event depends on the level being simulated. Some of the basic events are:

- **Species** — speciations and extinctions
- **Genomes** — duplications, transfers and losses
- **Sequences** — mutations
- **Traits** — phenotypic changes

Time is imposed by the species tree, and every rate is measured against that time scale. If your tree runs from 0 at the root to 1 at the tips, your simulation lasts one unit of time. Time is normally measured from the **crown** of the species tree, but you can instead set time zero at the **stem**. The difference is easiest to see in Figure 2.

![What `age` measures. With `age_type='crown'` (left) the age is the depth from the crown, the first speciation, to the present. With `age_type='stem'` (right) it is measured from the origin, so a stem branch precedes the crown.](figures/age_crown_print.png){width=92%}

A rate always has units of time⁻¹, on the scale imposed by the species tree. In a phylogenetic context, though, a single global rate rarely makes sense for most events. For example, a substitution happens at a **site**, so a mutation rate is counted per site (mutations × time⁻¹ × per site): each site is an independent chance to mutate. A speciation happens to a **lineage**, so the speciation rate is counted per lineage (speciations × time⁻¹ × per lineage): each branch alive is an independent chance for the tree to split. And a gene is lost one gene copy at a time, so gene loss is counted per copy (loss × time⁻¹ × per gene-copy).

Rates can also be modified, which makes ZOMBI2 a flexible platform for all sorts of scenarios. We might give a gene family a constant loss rate across the whole species tree, except in one clade that we know tends to shed genes, say a symbiotic bacterium, by multiplying the rate there by some number greater than one. Or we might let gene families evolve at different speeds: an antimicrobial-resistance family very prone to transfer, a ribosomal-protein family the opposite.

At the end of the day, the frequency at which an event fires depends on its **effective rate**:

\begin{center}
Effective rate = scope(base) × modifiers.
\end{center}

The **base** is the speed of a single event (how fast), in units of inverse time. The **scope** wraps that base to say how many independent chances the event has: per lineage, per copy, or per site. The **modifiers** are context multipliers, dimensionless, that let one lineage or one family run faster than another.

By default, this is the scope ZOMBI2 uses at each level:

| Level | Counted per | "How fast" is set by |
|---|---|---|
| Species | lineage | the diversification process |
| Genomes | copy (or lineage) | the duplication / transfer / loss rates |
| Sequences | site | the substitution rate (times a clock) |
| Traits | lineage | the trait model |

The scope is fixed by the level; the **modifiers** are where a rate gains its flexibility. A modifier reads some piece of context — the current time, the standing diversity, the lineage a branch sits on — and returns a dimensionless factor that multiplies the base. ZOMBI2 ships a small, shared set:

| Modifier | What it does to the rate |
|---|---|
| `OnTime` | Follows a **time schedule**: one factor up to a breakpoint, another after it, and so on — a skyline. |
| `OnTotalDiversity` | **Slows as the tree fills up**: the factor falls from 1 toward 0 as the number of lineages approaches a carrying capacity, and stays there. |
| `FromParent` | Is **inherited from the parent lineage and nudged at each split**, so the rate drifts gradually down the tree and close relatives keep similar rates. |
| `ByLineage` | Is an **independent draw for each lineage**, with no memory of its parent, so nearby branches are no more alike than distant ones. |

Two of these are **deterministic**: `OnTime` and `OnTotalDiversity` are fixed functions of the state of the world, so every lineage that meets the same time, or the same diversity, gets the same factor. The other two are **random and vary from lineage to lineage**, and they differ in *memory*: `FromParent` is passed down and drifts, so the rate is autocorrelated along the tree — a slowly wandering clock, or a clade that inherits a fast tempo — whereas `ByLineage` is drawn afresh on every branch, so the variation is scattered, an uncorrelated ("relaxed") clock. The random modifiers are **mean-corrected**, meaning their factors average to 1, so switching on heterogeneity spreads a rate around without secretly speeding the whole tree up.

Modifiers **stack by multiplication**, so they combine: `1.0 * mod.OnTime({0: 1, 5: 0.3}) * mod.FromParent(spread=0.3)` is a rate that both follows a schedule and drifts between lineages. And because a modifier attaches to *any* rate, the same handful reappears at every level — `OnTime` is a skyline for speciation and an early burst for a trait, `FromParent` is clade drift for diversification and the autocorrelated clock for sequences, `OnTotalDiversity` is diversity-dependence wherever a rate should ease off as lineages accumulate. Learn them once and you know them everywhere.

A more detailed introduction to rates is given in Appendix A (Gillespie).

## Going beyond the basic simulation: conditioning and joining levels

Some evolutionary scenarios cannot be simulated by the paradigm we have described so far. For example:

- A trait evolves along a tree and itself controls how fast that tree speciates.
- Gene content decides survival: lineages that acquire a key gene diversify faster than the rest.
- Two levels feed back on each other: a trait raises a gene family's loss rate, while carrying that family in turn pulls on the trait, so the two co-evolve.

ZOMBI2 adds two new connections between levels to reach scenarios like these. One is **conditioning**, the other is **joining**. Because they need some care, each has its own chapter; here we give only a brief overview.

When we **condition** one level on another, a parameter of the second level stops being a fixed number you set and instead reads the state of the first. In the run below, the gene-loss rate is no longer constant: it depends on the trait's value on each branch. For example, you can use this to simulate that aquatic lineages lose their olfactory genes faster. Because the trait can be simulated first and then held fixed, conditioning is still two ordinary runs in order: you simulate the driver level, write it out, and feed it to the second, exactly as you already feed a species tree to a genome run.

\begin{center}
P(Species) · P(Traits | Species) · P(Genomes | Species, Traits)
\end{center}

When we **join** two levels, neither can be simulated first, because each depends on the other as it unfolds. If a trait speeds up speciation, then faster-speciating lineages leave more descendants, so the tree's shape depends on the trait, but the trait is evolving along that very tree at the same time. No order works, so the two levels are grown together in a single run, step by step. When the coupling reaches back into the species tree, the tree itself becomes an output rather than an input.

\begin{center}
P(Species, Traits)
\end{center}

We **join** whenever a coupling would form a loop: when one level shapes another and is shaped back, directly or through the tree. If the influence runs only one way, we condition; if it runs in a loop, we join. In the directional case we can still name the variable the driver sets on its target: for a trait driving speciation, that variable is the speciation rate.
