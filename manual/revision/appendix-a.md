# Appendix A — Rates in detail, and the Gillespie algorithm

Chapter 2 introduced the shape every rate takes, `effective rate = scope(base) × modifiers`, and gave a
first feel for what the scope and the modifiers do. This appendix is the full reference: how a rate's
units work, the default scope at each level, the complete catalogue of modifiers, and the Gillespie
algorithm that turns these rates into the events of a simulation.

## How a rate is counted: the scope

A rate always has units of time⁻¹, on the scale imposed by the species tree. In a phylogenetic context,
though, a single global rate rarely makes sense for most events. A substitution happens at a **site**, so
a mutation rate is counted per site (mutations × time⁻¹ × per site): each site is an independent chance to
mutate. A speciation happens to a **lineage**, so the speciation rate is counted per lineage (speciations
× time⁻¹ × per lineage): each branch alive is an independent chance for the tree to split. And a gene is
lost one gene copy at a time, so gene loss is counted per copy (loss × time⁻¹ × per gene-copy). The unit a
rate is counted in — per lineage, per copy, per site — is what we call its **scope**.

By default, this is the scope ZOMBI2 uses at each level:

| Level | Counted per | "How fast" is set by |
|---|---|---|
| Species | lineage | the diversification process |
| Genomes (unordered) | copy (or lineage) | the duplication / transfer / loss rates |
| Sequences | site | the substitution rate (times a clock) |
| Traits | lineage | the trait model |

## Bending a rate: modifiers

Rates can also be altered through **modifiers**, which makes ZOMBI2 a flexible platform for all sorts of
scenarios. We might give a gene family a constant loss rate across the whole species tree, except in one
clade that we know tends to shed genes, say a symbiotic bacterium, by multiplying the rate there by some
number greater than one. Or we might let gene families evolve at different speeds: an antimicrobial-
resistance family very prone to transfer, a ribosomal-protein family the opposite.

The common modifiers you can reach are:

| Modifier | What it does to the rate |
|---|---|
| `OnTime` | Follows a **time schedule**: one factor up to a breakpoint, another after it. |
| `OnTotalDiversity` | **Slows as the tree fills up**: the factor falls from 1 toward 0 as the number of lineages approaches a carrying capacity, and stays there. |
| `FromParent` | Is **inherited from the parent lineage and nudged at each split**, so the rate drifts gradually down the tree and close relatives keep similar rates. |
| `ByLineage` | Is an **independent draw for each lineage**, with no memory of its parent, so nearby branches are no more alike than distant ones. |

Two of these are **deterministic**: `OnTime` and `OnTotalDiversity` are fixed functions of the state of
the world, so every lineage that meets the same time, or the same diversity, gets the same factor. The
other two are **random and vary from lineage to lineage**, and they differ in *memory*: `FromParent` is
passed down and drifts, so the rate is autocorrelated along the tree — a slowly wandering clock, or a
clade that inherits a fast tempo — whereas `ByLineage` is drawn afresh on every branch, so the variation
is scattered, an uncorrelated ("relaxed") clock. The random modifiers are **mean-corrected**, meaning
their factors average to 1, so switching on heterogeneity spreads a rate around without secretly speeding
the whole tree up.

Modifiers **stack by multiplication**, so they combine: `1.0 * mod.OnTime({0: 1, 5: 0.3}) *
mod.FromParent(spread=0.3)` is a rate that both follows a schedule and drifts between lineages. And because
a modifier attaches to *any* rate, the same handful reappears at every level — `OnTime` is a skyline for
speciation and an early burst for a trait, `FromParent` is clade drift for diversification and the
autocorrelated clock for sequences, `OnTotalDiversity` is diversity-dependence wherever a rate should ease
off as lineages accumulate. Learn them once and you know them everywhere.

## The Gillespie algorithm

*(To be carried, unchanged, from the current manual — how ZOMBI2 draws waiting times from these effective
rates and fires one event at a time.)*
