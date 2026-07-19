```{=latex}
\part{Coevolution}
```

# Coevolution: the framework

ZOMBI2 normally simulates in a **pipeline** — a species tree, then a trait along it, then gene
families along it, then sequences along the gene trees — because the couplings run one way and the
joint distribution factorises. This part is about what to do when they do *not*: when one level feeds
**back** into what produced another, and the levels must be grown **together**. A single idea covers
every such case; this chapter is that idea, and the chapters that follow work through the models,
grouped by *what each coupling shapes*.

## The pipeline, and when it breaks

The hierarchical pipeline works because the joint distribution factorises,

$$P(\text{tree}) \cdot P(\text{trait} \mid \text{tree}) \cdot P(\text{genes} \mid \text{tree}, \text{trait}),$$

so each stage can be drawn on the frozen output of the previous one. The moment a trait or a genome
feeds **back** into the process that generated the tree, this factorisation breaks: the levels must be
grown *together*. Coevolution is for these coupled scenarios.

## One grammar: driver, target, response

Every coupling in ZOMBI2 is one sentence:

> **driver → target-variable : response**

- **The driver** is whose state pushes — a *value along the tree* (a trait's value, a gene's presence
  or copy number) or an *event* (a speciation).
- **The target-variable** is the one quantity of the target that gets bent. Each level exposes a short,
  fixed menu, read straight off its process: for **species** it is `speciation` ($\lambda$) or
  `extinction` ($\mu$); for a **trait**, its OU `optimum` or its `value` at a split; for **gene
  content**, the `loss` / `gain` / `duplication` / `transfer` rates or `presence` at a split; for
  **sequences**, the `selection` ($\omega$ = dN/dS) or the substitution `speed`. Bending a *rate* is a
  **modulation** (a
  multiplier on it); bending a *state* is a **jump** (a shift at a node, or a moved optimum).
- **The response** is how the driver's value maps to the size of the effect: an exponential **scalar**
  link (`rate = base · exp(strength · driver)`, the default), a **table** of one value per discrete
  state (which recovers MuSSE), or a bounded **curve** (which recovers QuaSSE).

That one sentence spans everything from state-dependent diversification to trait-driven selection on
sequences. The named literature models — SSE, ClaSSE, key innovation, cladogenetic change — are
*instances* of it; the manual keeps them as searchable aliases, but the **structural** name
`driver:target` is primary, because it says exactly what is coupled to what.

## The diamond: four levels in three tiers

The couplings live on a diamond of four processes:

- **S** — species diversification (the birth–death process that grows the tree),
- **T** — a phenotypic trait (BM/OU/EB/Mk/threshold),
- **G** — genomes: the gene-family content (the DTL process),
- **Σ** — molecular sequences (substitution and selection).

They sit in three tiers: **S** is the substrate — the timeline; **T** and **G** are characters that
ride the species tree; and **Σ** rides the *gene* trees below them.

![The coevolution diamond. Each **directed** arrow driver $\to$ target is one model, selected with `--couple driver:target`; the two that point *into* S are drawn heavy, because an arrow into S makes the tree an output (grown jointly), while the others overlay a tree you supply. A straight **double-headed** arrow is a pair's *joint* model (both edges at once): ClaSSE, co-diversification, trait–gene feedback. Sequences (Σ) are a **target only** — a trait or gene content bends how they evolve (selection, dN/dS), but Σ drives nothing and rides the gene trees, so there is no species–sequence edge.](figures/coevolve_modes4.pdf){width=100%}

## The edges

Six edges couple species, traits and genomes, run from the `zombi2 coevolve` command:

| Edge (`coevolve --couple`) | Direction | Model |
|:------------------|:----------|:-----------------------------------------------|
| `traits:species` | T $\to$ S | state-dependent diversification (SSE / ClaSSE) |
| `species:traits` | S $\to$ T | cladogenetic trait jumps at speciation |
| `genomes:species` | G $\to$ S | key-innovation diversification |
| `species:genomes` | S $\to$ G | punctuational (cladogenetic) genome |
| `traits:genomes` | T $\to$ G | trait-linked gene families |
| `genomes:traits` | G $\to$ T | gene-conditioned trait |

: The six species/trait/gene coupling edges, each a `coevolve --couple driver:target` invocation. \label{tbl:edges}

Each of these six has a reverse, so switching **both** edges of a pair on together gives a **joint
(bidirectional) model** — three more: **ClaSSE** (traits $\leftrightarrow$ species), **co-diversification**
(species $\leftrightarrow$ genomes) and **trait–gene feedback** (traits $\leftrightarrow$ genomes). (The named
literature variants are choices *within* an edge: BiSSE, MuSSE and QuaSSE are three responses of the
single `traits:species` edge, picked with `--sse-model`.)

Three more edges point into **sequences**. Because a sequence rides its *gene* tree — downstream of the
genome layer — these live on the `zombi2 sequence` command rather than `coevolve`:

| Edge (`sequence --couple`) | Target-variable | Model |
|:------------------|:----------|:-----------------------------------------------|
| `traits:selection` | $\omega$ (dN/dS) | a trait sets each lineage's selection strength |
| `genomes:selection` | $\omega$ (dN/dS) | a gene event (e.g. duplication) relaxes selection |
| `traits:speed` | substitution rate | a trait scales the molecular clock |

Sequences are a **target only** — nothing is driven *by* a sequence — so these edges have no reverse
and no joint.

## Two rules

Two rules read any set of edges.

### Directional or bidirectional — does the tree grow?

The first: **does any active edge point into S?** If none does, the tree is fixed — read from
`-t/--tree`, with every coupling an *overlay* on it. If one does (`traits:species` or
`genomes:species`), the tree's shape depends on the coupled state and cannot be drawn first: the tree
becomes an **output**, grown jointly, and the run is forward-only (`--age`/`--tips`, no `-t`).

The asymmetry matters: an arrow pointing *out* of S (`species:traits`, `species:genomes`) does **not**
force a joint run — S drives the target but listens to nothing, so it can still be drawn first and the
target overlaid. "Touches S" is not the trigger; "points into S" is. (A bidirectional pair whose two
edges form a cycle is grown together — *fused* — for the same reason: neither side can be frozen first.)

### Adjacent tiers only

The second: **a coupling connects levels within one tier of each other.** Characters (T, G) couple to
their substrate (S) above and to sequences (Σ) below, and to each other; but there is **no
species–sequence edge** — a sequence rides a *gene* tree, not the species tree, so S and Σ are two
tiers apart. That diagonal is not a missing feature; it is ruled out by construction.

## Nulls: cut the arrow

Every coupling is a *claim* — this driver shapes that target — and the hard part is telling a real
coupling from the tree's own heterogeneity. So every edge ships a matched **decoupled null**: the same
process with the `driver → target` arrow **cut** but the target's variance kept, which is exactly the
grammar's *response set to zero*. Add `--null` to any run to generate it — `neutral` (the driver stops
setting the rates), `cid` (the variance comes from a *hidden*, uncorrelated driver — the honest
opponent), or `timing` (an at-speciation burst spread evenly along the branches). The
[null models](#coupling-that-shapes-traits-and-gene-content) section returns to these once the edges
are on the table.

## Where to next

The chapters that follow are organised by **what each coupling shapes**:

- [**Coupling that shapes the tree**](#coupling-that-shapes-the-tree) — the arrows into S that make the
  tree an output: state-dependent diversification (a trait), key innovations (gene content), and their
  joint models ClaSSE and co-diversification.
- [**Coupling that shapes traits and gene content**](#coupling-that-shapes-traits-and-gene-content) —
  the overlays on a given tree: cladogenetic and gene-conditioned traits, punctuational and trait-linked
  gene content, their feedback, and the null models in full.
- [**Coupling that shapes sequences**](#coupling-that-shapes-sequences) — the newest tier: a trait or a
  gene event bending selection (dN/dS) and the substitution rate on the gene trees.
