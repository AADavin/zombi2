```{=latex}
\part{Coevolution}
```

# Coevolution: the framework

ZOMBI2 normally simulates in a **pipeline** — a species tree, then a trait along it, then gene
families along it — because the couplings run one way and the joint distribution factorises. This
part is about what to do when they do *not*: when a trait or a genome feeds **back** into the
process that produced the tree, and the levels must be grown **together**. This chapter sets up the
one abstraction that covers every such case — a directed graph over the **four levels** — and the
chapters that follow work through the models.

## The pipeline, and when it breaks

The hierarchical pipeline works because the joint distribution factorises,

$$P(\text{tree}) \cdot P(\text{trait} \mid \text{tree}) \cdot P(\text{genes} \mid \text{tree}, \text{trait}),$$

so each stage can be drawn on the frozen output of the previous one. The moment a trait or a genome
feeds **back** into the process that generated the tree, this factorisation breaks: the tree, the
trait and the gene content must be grown *together*. The **`coevolve`** mode is for these coupled
scenarios.

## Four processes, directed edges

Coupling runs over four processes — the **diamond** of ZOMBI2's four levels:

- **S** — species diversification (the birth–death process that grows the tree),
- **T** — a phenotypic trait (BM/OU/EB/Mk/threshold),
- **G** — gene-family content (the DTL process),
- **Σ** — molecular sequences (substitution and selection, dN/dS).

They sit in three tiers: **S** is the substrate — the timeline; **T** and **G** are characters that
ride the species tree; and **Σ** rides the *gene* trees below them. A **coupling is a directed edge**
`driver -> target`: the driver's state modulates the target's rates. A coevolution scenario is a
**set of directed edges** on these nodes, and it makes the one thing that matters — *direction* —
explicit. Sequences are a **target only**: a trait or gene content can shape how a sequence evolves,
but a sequence drives nothing, so Σ has arrows coming in and none going out.

On the command line the coupling is one repeatable flag, `--couple driver:target`, where the order
reads as the arrow (driver first):

```bash
zombi2 coevolve --couple traits:species ...   # T->S: trait sets speciation (SSE)
zombi2 coevolve --couple species:traits ...   # S->T: speciation drives the trait
# both arrows at once = ClaSSE:
zombi2 coevolve --couple traits:species --couple species:traits ...
```

So `--couple species:traits` and `--couple traits:species` are deliberately **different models**, and
a bidirectional coupling is simply *both* edges. The `:` (rather than `->`) keeps the flag shell-safe.

## The edges and the joint models

The four processes give eight directed edges (Table \ref{tbl:edges}), each a distinct model — six
among species, traits and genes, and two more that point into sequences:

| Edge (`--couple`) | Direction | Model |
|:------------------|:----------|:-----------------------------------------------|
| `traits:species` | T $\to$ S | state-dependent diversification (SSE / ClaSSE) |
| `species:traits` | S $\to$ T | cladogenetic trait jumps at speciation |
| `genomes:species` | G $\to$ S | key-innovation diversification |
| `species:genomes` | S $\to$ G | punctuational (cladogenetic) genome |
| `traits:genomes` | T $\to$ G | trait-linked gene families |
| `genomes:traits` | G $\to$ T | gene-conditioned trait |
| `traits:sequences` | T $\to$ Σ | trait-driven selection (dN/dS) and substitution speed |
| `genomes:sequences` | G $\to$ Σ | post-duplication relaxed selection |

: The directed coupling edges — the `--couple driver:target` flag, the direction of the arrow, and the model each one selects. \label{tbl:edges}

The six S/T/G edges each have a reverse, so each node-pair's two edges can be switched on **together**,
giving a **joint (bidirectional) model** — three more: **ClaSSE** (traits $\leftrightarrow$ species),
**co-diversification** (species $\leftrightarrow$ genes) and **trait–gene feedback**
(traits $\leftrightarrow$ genes). The sequence edges have no reverse (Σ is a target only), so they
have no joint. (The named literature models are choices *within* an edge: BiSSE, MuSSE and QuaSSE are
three flavours of the single `traits:species` edge, picked with `--sse-model`.) One diagonal is
**forbidden**: there is no species–sequence edge, because a sequence rides its *gene* tree, not the
species tree.

![The coevolution diamond. Each **directed** arrow driver $\to$ target is one model, selected with `--couple driver:target`; the two that point *into* S are drawn heavy, because an arrow into S makes the tree an output (grown jointly), while the others overlay a tree you supply. A straight **double-headed** arrow is a pair's *joint* model (both edges at once): ClaSSE, co-diversification, trait–gene feedback. Sequences (Σ) are a **target only** — a trait or gene content bends how they evolve (selection, dN/dS), but Σ drives nothing and rides the gene trees, so there is no species–sequence edge.](figures/coevolve_modes4.pdf){width=100%}

## The one rule: does an edge point into S?

The one rule that governs their difficulty is: **does any active edge point into S?** If no edge
points into S, the tree is fixed — it is read from `-t/--tree` and every coupling is an *overlay* on a
frozen tree. If an edge does point into S (`traits:species` or `genomes:species`), the tree topology
depends on the coupled state and cannot be drawn first: the tree becomes an **output**, and those runs
are forward-only and take no `-t`.

Note the asymmetry: an arrow pointing *out* of S (`species:traits`, `species:genomes`) does **not**
trigger joint simulation. S drives the target but listens to nothing, so S can be drawn (or supplied
via `-t`) first and the target overlaid on it — the tree stays an input. It is only an arrow *into* S
that puts S downstream of its driver, breaks the pipeline factorisation, and forces the tree to be
grown jointly as an output. So "touches S" is not the trigger; "points into S" is.

## The map, and where to next

The chapters that follow take the couplings in turn.
[State-dependent diversification](#state-dependent-diversification) covers the species tree coupled
to a **trait** — the SSE family (the arrows that shape the tree), the reverse cladogenetic edge, and
their joint model ClaSSE. [Coupling gene content](#coupling-gene-content) then covers everything
involving **gene content** — genes driving diversification and being reshuffled at speciation, genes
and traits conditioning one another — and closes with the null models that let you tell a real
coupling from the tree's own heterogeneity. The two **sequence** edges (`traits:sequences`,
`genomes:sequences`) — trait- and gene-driven selection on molecular sequences — are the newest
additions to the diamond, and join the manual as the sequence-coupling code lands.

::: note
The one coevolution model still on the roadmap is the fully joint **`--all`** run: every edge active
at once, so all three pairs are bidirectional and the trait, the genome and the tree feed back on one
another with no single imposed direction (forward time resolves the mutual dependence). It *composes*
the existing edges rather than adding new science — think of it as `--couple` for all six arrows —
and is best treated as a stress-test showcase rather than a routine analysis mode.
:::
