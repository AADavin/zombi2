# Genomes I — Unordered

A genome, at its simplest, is a bag of genes: which gene families a lineage carries, and in how many copies. That is the **unordered** genome — a multiset of families, with no position and no neighbours. It is the resolution you want whenever families evolve independently of one another, which is exactly what a phylogenetic profile records, and it is where most gene-content studies live.

## The four events

Gene families change along the tree by four events:

- **Origination** creates a new gene family, one that was not there before.
- **Duplication** copies a gene within a genome, so a family gains a copy.
- **Transfer** copies a gene from one lineage into another one alive at the same time.
- **Loss** removes a copy.

Origination and transfer bring genes in, duplication multiplies them, and loss takes them away. A fifth, optional event, gene conversion, has its own section below; it is off unless you ask for it.

```python
from zombi2 import genomes
result = genomes.simulate_unordered(
    tree, duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,
    initial_families=20, seed=1)
```

## Rates

Each event has a rate, and there are two different things you can tune about it. Keeping them apart is the whole point.

**The scope: per what.** Duplication, transfer and loss are counted **per copy**: a family with ten copies loses genes ten times as fast as a family with one, so families grow and shrink multiplicatively. Origination is different, counted **per lineage**: a lineage seeds new families regardless of how many it already has. These are the sensible defaults, so you usually write nothing. To change one, wrap it:

```python
loss = scope.PerLineage(0.25)     # count loss per lineage, not per copy
origination = scope.Global(0.5)   # one constant total origination rate over the whole tree
```

**Per-family variation: shared or its own.** A bare number is a rate shared by every family. A `ByFamily` modifier makes each family draw its own, independently for each event, so a family that loses fast is not automatically duplicating fast:

```python
loss = 0.25 * mod.ByFamily(spread=0.5)   # each family gets its own loss rate around 0.25
```

If instead you want a family uniformly fast or slow at everything, that is a single per-family **speed** that scales all of its rates together — one factor per family, given as its own argument:

```python
family_speed = mod.Speed(spread=0.5)   # each family one speed, scaling all of its rates at once
```

And when you need exact control, you can hand specific families their own rates:

```python
families = [dict(duplication=0.5, transfer=0.8, loss=0.3), …]   # this family, exactly
```

The two never collide: the scope is a *wrapper* (per what), the variation is a *modifier* (how much each family differs).

## Transfers

Transfer is the one event with real mechanics, because a transferred gene has to land somewhere. When a transfer fires, a recipient is chosen among the lineages alive at that moment:

```python
transfer_to = "uniform"    # any contemporaneous lineage, equally likely
transfer_to = "distance"   # closer relatives more likely — the realistic case
```

and the gene either adds a new copy or overwrites one already there:

```python
replacement = False   # additive (default): the recipient gains a copy
replacement = True    # the gene replaces an existing copy in the recipient
```

By default a lineage does not transfer to itself (`self_transfer=False`). Finer control, where some lineages donate or receive more than others, is a **donor weight** (which is just a modifier on the transfer rate) and a **recipient weight** (a bias in the recipient rule).

## Gene conversion

*[Draft — the concept is settled; the conversion-mechanics API is still to be designed.]*

Gene conversion is a fifth event, off by default. Where transfer moves a gene *between* lineages, conversion acts *within* one: one copy of a family overwrites another copy of the **same** family, so the two become identical. It is how gene families undergo concerted evolution, and you turn it on with a conversion rate.

## Bounding growth

Under duplication a family can grow without limit, which is rarely what you want. Two caps hold it in check: a hard **maximum copy number** per family, and a soft **carrying capacity** that slows duplication as a family fills up, the family-level analogue of the diversity cap on the species tree.

## Gene trees and the event log

Everything above is recorded as an **event log**: for each family, the full sequence of originations, duplications, transfers, losses and speciations that produced it. The log is the genealogy, and everything else is read off it. From it ZOMBI2 reconstructs each family's **gene tree** and its **reconciliation**, the mapping of that gene tree into the species tree that labels every node with the event which made it. And the **profile matrix**, families by extant species, holds the copy numbers, which is the phylogenetic profile itself.

## The genome and profile objects

`simulate_unordered` returns a **`GenomesResult`** bundle. Every level returns the same shape of object — a `<Level>Result` sharing a common spine — so what you learn here carries over to sequences and traits. The spine is `.events` (the event log), `.tree` (the species tree it ran on), `.seed`, and `.write(dir, include=[...])` to materialise the outputs you choose to disk. The genome payload adds `.genomes` (per-node gene content), `.gene_trees` (each family's gene tree with its reconciliation), `.profiles` (the sparse families × species matrix), and `.transfers`. From it you can read a lineage's gene content, pull a single family's gene tree, or hand the whole thing to the sequence level.

The event log is the compact source of truth, and the rich views are read off it: `.gene_trees`, `.profiles` and the ancestral genomes are reconstructed lazily on access rather than all held in memory at once. When you only need one view at scale, declare it up front with `record=[...]` — for example `record=["profiles"]` makes the run a pure profile accumulator that never builds the gene-tree objects at all, and its footprint collapses to the sparse matrix.

## Usage from Python

```python
from zombi2 import genomes, modifiers as mod
from zombi2 import scope             # scope wrappers: Global, PerCopy, PerLineage, …

# the common case
genomes.simulate_unordered(tree, duplication=0.2, transfer=0.1, loss=0.25,
                           origination=0.5, initial_families=20, seed=1)

# per-family variation in loss, distance-weighted transfer, global origination
genomes.simulate_unordered(tree, duplication=0.2, loss=0.25 * mod.ByFamily(spread=0.5),
                           transfer=0.1, transfer_to="distance",
                           origination=scope.Global(0.5), seed=1)
```

## Usage from the CLI

*[Draft — the CLI re-fit to this API is still to be designed.]*

```bash
zombi2 genomes --unordered --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 \
    --initial-families 20 -t species_tree.nwk --seed 1 -o my_genomes
```

## Outputs

*[Draft — to finalise with Appendix B.]*

A run writes the gene trees, their reconciliations, the event log behind them, and the profile matrix (families × species). The full list of files lives in Appendix B.
