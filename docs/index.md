# ZOMBI2

**ZOMBI2** simulates the evolution of **species trees**, **genomes**, **sequences** and **traits** —
each on its own, conditioned on another, or jointly — and records the true history behind every
dataset. It is a ground-up redesign of [ZOMBI](https://github.com/AADavin/Zombi).

!!! note "Rebuild in progress"
    ZOMBI2 is being rebuilt as a **clean core, grown level by level from a single specification**.
    **Species trees** and **unordered genomes** are available now; sequences, traits, and the
    couplings between levels are being rebuilt and return here as they land. This guide grows a page
    per level.

## The four levels

Three levels form a chain, and traits branch off it:

```
Species → Genomes → Sequences     a genome lives on the species tree; a sequence lives inside a gene
Species → Traits                  a trait lives on the species tree
```

A genome, sequence, or trait always evolves **along a species tree**. You simulate each level on its
own, or let one **drive** another.

## Quickstart

Grow a species tree, then evolve gene families along it:

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_unordered

# a birth–death tree of 20 surviving species
sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)

# gene families along that tree — duplication, transfer, loss, origination
g = simulate_genomes_unordered(sp, duplication=0.2, transfer=0.1, loss=0.25,
                               origination=0.5, initial_families=20, seed=42)

# the genomes you observe are the extant tips
observed = {n.id: g.genomes[n.id] for n in sp.complete_tree.extant()}
```

Every rate is written the same way — a **scope** around a base, optionally times **modifiers**
(a rate that changes in time, saturates with diversity, or drifts along the tree):

```python
from zombi2.rates import scope, modifiers

sp = species.simulate_species_tree(
    birth = 1.0 * modifiers.Time({0: 1.0, 3: 0.5}),   # skyline: full rate, then half after time 3
    death = scope.Global(0.3),                        # one tree-wide death rate, not per lineage
    total_time = 8.0, seed = 1)
```

## Where next

- The [**Species trees**](guide/species-trees.md) guide — the birth–death process, the rate modifiers,
  sampling and fossils, and the objects a run returns.
- More level guides appear here as each level is rebuilt.
