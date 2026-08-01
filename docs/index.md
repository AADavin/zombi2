# ZOMBI2

**[Manual (pdf)](https://aadavin.github.io/zombi2/zombi2-manual.pdf)**

**ZOMBI2** simulates the evolution of **species trees**, **genomes**, **sequences** and **traits** —
each on its own, conditioned on another, or jointly — and records the true history behind every
dataset.

## The four levels

Species, genomes and sequences form a chain of ancestry — a genome lives on the species tree, a
sequence inside a gene. Traits are separate: a trait attaches directly to the species tree, so it
does not depend on genomes or sequences.

<figure markdown="span">
  ![The four levels of ZOMBI2](img/fig-2-1-four-levels.svg){ width="330" }
</figure>

A genome, sequence, or trait always evolves **along a species tree**. You simulate each level on its
own, or let one **drive** another.

## Quickstart

Grow a species tree, then evolve gene families along it:

```python
from zombi2 import species, genomes

# a birth–death tree of 20 surviving species
sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)

# gene families along that tree — duplication, transfer, loss, origination
g = genomes.simulate_genomes_family(sp, duplication=0.2, transfer=0.1, loss=0.25,
                                    origination=0.5, initial_families=20, seed=42)

# the genomes you observe are the extant tips
observed = {n.id: g.genomes[n.id] for n in sp.complete_tree.extant_leaves()}
```

ZOMBI2 lets you describe very specific evolutionary scenarios, and rates are how you do it.
Every rate is written the same way — a **scope** around a base, optionally times **modifiers**
(a rate that changes in time, saturates with diversity, or drifts along the tree):

```python
from zombi2.rates import scope, modifiers

sp = species.simulate_species_tree(
    birth = 1.0 * modifiers.OnTime({0: 1.0, 3: 0.5}),   # skyline: full rate, then half after time 3
    death = scope.Global(0.3),                        # one tree-wide death rate, not per lineage
    n_extant = 20, seed = 1)                           # grow until 20 species survive
```

## Where next

- New here? Start with the [**Introduction**](guide/introduction.md) and [**a tour of ZOMBI2**](guide/tour.md).
- The level guides — each grows one level and the objects a run returns:
  [**Species trees**](guide/species-trees.md),
  [**Genomes**](guide/genomes.md) (with [ordered](guide/genomes-ordered.md) and
  [nucleotide](guide/genomes-nucleotide.md) resolutions),
  [**Sequence evolution**](guide/sequences.md), and [**Trait evolution**](guide/traits.md).
- [**Conditioning**](guide/conditioning.md) and [**joining**](guide/joining.md) — letting one thing drive another.
- Reference: the [**API**](api/index.md), the [**output files**](output-files.md) each run writes, and the
  [**tools**](tools.md) that read a finished run.
