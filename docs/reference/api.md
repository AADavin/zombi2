# API reference

The Python API has **one canonical path per name**, reached through each level's package —
there are no top-level re-exports. A run always starts from a `simulate_*` entry point and
returns a `*Result` object that carries the true history behind the dataset.

```python
from zombi2 import species
result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)

from zombi2.genomes import simulate_genomes_unordered
from zombi2.rates import scope, modifiers
```

The reference below is generated from the source docstrings, level by level, in the same
order as the [user guide](../guide/introduction.md).

## Species trees

::: zombi2.species.simulate_species_tree

::: zombi2.species.SpeciesResult

::: zombi2.species.Tree

::: zombi2.species.Node

::: zombi2.species.Event

::: zombi2.species.prune

## Genomes

The genome level has three resolutions — unordered ⊂ ordered ⊂ nucleotide — one entry
point each.

::: zombi2.genomes.simulate_genomes_unordered

::: zombi2.genomes.simulate_genomes_ordered

::: zombi2.genomes.simulate_genomes_nucleotide

::: zombi2.genomes.GenomesResult

::: zombi2.genomes.OrderedGenomesResult

::: zombi2.genomes.NucleotideGenomesResult

::: zombi2.genomes.GeneTree

::: zombi2.genomes.GeneCopy

## Sequences

::: zombi2.sequences.simulate_sequences

::: zombi2.sequences.SequencesResult

## Traits

::: zombi2.traits.simulate_continuous

::: zombi2.traits.simulate_discrete

::: zombi2.traits.TraitsResult

::: zombi2.traits.Change

## Conditioning and joining

::: zombi2.joint.simulate_joint

::: zombi2.joint.JointResult

## Rates

Rates are written as expressions over **scopes** and **modifiers**; the same notation is
shared by the Python API, the CLI and a `--params` file.

::: zombi2.rates.scope

::: zombi2.rates.modifiers
