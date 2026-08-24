# API reference

The Python API has **one canonical path per name**, reached through each level's package —
there are no top-level re-exports. A run always starts from a `simulate_*` entry point and
returns a `*Result` object that carries the true history behind the dataset.

```python
from zombi2 import species, genomes
from zombi2.params import Global, PerCopy, PerLineage

result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
```

The reference is generated from the source docstrings, one page per level, in the same order as
the [user guide](../guide/introduction.md).

## The entry points

Every run starts here. The genome level has three **resolutions** — family ⊂ ordered ⊂
nucleotide — with one entry point each, so the resolution is chosen by which function you call.

| Level | Entry point | Returns | Guide |
|---|---|---|---|
| Species | [`simulate_species_tree`][zombi2.species.simulate_species_tree] | [`SpeciesResult`][zombi2.species.SpeciesResult] | [Species trees](../guide/species-trees.md) |
| Genomes · family | [`simulate_genomes_family`][zombi2.genomes.simulate_genomes_family] | [`FamilyGenomesResult`][zombi2.genomes.FamilyGenomesResult] | [Genomes I](../guide/genomes.md) |
| Genomes · ordered | [`simulate_genomes_ordered`][zombi2.genomes.simulate_genomes_ordered] | [`OrderedGenomesResult`][zombi2.genomes.OrderedGenomesResult] | [Genomes II](../guide/genomes-ordered.md) |
| Genomes · nucleotide | [`simulate_genomes_nucleotide`][zombi2.genomes.simulate_genomes_nucleotide] | [`NucleotideGenomesResult`][zombi2.genomes.NucleotideGenomesResult] | [Genomes III](../guide/genomes-nucleotide.md) |
| Sequences | [`simulate_sequences`][zombi2.sequences.simulate_sequences] | [`SequencesResult`][zombi2.sequences.SequencesResult] | [Sequences](../guide/sequences.md) |
| Traits · continuous | [`simulate_continuous`][zombi2.traits.simulate_continuous] | [`TraitsResult`][zombi2.traits.TraitsResult] | [Traits](../guide/traits.md) |
| Traits · discrete | [`simulate_discrete`][zombi2.traits.simulate_discrete] | [`TraitsResult`][zombi2.traits.TraitsResult] | [Traits](../guide/traits.md) |
| Traits · several at once | [`simulate_traits`][zombi2.traits.simulate_traits] | one `TraitsResult` per name | [Traits](../guide/traits.md) |
| Two levels at once | [`joint.simulate`][zombi2.joint.simulate] | [`JointResult`][zombi2.joint.JointResult] | [Joint runs](../guide/joining.md) |

Every result writes its outputs with `.write(directory)`; which files that leaves is catalogued
in [output files](../output-files.md).

## The supporting pieces

These are shared across levels rather than owned by one.

| Page | What it holds |
|---|---|
| [`zombi2.params`](rates.md) | the rate grammar — **scopes** (`PerCopy`, `PerLineage`, `Global`, …), the **verbs** chained onto them (`scaled_by`, `set_by`, `weighted_by`, and the two shortcuts `varying_among`, `changing_at`), the **drivers** they take (`Random`, `TotalDiversity`, `Time`, `Clade`, a filename) and the **mappings** the link carries (`Table`, `Curve`, `Scalar`, `Between`), the notation shared by the Python API, the CLI and a `--params` file |
| [`zombi2.tree`](species.md#trees) | trees — `Tree`, `Node`, `prune`, `read_newick`, and the tree-shape helpers. Its own module, documented on the species page because that is the level that grows one |

## The menus inside a level

Two levels offer a fixed set of choices. Each is owned by that level and documented as a section of
its page, not as a module of its own.

| Section | What it holds |
|---|---|
| [Who receives a transfer](genomes.md#who-receives-a-transfer) | the transfer-recipient rules — `Distance`, `Clades`, `Between` |
| [The substitution-model menu](sequences.md#the-substitution-model-menu) | the substitution models — `jc69` · `k80` · `hky85` · `gtr`, and the protein models `poisson` · `jtt` · `dayhoff` · `wag` · `lg` |

## Reading a run back

A run written to disk is read back with the `zombi2 tools` commands rather than through this API —
see [tools](../tools.md).
