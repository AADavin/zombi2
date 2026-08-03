# zombi2.species

Level 1: the species tree every other level lives on. One forward engine, shaped by
[scopes and modifiers](rates.md) rather than by a model zoo.

::: zombi2.species.simulate_species_tree

::: zombi2.species.SpeciesResult

::: zombi2.species.Event

## Trees

The tree object itself, and the readers and shape helpers that work on it. These live in
**`zombi2.tree`** — `from zombi2.tree import read_newick` — and are documented here because the
species level is the one that grows a tree; every other level takes one. `read_newick` reads a ZOMBI2
tree or an external one, so a genome run can start from a published phylogeny.

::: zombi2.tree.Tree

::: zombi2.tree.Node

::: zombi2.tree.prune

::: zombi2.tree.read_newick

::: zombi2.tree.as_tree
