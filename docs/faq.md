# FAQ

Short answers to the things that surprise newcomers first — the moments that look like a bug and
aren't. For everything else, see the guide and the reference.

## What do `n5`, `e14`, `g203` mean?

ZOMBI gives every node and gene copy a short id:

- **`n<number>`** — a species-tree node: an extant tip, or an internal ancestor.
- **`e<number>`** — a lineage that went **extinct** (a tip of the complete tree that died before the present).
- **`g<number>`** — a **gene copy**. Gene-tree tips read `n<node>_g<copy>` — which genome, which copy.

The ids are ZOMBI's own. If you supplied your own tree with real tip names, `names.tsv` maps each
`n<id>` back to your label.

## Why do the family (or file) numbers skip — `fam0, fam1, fam2, fam3, fam5, …`?

The gaps are families that **died out** before the present. `profiles.tsv` and `alignments/` list the
survivors, so their ids have holes. Nothing is lost, though: a family that died out still has its
complete gene tree (`gene_trees/gene_tree_fam<K>_complete.nwk`, written for every family ever born) and
its rows in `genome_events.tsv`, so you can see exactly what it did before dying. The report counts
them: `families N born · M surviving · K died out`.

## Why are there two of every tree — `_complete` and `_extant`?

- **`_complete`** is the whole tree — every node, including lineages that went extinct.
- **`_extant`** is that tree pruned to the sampled tips surviving to the present — the one most
  analyses want.

A gene tree is the same idea: `_complete` has every copy ever, `_extant` only the surviving copies.

## How do I actually look at a tree? The `.nwk` files are just text.

They are **Newick**, the standard tree format. Plot them with **Phylustrator**, ZOMBI2's companion
command-line tool for figures: `pip install phylustrator`, then its `phyl` command draws the trees —
and the genomes — a run produces, straight from the files above.

## `sequences` printed a warning about "saturated" alignments — did it fail?

No — it succeeded; that is a *warning*, not an error (the exit status is still 0). It means the default
substitution rate diverged the sequences so far that the alignments keep little history. Say how
diverged you want them instead: **`--divergence 0.2`** gives a readable alignment on any tree.
