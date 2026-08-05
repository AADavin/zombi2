# FAQ

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
plotting library: `pip install phylustrator`, then `phyl species_extant.nwk` draws any tree a run
produces. Its Python API draws the genomes too.

## My genomes came out empty. Is that a bug?

No. Loss is counted per copy and the last copy is a copy like any other, so a loss rate well above
duplication and origination strips a lineage of everything it has. The run tells you:
`genome_summary.json` reports `empty_genomes` and `zombi2 genomes` prints a line on standard error.
Lower `--loss`, raise `--origination`, or start with more `--initial-families`.

The chromosome-based resolutions do have a floor — a loss never takes a chromosome below its last
gene — but that is a floor on a chromosome, not on a genome.

## Why do my families stop at 10 copies?

`--max-family-size` defaults to 10 copies of one family **in one genome**, because a duplication rate
above the loss rate grows without bound. While the cap binds it discards events, so the realised
duplication and transfer rates fall below the ones you declared. Pass `--max-family-size none` when
you are measuring rates.

## I counted losses in the event log and got fewer than I expected.

A `transfer_replacing` row overwrites a copy on the recipient branch and writes **no separate `loss`
row** for it, so a script counting `loss` rows comes up short by exactly the number of replacing
transfers. `genome_summary.json` counts the biology instead: `transfer` as one number over both
kinds, each displaced copy under `loss`.

## Why does re-running a level refuse?

Because a later level was built from it, and re-running in place would leave that later output
mismatched. `--force` re-runs anyway and removes the now-stale downstream. The same guard covers
conditioning: a conditioned run writes a `conditioned_on` file naming the levels it reads, through a rate or `--transfer-to`, so
re-running the driver afterwards refuses rather than leaving the target silently stale. It works
**within one run directory** — a driver and a target written to two directories with `--from` are
not linked.

## `sequences` printed a warning about "saturated" alignments — did it fail?

No — it succeeded; that is a *warning*, not an error (the exit status is still 0). It compares the mean
identity of the alignments with the identity unrelated sequences already have under that model — 25%
for DNA with equal base frequencies. Close to that floor, the alignments have kept little history. The
substitution rate is per unit time, so at the default a tall tree saturates. Say how diverged you want
the sequences instead: **`--divergence 0.2`**.
