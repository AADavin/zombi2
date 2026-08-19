# FAQ

- [I ran `zombi2 genomes` and there is no FASTA. Where is it?](#i-ran-zombi2-genomes-and-there-is-no-fasta-where-is-it)
- [What do `n5`, `e14`, `g203` mean?](#what-do-n5-e14-g203-mean)
- [Why do the family (or file) numbers skip — `fam0, fam1, fam2, fam3, fam5, …`?](#why-do-the-family-or-file-numbers-skip-fam0-fam1-fam2-fam3-fam5)
- [Why are there two of every tree — `_complete` and `_extant`?](#why-are-there-two-of-every-tree-_complete-and-_extant)
- [How do I actually look at a tree? The `.nwk` files are just text.](#how-do-i-actually-look-at-a-tree-the-nwk-files-are-just-text)
- [How do I read a ZOMBI2 tree in ete3 or Biopython?](#how-do-i-read-a-zombi2-tree-in-ete3-or-biopython)
- [My genomes came out empty. Is that a bug?](#my-genomes-came-out-empty-is-that-a-bug)
- [Why do my families stop at 10 copies?](#why-do-my-families-stop-at-10-copies)
- [I counted losses in the event log and got fewer than I expected.](#i-counted-losses-in-the-event-log-and-got-fewer-than-i-expected)
- [Why does re-running a level refuse?](#why-does-re-running-a-level-refuse)
- [`sequences` printed a warning about "saturated" alignments — did it fail?](#sequences-printed-a-warning-about-saturated-alignments-did-it-fail)

## I ran `zombi2 genomes` and there is no FASTA. Where is it?

There isn't one, and that is the level working as intended. **Simulating a genome, in ZOMBI2, means
simulating its gene content and the history behind it** — which families exist in which lineage, how
many copies of each, and the duplication, transfer, loss and origination events that got them there.
A gene is an identity, not a string of bases. So the level writes `profiles.tsv` (families × species
copy counts), `genomes.tsv` (every node's genes, one row per copy), `genome_events.tsv` (every event)
and a gene tree per family. No sequence anywhere.

Sequences are the **next** level, and they are separate because they are a separate model: gene
content evolves by gaining and losing copies, sequence evolves by substitution along the gene trees
the genome level just produced. Run it on the genome output and you get `alignments/fam<f>.fasta`,
one alignment per family:

```bash
zombi2 sequences out --model jc69
```

If what you wanted was **whole genomes as FASTA** — a lineage's chromosomes, bases and all — that
needs `--resolution nucleotide` at the genome level, which lays genes out on a chromosome with
coordinates so there is something to assemble:

```bash
zombi2 genomes out --resolution nucleotide --root-length 3000 --genes 6
zombi2 sequences out --model jc69
```

The sequence level then writes `genomes/genome_<lineage>.fasta`, one file per node of the complete
tree, plus `genomes/genome_initial.fasta` for the genome the run started with, which is no node's.
The default resolution, `family`, has no coordinates and no order, so there is nothing to
concatenate.

Appendix B lists every file each level writes.

## What do `n5`, `e14`, `g203` mean?

ZOMBI2 gives every node and gene copy a short id:

- **`n<number>`** — a species-tree node: an extant tip, or an internal ancestor.
- **`e<number>`** — a lineage that went **extinct** (a tip of the complete tree that died before the present).
- **`g<number>`** — a **gene copy**. Gene-tree tips read `n<node>_g<copy>` — which genome, which copy.

The ids are ZOMBI2's own. If you supplied your own tree with real tip names, `names.tsv` maps each
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

## How do I read a ZOMBI2 tree in ete3 or Biopython?

ZOMBI2 labels **internal** nodes, not only the tips: a species tree carries each ancestor's `n<id>`,
and a gene tree carries the event that made the node — `speciation_n6`, `transfer_n2`. That is
ordinary Newick, but it is not what ete3 reads by default, so pass `format=1`:

<!-- doc-test: skip — ete3 is not a ZOMBI2 dependency -->
```python
from ete3 import Tree
t = Tree("out/genomes/gene_trees/gene_tree_fam0_extant.nwk", format=1)   # format=1: labelled internal nodes
```

Without it, ete3's default (`format=0`) expects an internal label to be a support value, and quotes
the one it could not read as a number: `NewickError: Unexpected newick format
'transfer_n20:0.2438901'`. This holds for every tree a run writes — the species trees, the gene
trees, the phylograms and the trait tree. With `format=1` the label is kept, so `node.name` on an
internal node of a gene tree is the event.

Biopython needs no flag and puts the same labels in `.name`:

<!-- doc-test: skip — Biopython is not a ZOMBI2 dependency -->
```python
from Bio import Phylo
t = Phylo.read("out/species/species_extant.nwk", "newick")
[c.name for c in t.get_nonterminals()]        # each ancestor's id — 'n0', 'n15', …
```

One exception: `trait_tree.nwk` annotates each node with its value, as `n19[&trait=-0.38831]`.
Biopython parses that into `.comment` and leaves `.name` as `n19`; ete3 with `format=1` reads the
whole thing as the name. For the numbers themselves, read `trait_values.tsv` and join on its `node`
column.

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
