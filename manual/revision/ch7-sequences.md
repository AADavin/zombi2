# Sequence evolution

The sequence level does two main things:

* It rescales the gene trees and the species tree from time into substitutions per site (**phylograms**).
* It evolves the nucleotides that sit inside every gene, so each family ends with an alignment.

The sequence level is always dependent on a genome-level run, and it takes that run's result directly:

```python
from zombi2 import genomes, sequences
from zombi2.sequences.substitution_models import hky85

my_genomes = genomes.simulate_genomes_unordered(tree, duplication=0.2, transfer=0.1,
                                                loss=0.25, origination=0.5, seed=1)
result = sequences.simulate_sequences(my_genomes, model=hky85(kappa=2.0),
                                      length=1000, seed=1)
```

Handing over the whole `GenomesResult` is the form to prefer. The sequence level needs each family's gene tree, which it takes from the result itself, but it also needs the *species* tree behind them if it is to report which lineages ran fast and which ran slow — and a bare collection of gene trees does not carry one. You can still pass a plain `{family: GeneTree}` mapping when that is all you have; you simply get no species phylogram out of it.

## Creating phylograms

A gene tree arrives as a **chronogram**, its branch lengths in time. What a sequence actually accumulates along a branch is not time but a number of *substitutions per site*, and that is time multiplied by an evolutionary rate. Turning one into the other is the whole job of the sequence level. Applying a rate to every branch rescales the tree from time into expected substitutions and yields a phylogram.

Two things therefore have to be chosen: *what* changes (the substitution model, the chemistry of which residue turns into which) and *how fast* it changes along each branch, which is the clock.

## The substitution models 

ZOMBI2 implements different standard models of sequence evolution:

```python
model = jc69()                    # equal rates, equal base frequencies — no free parameters
model = k80(kappa=2.0)            # a transition/transversion bias
model = hky85(kappa=2.0, freqs=(0.3, 0.2, 0.2, 0.3))            # bias + unequal frequencies
model = gtr(rates=(1,2,1,1,2,1), freqs=(0.25,0.25,0.25,0.25))   # six exchangeabilities + freqs
```

These are nucleotide models, so `length` counts sites and each site holds one of the four bases. They are a genuine menu rather than one model wearing four sets of parameters: each is a different rate matrix, a different claim about which substitutions are easy. They do nest, though, in the order written — `jc69` is `k80` with `kappa=1`, `k80` is `hky85` with equal frequencies — so moving down the list only ever adds free parameters.

## Relaxed molecular clocks

The rate itself is `substitution`, and it is counted **per site**: a gene-tree branch of Δ*t* time accrues `substitution · Δt` substitutions at every site. Leave it alone and it is `1.0` everywhere — the **strict clock**, one tempo for the whole tree.

Real lineages do not oblige, and a rate that varies from lineage to lineage is what the field calls a **relaxed clock**. In ZOMBI2 that is not a new kind of object but the ordinary rate grammar: multiply the rate by a modifier.

```python
# strict clock — one rate everywhere; the default, so write nothing
substitution = 1.0

# relaxed — each lineage draws its own rate, independently of its neighbours
substitution = 1.0 * mod.ByLineage(spread=0.3)                 # lognormal (the default)
substitution = 1.0 * mod.ByLineage(spread=0.3, dist="gamma")   # or gamma
```

**`ByLineage`** has *no memory*: each lineage is an independent draw, so a lineage's rate tells you nothing about its neighbours'. The distribution it draws from (`dist="lognormal"` or `"gamma"`) is a parameter of the modifier, not a modifier of its own.

One detail is worth stating, because it decides what the clock means. A clock is a property of a **lineage**: a whole species runs hot or cold, and every gene passing through that branch feels it. So the draw is made once per *species* branch and shared by all the families that pass through it, and each gene-tree branch reads the clock of the species branch it is reconciled to. That reconciliation is already known from the genome run, so nothing has to be wired up — but it is why the species tree has to come along, and why a lineage running hot shows up across all of its genes at once rather than in one family alone.

A reference table that can be handy to people who want to implement a specific model from the literature:

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| one rate everywhere | `substitution = 1.0` (default) | Strict / global clock |
| each lineage i.i.d. lognormal | `1.0 * mod.ByLineage(spread=…)` | Uncorrelated lognormal (UCLN) |
| each lineage i.i.d. gamma | `1.0 * mod.ByLineage(spread=…, dist="gamma")` | Uncorrelated gamma (UGAM) |
| each lineage i.i.d. — that *is* white-noise | `1.0 * mod.ByLineage(spread=…)` | White-noise clock |

Three of those four names are the same modifier. That is the table earning its keep: a reader who arrives wanting "a UCLN clock" or "a white-noise clock" finds that both are `mod.ByLineage`, differing in a parameter or in nothing at all.

## The objects

`simulate_sequences` returns a **`SequencesResult`**, which carries:

- `.alignments` — the observable data: for each family, the sequence at every **extant** gene copy. This is the alignment a phylogenetic method would be handed.
- `.ancestral` — the sequence at **every** node, internal and extinct alike. Because the run wrote a sequence at each node as it went, ancestral reconstruction is a byproduct rather than a separate step, and these are the exact ancestors, not estimates.
- `.phylograms` — for each family, its gene tree with branch lengths converted from time into substitutions per site: the tree the sequences were actually drawn along.
- `.species_phylogram` — the same conversion applied to the species tree, so the clock is visible as branch lengths. It is `None` when the run was given bare gene trees instead of a `GenomesResult`.

As with every level, the bundle also carries `.seed` and `.write(directory, outputs=[...])` to put the chosen outputs on disk.

## Usage from Python

An end-to-end run, from a species tree through genomes to alignments:

```python
from zombi2 import species, genomes, sequences
from zombi2.rates import modifiers as mod
from zombi2.sequences.substitution_models import hky85, gtr

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1).complete_tree
my_genomes = genomes.simulate_genomes_unordered(tree, duplication=0.2, transfer=0.1,
                                                loss=0.25, origination=0.5, seed=1)

# the common case: DNA under HKY85, a strict clock
result = sequences.simulate_sequences(my_genomes, model=hky85(kappa=2.0),
                                      length=1000, seed=1)
result.alignments          # {family: {gene copy: sequence}} — the observable data
result.ancestral           # the same, at every node
result.species_phylogram   # the species tree in substitutions per site

# GTR with unequal base frequencies, under a relaxed (uncorrelated) clock
result = sequences.simulate_sequences(my_genomes,
    model=gtr(rates=(1, 2, 1, 1, 2, 1), freqs=(0.3, 0.2, 0.2, 0.3)),
    substitution=1.0 * mod.ByLineage(spread=0.3),   # the relaxed clock
    length=500, seed=1)
```

## Usage from the CLI

On the command line the genome run is handed over as a **directory**: `zombi2 sequences --genomes DIR` reads that run's species tree and event log and replays the gene genealogy from them, so the two commands chain without anything else passing between them.

```bash
# 1. genomes along a species tree (from the previous chapters)
zombi2 genomes -t out/species_complete.nwk \
    --duplication 0.2 --transfer 0.1 --loss 0.25 --origination 0.5 --seed 1 -o out/

# 2. HKY85, 1000 sites, strict clock
zombi2 sequences --genomes out/ --model hky85 --kappa 2.0 \
    --length 1000 --seed 1 -o seqs/

# GTR with unequal frequencies under a relaxed clock, also writing the ancestral sequences
zombi2 sequences --genomes out/ --model gtr \
    --frequencies 0.3 0.2 0.2 0.3 \
    --substitution "1.0 * ByLineage(spread=0.3)" \
    --seed 1 -o seqs/ --write alignments phylograms ancestral species_phylogram
```

The clock keeps its written form on the command line, exactly as in Python — `"1.0 * ByLineage(spread=0.3)"` is the same expression either way, so a rate can be moved between a script, a flag and a `--params` file without being rewritten.

## Outputs

A run writes, by default, one **alignment** per gene family in FASTA, with the extant gene copies as the aligned rows, and the **phylograms** those sequences were drawn along — each family's gene tree in Newick, with branch lengths in substitutions per site rather than time, so the ground-truth tree behind every alignment is kept beside it.

Two more outputs are written on request. The **ancestral sequences** (`--write ancestral`) give the sequence at every internal node, the raw material for scoring an ancestral-reconstruction method against the truth. The **species phylogram** (`--write species_phylogram`) is the species tree under the same conversion, which is where the clock becomes visible: a lineage that ran hot is simply a longer branch there than its age would suggest. Every node is labelled `g<copy>`, so a phylogram's tips pair with its alignment and its internal nodes with the ancestral sequences. The full list of files lives in Appendix B.
