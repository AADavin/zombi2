# Sequence evolution

The sequence level does two main things:

* It rescales the gene trees and the species tree from time into substitutions per site (**phylograms**).
* It evolves the residues that sit inside every gene, so each family ends with an alignment.

The sequence level is always dependent on a genome-level run, and it takes that run's result directly:

```python
from zombi2 import genomes, sequences
from zombi2.sequences.substitution_models import hky85

my_genomes = genomes.simulate_genomes_unordered(tree, duplication=0.2, transfer=0.1,
                                                loss=0.25, origination=0.5, seed=1)
result = sequences.simulate_sequences(my_genomes, model=hky85(kappa=2.0),
                                      length=1000, seed=1)
```

## Creating phylograms

A gene tree in ZOMBI2 is by default a **chronogram**, its branch lengths measure time. What a sequence actually accumulates along a branch is not time but a number of *substitutions per site*, and that is time multiplied by an evolutionary rate. Turning one into the other is the whole job of the sequence level. Applying a rate to every branch rescales the tree from time into expected substitutions and yields a phylogram.

Two things therefore have to be chosen: *what* changes (the substitution model, the chemistry of which residue turns into which) and *how fast* it changes along each branch, which is the clock.

## The substitution models 

ZOMBI2 implements different standard models of sequence evolution:

```python
# --- nucleotide models (4 states, ACGT) ---
model = jc69()                    # equal rates, equal base frequencies — no free parameters
model = k80(kappa=2.0)            # a transition/transversion bias
model = hky85(kappa=2.0, freqs=(0.3, 0.2, 0.2, 0.3))            # bias + unequal frequencies
model = gtr(rates=(1,2,1,1,2,1), freqs=(0.25,0.25,0.25,0.25))   # six exchangeabilities + freqs

# --- protein models (20 states, amino acids) ---
model = poisson()                 # equal rates, equal frequencies — the JC69 of proteins
model = jtt()                     # Jones, Taylor & Thornton 1992
model = dayhoff()                 # Dayhoff, Schwartz & Orcutt 1978
model = wag()                     # Whelan & Goldman 2001
model = lg()                      # Le & Gascuel 2008
```

The model decides the alphabet, and `length` counts whatever that alphabet holds: bases for a nucleotide model, residues for a protein one.

The nucleotide models are four different rate matrices, not one model with four settings, but they do nest in the order written — `jc69` is `k80` with `kappa=1`, and `k80` is `hky85` with equal base frequencies — so each step down the list adds free parameters.

The protein models work differently. Their rate matrices are **empirical**: each was estimated once from a large set of real alignments and is then used as a fixed table. That is why they take no parameters. You choose one, you do not tune it. `poisson` is the exception that proves the rule — it gives every replacement the same rate, and is useful as a null rather than as a description of any real protein.

## Relaxed molecular clocks

The rate itself is `substitution`, and it is counted **per site**: a gene-tree branch of Δ*t* time accrues `substitution · Δt` substitutions at every site. Leave it alone and it is `1.0` everywhere — the **strict clock**, one tempo for the whole tree.

Real lineages, however, evolve at different paces. A substitution rate that changes from lineage to lineage is what the field calls a **relaxed clock**. In ZOMBI2 this is not a new kind of object: you multiply the rate by a modifier, exactly as you do at every other level.

```python
# strict clock — one rate everywhere; the default, so write nothing
substitution = 1.0

# relaxed — each lineage draws its own rate, independently of its neighbours
substitution = 1.0 * mod.ByLineage(spread=0.3)                 # lognormal (the default)
substitution = 1.0 * mod.ByLineage(spread=0.3, dist="gamma")   # or gamma
```

**`ByLineage`** has *no memory*: each lineage is an independent draw, so a lineage's rate tells you nothing about its neighbours'. The distribution it draws from (`dist="lognormal"` or `"gamma"`) is a parameter of the modifier.

One point decides what the clock actually means: **the clock belongs to the species, not to the gene.**

ZOMBI2 draws one rate for each species branch. Every gene that passes through that branch then evolves at that rate. Each gene-tree branch looks up the species branch it sits inside, which the genome run already recorded.

The consequence is that if a species evolves quickly, all of its genes evolve quickly together. A single gene cannot speed up on its own.

A reference table that can be handy to people who want to implement a specific model from the literature:

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| one rate everywhere | `substitution = 1.0` (default) | Strict / global clock |
| each lineage i.i.d. lognormal | `1.0 * mod.ByLineage(spread=…)` | Uncorrelated lognormal (UCLN) |
| each lineage i.i.d. gamma | `1.0 * mod.ByLineage(spread=…, dist="gamma")` | Uncorrelated gamma (UGAM) |
| each lineage i.i.d. — that *is* white-noise | `1.0 * mod.ByLineage(spread=…)` | White-noise clock |

Three of those four names are the same modifier. A reader who wants "a UCLN clock" or "a white-noise clock" will find that both are `mod.ByLineage`, differing only in a parameter, or in nothing at all.

## The objects

`simulate_sequences` returns a **`SequencesResult`**, which carries:

- `.alignments` — the observable data: for each family, the sequence at every **extant** gene copy. This is the alignment a phylogenetic method would be handed.
- `.ancestral` — the sequence at every node that is **not** an extant tip: internal nodes, and the tips where a copy was lost or its species died. The run wrote a sequence at each node as it went, so these are the exact ancestors, not estimates. With `.alignments` it accounts for every node of the tree exactly once, so every label in a complete phylogram names a sequence.
- `.founding` — for each family, the sequence it began with, at its origination.
- `.phylograms` — for each family, its gene tree with branch lengths converted from time into substitutions per site: the tree the sequences were drawn along.
- `.species_phylogram` — the same conversion applied to the species tree, so the clock is visible as branch lengths.
- `.genomes`, `.ancestral_genomes`, `.initial_genome` — the assembled genome of each lineage, and of the run's starting point, present only when the run came from a **nucleotide** genome. See below.

As with every level, the bundle also carries `.seed` and `.write(directory, outputs=[...])` to put the chosen outputs on disk.

### Where a sequence starts

A family does not begin at the first branching of its gene tree. It begins when it originates, and the founding gene then lives for a while — its **stem** — before anything splits it. So that is where the sequence starts: one draw from the model's stationary frequencies at the origination, which then evolves across the stem in the ordinary way and arrives at the root gene as a sequence that has already changed. `.founding` is that first draw; `.ancestral` holds what the root gene ended up with, and the two differ by however much the stem allowed.

This is why a phylogram's root carries a branch length. It is the stem in substitutions per site, exactly as every other branch is its own stretch of time converted by the rate. Under a strict clock of rate 1 a phylogram is its chronogram, root branch included.

## Usage from Python

An end-to-end run, from a species tree through genomes to alignments:

```python
from zombi2 import species, genomes, sequences
from zombi2.rates import modifiers as mod
from zombi2.sequences.substitution_models import hky85, gtr, lg

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

# proteins under LG — 300 residues per gene
result = sequences.simulate_sequences(my_genomes, model=lg(), length=300, seed=1)
```

## Running on a nucleotide genome

Hand it a **nucleotide** genome run (Chapter 6) instead and the level does more, because that run knows where everything is. Every root block evolves — the spacer as well as the genes — each at its own length in base pairs, read from the genome. So `length` does not apply and is refused: one number would contradict the coordinates the genomes run already wrote.

Genes and spacer get their own models. `model` evolves the genes; `intergene_model` evolves the spacer, at `intergene_speed` times the rate — 3× by default, and `jc69` by default, which is flat and has no free parameters.

```python
from zombi2 import species, genomes, sequences
from zombi2.sequences.substitution_models import hky85

tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=5, seed=1).complete_tree
my_genomes = genomes.simulate_genomes_nucleotide(
    tree, gff="ecoli.gff", inversion=1.0, inversion_length=5000,
    duplication=0.3, loss=0.3, seed=1)

result = sequences.simulate_sequences(my_genomes, model=hky85(kappa=3.0),
                                      intergene_speed=3.0, substitution=0.05, seed=1)

result.genomes["n5"]             # {chromosome: sequence} — a whole assembled genome
result.ancestral_genomes["n0"]   # the same, at an internal node — a reconstructed ancestor
result.initial_genome            # the genome the run started with, before anything happened
```

Because the whole genome is covered, the run can put the genomes back together. `.genomes` holds one entry per extant lineage: its chromosomes, each one its blocks concatenated in physical order, reverse-complemented wherever the genome carries them inverted. This is the only place ZOMBI2 emits a genome as sequence rather than as coordinates, and it is a genome with a known history — every base of it traces to a block, a tree and an event.

`.ancestral_genomes` is the same for every other node — the ancestors, and the lineages that went extinct — and pairs with `.genomes` exactly as `.ancestral` pairs with `.alignments`: the extant tips are the observable half, everything else is the history behind them. So a run gives you the whole tree's genomes, reconstructed rather than estimated.

No node is left out. That is worth stating plainly, because it depends entirely on one decision: **every node votes on where the blocks are cut**, not just the surviving ones. A partition cut from the survivors alone leaves two kinds of hole — material no survivor kept has no block at all, and an ancestor can be left holding a *fragment* of a block whose genealogy, being the survivors', has no lineage for it. Counting every node closes both, because each node's own breakpoints are then in the partition, so each of its blocks is a whole number of blocks.

One genome is not a node's, and it is the one you are most likely to want: `.initial_genome`, the genome the run **started** with. The root *node* is not it. A node sits at the end of its branch, and the root branch is real simulated time — hand in an E. coli annotation at these rates and eight inversions, a transposition, a duplication and a loss can land on the stem before anything has speciated. So `.initial_genome` is the genome you handed in, and `.ancestral_genomes["n0"]` is that genome after the stem. This is the same distinction as `.founding` against `.ancestral` for a single gene, and it is why the initial genome gets its own field and its own file.

In a deliberately punishing 91-node run — 39 kb over three replicons, every event kind on, 128 losses, 80 transfers, chromosomes splitting and merging — all 91 genomes came back exactly: 25 extant lineages, 45 ancestors, 21 extinct lineages, 3.26 Mb compared base by base against what the run actually held. Cutting at every node instead of the survivors cost 1.34× as many blocks.

The model has to be a nucleotide one, since a genome is measured in base pairs and its blocks are read on either strand. A protein model is refused.

From the command line it is the same two commands as any other run — the handoff says which resolution wrote it, so nothing has to be repeated:

```bash
zombi2 genomes out/ --resolution nucleotide --gff ecoli.gff --trim-overlaps \
  --inversion 5.0 --inversion-length 50000 --loss 2.0 --loss-length 8000 --seed 7

zombi2 sequences out/ --model hky85 --kappa 3.0 --substitution 0.02 \
  --intergene-speed 3.0 --write genomes initial_genome --seed 7
```

On the real *E. coli* annotation that is 4 477 genes over 8 272 blocks and 23 Mb of assembled genome, in about nine seconds. `--length` is refused here — the genome sets the lengths — and `--intergene-model` / `--intergene-speed` apply only here.

### Blocks are numbered, not families

This is the one thing to get right when reading a nucleotide run's output. The integer key of `.alignments`, `.ancestral`, `.founding` and `.phylograms` is a **block index** — a position in the genome run's `.root_blocks` — because every block evolves and spacer belongs to no family. Meanwhile `.gene_spans` and `.gene_trees` are keyed by **gene family id**. Two numbering schemes, both plain integers over overlapping ranges, so mixing them up returns a different locus without complaining.

`block_of` is the join, and it is the only thing you should use:

```python
block = my_genomes.block_of(my_genomes.gene_names["dnaA"])
result.alignments[block]     # dnaA's alignment
result.phylograms[block]     # the tree it was drawn along
```

For the run above, `dnaA` is gene family 1 and root block 6. Asking for `phylograms[1]` gets you a stretch of spacer evolving at three times the rate — a plausible answer to the wrong question. `result.unit` says which scheme a result uses (`"family"` or `"block"`), and the filenames follow it: a nucleotide run writes `block6.fasta` and `phylogram_block6_complete.nwk`, never `fam`.

### The round trip

The check worth doing once, because it is the claim the whole level rests on: declare a genome in a GFF, evolve it down a tree until the leaves are thoroughly rearranged, then rebuild it from the descendants and compare with what you declared. At a substitution rate of zero nothing mutates, so the comparison is exact.

On real E. coli — NC_000913.3, 4,641,652 bp, 4,498 genes, a 25-node tree with 1,236 events — `.initial_genome` comes back identical, base for base, with all 4,477 annotated genes at the coordinates and on the strands the GFF gave them. Compare `.ancestral_genomes["n0"]` instead and it will usually *not* match, and correctly so: that is the genome after the stem's events, which is what was actually there. Both are exact; they are answers to different questions.

## Usage from the CLI

On the command line the genome run is handed over as a **directory** — the run directory itself, which by then holds the genomes. `zombi2 sequences out/` reads that run's species tree and event log and replays the gene genealogy from them, so the two commands chain without anything else passing between them. Point `--from` at another run to read one and write somewhere else.

```bash
# 1. genomes along a species tree (from the previous chapters)
zombi2 genomes out/ \
    --duplication 0.2 --transfer 0.1 --loss 0.25 --origination 0.5 --seed 1

# 2. HKY85, 1000 sites, strict clock
zombi2 sequences seqs/ --from out/ --model hky85 --kappa 2.0 \
    --length 1000 --seed 1

# GTR with unequal frequencies under a relaxed clock, also writing the ancestral sequences
zombi2 sequences seqs/ --from out/ --model gtr \
    --frequencies 0.3 0.2 0.2 0.3 \
    --substitution "1.0 * ByLineage(spread=0.3)" \
    --seed 1 --write alignments phylograms ancestral species_phylogram
```

A protein model is the same command with a different `--model`:

```bash
# proteins under LG, 300 residues per gene
zombi2 sequences seqs/ --from out/ --model lg --length 300 --seed 1
```

Because a protein model has no parameters, passing one is an error rather than a flag that gets quietly ignored: `--model lg --kappa 2.0` stops with *"these options don't apply to --model lg: --kappa"*.

The clock keeps its written form on the command line, exactly as in Python — `"1.0 * ByLineage(spread=0.3)"` is the same expression either way, so a rate can be moved between a script, a flag and a `--params` file without being rewritten.

## Outputs

A run writes, by default, one **alignment** per gene family in FASTA, with the extant gene copies as the aligned rows, and the **phylograms** those sequences were drawn along — each family's gene tree in Newick, with branch lengths in substitutions per site rather than time, so the ground-truth tree behind every alignment is kept beside it.

The **clock species tree** (`clock_species_tree_complete.nwk`, and `…_extant.nwk`) comes with them: the species tree under the same conversion, and where the clock becomes visible — a fast-evolving lineage has a longer branch there than its age alone would give it. One output is written only on request: the **ancestral sequences** (`--write ancestral`) give the sequence at every internal node, which is what you need to score an ancestral-reconstruction method against the truth.

Every node is labelled `g<copy>`, so a phylogram's tips match its alignment and its internal nodes match the ancestral sequences. A run on a nucleotide genome adds `genome_<lineage>.fasta`, one file per extant lineage with one record per chromosome — the assembled genome — and, on request, `genome_ancestral_<lineage>.fasta` for the internal nodes. The full list of files lives in Appendix B.
