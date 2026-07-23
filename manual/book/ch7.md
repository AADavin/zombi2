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

The model decides the alphabet, and `length` counts whatever that alphabet holds: bases for a nucleotide model, residues for a protein one. The nucleotide models are four different rate matrices, not one model with four settings, but they do nest in the order written — `jc69` is `k80` with `kappa=1`, and `k80` is `hky85` with equal base frequencies — so each step down the list adds free parameters. The protein models work differently. Their rate matrices are **empirical**: each was estimated once from a large set of real alignments and is then used as a fixed table. That is why they take no parameters.

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

One important point: **the clock belongs to the species tree, not to the gene trees.**

ZOMBI2 draws one rate for each species branch. Every gene that passes through that branch then evolves at that rate. Each gene-tree branch looks up the species branch it sits inside, which the genome run already recorded. The consequence is that if a species evolves quickly, all of its genes evolve quickly together. 

A reference table that can be handy to people who want to implement a specific model from the literature:

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| one rate everywhere | `substitution = 1.0` (default) | Strict / global clock |
| each lineage i.i.d. lognormal | `1.0 * mod.ByLineage(spread=…)` | Uncorrelated lognormal (UCLN) |
| each lineage i.i.d. gamma | `1.0 * mod.ByLineage(spread=…, dist="gamma")` | Uncorrelated gamma (UGAM) |
| each lineage i.i.d. — that *is* white-noise | `1.0 * mod.ByLineage(spread=…)` | White-noise clock |

## The objects

`simulate_sequences` returns a **`SequencesResult`**, which carries:

- `.alignments` — the observable data: for each family, the sequence at every **extant** gene copy. This is the alignment a phylogenetic method would be handed.
- `.ancestral` — the sequence at every node that is **not** an extant tip: internal nodes, and the tips where a copy was lost or its species died. The run wrote a sequence at each node as it went, so these are the exact ancestors, not estimates. With `.alignments` it accounts for every node of the tree exactly once, so every label in a complete phylogram names a sequence.
- `.founding` — for each family, the sequence it began with, at its origination.
- `.phylograms` — for each family, its gene tree with branch lengths converted from time into substitutions per site: the tree the sequences were drawn along.
- `.species_phylogram` — the same conversion applied to the species tree, so the clock is visible as branch lengths.
- `.genomes`, `.initial_genome` — the assembled genome of every node, and of the run's starting point, present only when the run came from a **nucleotide** genome. See below.

As with every level, the bundle also carries `.seed` and `.write(directory, outputs=[...])` to put the chosen outputs on disk.

### Where a sequence starts

A family does not begin at the first branching of its gene tree. It begins when it originates, and the founding gene then lives for a while — its **stem** — before anything splits it. So that is where the sequence starts: one draw from the model's stationary frequencies at the origination, which then evolves across the stem in the ordinary way and arrives at the root gene as a sequence that has already changed. `.founding` is that first draw; `.ancestral` holds what the root gene ended up with, and the two differ by however much the stem allowed.

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

Hand it a **nucleotide** genome run instead you get the full fasta genomes. Genes and spacer get their own models. `model` evolves the genes; `intergene_model` evolves the spacer, at `intergene_speed` times the rate — 3× by default, and `jc69` by default, which is flat and has no free parameters.

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
result.genomes["n0"]             # the same at an ancestor — reconstructed, not estimated
result.initial_genome            # the genome the run started with, before anything happened
```

From the command line it is the same two commands as any other run:

```bash
zombi2 genomes out/ --resolution nucleotide --gff ecoli.gff --trim-overlaps \
  --inversion 5.0 --inversion-length 50000 --loss 2.0 --loss-length 8000 --seed 7

zombi2 sequences out/ --model hky85 --kappa 3.0 --substitution 0.02 \
  --intergene-speed 3.0 --seed 7
```

### Starting from a real sequence

So far the founding sequence of each block is *drawn* — from the model's frequencies, random ACGT. Hand the genomes run a **FASTA** alongside the GFF and it starts from the sequence you supply instead:

```python
my_genomes = genomes.simulate_genomes_nucleotide(
    tree, gff="ecoli.gff", fasta="ecoli.fasta",     # layout AND letters
    inversion=1.0, loss=0.3, seed=1)
result = sequences.simulate_sequences(my_genomes, model=hky85(kappa=3.0), substitution=0.05, seed=1)
```

The FASTA has one `>seqid` record per GFF `##sequence-region`, each exactly its declared length. Every block is then founded from the real DNA at its own initial coordinates, so an assembled genome descends from exactly what you gave. A gene that origination invents mid-run has no supplied DNA (it did not exist initially), so its block still draws from the model.

## Usage from the CLI

On the command line the genome run is handed over as a **directory** — the run directory itself, which by then holds the genomes. `zombi2 sequences out/` reads that run's species tree and event log and replays the gene genealogy from them, so the two commands chain without anything else passing between them. A nucleotide run given a `--fasta` also hands its initial DNA across (in `initial_sequence.fasta`), so the sequences descend from your real sequence without you naming it twice. Point `--from` at another run to read one and write somewhere else.

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

## Outputs

A run writes, by default, one **alignment** per gene family in FASTA, with the extant gene copies as the aligned rows, and the **phylograms** those sequences were drawn along — each family's gene tree in Newick, with branch lengths in substitutions per site rather than time, so the ground-truth tree behind every alignment is kept beside it.

The **clock species tree** (`clock_species_tree_complete.nwk`, and `…_extant.nwk`) comes with them: the species tree under the same conversion, and where the clock becomes visible — a fast-evolving lineage has a longer branch there than its age alone would give it. One output is written only on request: the **ancestral sequences** (`--write ancestral`) give the sequence at every internal node, which is what you need to score an ancestral-reconstruction method against the truth.
