# Sequence evolution

The sequence level does two main things:

* It produces gene trees and species tree scaled measured in substitutions per branch (**phylograms**)
* It simulates nucleotides, amino acids, or codons that sit inside every gene. 

The sequence level is always dependent on a genome-level run:

```python
from zombi2 import sequences
from zombi2.sequences.substitution_models import hky85
result = sequences.simulate_sequences(gene_trees, model=hky85(kappa=2.0),
                                      length=1000, seed=1)
```

CLAUDE -> DISCUSS WITH ME: gene_trees sould be in fact a genome level object

## Creating phylograms

A gene tree arrives as a **chronogram**, its branch lengths in time. What a sequence actually accumulates along a branch is not time but a number of *substitutions per site*, and that is time multiplied by an evolutionary rate. Turning one into the other is the whole job of the sequence level. Applying a rate to every branch rescales the tree from time into expected substitutions and yields a phylogram.

Two things therefore have to be chosen: *what* changes (the substitution model, the chemistry of which residue turns into which) and *how fast* it changes along each branch, which is the clock.

## The substitution models 

ZOMBI2 implements different standard models of sequence evolution:

```python
model = jc69()                    # equal rates, equal base frequencies — no free parameters
model = k80(kappa=2.0)            # a transition/transversion bias
model = hky85(kappa=2.0, freqs=(0.3, 0.2, 0.2, 0.3))
model = gtr(rates=(1,2,1,1,2,1), freqs=(0.25,0.25,0.25,0.25))   # six exchangeabilities + freqs
model = lg()                      # an empirical 20-state amino-acid matrix
model = gy94(omega=0.2)           # a codon model: dN/dS below 1 is purifying selection
```

The alphabet follows from the model: the DNA models (`jc69`, `k80`, `hky85`, `gtr`) evolve nucleotides, the protein models (`lg`, `wag`, `jtt`, `dayhoff`, `poisson`) evolve amino acids, and the codon models (`gy94`, `mg94`) evolve in-frame triplets and read `length` as a number of codons. 

## Rate variation

A substitution rate can vary along the tree in two ways:

- **Across lineages** — some lineages evolve faster than others. 
- **Across sites** — within one gene, some positions evolve faster than others. This is the classic **+Γ** of phylogenetics: each site draws a multiplier from a Gamma distribution of shape α. 

```python
sequences.simulate_sequences(gene_trees,
    model=gtr(rates=..., freqs=...),
    substitution=1.0 * mod.ByLineage(spread=0.3),   # the clock: across-branch variation
    gamma=0.5,                                      # +Γ: across-site variation (the shape α)
    length=1000, seed=1)
```

## Relaxed molecular clocks

By using the modifiers of ZOMBI2, it is easy to reproduce some of the most well known relaxed molecular clocks in the literature. 

```python
# strict clock — one rate everywhere; the default, so write nothing
substitution = 1.0

# uncorrelated / relaxed — each lineage draws its own rate, independently (no memory)
substitution = 1.0 * mod.ByLineage(spread=0.3)                  # lognormal (the default)
substitution = 1.0 * mod.ByLineage(spread=0.3, dist="gamma")   # or gamma

# autocorrelated — the rate drifts continuously down the tree (continuous memory)
substitution = 1.0 * mod.FromParent(spread=0.3)

# CIR — the same continuous drift, but pulled back toward a mean (mean-reversion)
substitution = 1.0 * mod.FromParent(spread=0.3, reverts_to=1.0, pull=0.5)

# the Markov clock — the rate hops between a few discrete categories (discrete memory)
substitution = 1.0 * mod.Markov(rates=[0.5, 1.0, 2.0], switch=0.1)
```

The important modifiers to remember are:

- **`ByLineage`** — *no memory*: each lineage is an independent draw, so a lineage's rate tells you nothing about its neighbours'. This is the uncorrelated or relaxed family; the distribution it draws from (`dist="lognormal"` or `"gamma"`) is a parameter, not a new modifier. An i.i.d. per-lineage draw is itself the "white-noise" clock.
- **`FromParent`** — *continuous memory*: a lineage starts from its parent's rate and drifts, so close relatives resemble each other. This is the autocorrelated clock; adding `reverts_to` (a mean) and `pull` (how hard it is drawn back) gives the CIR clock.
- **`Markov`** — *discrete memory*: the rate is one of a few fixed categories and hops between them along the tree at a switching rate.

A reference table that can be handy to people who want to implement a specific model from the literature:

| Literature name | What it does | ZOMBI2 |
|---|---|---|
| Strict / global clock | one rate everywhere | `substitution = 1.0` (default) |
| Uncorrelated lognormal (UCLN) | each lineage i.i.d. lognormal | `1.0 * mod.ByLineage(spread=…)` |
| Uncorrelated gamma (UGAM) | each lineage i.i.d. gamma | `1.0 * mod.ByLineage(spread=…, dist="gamma")` |
| White-noise clock | each lineage i.i.d. — that *is* white-noise | `1.0 * mod.ByLineage(spread=…)` |
| Autocorrelated lognormal (Thorne–Kishino) | rate drifts along the tree | `1.0 * mod.FromParent(spread=…)` |
| CIR clock | drift with mean-reversion | `1.0 * mod.FromParent(spread=…, reverts_to=…, pull=…)` |
| Discrete-category / random local clock | rate hops between categories | `1.0 * mod.Markov(rates=[…], switch=…)` |
| +Γ rate heterogeneity | variation across sites | `gamma=α` (not a clock) |
| GY94 / MG94 | codon model, one dN/dS | `model=gy94(omega=…)` |
| M-series site models (M1a…M8) | dN/dS varies across codon sites | `model=m8(...)` |

## Families evolving at different speeds

There is one more kind of rate variation: a **per-family speed**. Independently of any lineage or site effect, one whole gene family can evolve faster than another. This is `mod.ByFamily`, exactly the sequence-level reading of the genome level's modifier of the same name: each family draws one constant speed multiplier. 

```python
substitution = 1.0 * mod.ByFamily(spread=0.4)          # each gene family its own speed
substitution = 1.0 * mod.ByLineage(spread=0.3) * mod.ByFamily(spread=0.4)  # lineages and families
```

## The objects

(write briefly following what is written in Ch4)

## Usage from Python

An end-to-end run, from the gene trees a genome simulation produced through to alignments:

```python
from zombi2 import sequences
from zombi2.rates import modifiers as mod
from zombi2.sequences.substitution_models import hky85

# the common case: DNA under HKY, a strict clock, one rate per site
result = sequences.simulate_sequences(gene_trees, model=hky85(kappa=2.0),
                                      length=1000, seed=1)

# proteins under LG, a relaxed (uncorrelated) clock, and +Γ across sites
result = sequences.simulate_sequences(gene_trees,
    model=lg(),
    substitution=1.0 * mod.ByLineage(spread=0.3),   # the relaxed clock
    gamma=0.5,                                      # across-site heterogeneity
    length=500, seed=1)

# codons under GY94 with purifying selection and an autocorrelated clock
result = sequences.simulate_sequences(gene_trees,
    model=gy94(omega=0.2),
    substitution=1.0 * mod.FromParent(spread=0.3),  # ClaDS's modifier, one level down
    length=300, seed=1)
```

## Usage from the CLI

*[Draft — the sequence CLI is not part of the clean core yet, and its re-fit to this API is still to be designed; the command below is a provisional sketch.]*

```bash
zombi2 sequences --genomes my_genomes/ --model hky85 --kappa 2.0 \
    --length 1000 --clock relaxed --spread 0.3 --gamma 0.5 --seed 1 -o my_sequences
```

## Outputs

A run writes one **alignment** per gene family, in FASTA, with the extant genes as the aligned rows. And because every internal node was recorded, the **ancestral sequences** at each node can be written too, the raw material for scoring any ancestral-reconstruction method against the truth. It also writes the **phylograms** the sequences were drawn along — each gene tree, and the species tree, in Newick with branch lengths in substitutions per site rather than time — so the ground-truth tree behind each alignment is kept, and the molecular clock (which lineages ran fast or slow) is visible on the species tree. The full list of files lives in Appendix B.
