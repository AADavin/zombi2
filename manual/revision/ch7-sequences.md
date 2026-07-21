# Sequence evolution

The sequence level is where the model finally reaches the letters: the nucleotides, amino acids, or codons that sit inside every gene. A genome run tells you which gene families a lineage carries and hands you, for each family, its **gene tree**, a timetree whose branch lengths are measured in time. This chapter takes those gene trees and evolves a sequence down each one, so that every gene ends with an aligned column of residues at its tips. The entry point is one function:

```python
from zombi2 import sequences
from zombi2.sequences.substitution_models import hky85
result = sequences.simulate_sequences(gene_trees, model=hky85(kappa=2.0),
                                      length=1000, seed=1)
```

*[Draft — `simulate_sequences` is built as a first slice (the nucleotide models `jc69`/`k80`/`hky85`/`gtr` under a strict or `ByLineage` clock, following `docs/design/sequence-api.md`); the rest of the menu — the other clocks, +Γ across sites, and the protein and codon models — is not yet wired. The chapter documents the target; the divergences are noted as they arise.]*

## From a timetree to a phylogram

A gene tree arrives as a **chronogram**, its branch lengths in time. What a sequence actually accumulates along a branch is not time but a number of *substitutions per site*, and that is time multiplied by an evolutionary rate. Turning one into the other is the whole job of the sequence level. Applying a rate to every branch rescales the tree from time into expected substitutions and yields a **phylogram**, and the phylogram is the real deliverable: one branch length there is one expected substitution per site. ZOMBI2 then draws an actual sequence down the phylogram, letter by letter, under a substitution model.

Two things therefore have to be chosen: *what* changes (the substitution model, the chemistry of which residue turns into which) and *how fast* it changes along each branch, which is the clock. They are independent choices, and the chapter takes them in that order.

## The substitution model is a menu, not a zoo

Everywhere else in ZOMBI2 a family of models collapses into one shape: the seven species processes became one birth–death call, the eight clocks (below) become three modifiers. The substitution model does **not** collapse, and that is correct. JC69, HKY, GTR, and LG are not one model wearing different parameters; they are genuinely different rate matrices, different physical claims about how substitutions happen. Forcing a grammar over them would be a worse description than the honest one, which is a **menu of constructors**, each taking the parameters its chemistry actually has:

```python
model = jc69()                    # equal rates, equal base frequencies — no free parameters
model = k80(kappa=2.0)            # a transition/transversion bias
model = hky85(kappa=2.0, freqs=(0.3, 0.2, 0.2, 0.3))
model = gtr(rates=(1,2,1,1,2,1), freqs=(0.25,0.25,0.25,0.25))   # six exchangeabilities + freqs
model = lg()                      # an empirical 20-state amino-acid matrix
model = gy94(omega=0.2)           # a codon model: dN/dS below 1 is purifying selection
```

The alphabet follows from the model: the DNA models (`jc69`, `k80`, `hky85`, `gtr`) evolve nucleotides, the protein models (`lg`, `wag`, `jtt`, `dayhoff`, `poisson`) evolve amino acids, and the codon models (`gy94`, `mg94`) evolve in-frame triplets and read `length` as a number of codons. This is the theme worth stating outright: ZOMBI2 is elegant where the concepts really are one thing, and a plain menu where they are not. The menu is not a failure to unify; it is the truthful shape.

Two codon constructors need a word of honesty. `gy94` and `mg94` take an `omega` (the dN/dS ratio) that is constant across sites. The M-series codon models (`m1a`, `m2a`, `m3`, `m7`, `m8`) are also menu constructors, but each one bakes in a *distribution* of dN/dS across sites: some conserved, some free, a few under positive selection. They are the codon world's own way of writing across-site rate variation, and they stand in place of +Γ rather than beside it.

## Two axes of rate variation, and they are different

A substitution rate can vary along the tree in two ways that are constantly confused in the literature and must be kept apart here, because ZOMBI2 gives each its own argument.

- **Across branches** — some lineages evolve faster than others. This is the **clock**, and it is a modifier on the substitution rate. It is the one and only place in ZOMBI2 the word *clock* is allowed.
- **Across sites** — within one gene, some positions evolve faster than others. This is the classic **+Γ** of phylogenetics: each site draws a multiplier from a Gamma distribution of shape α. It is its own argument, `gamma=α`, and it is not a clock.

```python
sequences.simulate_sequences(gene_trees,
    model=gtr(rates=..., freqs=...),
    substitution=1.0 * mod.ByLineage(spread=0.3),   # the clock: across-branch variation
    gamma=0.5,                                      # +Γ: across-site variation (the shape α)
    length=1000, seed=1)
```

*[Figure 7.x — the two axes: the clock stretches whole branches, +Γ speeds and slows individual columns within a gene. To draw.]*

The scope never enters this level as a choice you make. Species are counted per lineage and genomes per copy, and you pick; a sequence is always counted **per site**, and there is nothing to tune. All of the interesting variation at this level therefore lives in the modifiers and in +Γ, which is exactly why the clock and +Γ carry the weight here.

## The clock collapses to three modifiers

The old code's `legacy/sequences/clocks.py` carries a zoo of clock classes (`StrictClock`, `UncorrelatedLogNormalClock`, `UncorrelatedGammaClock`, `WhiteNoiseClock`, `AutocorrelatedLogNormalClock`, `CIRClock`, `RateVariation`), exactly the kind the rest of the manual has been dismantling. A clock is not a special kind of object. It is a **modifier on the substitution rate**, written in the same `base × modifiers` grammar as every other rate in the book. Once you see that, the zoo becomes three modifiers plus the default, sorted by one question: *what does a lineage's rate remember?*

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

- **`ByLineage`** — *no memory*: each lineage is an independent draw, so a lineage's rate tells you nothing about its neighbours'. This is the uncorrelated or relaxed family; the distribution it draws from (`dist="lognormal"` or `"gamma"`) is a parameter, not a new modifier. An i.i.d. per-lineage draw is itself the "white-noise" clock — it needs no label of its own.
- **`FromParent`** — *continuous memory*: a lineage starts from its parent's rate and drifts, so close relatives resemble each other. This is the autocorrelated clock; adding `reverts_to` (a mean) and `pull` (how hard it is drawn back) gives the CIR clock — the same two knobs as the OU trait, one level down.
- **`Markov`** — *discrete memory*: the rate is one of a few fixed categories and hops between them along the tree at a switching rate.

The payoff is that two of these three are not new. **`mod.ByLineage` is the per-lineage twin of the genome level's `mod.ByFamily`**: the same idea of independent, this-unit-gets-its-own-rate heterogeneity, applied to lineages instead of families. And **`mod.FromParent` is literally the same modifier you met on the species tree**, where it was ClaDS: a rate that drifts as lineages split. The autocorrelated molecular clock and ClaDS diversification are one modifier, spelled the same way, at two different levels. That is the whole reason for a shared vocabulary, and this is where it pays off most visibly.

## The literature → command bridge

The old model names are useful (a reader arrives knowing they want "a CIR clock" or "+Γ"), but they do not organise anything. They live here, in one table, and nowhere else.

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| one rate everywhere | `substitution = 1.0` (default) | Strict / global clock |
| each lineage i.i.d. lognormal | `1.0 * mod.ByLineage(spread=…)` | Uncorrelated lognormal (UCLN) |
| each lineage i.i.d. gamma | `1.0 * mod.ByLineage(spread=…, dist="gamma")` | Uncorrelated gamma (UGAM) |
| each lineage i.i.d. — that *is* white-noise | `1.0 * mod.ByLineage(spread=…)` | White-noise clock |
| rate drifts along the tree | `1.0 * mod.FromParent(spread=…)` | Autocorrelated lognormal (Thorne–Kishino) |
| drift with mean-reversion | `1.0 * mod.FromParent(spread=…, reverts_to=…, pull=…)` | CIR clock |
| rate hops between categories | `1.0 * mod.Markov(rates=[…], switch=…)` | Discrete-category / random local clock |
| variation across sites | `gamma=α` (not a clock) | +Γ rate heterogeneity |
| codon model, one dN/dS | `model=gy94(omega=…)` | GY94 / MG94 |
| dN/dS varies across codon sites | `model=m8(...)` | M-series site models (M1a…M8) |

## A third source, and which tree the clock rides

There is one more kind of rate variation the two-axis picture above does not name: a **per-family speed**. Independently of any lineage or site effect, one whole gene family can evolve faster than another. This is `mod.ByFamily`, exactly the sequence-level reading of the genome level's modifier of the same name: each family draws one constant speed multiplier. Under the shared grammar it is written the same way:

```python
substitution = 1.0 * mod.ByFamily(spread=0.4)          # each gene family its own speed
substitution = 1.0 * mod.ByLineage(spread=0.3) * mod.ByFamily(spread=0.4)  # lineages and families
```

The composition is settled, and it turns on a point about *which tree the clock rides*. A clock is a property of a **lineage**: a whole species runs hot or cold, and every gene passing through that branch feels it. So the across-branch modifiers — `ByLineage`, `FromParent` — ride the **species** tree and give one clock value per species branch, shared by all its genes; each gene-tree branch simply reads the clock of the species branch it is reconciled to, which ZOMBI2 knows exactly, so nothing has to be wired up. The per-family speed, by contrast, is `mod.ByFamily`, the same modifier as the genome level: some families are fast whatever the lineage. The two axes are orthogonal and **compose**, reproducing the lineage-clock × per-family-speed pairing. (A fully per-gene-tree-branch clock — a single family fluctuating branch by branch *within* one lineage — is deferred as an exotic case; `ByLineage` is the species-branch clock.)

## The objects

*[Draft — `SequencesResult` is built with the payload below; the fuller shared spine (`.events`, `.tree`, and the lazy `record=` memory dial of `docs/design/result-api.md`) is the deferred target, noted where it arises.]*

`simulate_sequences` returns a **`SequencesResult`**, the sequence-level member of the one-bundle-per-level family (`SpeciesResult`, `GenomesResult`, `SequencesResult`, `TraitsResult`), all sharing a common shape and differing only in payload. Its payload is `.alignments` (the per-family alignments, one row per extant gene-tree tip), `.ancestral` (the reconstructed sequence at every internal node), and `.seed`, so the run reproduces. The fuller spine the other levels are growing — `.events`, the substitution-scaled phylograms reached through `.tree`, and the `record=` memory dial — is the deferred target.

A run records a sequence at **every** node, internal as well as tip, so ancestral sequence reconstruction is a byproduct, not an extra step: the tip sequences are the observable alignment (`.alignments`), and the internal-node sequences (`.ancestral`) are the ancestors at every split, which is also how the nucleotide genome model reconstructs the DNA of ancestral genomes. Today `.write(directory, outputs=[...])` streams the chosen outputs — the alignments, the ancestral sequences — to FASTA on disk; the lazy event-log replay and the `record=[...]` dial that scopes what a run keeps are the deferred target.

Today `simulate_sequences` applies the substitution scaling per gene-tree branch internally; surfacing those scaled phylograms on the result (as `.tree`) arrives with the deferred spine above.

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
