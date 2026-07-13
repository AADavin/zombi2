# Sequences

Species trees and gene trees in ZOMBI2 are **timetrees** — branch lengths are time. To get the
branch lengths you would infer from *sequence data*, ZOMBI2 turns each gene tree into evolving
**molecular sequences** in two stages: first a **relaxed clock** rescales the tree from time into
expected **substitutions per site** (a chronogram into a phylogram), then a **substitution process**
evolves an actual alignment down it. Both run from the [`zombi2 sequence`](../cli.md) command (or
`SequenceEvolution` / `evolve_on_tree` in Python).

<figure markdown="span">
  ![DNA substitution models evolving alignments along the gene trees.](../img/seq_subst_models.svg){ width="560" }
</figure>

## Relaxed molecular clocks: time into substitutions

Sequence evolution accumulates *substitutions per site*, which is time multiplied by an evolutionary
**rate** that varies across the tree. A relaxed molecular clock is the model of that variation: it
draws a rate for each branch and rescales the tree from time into expected substitutions, turning a
**chronogram** into a **phylogram**. The *strict* clock is the baseline — one global rate, no
variation. A *relaxed* clock lets the rate change branch to branch, either **uncorrelated** (i.i.d.
draws, a branch's rate says nothing about its neighbours') or **autocorrelated** (a branch's rate is
anchored to its parent's, so nearby lineages evolve at similar rates). ZOMBI2 ships a whole family of
clocks in the `zombi2.sequences` namespace, all sharing one `Clock` interface — `scale(tree, seed=...)`
returns the phylogram — differing only in how the per-branch rate is drawn. In the language of
[Rates: a primer](rates.md), a relaxed clock is a **per-branch modifier** on the substitution rate,
and among-site variation (gamma rates) is the matching **per-site modifier** — the same idea applied
to two different contexts.

| Model | Family | Reach for it when |
| --- | --- | --- |
| **StrictClock** | baseline | you want a null: one constant `rate` on every branch |
| **UCLN** (`UncorrelatedLogNormalClock`) | uncorrelated | lineage-independent lognormal rate noise, `E[rate] = mean` |
| **UGAM** (`UncorrelatedGammaClock`) | uncorrelated | the same idea with a gamma spread tuned by a single `shape` |
| **WhiteNoise** (`WhiteNoiseClock`) | uncorrelated | short branches should vary most (variance ∝ 1/duration) |
| **AutocorrelatedLogNormal** (`AutocorrelatedLogNormalClock`) | autocorrelated | rate heritability matters — a geometric walk down lineages |
| **CIR** (`CIRClock`) | autocorrelated | mean-reverting rates that also vary *within* a branch |
| **RateVariation** | autocorrelated | discrete ordered rate **bins** with a nearest-bin walk (GTDB) |

### The clock models

#### StrictClock

A single rate `rate` on every branch (no rate variation): the phylogram is the chronogram uniformly
stretched by `rate`, so relative branch proportions are unchanged. This is the baseline every relaxed
clock relaxes; reach for it as a null.

#### Uncorrelated lognormal (UCLN)

Each branch draws an **independent** lognormal multiplier with `E[rate] = mean` for any spread:
`rate = mean · exp(𝒩(−σ²/2, σ))`. Larger `sigma` means more rate heterogeneity; `sigma = 0` is the
strict clock. Because the draws are i.i.d., a branch's rate is uninformative about its neighbours'
(Drummond et al. 2006).

#### Uncorrelated gamma (UGAM)

Each branch draws an **independent** gamma rate with mean `mean` and variance `mean²/shape`. The
single `shape` knob controls dispersion: large `shape` concentrates rates near `mean` (→ strict),
small `shape` spreads them widely (Drummond et al. 2006; PhyloBayes `-ugam`).

#### WhiteNoise

An uncorrelated clock whose branch multiplier is the integral of a white-noise rate over the branch:
gamma-distributed with mean `mean` and variance `mean²·σ²/Δt`, inversely proportional to branch
duration `Δt`. Long branches average the noise away (rate → `mean`); short branches are highly
variable. That branch-length dependence is what distinguishes it from UGAM; `sigma = 0` is strict
(PhyloBayes `-wn`).

#### AutocorrelatedLogNormal

The rate evolves down the tree as a geometric random walk anchored to the parent,
`R_child = R_parent · exp(𝒩(0, σ·√ℓ))` with `ℓ` the branch length in time. A child's rate is centred
on its parent's, so nearby lineages have similar rates. `sigma = 0` freezes the walk into a strict
clock at `root_rate`. This is the shorthand `--branch-speed` clock (Thorne, Kishino & Painter 1998).

#### CIR

The instantaneous rate follows a mean-reverting Cox–Ingersoll–Ross diffusion,
`dr = θ(μ − r)dt + σ√r dW`, which stays strictly positive and pulls back toward the long-run mean
`mean` (μ) at speed `theta` (θ), with volatility `sigma` (σ). The path is simulated *within* each
branch by Euler–Maruyama on sub-steps of at most `max_step`, so — unlike the lognormal walk — the rate
also varies within a branch, while a child still starts where its parent ended (Lepage et al. 2007).

#### RateVariation

The discrete-bin, within-branch clock used in the GTDB archaea study. An **ordered** set of positive
rate `bins` is laid down and a continuous-time Markov process runs along the phylogeny, stepping only
to an **adjacent** bin (index ± 1) at `switch_rate`, with `up_bias` the probability a step goes to the
faster neighbour. Because the rate changes gradually, a single branch may split into several segments
in neighbouring bins; `switch_rate = 0` freezes it in its `start` bin (a strict clock).

The [`RateVariation`](#the-discrete-bin-gtdb-clock) worked example below covers the discrete-bin clock
in detail; the others follow the same `scale(tree, seed=...)` pattern.

### The discrete-bin (GTDB) clock

`RateVariation` — the discrete-bin, Markov-switching clock from the GTDB archaea study — works like
this:

- An **ordered** set of rate **bins** — multipliers, some fast (> 1), some slow (< 1).
- A continuous-time Markov process runs **along the phylogeny** and can only step to an
  **adjacent bin** (index ± 1). Because the rate changes gradually, nearby lineages have
  similar rates — the process is *autocorrelated*. The current bin is inherited by both
  descendants at each node.
- A branch may thus be split into several **segments** in adjacent bins; its substitution
  length is `Σ (segment_duration × bin_rate)`.

```python
from zombi2.species import simulate_species_tree, BirthDeath
from zombi2.sequences import RateVariation

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=30, age=8.0, seed=1)

rv = RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0],  # ordered slow -> fast
                   switch_rate=1.0,   # rate of stepping to a neighbouring bin
                   up_bias=0.5,       # P(step up); 0.5 = symmetric walk
                   start=2)           # initial bin index (default: the middle bin)
scaled = rv.scale(tree, seed=1)

print(scaled.to_newick())          # the phylogram (substitution lengths)
scaled.branch_lengths[node]        # substitution length of one branch
scaled.segments[node]              # [(bin_index, duration), ...] pieces of that branch
```

- Bin **order matters**: the walk only moves between neighbours, so put the bins in
  increasing (or decreasing) order of rate.
- `switch_rate=0` gives a strict clock (frozen in `start`).
- A symmetric walk (`up_bias=0.5`) has a uniform stationary distribution, so over a large
  tree the average `substitution / time` ratio approaches `mean(bins)`; `up_bias` shifts
  the chain toward faster (>0.5) or slower (<0.5) bins.

#### Works on gene trees too

`RateVariation` operates on any `Tree`. Gene trees come out of reconstruction as Newick, so
load them with `read_newick` first:

```python
from zombi2.genomes import simulate_genomes
from zombi2.sequences import RateVariation
from zombi2 import read_newick

genomes = simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5, seed=1)
_, extant = genomes.gene_trees()["1"]
gene_tree = read_newick(extant)          # Newick -> Tree
phylogram = RateVariation(bins=[0.5, 2.0], switch_rate=1.0).scale(gene_tree, seed=1)
```

The same rate process can be applied independently to the species tree and to each gene
tree, or you could drive them from a shared process — a natural building block for
simulating realistic, non-clocklike branch lengths. That shared-process version is
`SequenceEvolution`, next.

### Sequence evolution: the gene × lineage clock

`RateVariation` (and any clock) scales **one tree in isolation**. But in a real dataset every gene
shares the same underlying lineage history: a fast-evolving clade is fast for *all* the genes passing
through it. `SequenceEvolution` builds that in. It converts a whole
[`simulate_genomes`](genomes.md) result — all the reconciled gene trees at once —
from time into substitutions/site under a **gene × lineage** clock:

```
rate(family g, species branch b) = R_b · s_g
```

- **`R_b`** — a lineage rate drawn **once on the species tree** and **shared by every
  family**. This is any relaxed clock from the family: the autocorrelated **lognormal** walk
  `R_child = R_parent · exp(N(0, σ·√branch_length))` via `branch_sigma` (`0` = strict clock), or
  **any `Clock`** passed as `lineage=` — the **discrete-bin** `RateVariation` (the GTDB model), an
  uncorrelated clock, `CIRClock`, and so on.
- **`s_g`** — each family's intrinsic speed, one constant drawn from a distribution
  (`family_speed`), so some families are globally fast, others slow.

A gene-tree branch on species branch `b` over `[t0, t1]` gets substitution length
`s_g · R_b · (t1 − t0)`; because reconciliation is exact, a branch spanning several species
branches (after pruning) just sums the pieces.

```python
from zombi2.species import simulate_species_tree, BirthDeath
from zombi2.genomes import simulate_genomes
from zombi2.sequences import SequenceEvolution
from zombi2 import LogNormal

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=30, age=8.0, seed=1)
genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.2,
                           origination=0.5, seed=1)

se = SequenceEvolution(branch_sigma=0.5,                    # shared lineage clock
                       family_speed=LogNormal(0.0, 0.4))    # per-family speed
phylo = se.scale(genomes, seed=2)

phylo.extant["1"]        # family 1's extant tree, branch lengths in substitutions/site
phylo.complete["1"]      # the complete tree (losses included)
phylo.family_speed["1"]  # the drawn s_g
phylo.branch_rate["i5"]  # the shared R_b for species branch i5
```

With `branch_sigma=0` and unit family speeds the phylogram is identical to the input
chronogram — the two rate sources are clean multiplicative overlays. `branch_sigma=0` alone
gives *gene-family speed only* (each family a constant multiple of time); `family_speed=1.0`
alone gives a *shared lineage clock only*.

To use the discrete-bin GTDB clock for the lineage part instead of the lognormal one, pass a
`RateVariation` as `lineage=` (mutually exclusive with `branch_sigma`) — a branch may then
carry several rate segments, integrated exactly under each gene branch:

```python
from zombi2.sequences import RateVariation, SequenceEvolution
from zombi2 import LogNormal

rv = RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=1.0)  # slow -> fast
phylo = SequenceEvolution(lineage=rv, family_speed=LogNormal(0.0, 0.4)).scale(genomes, seed=2)
```

### Python: any clock in the family

Import any clock from `zombi2.sequences` (all are re-exported at the `zombi2` top level too, so
`zombi2.UncorrelatedLogNormalClock` also works):

```python
import zombi2 as z
from zombi2.sequences import (
    StrictClock, UncorrelatedLogNormalClock, UncorrelatedGammaClock,
    WhiteNoiseClock, AutocorrelatedLogNormalClock, CIRClock, RateVariation,
)

# A timetree (chronogram): branch lengths are time.
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=30, age=8.0, seed=1)

# Pick any clock and rescale time -> expected substitutions/site (a phylogram).
clock = UncorrelatedLogNormalClock(sigma=0.5, mean=1.0)   # i.i.d. lognormal per branch
phylo = clock.scale(tree, seed=2)

phylo.to_newick()             # the phylogram in Newick (substitution branch lengths)
node = next(n for n in tree.nodes_preorder() if n.parent is not None)
phylo.branch_lengths[node]    # substitution length of that branch
phylo.branch_rate[node]       # the (time-averaged) rate applied to it

# The other members — same scale(...) call, different rate rule:
UncorrelatedGammaClock(shape=3.0, mean=1.0).scale(tree, seed=2)
WhiteNoiseClock(sigma=0.5).scale(tree, seed=2)
AutocorrelatedLogNormalClock(sigma=0.4).scale(tree, seed=2)     # rate anchored to the parent's
CIRClock(theta=1.0, sigma=0.4, mean=1.0).scale(tree, seed=2)    # mean-reverting, varies within a branch
RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=1.0).scale(tree, seed=2)
StrictClock(1.0).scale(tree, seed=2)                            # the sigma-free baseline
```

A clock scales any `Tree`, so the same call rescales a gene tree loaded with `z.read_newick(...)`. To
drive one shared lineage clock across every gene family at once, pass the clock to
`zombi2.sequences.SequenceEvolution`.

### From the CLI

Sequence evolution is a **separate command**, `zombi2 sequence`, run on a prior `genomes`
result — so you can retune the clock without re-simulating gene content. `genomes` just has to
have written the event trace (`trace` in `--write`):

```bash
zombi2 genomes  -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace profiles -o run/

# any clock in the family via --clock (here an uncorrelated lognormal)
zombi2 sequence --genomes run/ --clock uncorrelated-lognormal --clock-sigma 0.5 \
    --family-speed 0.5 -o run/

# the autocorrelated lognormal and the discrete-bin GTDB clock also have shorthands
zombi2 sequence --genomes run/ --branch-speed 0.4 --family-speed 0.5 -o run/
zombi2 sequence --genomes run/ --branch-bins 0.25,0.5,1,2,4 --branch-switch-rate 1.0 \
    --family-speed 0.5 -o run/
```

`--clock` selects any model in the family —
`{strict, autocorrelated-lognormal, uncorrelated-lognormal, uncorrelated-gamma, white-noise, cir,
discrete-bin}` — with its parameter given by `--clock-sigma` (lognormal / white-noise / CIR spread),
`--clock-shape` (gamma), `--clock-theta` (CIR mean-reversion), or `--clock-mean` (the target/strict
rate). `--branch-speed` and `--branch-bins` remain as shorthands for the autocorrelated-lognormal and
discrete-bin clocks. Parameter defaults are `--clock-mean 1.0`, `--clock-sigma 0.5`,
`--clock-shape 3.0`, and `--clock-theta 1.0`. The discrete-bin shorthand also takes
`--branch-switch-rate` and `--branch-up-bias`.

`sequence` replays `run/Events_trace.tsv`, rebuilds the reconciled gene trees, and writes
`run/gene_trees/<family>_extant_subst.nwk` (and `_complete_subst.nwk`), plus
`gene_family_speeds.tsv` and `branch_rates.tsv` recording the drawn `s_g` and `R_b` for
reproducibility. A lineage clock alone (any `--clock`, `--branch-speed`, or `--branch-bins`) is a
shared lineage clock; `--family-speed` alone is per-family speed; together they are the full
gene × lineage model. `--family-speed SIGMA` draws each family's constant multiplier ~ LogNormal(0, SIGMA).

### Clock output

`zombi2 sequence` writes, into the output directory: `gene_trees/` — the rescaled substitution-unit
gene trees, `<family>_complete_subst.nwk` and `<family>_extant_subst.nwk` per family;
`branch_rates.tsv` — the per-species-branch clock rate applied to the shared lineage clock;
`gene_family_speeds.tsv` — each family's intrinsic `--family-speed` multiplier; and `sequence.log`,
the run manifest. Adding `--subst-model` also writes an `alignments/` directory of simulated
DNA/protein alignments, one per family.

### Clock validation

Each clock is validated by its defining property; every check below is a real test in `tests/`.

- **StrictClock** — a unit strict clock reproduces the chronogram exactly
  (`test_clocks.py::test_strict_clock_of_rate_one_is_the_chronogram`).
- **UCLN** — at `sigma = 0` it equals the strict clock
  (`test_clocks.py::test_uncorrelated_lognormal_sigma_zero_is_strict`), and unit-mean it leaves total
  tree length ≈ unchanged (`test_clocks.py::test_unit_mean_clocks_average_to_one`).
- **UGAM** — as an uncorrelated clock it shows near-zero parent↔child rate correlation, alongside UCLN
  and WhiteNoise (`test_clocks.py::test_uncorrelated_clocks_have_near_zero_parent_child_correlation`).
- **WhiteNoise** — the same near-zero correlation check, and unit-mean averaging to one
  (`test_clocks.py::test_uncorrelated_clocks_have_near_zero_parent_child_correlation`,
  `::test_unit_mean_clocks_average_to_one`).
- **AutocorrelatedLogNormal** — clearly positive parent↔child rate correlation, the split that names
  the family (`test_clocks.py::test_autocorrelated_clocks_have_positive_parent_child_correlation`).
- **CIR** — positive parent↔child correlation
  (`test_clocks.py::test_autocorrelated_clocks_have_positive_parent_child_correlation`) and reversion
  to its long-run mean (`test_clocks.py::test_cir_reverts_to_its_long_run_mean`).
- **RateVariation** — a symmetric discrete-bin walk averages to `mean(bins)`
  (`test_rate_variation.py::test_symmetric_walk_mean_rate_matches_bin_mean`).

These are strong oracle and statistical checks — the σ = 0 reductions are exact equalities, and the
uncorrelated/autocorrelated correlation contrast and mean-rate properties are pinned within tolerances
on large simulated trees. They confirm each clock's defining behaviour, not that inference recovers the
parameters.

## Simulating alignments (DNA and protein)

Add `--subst-model` and `sequence` goes one step further: it evolves an **actual sequence
alignment** down each rescaled gene tree (the rescaled branch lengths *are* the expected
substitutions per site). One FASTA per family is written to `run/alignments/<family>.fasta`,
its records headed by the same `<species>_<gene-id>` labels the leaves carry in the
`_extant_subst.nwk` tree. DNA versus protein is auto-detected from the model name.

```bash
# DNA alignments (HKY85), 600 bp
zombi2 sequence --genomes run/ --branch-speed 0.4 --family-speed 0.5 \
    --subst-model hky85 --seq-length 600 -o run/

# protein alignments under LG, 300 aa, with +Γ across-site rate heterogeneity
zombi2 sequence --genomes run/ --branch-speed 0.4 \
    --subst-model lg --seq-length 300 --gamma-shape 0.5 -o run/
```

Available models:

- **DNA (4 states, `ACGT`):** `jc69`, `k80` (`--kappa`), `hky85` (`--kappa`, `--base-freqs`),
  `gtr` (`--gtr-rates`, `--base-freqs`).
- **Amino acid (20 states):** `poisson` (equal rates, uniform frequencies) and the empirical
  matrices `lg` (Le & Gascuel 2008, the modern default), `wag` (Whelan & Goldman 2001), `jtt`
  (Jones et al. 1992), and `dayhoff` (Dayhoff et al. 1978) — exchangeabilities and frequencies
  transcribed byte-for-byte from the reference PAML data files. The protein models take no
  parameters.

`--gamma-shape ALPHA` adds discrete-Gamma (+Γ) across-site rate heterogeneity for any model;
`--seq-length N` sets the alignment length; `--root-fasta FILE` seeds each family's root from a
FASTA keyed by family id instead of a random draw from the stationary distribution. Omit
`--subst-model` and `sequence` only rescales the trees, exactly as before — no alignments are
written.

## DNA substitution models

Once `genomes` has produced gene trees and `sequence` has rescaled them into **substitutions per
site** (a phylogram), ZOMBI2 evolves an actual DNA alignment down each tree under a **continuous-time
Markov substitution model** over the four bases `ACGT`. Every model here is time-reversible and
normalised to one expected substitution per site per unit branch length, so the branch lengths mean
exactly what they say. The transition matrix over a branch of length `t` is `P(t) = exp(Q·t)`,
computed numpy-only by a reversible eigendecomposition. The models below are the classic nested
family — JC → K80 → HKY85 → GTR, each relaxing one assumption of the last — and any of them can carry
**+Gamma** among-site rate variation.

| Model | Parameters | Reach for it when |
| --- | --- | --- |
| **JC69** (JC) | none — equal rates, equal base freqs | the simplest baseline, no free parameters |
| **K80** (K2P) | `kappa` (transition/transversion ratio) | transitions and transversions differ but bases are even |
| **HKY85** (HKY) | `kappa` + unequal base freqs | a ti/tv bias *and* skewed base composition |
| **GTR** | 6 exchangeabilities + base freqs | the most general reversible DNA model |
| **+Gamma** (GammaRates) | `shape` (α), overlaid on any of the above | rates vary from site to site |

### The DNA models

#### Jukes–Cantor (JC69)

`jc69()` — the one-parameter-free baseline: all six exchangeabilities equal, all four base
frequencies `0.25`. Every substitution is equally likely, so a branch of length `t` gives the closed
form `P_ii = 1/4 + 3/4·e^{-4t/3}`. Reach for it as the null model or a sanity check.

#### Kimura 2-parameter (K80)

`k80(kappa=2.0)` — separates **transitions** (A↔G, C↔T) from **transversions** by a single ratio
`kappa`, keeping base frequencies equal at `0.25`. `kappa 1` reduces to JC69. Use it when the ti/tv
bias matters but base composition is even.

#### HKY85 (HKY)

`hky85(kappa=2.0, freqs=(0.25,0.25,0.25,0.25))` — K80's transition bias `kappa` plus **unequal
equilibrium base frequencies** `(A,C,G,T)`. Use it when the sequences have both a ti/tv bias and a
skewed base composition; the long-branch stationary distribution recovers `freqs`.

#### GTR

`gtr(rates=(1,1,1,1,1,1), freqs=(0.25,0.25,0.25,0.25))` — the **general time-reversible** model: six
free exchangeabilities `[AC,AG,AT,CG,CT,GT]` and arbitrary base frequencies. It is the most general
reversible DNA model and the superclass of all the others (`rates` all `1`, `freqs` all `0.25` gives
JC69). Reach for it when you want no built-in symmetry assumptions.

#### +Gamma among-site rate variation (GammaRates)

`GammaRates(shape, k=4)` — not a substitution model on its own but an **overlay** on any of the four
above. Sites are binned into `k` equal-probability discrete-Gamma categories (Yang 1994) with mean
rate 1; a small `shape` (α) makes rates highly heterogeneous across sites, a large `shape` makes them
nearly uniform. Use it whenever real substitution rates vary from one site to the next.

### DNA from the command line

Sequence simulation is its own step: run [`zombi2 sequence`](../cli.md#sequence) on a `genomes` run
that was written with `trace` in `--write`. `--subst-model` picks the model (`jc69`/`k80`/`hky85`/`gtr`
for DNA; the name auto-detects DNA vs protein), and the model-specific knobs are `--kappa`,
`--base-freqs`, `--gtr-rates`, and `--gamma-shape`. Omit `--subst-model` to only rescale the trees
without simulating sequences.

```bash
# a genomes run written with the event trace (prerequisite)
zombi2 species --birth 1 --death 0.3 --tips 8 --age 3 --seed 1 -o run/
zombi2 genomes -t run/species_tree.nwk --dup 0.3 --trans 0.1 --loss 0.3 \
    --orig 0.5 --write trace --seed 1 -o run/

# JC69 — no free parameters
zombi2 sequence --genomes run/ --subst-model jc69 --branch-speed 0.4 --seed 7 -o jc/

# K80 with a transition/transversion ratio
zombi2 sequence --genomes run/ --subst-model k80 --kappa 4 --seed 7 -o k80/

# HKY85: ti/tv bias plus unequal base frequencies (A C G T)
zombi2 sequence --genomes run/ --subst-model hky85 --kappa 4 \
    --base-freqs 0.4 0.1 0.1 0.4 --seed 7 -o hky/

# GTR: 6 exchangeabilities [AC AG AT CG CT GT] + base freqs, with +Gamma across sites
zombi2 sequence --genomes run/ --subst-model gtr \
    --gtr-rates 1 2.5 1 1 2.5 1 --base-freqs 0.3 0.2 0.2 0.3 \
    --gamma-shape 0.5 --branch-speed 0.4 --seed 7 -o gtr/
```

`--seq-length N` sets the alignment length (default 300); `--root-fasta` seeds each family's root
from a FASTA instead of a random draw. The `--branch-speed`/`--family-speed`/`--clock` knobs govern
the relaxed clock that turns the time tree into a phylogram — see the
[substitution branch-lengths](../cli.md#sequence) section.

### DNA from Python

The models live in `zombi2.sequences` as factory functions that build a `SubstitutionModel`, plus
`GammaRates` and the `evolve_on_tree` simulator (the same objects also re-export at the top level, so
`zombi2.hky85` works too):

```python
import numpy as np
from zombi2.sequences import jc69, k80, hky85, gtr, GammaRates, evolve_on_tree

# a minimal tree node has a .gid and .children (real trees come from a genomes run)
class Node:
    def __init__(self, gid, children=()):
        self.gid, self.children = gid, list(children)

a, b = Node("a"), Node("b")
root = Node("r", [a, b])
subst = {root: 0.0, a: 0.2, b: 0.2}          # branch lengths in substitutions/site

# HKY85 (kappa + unequal freqs) with +Gamma across-site rate variation
model = hky85(kappa=4.0, freqs=(0.4, 0.1, 0.1, 0.4))
seqs = evolve_on_tree(root, subst, model, np.random.default_rng(0),
                      length=1000, gamma=GammaRates(shape=0.5))
# seqs maps each node's gid -> its DNA string (internal nodes and tips)

# GTR with six explicit exchangeabilities [AC AG AT CG CT GT]
g = gtr(rates=(1, 2.5, 1, 1, 2.5, 1), freqs=(0.3, 0.2, 0.2, 0.3))
```

In practice you drive this through the `zombi2 sequence` command (or `SequenceEvolution`), which
supplies the rescaled gene trees; `evolve_on_tree` is the low-level engine underneath.

### DNA output

`zombi2 sequence --subst-model MODEL` writes, under `-o`:
`gene_trees/` with the substitution-unit phylograms (`<fam>_complete_subst.nwk` and
`<fam>_extant_subst.nwk` per family), `alignments/<fam>.fasta` with the simulated per-family DNA
alignment (one record per surviving gene copy), `branch_rates.tsv` and `gene_family_speeds.tsv`
recording the clock rate on each branch and the per-family speed multiplier, and `sequence.log` the
run manifest (every parameter, including the substitution model and its `--kappa`/`--base-freqs`/
`--gtr-rates`/`--gamma-shape` settings).

### DNA validation

- **JC69.** The numpy-only `exp(Qt)` matches the JC69 closed form `P_ii = 1/4 + 3/4·e^{-4t/3}` at
  several branch lengths (`test_sequence_sim.py::test_p_matrix_matches_jc_closed_form`), and a pair of
  tips `0.2` from the root recovers the true `0.4` JC-corrected distance to within `0.02`
  (`test_sequence_sim.py::test_jc_distance_recovered`).
- **K80.** The observed transition/transversion structure across an evolved branch matches the
  Kimura-1980 closed form implied by `kappa` — the per-site transition and transversion difference
  fractions, and their ratio, agree with the closed-form `p_ti`/`p_tv` derived from `kappa` (an
  oracle, not just matrix invariants)
  (`test_sequence_sim.py::test_k80_transition_transversion_matches_kappa`).
- **HKY85.** A single long-branch star run recovers the specified unequal base frequencies
  `(0.4,0.1,0.1,0.4)` to within `0.01`
  (`test_sequence_sim.py::test_stationary_frequencies_recovered`).
- **GTR.** A long GTR run, started from a uniform (25%-each) root that is *not* the target, recovers
  the asymmetric stationary base frequencies `pi=(0.1,0.2,0.3,0.4)` to within `0.01`
  (`test_sequence_sim.py::test_gtr_stationary_frequencies_recovered`).
- **+Gamma.** Overlaying the Gamma over-disperses the per-site divergence versus no-Gamma: the mean
  and variance of the per-site count of differing tips match the closed-form Binomial-vs-Gamma-mixture
  moments (law of total variance over the category mixture), and the +Gamma variance is strictly
  larger (`test_sequence_sim.py::test_gamma_overdisperses_per_site_divergence`). The discrete-Gamma
  category rates themselves average to exactly 1 with increasing, distinct categories
  (`test_sequence_sim.py::test_gamma_rates_numpy_only_mean_one`).

### DNA references

- Jukes, T. H. & Cantor, C. R. (1969). Evolution of protein molecules. In *Mammalian Protein
  Metabolism*, 21–132.
- Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through
  comparative studies of nucleotide sequences. *Journal of Molecular Evolution* 16: 111–120.
- Hasegawa, M., Kishino, H. & Yano, T. (1985). Dating of the human–ape splitting by a molecular
  clock of mitochondrial DNA. *Journal of Molecular Evolution* 22: 160–174.
- Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis of DNA sequences.
  *Lectures on Mathematics in the Life Sciences* 17: 57–86.
- Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates
  over sites: approximate methods. *Journal of Molecular Evolution* 39: 306–314.

## Protein substitution models

Once a `genomes` run's gene trees have been rescaled from time into substitutions per site by the
relaxed molecular clock, ZOMBI2 can evolve an **amino-acid alignment** along each tree. A protein
substitution model is a 20-state **time-reversible** rate matrix `Q_ij = S_ij · π_j`,
built from a symmetric exchangeability matrix `S` and the stationary amino-acid frequencies `π`, and
normalised to one expected substitution per site. All the empirical models below are transcribed
byte-for-byte from the reference PAML data files (Ziheng Yang, `abacus-gene/paml`); ZOMBI2 differs
from a nucleotide model only in the alphabet (`K = 20`) — the same `evolve_on_tree` engine runs both.

| Model | Exchangeabilities `S` / frequencies `π` | Reach for it when |
| --- | --- | --- |
| **Poisson** | all equal / uniform (1/20) | you want a neutral, parameter-free baseline (F81-for-proteins) |
| **LG** | Le & Gascuel 2008 empirical | the modern default for most protein data |
| **WAG** | Whelan & Goldman 2001 empirical | a widely used, well-established general matrix |
| **JTT** | Jones, Taylor & Thornton 1992 empirical | matching an older analysis or JTT-based pipeline |
| **Dayhoff** | Dayhoff, Schwartz & Orcutt 1978 empirical | reproducing classic Dayhoff-PAM results |

### The protein models

#### Poisson

Equal off-diagonal exchangeabilities and uniform frequencies (`π_i = 1/20`) — the protein analogue of
Felsenstein-81. No parameters and no amino-acid preference: every residue is equally likely and every
substitution equally rated. Use it as a neutral baseline or a sanity check, not as a realistic model
of protein evolution.

#### LG

The Le & Gascuel (2008) exchangeabilities and frequencies, estimated across a large curated alignment
database with across-site rate variation accounted for during estimation. The current default choice
for most empirical protein work; reach for it unless a specific pipeline dictates otherwise.

#### WAG

The Whelan & Goldman (2001) matrix, estimated by maximum likelihood from a broad set of globular
protein families. A long-established, general-purpose model; a reasonable choice when comparability
with the extensive WAG literature matters.

#### JTT

The Jones, Taylor & Thornton (1992) matrix (PAML's `jones.dat`), derived from a large survey of
protein sequences. Use it to match an older analysis or a JTT-based reconstruction pipeline.

#### Dayhoff

The original Dayhoff, Schwartz & Orcutt (1978) exchangeabilities (PAML's `dayhoff.dat`), the ancestor
of the PAM family. Mostly of historical and reproduction value; reach for it to recreate classic
Dayhoff-PAM results.

### Protein from the command line

Protein evolution runs through the `sequence` command: pass a prior `genomes` output directory and
select a protein model with `--subst-model` (DNA vs protein is auto-detected from the name). The
`genomes` run must have been done with `trace` in `--write`. Alignment length is set with
`--seq-length` (amino acids; default 300), and `--gamma-shape` adds discrete-Gamma across-site rate
heterogeneity. The empirical models take no further parameters — `S` and `π` are fixed. The lineage
clock knobs (e.g. `--branch-speed`, `--family-speed`) rescale the trees exactly as for DNA.

```bash
# 1) a genomes run recorded with the event trace
zombi2 species --birth 1 --death 0.3 --tips 8 --age 3 --seed 1 -o run/
zombi2 genomes -t run/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace --seed 1 -o run/

# 2) evolve LG protein alignments along the rescaled gene trees
zombi2 sequence --genomes run/ --subst-model lg --seq-length 200 \
    --branch-speed 0.4 --seed 7 -o seqs/

# WAG with across-site rate heterogeneity
zombi2 sequence --genomes run/ --subst-model wag --gamma-shape 0.5 --seed 7 -o seqs/
```

### Protein from Python

The model builders and the sequence engine live in `zombi2.sequences` (each name is re-exported at
the top level too, so `zombi2.lg` also works):

```python
import numpy as np
from zombi2.sequences import lg, poisson, make_model, evolve_on_tree, AMINO_ACIDS, PROTEIN_MODELS

model = lg()                       # or wag(), jtt(), dayhoff(), poisson()
model = make_model("wag")          # ...or build one by name
print(PROTEIN_MODELS)              # ('poisson', 'lg', 'wag', 'jtt', 'dayhoff')

# evolve amino-acid sequences over a gene tree (nodes carry .gid and .children;
# `subst` maps each node to the substitution length of the branch ending at it)
class Node:
    def __init__(self, gid, children=()):
        self.gid, self.children = gid, list(children)

tips = [Node(f"t{i}") for i in range(4)]
root = Node("r", tips)
seqs = evolve_on_tree(root, {t: 0.3 for t in tips}, model,
                      np.random.default_rng(0), length=60)
assert set("".join(seqs[t.gid] for t in tips)) <= set(AMINO_ACIDS)   # genuine 20-AA protein
```

### Protein output

`sequence` writes the rescaled trees to `gene_trees/` (`<family>_complete_subst.nwk` and
`<family>_extant_subst.nwk`, branch lengths in substitutions per site) and, with `--subst-model`, one
protein alignment per family to `alignments/<family>.fasta` (FASTA headers are the standard
[node names](../contributing/conventions.md#naming); each record is a 20-letter amino-acid sequence).
`branch_rates.tsv` and `gene_family_speeds.tsv` record the per-lineage clock rate and per-family speed
multiplier, and `sequence.log` is the run manifest.

### Protein validation

- **Empirical models (LG / WAG / JTT / Dayhoff).** The stored stationary frequencies match the
  published `π` vectors to 1e-4, catching any transcription error a mere reversibility check would
  miss (`test_sequence_sim.py::test_empirical_aa_frequencies_match_published`).
- **Poisson.** Uniform stationary frequencies and equal off-diagonal rates hold exactly by
  construction (`test_sequence_sim.py::test_poisson_is_exact`).
- **Stationarity (LG).** A star tree with long branches recovers the model's stationary amino-acid
  frequencies to within 0.02 (`test_sequence_sim.py::test_protein_stationary_recovered_on_long_branch`).

### Protein references

- Le, S. Q. & Gascuel, O. (2008). An improved general amino acid replacement matrix.
  *Molecular Biology and Evolution* 25(7): 1307–1320.
- Whelan, S. & Goldman, N. (2001). A general empirical model of protein evolution derived from
  multiple protein families using a maximum-likelihood approach.
  *Molecular Biology and Evolution* 18(5): 691–699.
- Jones, D. T., Taylor, W. R. & Thornton, J. M. (1992). The rapid generation of mutation data
  matrices from protein sequences. *Computer Applications in the Biosciences* 8(3): 275–282.
- Dayhoff, M. O., Schwartz, R. M. & Orcutt, B. C. (1978). A model of evolutionary change in proteins.
  In *Atlas of Protein Sequence and Structure*, vol. 5, suppl. 3, pp. 345–352.

## Clock references

- Drummond, A. J., Ho, S. Y. W., Phillips, M. J. & Rambaut, A. (2006). Relaxed phylogenetics and
  dating with confidence. *PLoS Biology* 4(5): e88. *(uncorrelated relaxed clocks)*
- Thorne, J. L., Kishino, H. & Painter, I. S. (1998). Estimating the rate of evolution of the rate of
  molecular evolution. *Molecular Biology and Evolution* 15(12): 1647–1657. *(autocorrelated clock)*
- Lepage, T., Bryant, D., Philippe, H. & Lartillot, N. (2007). A general comparison of relaxed
  molecular clock models. *Molecular Biology and Evolution* 24(12): 2669–2680. *(CIR clock)*
