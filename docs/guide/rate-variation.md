# Rate variation (relaxed clock)

> **Reference:** see the [relaxed-clock](../models/relaxed-clocks.md) and
> [DNA-substitution](../models/dna-substitution.md) catalog pages.

Species trees and gene trees in ZOMBI2 are **timetrees** — branch lengths are time. To get
the branch lengths you would infer from *sequence evolution*, overlay a substitution rate
that varies across the tree — a **relaxed molecular clock**, which turns a chronogram into a
**phylogram**.

ZOMBI2 provides a whole family of these clocks in the `zombi2.sequences` namespace, all sharing one
`Clock` interface (`scale(tree, seed=...)` returns the phylogram): the strict clock (`StrictClock`),
the *uncorrelated* clocks (`UncorrelatedLogNormalClock`, `UncorrelatedGammaClock`, `WhiteNoiseClock`
— each branch draws an i.i.d. rate), and the *autocorrelated* clocks (`AutocorrelatedLogNormalClock`,
`CIRClock`, and `RateVariation` — a branch's rate is anchored to its parent's). This page focuses on
**`RateVariation`**, the discrete-bin, Markov-switching model from the GTDB archaea study; the rest of
the family is covered in the manual's *Relaxed molecular clocks* chapter and used exactly the same
way.

## The model

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

## Works on gene trees too

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
`SequenceEvolution`, below.

## Family sequence evolution

The gene × lineage clock. `RateVariation` scales **one tree in isolation**. But in a real dataset every gene shares
the same underlying lineage history: a fast-evolving clade is fast for *all* the genes
passing through it. `SequenceEvolution` builds that in. It converts a whole
[`simulate_genomes`](gene-families.md) result — all the reconciled gene trees at once —
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
discrete-bin clocks.

`sequence` replays `run/Events_trace.tsv`, rebuilds the reconciled gene trees, and writes
`run/gene_trees/<family>_extant_subst.nwk` (and `_complete_subst.nwk`), plus
`gene_family_speeds.tsv` and `branch_rates.tsv` recording the drawn `s_g` and `R_b` for
reproducibility. A lineage clock alone (any `--clock`, `--branch-speed`, or `--branch-bins`) is a
shared lineage clock; `--family-speed` alone is per-family speed; together they are the full
gene × lineage model.

### Simulating alignments (DNA and protein)

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
