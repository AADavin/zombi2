# Rate variation (relaxed clock)

Species trees and gene trees in ZOMBI2 are **timetrees** — branch lengths are time. To get
the branch lengths you would infer from *sequence evolution*, overlay a substitution rate
that varies across the tree. `RateVariation` implements the discrete-bin, Markov-switching
model from the GTDB archaea study, turning a chronogram into a **phylogram**.

## The model

- An **ordered** set of rate **bins** — multipliers, some fast (> 1), some slow (< 1).
- A continuous-time Markov process runs **along the phylogeny** and can only step to an
  **adjacent bin** (index ± 1). Because the rate changes gradually, nearby lineages have
  similar rates — the process is *autocorrelated*. The current bin is inherited by both
  descendants at each node.
- A branch may thus be split into several **segments** in adjacent bins; its substitution
  length is `Σ (segment_duration × bin_rate)`.

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=30, age=8.0, seed=1)

rv = z.RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0],  # ordered slow -> fast
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
genomes = z.simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5, seed=1)
_, extant = genomes.gene_trees()["1"]
gene_tree = z.read_newick(extant)          # Newick -> Tree
phylogram = z.RateVariation(bins=[0.5, 2.0], switch_rate=1.0).scale(gene_tree, seed=1)
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
  family**. Two clocks to choose from: an autocorrelated **lognormal** relaxed clock,
  `R_child = R_parent · exp(N(0, σ·√branch_length))` (`branch_sigma`; `0` = strict clock), or
  the **discrete-bin** within-branch `RateVariation` (the GTDB model; `lineage=`), which can
  vary the rate *within* a branch.
- **`s_g`** — each family's intrinsic speed, one constant drawn from a distribution
  (`family_speed`), so some families are globally fast, others slow.

A gene-tree branch on species branch `b` over `[t0, t1]` gets substitution length
`s_g · R_b · (t1 − t0)`; because reconciliation is exact, a branch spanning several species
branches (after pruning) just sums the pieces.

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=30, age=8.0, seed=1)
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.2,
                             origination=0.5, seed=1)

se = z.SequenceEvolution(branch_sigma=0.5,                    # shared lineage clock
                         family_speed=z.LogNormal(0.0, 0.4))  # per-family speed
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
rv = z.RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=1.0)  # slow -> fast
phylo = z.SequenceEvolution(lineage=rv, family_speed=z.LogNormal(0.0, 0.4)).scale(genomes, seed=2)
```

### From the CLI

Sequence evolution is a **separate command**, `zombi2 sequence`, run on a prior `genomes`
result — so you can retune the clock without re-simulating gene content. `genomes` just has to
have written the event trace (`trace` in `--write`):

```bash
zombi2 genomes  -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace profiles -o run/

# lognormal lineage clock + per-family speed
zombi2 sequence --genomes run/ --branch-speed 0.4 --family-speed 0.5 -o run/

# discrete-bin (GTDB) lineage clock instead of --branch-speed
zombi2 sequence --genomes run/ --branch-bins 0.25,0.5,1,2,4 --branch-switch-rate 1.0 \
    --family-speed 0.5 -o run/
```

`sequence` replays `run/Events_trace.tsv`, rebuilds the reconciled gene trees, and writes
`run/gene_trees/<family>_extant_subst.nwk` (and `_complete_subst.nwk`), plus
`gene_family_speeds.tsv` and `branch_rates.tsv` recording the drawn `s_g` and `R_b` for
reproducibility. `--branch-speed` **or** `--branch-bins` alone is a shared lineage clock;
`--family-speed` alone is per-family speed; together they are the full gene × lineage model.

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
