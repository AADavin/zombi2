# Sequence evolution

Species trees and gene trees in ZOMBI2 are **timetrees**: every branch length is an amount of *time*.
Real sequence data does not measure time. What an aligner and a phylogenetic reconstruction see is
the number of substitutions accumulated along each branch, and that number depends not only on how
long a branch lasted but on how fast it was evolving. A slowly evolving genome and a fast one that
lived for a quarter as long can look identical in an alignment.

Sequence evolution is the step that closes this gap. It overlays a substitution rate that varies
across the tree and rescales every branch from time into substitutions per site, turning a
**chronogram** (branch lengths in time) into a **phylogram** (branch lengths in substitutions). This
chapter covers the shared gene $\times$ lineage clock and the `zombi2 sequences` command that applies
it to a completed gene-family simulation.

## The problem

A gene tree that comes out of a `genomes` run has branch lengths in time, inherited from the species
timetree it was reconciled onto. To obtain the branch lengths you would actually *infer from
sequences*, each branch must be multiplied by a rate. If that rate were a single constant everywhere,
the tree would just be uniformly stretched: a **strict molecular clock**, and the relative branch
proportions would be unchanged. Real lineages violate the strict clock. Rates drift over time, some
clades evolve faster than others, and different gene families evolve at systematically different
speeds even within the same lineage.

The model therefore needs two things at once: rate variation *across the tree* that is shared by all
the genes living in a given lineage, and rate variation *across families* that is intrinsic to each
gene. Both are multiplicative overlays on the timetree.

![From time to substitutions. The same tree drawn as a **chronogram** (branch lengths in time — ultrametric, so every tip lines up at the present) and as a **phylogram** (branch lengths in substitutions per site, after multiplying each branch by an evolutionary rate). Both trees are painted by the same per-lineage rate: fast (yellow) branches that were short in time stretch out, slow (purple) ones collapse, and the tips no longer line up. Rescaling a chronogram into a phylogram this way is exactly what sequence evolution does.](figures/seq_chrono_phylo.pdf){width=100%}

## The model

Every gene family shares the same underlying lineage history. If a clade of the species tree is
fast-evolving, it is fast for *all* the families whose gene branches pass through it.
`SequenceEvolution` builds this in. It converts a whole `simulate_genomes` result — all the reconciled
gene trees at once — from time into substitutions per site under a **gene $\times$ lineage** clock:

$$\text{rate}(\text{family } g,\ \text{species branch } b) = R_b \cdot s_g$$

The rate factorises into two independent pieces.

- **$R_b$ — the shared lineage clock.** A rate is drawn **once on the species tree** and **shared by
  every family**. It captures the fact that a whole region of the tree can be fast or slow regardless
  of which gene you look at.
- **$s_g$ — the per-family speed.** Each family draws one intrinsic constant from a distribution, so
  some families are globally fast and others slow, independently of where in the tree they live.

A gene-tree branch that sits on species branch $b$ over the time interval $[t_0, t_1]$ receives
substitution length $s_g \cdot R_b \cdot (t_1 - t_0)$. Because reconciliation is exact, a gene branch
that spans several species branches (after pruning away losses) simply sums the contributions of the
pieces it crosses.

![The gene $\times$ lineage clock, $\text{rate} = R_b \cdot s_g$. **Left**, the shared lineage clock $R_b$ painted on the species tree — drawn once and shared by every family. **Right**, three gene families, each with its own speed $s_g$. Because $R_b$ is shared, the three phylograms have an identical colour pattern (the same clades run fast or slow); the family speed only sets each tree's overall length. Length carries $s_g$, colour carries $R_b$.](figures/seq_gene_lineage.pdf){width=100%}

### The shared lineage clock

The lineage clock can be **any** of ZOMBI2's relaxed molecular clocks — the full family (strict,
uncorrelated lognormal, gamma and white noise, autocorrelated lognormal, and Cox–Ingersoll–Ross) is
the subject of Chapter 15. Two are shown here; they are selected differently, and mutually exclusive:
the autocorrelated lognormal via `branch_sigma`, and any other clock (including the discrete-bin model
below) via `lineage=`.

**Autocorrelated lognormal relaxed clock.** The rate evolves down the species tree as a geometric
random walk,
$$R_{\text{child}} = R_{\text{parent}} \cdot \exp\!\big(\mathcal{N}(0,\ \sigma\sqrt{\ell})\big),$$
where $\ell$ is the branch length in time and $\sigma$ is the `branch_sigma` parameter. Because a
child's rate is anchored to its parent's, nearby lineages have similar rates — the clock is
*autocorrelated* [@thorne1998autocorrelated]. Setting $\sigma = 0$ freezes the walk and recovers a
strict clock, $R_b = 1$ everywhere.

**Discrete-bin within-branch clock.** The alternative is the discrete-bin, Markov-switching model
implemented by `RateVariation`. An **ordered** set of rate **bins** (multipliers, some faster than 1,
some slower) is laid down, and a continuous-time Markov process runs *along the phylogeny*, stepping
only to an **adjacent bin**. Because the rate can change gradually within a branch, a single branch
may be split into several **segments** in neighbouring bins; the substitution length of that branch is
the sum of `segment_duration` $\times$ `bin_rate` over its segments. Unlike the lognormal clock it can
vary the rate *within* a branch, not only at nodes. You pass a configured `RateVariation` — or any
other clock from the family (Chapter 15) — as the `lineage` argument.

## Usage from Python

`SequenceEvolution` takes a completed `simulate_genomes` result and returns a scaled object holding
one phylogram per family, together with the rates it drew.

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

To use the discrete-bin clock for the lineage part instead of the lognormal one, pass a
`RateVariation` as `lineage=` (mutually exclusive with `branch_sigma`); any other clock from the
family (Chapter 15) is passed the same way:

```python
rv = z.RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=1.0)  # slow -> fast
phylo = z.SequenceEvolution(lineage=rv,
                            family_speed=z.LogNormal(0.0, 0.4)).scale(genomes, seed=2)
```

Bin order matters: the walk only moves between neighbours, so the bins must be given in increasing (or
decreasing) order of rate.

## Usage from the CLI

Sequence evolution is a **separate command**, `zombi2 sequences`, run on a prior `genomes` result.
Decoupling the two means you can retune the clock without re-simulating gene content. The only
requirement is that the `genomes` run wrote its event trace — pass `trace` to `--write`:

```bash
zombi2 genomes -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace profiles -o run/

# lognormal lineage clock + per-family speed
zombi2 sequences --genomes run/ --branch-speed 0.4 --family-speed 0.5 -o run/

# discrete-bin lineage clock instead of --branch-speed
zombi2 sequences --genomes run/ --branch-bins 0.25,0.5,1,2,4 --branch-switch-rate 1.0 \
    --family-speed 0.5 -o run/
```

Here `--branch-speed` sets the lognormal $\sigma$ and `--family-speed` the spread of the per-family
speed distribution. The `--branch-bins` form supplies the discrete-bin clock's ordered multipliers,
with `--branch-switch-rate` the rate of stepping to a neighbouring bin.

`sequence` replays `run/events_trace.tsv`, rebuilds the reconciled gene trees from the recorded
events, and writes one phylogram per family to `run/gene_trees/<family>_extant_subst.nwk` (and
`_complete_subst.nwk` for the trees that include losses). It also writes `gene_family_speeds.tsv` and
`branch_rates.tsv`, recording the drawn $s_g$ and $R_b$ so a run is fully reproducible.

::: note
`--branch-speed` **or** `--branch-bins` alone is a shared lineage clock; `--family-speed` alone is
per-family speed; together they are the full gene $\times$ lineage model.
:::

::: tip
Because `sequence` reads only the event trace on disk, a single expensive `genomes` run can feed many
sequence-evolution runs. Sweep the clock parameters — different `--branch-speed`, different bins,
different `--family-speed` — reusing the same gene content each time.
:::

## Simulating sequence alignments

Rescaling turns each gene tree into a phylogram whose branch lengths *are* the expected number of
substitutions per site. The natural next step is to draw an **actual alignment** down that tree — a
sequence at every leaf, evolved base by base (or residue by residue) under a substitution model. Add
`--subst-model` to the `sequence` command and it does exactly that, after the rescale step:

```bash
# DNA alignments (HKY85) along the rescaled gene trees
zombi2 sequences --genomes run/ --branch-speed 0.4 --family-speed 0.5 \
    --subst-model hky85 --seq-length 600 -o run/

# protein alignments under the LG model, with +Γ rate heterogeneity across sites
zombi2 sequences --genomes run/ --branch-speed 0.4 \
    --subst-model lg --seq-length 300 --gamma-shape 0.5 -o run/
```

DNA versus protein is **auto-detected** from the model name. One FASTA per gene family is written to
`run/alignments/<family>.fasta`, containing the leaf sequences of that family's extant tree; each
record is headed by the same `<species>_<gene-id>` label the leaf carries in
`gene_trees/<family>_extant_subst.nwk`, so alignment and tree line up one-to-one. Omit `--subst-model`
and `sequence` behaves exactly as before — it only rescales the trees and writes no alignments.

![Substitution models over the four bases. **Top**, the nucleotide models as exchange graphs (purines A, G above; pyrimidines C, T below), with edge width the exchange rate and node area the stationary base frequency: across JC69, K80, HKY85 and GTR the model adds a transition/transversion bias, then unequal frequencies, then six free exchange rates. **Bottom**, one HKY85 alignment simulated down a small phylogram — sister tips share near-identical sequences and divergent clades differ, so the alignment tracks the tree.](figures/seq_subst_models.pdf){width=100%}

### Available models

A continuous-time Markov substitution model is a normalised rate matrix $Q$ (one expected
substitution per site per unit branch length) plus its stationary frequencies. The transition matrix
over a branch of length $t$ is $P(t) = e^{Qt}$; every model here is time-reversible, so $e^{Qt}$ is
computed by eigendecomposition of the symmetric matrix $B = \operatorname{diag}(\sqrt{\pi})\, Q\,
\operatorname{diag}(1/\sqrt{\pi})$.

**Nucleotide (4 states, `ACGT`).** The four nucleotide models are listed in Table \ref{tbl:ntmodels}.

| model    | description                                                 | parameters |
|:---------|:------------------------------------------------------------|:-----------------------------|
| `jc69`   | Jukes–Cantor: equal rates, equal base frequencies           | —          |
| `k80`    | Kimura 2-parameter: transition/transversion ratio `--kappa` | `--kappa`  |
| `hky85`  | HKY85: transition bias with unequal base frequencies        | `--kappa`, `--base-freqs` |
| `gtr`    | general time-reversible: 6 exchangeabilities + frequencies  | `--gtr-rates`, `--base-freqs` |

: The nucleotide substitution models (`--seq-model`), their character, and the parameters each takes. \label{tbl:ntmodels}

**Amino acid (20 states).** Empirical exchangeability matrices, transcribed byte-for-byte from the
reference PAML data files (Ziheng Yang), plus a parameter-free Poisson model (Table \ref{tbl:aamodels}).

| model     | description                                                |
|:----------|:-----------------------------------------------------------|
| `poisson` | equal exchangeabilities, uniform frequencies (protein F81) |
| `lg`      | LG [@le2008improved] — the modern default                  |
| `wag`     | WAG [@whelan2001general]                                   |
| `jtt`     | JTT [@jones1992rapid]                                      |
| `dayhoff` | Dayhoff [@dayhoff1978model] (PAML values)                  |

: The empirical amino-acid substitution models (`--seq-model`); all are parameter-free. \label{tbl:aamodels}

The protein models are empirical and take no parameters. `--gamma-shape ALPHA` adds discrete-Gamma
across-site rate heterogeneity (+Γ) for any model; `--seq-length N` sets the alignment length. To
start each family from a chosen root sequence rather than a random draw from the stationary
distribution, pass `--root-fasta FILE` — a FASTA keyed by family id, whose per-record length overrides
`--seq-length` for that family.

::: note
The rescaled branch lengths already carry the gene $\times$ lineage clock, so `--subst-model` needs no
extra rate: one unit of branch length is one expected substitution per site by construction.
:::
