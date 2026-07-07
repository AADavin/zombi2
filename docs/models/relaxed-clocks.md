# Relaxed clocks

> **Tutorial:** see the [Rate variation](../guide/rate-variation.md) guide.

ZOMBI2 species trees and gene trees are **timetrees** — every branch length is an amount of *time*.
Sequence evolution accumulates *substitutions per site*, which is time multiplied by an evolutionary
**rate** that varies across the tree. A relaxed molecular clock is the model of that variation: it
draws a rate for each branch and rescales the tree from time into expected substitutions, turning a
**chronogram** into a **phylogram**. The *strict* clock is the baseline — one global rate, no
variation. A *relaxed* clock lets the rate change branch to branch, either **uncorrelated** (i.i.d.
draws, a branch's rate says nothing about its neighbours') or **autocorrelated** (a branch's rate is
anchored to its parent's, so nearby lineages evolve at similar rates). All seven share one `Clock`
interface — `scale(tree, seed=...)` → a phylogram — differing only in how the per-branch rate is drawn.

| Model | Family | Reach for it when |
| --- | --- | --- |
| **StrictClock** | baseline | you want a null: one constant `rate` on every branch |
| **UCLN** | uncorrelated | lineage-independent lognormal rate noise, `E[rate] = mean` |
| **UGAM** | uncorrelated | the same idea with a gamma spread tuned by a single `shape` |
| **WhiteNoise** | uncorrelated | short branches should vary most (variance ∝ 1/duration) |
| **AutocorrelatedLogNormal** | autocorrelated | rate heritability matters — a geometric walk down lineages |
| **CIR** | autocorrelated | mean-reverting rates that also vary *within* a branch |
| **RateVariation** | autocorrelated | discrete ordered rate **bins** with a nearest-bin walk (GTDB) |

## The models

### StrictClock

A single rate `rate` on every branch (no rate variation): the phylogram is the chronogram uniformly
stretched by `rate`, so relative branch proportions are unchanged. This is the baseline every relaxed
clock relaxes; reach for it as a null.

### Uncorrelated lognormal (UCLN)

Each branch draws an **independent** lognormal multiplier with `E[rate] = mean` for any spread:
`rate = mean · exp(𝒩(−σ²/2, σ))`. Larger `sigma` means more rate heterogeneity; `sigma = 0` is the
strict clock. Because the draws are i.i.d., a branch's rate is uninformative about its neighbours'
(Drummond et al. 2006).

### Uncorrelated gamma (UGAM)

Each branch draws an **independent** gamma rate with mean `mean` and variance `mean²/shape`. The
single `shape` knob controls dispersion: large `shape` concentrates rates near `mean` (→ strict),
small `shape` spreads them widely (Drummond et al. 2006; PhyloBayes `-ugam`).

### WhiteNoise

An uncorrelated clock whose branch multiplier is the integral of a white-noise rate over the branch:
gamma-distributed with mean `mean` and variance `mean²·σ²/Δt`, inversely proportional to branch
duration `Δt`. Long branches average the noise away (rate → `mean`); short branches are highly
variable. That branch-length dependence is what distinguishes it from UGAM; `sigma = 0` is strict
(PhyloBayes `-wn`).

### AutocorrelatedLogNormal

The rate evolves down the tree as a geometric random walk anchored to the parent,
`R_child = R_parent · exp(𝒩(0, σ·√ℓ))` with `ℓ` the branch length in time. A child's rate is centred
on its parent's, so nearby lineages have similar rates. `sigma = 0` freezes the walk into a strict
clock at `root_rate`. This is the shorthand `--branch-speed` clock (Thorne, Kishino & Painter 1998).

### CIR

The instantaneous rate follows a mean-reverting Cox–Ingersoll–Ross diffusion,
`dr = θ(μ − r)dt + σ√r dW`, which stays strictly positive and pulls back toward the long-run mean
`mean` (μ) at speed `theta` (θ), with volatility `sigma` (σ). The path is simulated *within* each
branch by Euler–Maruyama on sub-steps of at most `max_step`, so — unlike the lognormal walk — the rate
also varies within a branch, while a child still starts where its parent ended (Lepage et al. 2007).

### RateVariation

The discrete-bin, within-branch clock used in the GTDB archaea study. An **ordered** set of positive
rate `bins` is laid down and a continuous-time Markov process runs along the phylogeny, stepping only
to an **adjacent** bin (index ± 1) at `switch_rate`, with `up_bias` the probability a step goes to the
faster neighbour. Because the rate changes gradually, a single branch may split into several segments
in neighbouring bins; `switch_rate = 0` freezes it in its `start` bin (a strict clock).

## Command line

Sequence evolution is a **separate command**, `zombi2 sequence`, run on a prior `genomes` output that
was written with the event trace (`trace` in `--write`), so you can retune the clock without
re-simulating gene content. Pick the clock with `--clock`:

```bash
# a 'genomes' run carrying the trace 'sequence' replays
zombi2 genomes -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace profiles -o run/

# any clock in the family via --clock (here an uncorrelated lognormal)
zombi2 sequence --genomes run/ --clock uncorrelated-lognormal --clock-sigma 0.5 \
    --family-speed 0.5 --seed 7 -o run/
```

`--clock` selects any of `{strict, autocorrelated-lognormal, uncorrelated-lognormal,
uncorrelated-gamma, white-noise, cir, discrete-bin}`, parameterised by `--clock-mean` (the
mean/strict/root rate, default `1.0`), `--clock-sigma` (lognormal / white-noise / CIR spread, default
`0.5`), `--clock-shape` (gamma, default `3.0`), or `--clock-theta` (CIR mean-reversion, default
`1.0`). Two clocks keep shorthands, used when `--clock` is omitted: `--branch-speed SIGMA` for the
autocorrelated lognormal, and `--branch-bins R1,R2,…` (with `--branch-switch-rate`, `--branch-up-bias`)
for the discrete-bin GTDB clock.

```bash
zombi2 sequence --genomes run/ --branch-speed 0.4 --family-speed 0.5 -o run/          # autocorrelated lognormal
zombi2 sequence --genomes run/ --branch-bins 0.25,0.5,1,2,4 --branch-switch-rate 1.0 -o run/  # discrete-bin
```

`--family-speed SIGMA` overlays a per-family intrinsic speed on top of the shared lineage clock (each
family draws a constant multiplier ~ LogNormal(0, SIGMA)); give `--subst-model` to also evolve a DNA
or protein alignment down each rescaled tree.

## Python

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

## Output

`zombi2 sequence` writes, into the output directory: `gene_trees/` — the rescaled substitution-unit
gene trees, `<family>_complete_subst.nwk` and `<family>_extant_subst.nwk` per family;
`branch_rates.tsv` — the per-species-branch clock rate applied to the shared lineage clock;
`gene_family_speeds.tsv` — each family's intrinsic `--family-speed` multiplier; and `sequence.log`,
the run manifest. Adding `--subst-model` also writes an `alignments/` directory of simulated
DNA/protein alignments, one per family.

## Validation

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

## References

- Drummond, A. J., Ho, S. Y. W., Phillips, M. J. & Rambaut, A. (2006). Relaxed phylogenetics and
  dating with confidence. *PLoS Biology* 4(5): e88. *(uncorrelated relaxed clocks)*
- Thorne, J. L., Kishino, H. & Painter, I. S. (1998). Estimating the rate of evolution of the rate of
  molecular evolution. *Molecular Biology and Evolution* 15(12): 1647–1657. *(autocorrelated clock)*
- Lepage, T., Bryant, D., Philippe, H. & Lartillot, N. (2007). A general comparison of relaxed
  molecular clock models. *Molecular Biology and Evolution* 24(12): 2669–2680. *(CIR clock)*
