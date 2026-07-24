# Can you trust RED? Validating a tree-rescaling method

**What we test:** whether **Relative Evolutionary Divergence (RED)** — the measure GTDB uses to turn
a phylogram into a relative divergence scale — recovers a tree's divergence times once uneven
molecular rates are in play. **How:** RED cannot be checked on a real tree, because a real tree's
true divergence times are unknown. So we read how rate-variable real trees are off one model-free
number on a real phylogeny (the GTDB archaeal tree), simulate trees that variable — where the true
times are known by construction — and grade RED against them.

## The question

RED (Parks et al. 2018) gives every node a number between 0 (the root) and 1 (the tips) from the
branch lengths alone: walking root to tip, it places each node along its path in proportion to
accumulated branch length. GTDB uses it to normalise taxonomic ranks across the tree of life, so a
phylum sits at a comparable RED whatever lineage it belongs to. That works only if RED is a faithful
stand-in for relative divergence time. Under a strict clock it is exact: branch length is
proportional to time, so RED equals node age over root age. Under rate variation branch length is no
longer time. The question is whether RED still recovers the ages.

We check the exact case first. On a dated (ultrametric) tree, RED equals `node age / root age` to
machine precision — in this port, the largest deviation over a 400-tip tree is 2×10⁻¹⁶. So RED of
the dated tree is the ground truth, and RED of the rate-distorted phylogram is the estimate we grade
against it.

## Why we cannot check on real data

On a real tree the answer is hidden. A phylogram measures substitutions, which are rate × time; from
substitutions alone rate and time are jointly unidentifiable, so the true node ages we would grade
RED against cannot be read off. Dating the tree first would assume a rate model — the thing in
question — and make the test circular.

Simulation avoids this. In a simulated tree the true node ages are known, so RED can be graded
exactly. The only thing we take from the real world is how much rate variation to put in. That
single number is all we borrow from GTDB; no real branch length enters the forward model.

## The observable — how ragged are real trees?

Every genome in the GTDB archaeal tree — 10,122 genomes, branch lengths in substitutions per site
(Parks et al. 2018; Rinke et al. 2021) — is extant: all tips sit at the present. Every tip is
therefore the same amount of time from the root. Any spread in root-to-tip substitutions can only
come from rate variation. We summarise it with the coefficient of variation (CV) of root-to-tip
substitution distances. A strict clock gives CV = 0; heterogeneity spreads the tips out.

![Figure 1](figures/observable.png)

*Figure 1. The distribution of root-to-tip substitution distances across the 10,122 genomes of the
GTDB archaeal tree. The dashed line is the mean; the CV is 0.232 — the fastest lineage has
accumulated roughly four times the substitutions of the slowest. This is the raggedness the test
must reproduce to be realistic. It is read straight off the phylogram — no dating, no
ultrametricising, no rate model assumed; its one assumption is that the tree is correctly rooted,
which GTDB provides.*

## Calibrating the clock

We build the simulated trees the way ZOMBI2 builds any tree, and calibrate their raggedness to the
real number. We simulate 8 species trees under the Yule process (400 tips, known node ages), evolve
each under ZOMBI2's relaxed lineage clock — `substitution = 1.0 * ByLineage(spread=σ)` (Drummond et
al. 2006) — and read the `species_phylogram`, whose branch lengths are substitutions. Sweeping σ and
finding where the mean root-to-tip CV crosses 0.232 gives the σ that makes a simulated tree as ragged
as real archaea (Figure 2): σ ≈ 0.54 for the lognormal tail, σ ≈ 0.59 for the gamma tail.

![Figure 2](figures/clock_recovery.png)

*Figure 2. Root-to-tip substitution CV as the clock heterogeneity σ grows, for the two tails ZOMBI2
ships for its uncorrelated lineage clock (`ByLineage`, `dist="lognormal"` and `dist="gamma"`). Where
a curve crosses the GTDB target (CV = 0.232, dashed) is the σ that reproduces real raggedness, marked
with a filled circle. Shaded bands are ±1 s.d. across the 8 Yule trees. The amount of variation (CV)
is identifiable; the σ read off is tail-dependent, so both tails are carried into the test.*

The forward model contains no real branch lengths: the GTDB tree sets only the target CV.

## Testing RED

The test itself: on the same simulated trees, apply the relaxed clock, compute RED from the resulting
substitution branch lengths, and compare RED's recovered node ages to the truth. Plotting RED's
accuracy against root-to-tip CV — the quantity measured on GTDB — puts the whole recipe on one axis:
read up from the real value (CV = 0.232) to find how well RED does at realistic raggedness (Figure 3).
Figure 4 is the per-tree picture.

![Figure 3](figures/red_bridge.png)

*Figure 3. RED accuracy (left, Pearson r between RED and true relative age) and error (right, nRMSE as
a percentage of tree depth) against how rate-variable the tree is, one curve per clock tail. The
dashed line is real archaea (CV = 0.232); filled markers read off RED's accuracy there. RED is
near-exact for mild variation and degrades as trees get raggeder. Shaded bands are ±1 s.d. across the
8 Yule trees.*

![Figure 4](figures/red_scatter.png)

*Figure 4. RED-recovered versus true node ages on one 500-tip tree under the lognormal clock, at three
raggedness levels: below, at, and above the real archaeal value. At CV = 0.22 (centre) RED still
tracks the diagonal (r = 0.92); the scatter opens up well beyond real heterogeneity.*

## What it means

**At the raggedness real archaea show, RED holds up.** Reading off CV = 0.232, RED recovers node ages
with Pearson **r = 0.95 (lognormal tail) and r = 0.94 (gamma tail)**, and **nRMSE ≈ 6% of tree
depth**. The amount of rate variation in real archaea is not enough to break RED — a quantitative
version of the assumption GTDB relies on (Rinke et al. 2021).

Three things are worth being precise about.

- **RED is an ordinal proxy, not exact ages.** Even at the best fit there is a few-percent age error.
  Use RED to order divergences and normalise ranks — its designed job — not to read absolute times off.
- **RED only breaks down past real data.** Beyond CV ≈ 0.23 the accuracy falls away: by CV ≈ 0.46,
  r drops to ≈ 0.79 and nRMSE rises to ≈ 9%. Real archaea sit on the safe side of that.
- **The CV pins the amount of variation, not its structure.** The identifiable quantity is how much
  rate variation there is, not whether it is autocorrelated (neighbouring lineages evolving at similar
  rates). See the limitations below: the clean core's sequence-level clock is uncorrelated only, so
  this port cannot vary that structure — the two tails it does sweep agree closely (r within 0.01),
  which is the residual this single observable leaves open.

> **Calibrate realism, then test against known truth.** A method like RED cannot be graded on the data
> it is meant for, because that data hides the answer. But one honest number can say how demanding the
> real case is (here, CV = 0.232), that number can be reproduced in a simulation where the answer is
> known, and the method graded there. This is the same move as the [synteny recipe](../synteny_inversions/REPORT.md),
> one level up: a summary pins one thing cleanly and hands what it cannot pin to a stated modelling choice.

## Assumptions and limitations

- **Model-free observable, rooted tree.** The root-to-tip CV assumes the tree is correctly rooted.
- **Uncorrelated clock only.** The clean core wires a single sequence-level clock, the uncorrelated
  `ByLineage` (`sequences/__init__.py`, `WIRED_MODIFIERS = (ByLineage,)`), with a lognormal or gamma
  tail. The autocorrelated clock (`FromParent`) is a species-level modifier and is rejected at the
  sequence level, so — unlike the retired six-clock version of this recipe — this port cannot test how
  RED responds to autocorrelated rate variation. It reports the uncorrelated case, which is the harder
  one for RED: autocorrelation preserves local order and would only help.
- **One domain.** This is archaea (GTDB). Bacteria, or a dated eukaryote phylogeny, would each set
  their own CV and could land on a different part of the RED curve.

## Reproducing this recipe

```bash
cd analyses/red
python observable.py     # the observable: GTDB root-to-tip substitution CV (= 0.232)
python experiment.py     # calibrate the clock, then grade RED vs raggedness -> results.json
python figures.py        # Figures 1-4 from results.json
```

The GTDB archaeal reference tree ships with the recipe (`data/ar53.tree`); refresh it with
`curl -fsSL -o data/ar53.tree https://data.gtdb.ecogenomic.org/releases/latest/ar53.tree`. The relaxed
clock is ZOMBI2's `zombi2.sequences` `ByLineage` modifier. The RED estimator is the shipped
`zombi2.tree.relative_evolutionary_divergence` (also on the CLI as `zombi2 tools tree --red`), which is
exact on an ultrametric tree to machine precision. Every number here regenerates from fixed seeds.

## References

- Drummond AJ, Ho SYW, Phillips MJ, Rambaut A (2006). *Relaxed phylogenetics and dating with
  confidence.* PLoS Biology 4:e88.
- Parks DH, Chuvochina M, Waite DW, et al. (2018). *A standardized bacterial taxonomy based on genome
  phylogeny substantially revises the tree of life.* Nature Biotechnology 36:996–1004.
- Rinke C, Chuvochina M, Mussig AJ, et al. (2021). *A standardized archaeal taxonomy for the Genome
  Taxonomy Database.* Nature Microbiology 6:946–959.
- Thorne JL, Kishino H, Painter IS (1998). *Estimating the rate of evolution of the rate of molecular
  evolution.* Molecular Biology and Evolution 15(12):1647–1657.
