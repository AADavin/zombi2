# What can gene order constrain?

**What we test:** whether the parameters of a genome rearrangement model can be recovered from gene
order alone. Genomes are simulated at known rates down a dated tree, the parameters are then
inferred back from the simulated genomes, and the answer is compared to the values that generated
the data. Everything here is simulated, so every claim is checked against a known truth.

**Verdict.** Rates are recoverable, event extent is not. In a mixed model of inversions and
translocations, both rates are recovered exactly, and so is the share of events that are
translocations. The mean extent of an inversion is not recovered at all: its credible interval
covers the whole grid, from 1 gene to 64 genes, in every replicate. Fixing the extent at the wrong
value is not free. Setting it to twice the truth pulls the inversion rate 29% low and raises the
apparent translocation share from 0.20 to 0.26.

## The question

An inversion reverses a segment of a chromosome. A translocation moves a segment from one chromosome
to another. Both break the adjacencies between neighbouring genes, so gene order slowly diverges
between species, and the amount of divergence carries information about how often the events happen.
A dated tree turns that into a rate. Two things about each event type could in principle be read
from the data: how often it happens, and how large it is. This study asks which of the two the data
actually pin, and whether two event types acting together can be told apart.

The question is asked on simulated data on purpose. Given a real pair of genomes there is no way to
check an inferred rate, because the true rate is unknown. Simulating the data first makes the answer
available, so the method can be graded rather than trusted.

## The design

**The tree.** A birth-death tree (birth 1.0, death 0.3) grown to 40 extant tips, then rescaled to a
crown depth of 100 time units. The tree is ultrametric, so the patristic distance between two tips
runs from 0 to 200 time units. It gives 780 tip pairs. Pairs are grouped into three divergence bins,
0 to 60, 60 to 120 and 120 to 200 time units, holding 138, 199 and 443 pairs, with mean divergence
30.6, 80.6 and 167.1 time units.

**The genomes.** 2000 genes on 8 linear chromosomes, each gene 1000 bp, with 1.4 bp of chromosome
per bp of coding, so each chromosome is 350,000 bp long. All genomes are simulated, down the whole
tree, which yields one genome at every extant tip. That is what the statistics need, because they
are computed over pairs of tips. No DNA is simulated. Only the gene layout descends the tree.

**The truth.** The inversion rate is 1.0e-3 per gene per time unit and the mean inversion extent is
4 genes. In the two mixed arms the translocation rate is 2.5e-4 per gene per time unit, which makes
the true translocation share 0.20.

**The arms.**

| Arm | Data generated with | Grid searched | Question |
|---|---|---|---|
| A | inversions only | inversion rate x inversion extent, 13 x 7 = 91 cells | which inversion parameter is pinned |
| B | inversions and translocations | inversion rate x translocation rate, 9 x 6 = 54 cells, extent fixed at the truth | can two event types be told apart |
| C | the same data as arm B | the same grid, extent fixed at twice the truth | what does a wrong extent cost |

The inversion rate grid steps by a factor of the square root of 2, from 1.25e-4 to 8.0e-3 in arm A
and from 2.5e-4 to 4.0e-3 in arms B and C. The extent grid is 1, 2, 4, 8, 16, 32 and 64 genes. The
translocation rate grid is 0, 6.25e-5, 1.25e-4, 2.5e-4, 5.0e-4 and 1.0e-3. Every grid contains the
true value, so exact recovery is possible.

**Replication.** Each grid cell is simulated 3 times and the statistics are averaged. Each arm has 3
independent observed datasets, and each is fitted separately, so the three answers per arm are three
independent attempts at the same recovery. All seeds are fixed by the arm, the cell and the
replicate, and are recorded in `results.json`.

## The three summary statistics

**Gene-order conservation.** Consecutive genes on one chromosome form an adjacency. The breakpoint
distance between two genomes is the number of adjacencies present in one and absent in the other.
Gene-order conservation is its normalised complement, the fraction of adjacencies the two genomes
share. The two carry the same information; the fraction is used so that genomes with different
adjacency counts stay comparable.

**Conserved segment length.** The maximal runs of genes whose consecutive adjacencies all survive in
the other genome, reported as a mean run length in genes. This is the spacing between breakpoints,
which is the same set of breaks counted from the other side.

**Cross-chromosome breaks.** Of the adjacencies that are broken, the fraction whose two genes lie on
different chromosomes in the other genome. An inversion acts inside one chromosome, so it breaks
adjacencies but leaves both genes on the same chromosome. A translocation separates them. This is
the statistic that tells the two event types apart, and it is used only in the mixed arms.

Each statistic is averaged within each divergence bin, so the fitted vector is a statistic against
divergence time. What the observed data look like:

| Arm | Statistic | bin 0 to 60 | bin 60 to 120 | bin 120 to 200 |
|---|---|---|---|---|
| A | gene-order conservation | 0.960 | 0.899 | 0.802 |
| A | conserved segment length | 34.3 | 9.6 | 5.0 |
| A | cross-chromosome breaks | 0.000 | 0.000 | 0.000 |
| B | gene-order conservation | 0.940 | 0.847 | 0.726 |
| B | conserved segment length | 25.1 | 6.4 | 3.6 |
| B | cross-chromosome breaks | 0.241 | 0.282 | 0.305 |

Cross-chromosome breaks are exactly zero in arm A across all three replicates, which is the check
that the statistic reads translocations and nothing else. Inversions cannot produce a single one.

## Why ABC, and what else is available

For an inversion-only model the likelihood is not out of reach, and this study does not claim
otherwise. Caprara and Lancia (2000) analysed sorting by reversals experimentally and statistically.
York, Durrett and Nielsen (2002) built a Bayesian sampler that estimates the number of inversions
separating two chromosomes, and extended it in 2007 to inversions whose tract length varies. If the
question here were only about inversions, one of those would be the sharper tool, and approximate
Bayesian computation would need an excuse.

The case for approximate Bayesian computation is the mixed arm. Once translocations act alongside
inversions, the event model is mixed, and the distance and likelihood machinery developed for pure
reversal models no longer applies to it. Approximate Bayesian computation needs only two things: a
simulator that produces both event types, and a summary statistic that keeps its meaning under both.
ZOMBI2 supplies the first and the breakpoint distance supplies the second. That is why arm B is the
headline of this study, and why arm A is best read as a control that shows the method works in the
regime where a sharper tool already exists.

## Why the breakpoint distance

Not because it is cheaper, and not because it holds out longer. The breakpoint distance and the
inversion distance are both monotone in the number of events, both can be computed in time linear in
the number of genes, and of the two it is the inversion distance that saturates later, so it keeps
information for longer. The reason to use the breakpoint distance here is that it counts broken
adjacencies whatever broke them. It therefore means the same thing when translocations act alongside
inversions, which is this study's target regime, while the inversion distance is defined by a pure
reversal model.

## Results

![Figure 1](figures/rearrangements.png)

*Figure 1. Panel A, the misfit over arm A's grid of inversion rate against inversion extent. The
band of good fits is narrow in rate and runs the full height of the grid. Panel B, the misfit
against each parameter as a multiple of its true value, profiled over the other parameter. The rate
has a minimum at the truth and the extent is flat. Panel C, the misfit over arm B's grid of the two
rates, with a single minimum on the truth. Panel D, the recovered value as a multiple of the truth,
for all three observed replicates of every arm.*

**Arm A: the rate is pinned, the extent is not.** The recovered inversion rate is 7.07e-4 in all
three replicates, one grid step below the true 1.0e-3. The recovered extent is 32, 64 and 64 genes
against a true 4 genes, and its credible interval covers the entire grid, 1 to 64 genes, in all
three replicates. The clearest way to see the difference is to profile the misfit along each axis,
taking the best value over the other axis. Along the rate axis the misfit runs from 0.022 to 1.051, a
range of 1.03. Along the extent axis it runs from 0.022 to 0.055, a range of 0.033. The data respond
about 31 times more strongly to the rate than to the extent.

The reason is that an inversion breaks two adjacencies whatever its extent, so the number of events
sets how much gene order erodes and the size of each event barely matters. That produces the ridge in
panel A, and it is also why the rate lands one step low here: with the extent free, the best cell
sits somewhere along the ridge rather than on the truth. The bottom row of panel A is the exception.
A one-gene inversion often reverses a single gene without separating it from its neighbours, so it
erodes gene order less, and a higher rate is needed to match the data.

**Arm B: a mixed model's two rates can be told apart.** With the extent fixed at the truth, the
recovered inversion rate is 1.0e-3 and the recovered translocation rate is 2.5e-4 in all three
replicates, both exactly the true values. The recovered translocation share is 0.20 against a true
0.20. The credible interval is 5.0e-4 to 2.0e-3 for the inversion rate and 6.25e-5 to 5.0e-4 for the
translocation rate. The best inversion-only cell, meaning the best cell of the whole translocation
rate of zero row, reaches a misfit of 0.809, against 0.098 for the best mixed cell, a factor of 8.2.
So the data do not merely allow a translocation rate to be estimated, they reject the inversion-only
model outright.

**Arm C: a wrong extent has a price.** Fitting the same observed data with the extent fixed at 8
genes, twice the truth, moves the recovered inversion rate to 7.07e-4 in all three replicates, 29%
below the truth. The recovered translocation rate is unchanged at 2.5e-4, still exact. The recovered
translocation share therefore rises to 0.26 against a true 0.20. The practical reading is that the
extent has to be fixed from outside the data, because the data will not supply it, and that the
choice is visible in the inversion rate and in the translocation share but not in the translocation
rate.

## Assumptions and limitations

- **The credible interval is a rough one.** Grid cells are weighted by how closely they match, with
  the tolerance set at the tenth percentile of the misfits over the grid. That is a working choice,
  not a calibrated posterior, and the width of the interval depends partly on it and partly on the
  span of the grid. The recovery claims rest on the best cell and on the shape of the profiled
  misfit, not on the interval.
- **The point estimate cannot be finer than the grid.** The rate grid steps by a factor of the square
  root of 2, so the smallest deviation this study can resolve is 41% up or 29% down. Arm A's and arm
  C's rate estimates are one step low, which is the finest error the design can report.
- **Breakpoints never fall inside a gene.** This is a rule of ZOMBI2's nucleotide model, so genes are
  never split and every genome keeps all 2000 of them, one to one. The core gene set is therefore the
  whole gene set. No gene gain, loss or duplication is simulated.
- **The requested extent is not exactly the realised extent.** The engine snaps a requested arc to
  the nearest legal breakpoint, and the chromosome ends truncate long draws. Measured over one run,
  requests of 1, 2, 4 and 8 genes are realised at 1.0, 2.0, 4.1 and 8.0 genes, while requests of 32
  and 64 genes are realised at 27.4 and 42.0. The extent axis therefore compresses at its top end.
  This does not weaken the finding, since the axis still spans a 42-fold range of realised extent and
  the misfit barely responds to it.
- **Two event types, not the full set.** Only inversions and translocations act here.
  Transpositions, and rearrangements that change the number of chromosomes, are left out. The
  statistics used would extend to them, since the breakpoint distance is indifferent to what broke an
  adjacency, but this study does not test that.
- **One tree.** All arms use a single 40-tip tree, so the results do not say how the answer varies
  with tree shape, depth or tip count.

## Reproducing this study

```bash
cd analyses/rearrangements
python experiment.py    # 597 grid simulations + 6 observed, writes results.json and data/tree.nwk
python figures.py       # Figure 1 from results.json
```

The full run took 981 s on 9 worker processes, split as 478 s for arm A, 249 s for arm B and 246 s
for arm C. It needs no input data and no network access. The dated tree is regenerated from
`TREE_SEED` on every run and written to `data/tree.nwk`; every simulation seed is derived from
`MASTER_SEED` together with the arm, the cell and the replicate, so the run does not depend on how
many worker processes are used. `results.json` is committed and holds the provenance, the design, all
grid and observed summary statistics, the full misfit surfaces and the verdicts. Figures are
regenerated rather than committed.

Verified against ZOMBI2 0.42.2 at commit 98cb4f9.

## References

- Caprara A, Lancia G (2000). *Experimental and statistical analysis of sorting by reversals.* In:
  Sankoff D, Nadeau JH (eds), *Comparative Genomics: Empirical and Analytical Approaches to Gene
  Order Dynamics, Map Alignment and the Evolution of Gene Families.* Kluwer Academic Publishers,
  Dordrecht, 171-183. doi:10.1007/978-94-011-4309-7_16
- York TL, Durrett R, Nielsen R (2002). *Bayesian estimation of the number of inversions in the
  history of two chromosomes.* Journal of Computational Biology 9(6):805-818.
  doi:10.1089/10665270260518281
- York TL, Durrett R, Nielsen R (2007). *Dependence of paracentric inversion rate on tract length.*
  BMC Bioinformatics 8:115. doi:10.1186/1471-2105-8-115
