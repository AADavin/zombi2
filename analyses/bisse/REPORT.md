# Can BiSSE find the gene that drives speciation?

**What we test:** a named gene family multiplies a lineage's speciation rate while it is
present, in a joint ZOMBI2 run. Two questions. First, what does that dependency do to the
data an analysis would be handed? Second, does the standard test for state-dependent
diversification (BiSSE, fit with `diversitree` exactly as its documentation recommends)
detect the family that drives, and stay quiet on a family that does not?

## The design

One joint run grows the species tree and the genome together: a `driver` family, present
in the root genome and lost per copy at rate 0.12, multiplies the speciation rate
(base 1.0, extinction 0.3) by a factor *f* while present. We grew 200 replicate clades to
150 extant tips at each of six factors, *f* = 1, 1.25, 1.5, 2, 3, 5, on matched seeds.
At *f* = 1 the family multiplies speciation by one, which is to say not at all: those 200
runs are the null. Every genome also carries a second family, `control`, which appears at
the root and is lost at the same rate as the driver but which no rate reads, so it cannot
influence the tree; any rise in its prevalence can only come from the shape of the tree.
A separate duration-matched null arm reruns the *f* = 3 replicates undriven to the same
total time.

The written form of the dependency, in full:

```python
joint.simulate(
    species.birth_death(
        birth=PerLineage(1.0).scaled_by("genomes:driver",
                                        {"present": f, "absent": 1.0}),
        death=0.3, n_extant=150),
    genomes.genome(duplication=0.0, loss=0.12, origination=0.0,
                   families=[genomes.family("driver"), genomes.family("control")]),
    seed=seed)
```

## What the simulation shows (panel A)

The driver's prevalence among extant tips rises with the factor, from 0.41 in the null to
0.95 at five-fold. The control says most of that rise is not selection: carried through
the very same runs, it rises from 0.42 to 0.83 with nothing driving it, because a driven
clade reaches 150 tips sooner (mean tree height falls from 7.6 to 1.4 time units) and a
younger tree has lost less of everything. At *f* = 3 the decomposition is exact: of the
raw +0.50 gap between the driven run and its matched null (95% CI 0.46–0.55), +0.40 is
the shorter tree, shared with the control, and +0.11 (CI 0.08–0.14, sign-flip
p < 5×10⁻⁶) is carriers proliferating differentially. The duration-matched null agrees
with the in-run control (0.80 against 0.81), two independent estimates of the same
tree-age effect. The control also lands on its closed form exp(−loss × height) at every
factor, which checks the joint engine's genome half against something it does not know.

## What BiSSE reports (panel B)

For every replicate we fit the six-parameter BiSSE model and the state-independent
constraint (λ₁ = λ₀, μ₁ = μ₀) on each family's tip presence, a likelihood-ratio test on
two degrees of freedom at α = 0.05: 2,165 fits, 235 fits skipped because an
invariant character cannot be fit, no fit errors.

- **Calibration is clean.** The driver at the null is rejected in 6/171 fits (3.5%,
  Wilson 95% CI 1.6–7.4%). The control is rejected in 5.2% pooled over the driven trees
  (47/908, CI 3.9–6.8%), even though those trees carry real rate heterogeneity; it
  belongs to the other family, and BiSSE does not misattribute it here.
- **Direction is right.** Among significant fits at *f* > 1, 255 of 256 put the higher
  speciation rate on the carrier state.
- **Power is the limit, and it is not monotone.** A three-fold effect on speciation,
  enough to move prevalence by +0.50, is found in 36% of 150-tip clades. A five-fold
  effect is found *less* often, 30%, because a stronger driver pushes the family toward
  fixation: prevalence 0.95, trees of height 1.4, and few losses left to inform the
  likelihood. On data like these the risk is not the false positive; it is reading a
  non-significant result as the absence of the effect.

![The two panels](figures/bisse.png)

## Reproducing it

The Python side needs only ZOMBI2 (this run: 0.42.2, commit `c70d8dd`) and numpy;
matplotlib for the figure. The fits need R with `diversitree` (this run: 0.10.1). A
working recipe on macOS/arm64, where conda-forge lacks `r-subplex`:

```bash
conda create -y -n bisse -c conda-forge r-base r-ape r-desolve compilers make
conda run -n bisse Rscript -e 'install.packages(c("subplex","diversitree"),
    repos="https://cloud.r-project.org", type="source")'
```

Then, in this directory:

```bash
python experiment.py            # ~5 min: 1,400 joint runs -> results.json + data/
conda run -n bisse Rscript fit_bisse.R 8   # ~7 min on 8 cores -> fits.tsv
python aggregate.py             # -> verdicts.json + the table on stdout
python figures.py               # -> figures/bisse.{png,pdf}
```

`data/` (1,200 trees and tip-state tables, ~14 MB) is not committed: `experiment.py`
regenerates it byte-identically from the master seed recorded in `results.json`
(20260802), by ZOMBI2's reproducibility contract. `results.json`, `fits.tsv` and
`verdicts.json`, the run's outcomes, are committed.

## Relation to the paper

This analysis is the worked example "Can BiSSE find the gene that drives speciation?" in
the ZOMBI2 manuscript. The manuscript's figure is rendered from the same two JSON files
by the paper workspace's own script; `figures.py` here is a self-contained equivalent.

## Does a bigger tree fix the power? (panel C)

The 150-tip ceiling asks whether the limit is the method or the information. Rerunning
the *f* = 3 and *f* = 5 arms at 500 and 1,000 extant tips answers it (100 replicates
each, the same first hundred seeds; `size_experiment.py`, `fit_bisse_size.R`,
`aggregate_size.py`): power at *f* = 3 climbs from 36% through 96% to 99%, and at
*f* = 5 from 30% through 86% to 100%. Two details matter. The *f* = 5 deficit is not a
small-tree artifact: it persists at 500 tips (86% against 96%) and closes only at 1,000.
And calibration is not free at scale: the control, at the nominal level through 500
tips, is rejected in 8.2% of the 1,000-tip fits pooled over both factors (16 of 196,
Wilson CI 5.1-12.8%), consistent with the false-positive inflation reported for
state-dependent diversification on large trees. On these data, small trees hide true
drivers, and on the largest trees false positives begin to appear.

The size arms add 400 joint runs (about 30 minutes) and 800 fits (15 minutes at 500
tips and 28 at 1,000, on 8 cores). The outcome files are committed
(`results_size.json`, `fits_n500.tsv`, `fits_n1000.tsv`, `verdicts_size.json`); the
trees regenerate from the master seed as before.
