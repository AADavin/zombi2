# Can Pagel's test detect a feedback?

**What we test:** a genome and a trait shape each other in one joint ZOMBI2 run: a binary
habitat multiplies the loss rate of every gene family, and the eye family's absence
multiplies the rate of switching into the cave. The standard test for correlated
evolution of two binary characters (Pagel 1994, fit with `phytools::fitPagel`) is applied
to the two tip characters the run produces. Does it detect the feedback? Does it detect
each direction alone? And does it reject at the nominal rate when there is
nothing to find?

## The design

Four arms on matched seeds, 150 replicates each. Every replicate simulates a dated
species tree to 150 extant tips (birth 1.0, seed = the replicate number); the tree is
shared across the arms. On that tree one joint run simulates the genome and the habitat
together:

```python
r = joint.simulate(
    genome(duplication=0.05, origination=8.0, initial_families=40,
           loss=PerCopy(0.25).scaled_by("trait", {"cave": L, "surface": 1.0}),
           families=[family("eye")]),
    traits.discrete(states=["surface", "cave"], start="surface",
                    switch={"surface->cave": PerLineage(0.08).scaled_by(
                                "genomes:eye", {"present": 1.0, "absent": S}),
                            "cave->surface": 0.10}),
    tree=ct, seed=seed)
```

The arms set the two factors: `feedback` has L = 5 and S = 12, `trait2gen` only L = 5,
`gen2trait` only S = 12, and `null` has both at 1, so the connections are written but
multiply by one. Every replicate also carries a **control** character, connected
to nothing: one family from a separate genome run on the same tree (loss 0.2,
transfer 0.25, seed + 1000). Any signal the test finds on the control is an artifact of
the tree or the method.

## What the test reports

For every replicate we fit `fitPagel` twice: the habitat against the eye family's tip
presence, and the habitat against the control. The test compares the dependent
eight-parameter Markov model against the independent four-parameter one with a
likelihood-ratio test at α = 0.05. That is 1,200 fits; 341 were skipped because a
character invariant at the tips cannot be fit.

| arm | habitat × eye | habitat × control |
|---|---|---|
| feedback (both directions) | **90.6%** (87/96, CI 83.1–95.0%) | 6.0% (7/116) |
| gen2trait (family drives the switch) | **87.3%** (89/102, CI 79.4–92.4%) | 7.8% (9/116) |
| trait2gen (habitat drives loss) | **19.8%** (19/96, CI 13.1–28.9%) | 6.0% (7/116) |
| null (no dependency) | 5.0% (5/101, CI 2.1–11.1%) | 4.3% (5/116) |

- **Calibration is correct.** The null rejects at the nominal rate, and so does the
  control in every arm, including on trees whose gene content and habitat are genuinely
  dependent.
- **Power is asymmetric.** In the feedback arm and in the switch-rate arm the test
  detects the dependency in about nine runs of ten. In the loss arm it detects the
  dependency in one run of five, even though that dependency is a five-fold change in
  the loss rate of every family in the genome.
- **The reason is mechanical.** A connection produces extra events only where the
  driving state is present. Lineages missing the eye family are common, because copies
  are steadily lost, so the twelve-fold switch rate applies across much of the tree;
  cave lineages are rare at the base switch rates, so the five-fold loss applies to few
  branches. Tip presence is also a coarse measure of the loss rate: a family present in
  several copies must lose them all before the character changes.
- **The test reports dependence, not a direction.** The feedback arm and the
  switch-rate arm cannot be told apart from the test result.

![Rejection rates by arm](figures/pagel.png)

## Reproducing it

The Python side needs only ZOMBI2 (this run: 0.43.2) and numpy; matplotlib for the
figure. The fits need R with `phytools` (the `bisse` conda env from `analyses/bisse`
plus `install.packages("phytools")`). Then, in this directory:

```bash
python experiment.py        # ~15 min: 600 joint runs + 600 control runs -> data/
Rscript fit_pagel.R         # ~30 min on 8 cores -> fits.tsv
python aggregate.py         # -> results.json + the table on stdout
python figures.py           # -> figures/pagel.{png,pdf}
```

`data/` (600 trees and tip-state tables) is not committed: `experiment.py` regenerates
it byte-identically, because a ZOMBI2 run is determined by its seed and the seeds here
are the replicate numbers. The run's outcomes are committed: `fits.tsv`, `results.json` and
`manifest.json`.

## Relation to the paper

This analysis is the worked example "Can Pagel's test detect a feedback?" in the ZOMBI2
manuscript, and `figures/pagel.pdf` is the manuscript's figure.
