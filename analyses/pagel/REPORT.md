# Can Pagel's test detect a feedback?

**What we test:** a genome and a trait shape each other in one joint ZOMBI2 run: a binary
habitat multiplies the loss rate of every gene family, and one named family's absence
multiplies the rate of switching into the cave. The standard test for correlated
evolution of two binary characters (Pagel 1994, fit with `phytools::fitPagel`) is applied
to the two tip characters the run produces. Does it detect the feedback? Does it detect
each direction alone? And does it stay quiet when there is nothing to find?

## The design

Four arms on matched seeds, 150 replicates each. Every replicate simulates a dated
species tree to 150 extant tips (birth 1.0, seed = the replicate number), shared across
the arms. On that tree one joint run grows the genome and the habitat together:

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
`gen2trait` only S = 12, and `null` has both at 1, which is the same machinery with the
dependencies multiplying by one. Every replicate also carries a **control** character
from an independent genome run on the same tree (one family, loss 0.2, transfer 0.25,
seed + 1000), connected to nothing, so any signal the test finds on it is an artifact of
the tree or the method.

## What the test reports

For every replicate we fit `fitPagel` (the dependent eight-parameter model against the
independent four-parameter one, a likelihood-ratio test at α = 0.05) on the habitat
paired with the eye family's tip presence, and on the habitat paired with the control:
1,200 fits, of which 341 were skipped because a character invariant at the tips cannot
be fit.

| arm | habitat × eye | habitat × control |
|---|---|---|
| feedback (both directions) | **90.6%** (87/96, CI 83.1–95.0%) | 6.0% (7/116) |
| gen2trait (family drives the switch) | **87.3%** (89/102, CI 79.4–92.4%) | 7.8% (9/116) |
| trait2gen (habitat drives loss) | **19.8%** (19/96, CI 13.1–28.9%) | 6.0% (7/116) |
| null (no dependency) | 5.0% (5/101, CI 2.1–11.1%) | 4.3% (5/116) |

- **Calibration is clean.** The null sits at the nominal level, and the control does too
  in every arm, including on trees whose gene content and habitat are genuinely
  dependent.
- **Power is asymmetric.** The feedback and the arm in which the family drives the
  habitat's switch rate are detected in about nine runs of ten. The arm in which the
  habitat drives loss, a five-fold change in the loss rate of every family in the
  genome, is detected in one run of five.
- **The reason is visible in the generator.** A connection produces extra events only
  where the driving state is occupied. Lineages missing the eye family are common,
  because copies are steadily lost, so the twelve-fold switch rate acts across much of
  the tree; cave lineages are rare at the base switch rates, so the five-fold loss rate
  acts on little of it. Tip presence is also a coarse readout of loss: a family present
  in several copies must lose them all before the character changes.
- **Detection is not direction.** The test's verdict is dependence or independence; the
  feedback and the one-way switch arm look alike in it.

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
it byte-identically, the seeds being the replicate numbers, by ZOMBI2's reproducibility
contract. `fits.tsv`, `results.json` and `manifest.json`, the run's outcomes, are
committed.

## Relation to the paper

This analysis is the worked example "Can Pagel's test detect a feedback?" in the ZOMBI2
manuscript, and `figures/pagel.pdf` is the manuscript's figure.
