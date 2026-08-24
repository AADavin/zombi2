# Who trades genes with whom?

**What we test:** gene transfers happen preferentially between lineages that share a
habitat, so an inferred transfer network also records which lineages lived in the same
habitat. Can the habitats of ancestral lineages be recovered from that record, using a
standard reconciliation program on realistic data? Only a simulation can answer it,
because only a simulation knows the true ancestral habitats.

## The design

One conditioned run. A two-state habitat is simulated first, on a complete tree of 100
extant and 75 extinct lineages (birth 1.0, death 0.5, tree seed 11; switch rate 0.05
each way, trait seed 1). The switches are rare, so each habitat covers whole clades:
14 switches, and 34 of the 100 extant tips end in habitat B. One thousand gene families
then evolve along the same tree under constant rates (duplication 0.02, transfer 0.05,
loss 0.15 per copy, seed 7), with one connection:

```python
transfer_to = Recipients().weighted_by(habitat,
                  Between({("A", "A"): 10.0, ("B", "B"): 10.0}, default=1.0))
```

The recipient of every transfer is drawn with a tenfold preference for lineages sharing
the donor's habitat at that moment. With this preference, 90.1% of the 4,750 transfers
connect lineages in the same habitat.

## The ceiling, before any tool runs

The observed data are the extant species tree and the extant gene trees. Because the run
records the complete history, the limit of any inference can be computed first:

| what | share |
|---|---|
| true transfers | 4,750 |
| leave any extant trace | 67.2% |
| detectable in principle (trace and donor context in the gene tree) | 66.5% |
| dead donor among detectable | 22.5% |
| dead recipient among detectable | 0.7% |

A third of the transfers left no trace, because every descendant of the transferred
copy died: no method can find them. A transfer from a donor with no surviving
descendants can at best be assigned to the extant branch from which the donor's lineage
diverged. After extinction and this reassignment, the same-habitat share of the
detectable transfers is 86.6%. That share is the ceiling: no method can exceed it.

## What ALE recovers

The true gene tree of every family in at least four species (718 families) was
reconciled against the extant species tree with the ALE program `ALEml_undated`. Feeding
ALE the true gene trees is deliberate: it isolates reconciliation error from gene tree
error.

- Of the transfers inferred with frequency above 0.5, 88% are true, and they recover
  53% of the detectable transfers. Lowering the threshold to 0.1 raises recall to 72%
  at 45% precision.
- The total inferred count is nearly unbiased: 3,167 frequency-weighted transfers
  against 3,160 detectable.
- The same-habitat share of the inferred network is 76.1%, against the 86.6% ceiling
  and a random-mixing baseline near 55%.

## Classifying the ancestral habitats

Each tip keeps its observed habitat, and each inferred transfer contributes one vote to
each of the two branches it connects; the vote is for the partner's habitat, and it is
weighted by the transfer's frequency. The votes assign a habitat to all 99 ancestral
branches, and 92% of the assignments are correct. The same procedure on the error-free
network built from the truth is also 92% correct, so the reconciliation's errors do not
lower the accuracy.

Two honest limits:

- **Depth.** The errors concentrate where lineages are few: the share of contemporaneous
  lineages correctly labelled rises from 50% near the root to nearly 100% at the
  present. The deepest branches are misclassified even on the error-free network.
- **Switches inside branches.** An undated reconciliation cannot place a habitat switch
  inside a branch. Seven branches switched habitat within their span; on them, part of
  the branch is mislabelled whichever habitat is assigned. The vote shares do not flag
  these branches either: switched branches average a minority-vote share of 0.32
  against 0.25 for unswitched ones, which does not separate them.

One further observation for practitioners: letting classified branches vote in later
rounds (label propagation) helps on the error-free network but lowers the accuracy on
ALE's network from 92% to 86%, because wrong calls then propagate. On inferred
networks, the tip-anchored votes are the ones to trust.

## The figure

`figures/transfers.png`, the manuscript's figure, from `figures.py`: **A** the complete
tree, coloured by the true habitat, extinct lineages included; **B** the extant tree
with the same true history, changing colour within a branch at the instant of each
switch; **C** each branch coloured by the habitat of its inferred transfer partners,
with pale branches where the evidence is little or mixed; **D** the transfer counts by
habitat, realized in the run and inferred by ALE; **E** each branch's assigned habitat
against the true history, cut at the true switch points, with the share of
contemporaneous lineages correctly labelled on the same time axis.

## Reproducing it

The Python side needs ZOMBI2 (this run: 0.43.2), numpy, and, for the figure, matplotlib,
Pillow and Phylustrator. The reconciliations need ALE (the runs here used the `ale.sif`
container on ETH Euler; `cluster/README.md` has the folder layout and the two batch
scripts). Then, in this directory:

```bash
python experiment.py            # 1000 families -> data/: gene trees + truth tables
sbatch cluster/01_observe.slurm # on the cluster: 718 .ale files
sbatch cluster/02_undated.slurm # on the cluster: 718 reconciliations -> recs/
python oracle.py                # the visibility and classification ceilings
python score.py                 # ALE against the truth
python figures.py               # -> figures/transfers.png
```

`data/` and `recs/` are not committed: `experiment.py` regenerates the data
byte-identically from the seeds, and the reconciliations regenerate from the two batch
scripts. The numbers above are from a single simulated world; the pipeline reruns
end-to-end from the seeds.

## Relation to the paper

This analysis is the worked example "Who trades genes with whom?" in the ZOMBI2
manuscript, and `figures/transfers.png` is the manuscript's figure.
