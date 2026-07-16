# Coupling that shapes the tree

The couplings in this chapter all point **into S**: a trait or the gene content sets each lineage's
speciation and extinction rate, so the shape of the tree *depends* on the coupled state. That breaks
the pipeline — the tree can no longer be drawn first — and the run becomes forward-only, growing the
tree jointly with its driver as an **output** (give `--age`/`--tips`, no `-t`). Two drivers can shape
the tree — a trait (state-dependent diversification) or gene content (key innovation) — and each has a
joint model in which speciation feeds back on the driver.

## Trait-driven diversification (`traits:species`)

A discrete or continuous trait drives speciation and extinction, and the tree is grown *jointly* with
the trait — the SSE family. Because the trait shapes the topology, this edge **produces** the tree (no
`-t`) and takes a stopping condition (`n_tips` or `age`). The Python driver is `simulate_sse`:

```python
from zombi2.coevolve import simulate_sse, BiSSE

# BiSSE: state 1 speciates 3x faster, so it dominates the standing tips
res = simulate_sse(
    BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.2, q01=0.1, q10=0.1),
    n_tips=200, seed=1)
res.tree               # complete tree; prune() for the reconstructed one
res.labeled_values()   # the trait at the extant tips
```

![State-dependent diversification (BiSSE). **A**, the model: two states, the anagenetic transitions (arrow width = rate), and each state's speciation rate drawn as a fork (width = $\lambda$) — state 1 branches three times faster. **B**, one `simulate_sse` realization: branches are heavy where the lineage is in state 1 and light in state 0, lineages that go extinct end in a cross, and the extant tips carry chips. The fast (state-1) lineages proliferate and fill most of the standing tips — the diversification signal, written into the shape of the tree itself.](figures/sse.pdf){width=100%}

The four SSE variants are four **responses** of the single `traits:species` edge — they differ only
in the kind of trait doing the driving. `BiSSE` is the **binary** state-dependent birth–death process
[@maddison2007bisse]; `MuSSE` the **k-state** variant [@fitzjohn2012diversitree]; `QuaSSE` the
**continuous-trait** variant, whose speciation and extinction are functions of the trait value
[@fitzjohn2010quasse]; and `HiSSE` the **hidden-state** model [@beaulieu2016hisse], where
diversification is driven by an unobserved class rather than the focal trait — the honest null against
which a real association must be judged.

```python
import numpy as np
from zombi2.coevolve import MuSSE, QuaSSE, HiSSE

# k-state:
MuSSE(birth=[1, 3], death=[0.2, 0.2], Q=np.array([[-0.1, 0.1], [0.1, -0.1]]))
# continuous trait:
QuaSSE(speciation=lambda x: 1 + 2 / (1 + np.exp(-x)),
       extinction=lambda x: 0.2, sigma2=0.5, rate_bound=5.0, x0=0.0)
# hidden classes:
HiSSE(classes=[BiSSE(0.5, 0.7, 0.2, 0.2, 0.1, 0.1),
               BiSSE(2.0, 3.0, 0.2, 0.2, 0.1, 0.1)], hidden_transition=0.1)
```

From the command line, `--sse-model` picks the variant:

```bash
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/
```

It writes `species_tree.nwk` (the tree the trait's rates shaped), `traits.tsv` (every node — tips
*and* ancestral states) and `trait_tree.nwk`. `--sse-model musse` takes `--birth`/`--death` vectors
and a `--q-matrix` file; `--sse-model quasse` takes a sigmoidal speciation (`--spec-low/high/center/slope`)
and Brownian `--diffusion`; and `--sse-model hisse` adds `--hidden-classes` diversification regimes
spanning the base rates up to `--hidden-scale`× faster (`--hidden-switch` between them). **What it
recovers:** the fast-speciating state accumulates lineages, so it dominates the standing tips — the
diversification signal is written into the tree shape itself.

![The continuous variant, QuaSSE. **A**, the model: the speciation rate is a rising function $\lambda(x)$ of the trait while extinction is flat, and the trait itself diffuses by Brownian motion (the axis is tinted with the viridis ramp used to paint the tree). **B**, one realization, each branch painted by its trait value: the high-value (yellow) lineages branch faster and proliferate, while the low-value (blue) lineages stay sparse and go extinct — the same "fast state fills the tips" signal as BiSSE, now on a continuous character.](figures/sse_quasse.pdf){width=100%}

## Key-innovation diversification (`genomes:species`)

Gene content can drive the tree just as a trait can: a small panel of binary **driver** ("key
innovation") families whose *presence* sets each lineage's speciation and extinction rate. Drivers are
gained de novo (origination) and — the interesting part — by **transfer**, which is
*frequency-dependent*: a driver carried by more of the live population is donated more often, so it
spreads as the tree grows. That feedback is why the tree and the gene content must grow together; this
edge produces the tree. The neutral bulk genome does not affect diversification, so it is overlaid
afterward on the finished tree with the ordinary `genomes` machinery (exact under independent families,
not an approximation).

```python
from zombi2.coevolve import simulate_gene_diversification, GeneDiversification

res = simulate_gene_diversification(
    GeneDiversification(n_drivers=2, lambda0=1.0, mu0=0.2, driver_speciation=1.2,
                        transfer=0.8, loss=0.3, root_drivers=1),
    n_tips=200, seed=1)
res.tree               # the tree the drivers shaped
res.driver_names()     # ['D0', 'D1']
res.tip_prevalence()   # fraction of extant tips carrying each driver
```

From the command line, the neutral genome overlay is a second, ordinary `genomes` call on the tree the
first command wrote:

```bash
zombi2 coevolve --couple genomes:species --drivers 2 --root-drivers 1 \
    --lambda0 1 --mu0 0.2 --driver-speciation 1.2 \
    --driver-transfer 0.8 --driver-loss 0.3 \
    --tips 200 --seed 1 -o out/

# overlay the neutral genome on the resulting tree (factorisation made explicit):
zombi2 genomes -t out/species_tree.nwk --trans 1 --loss 0.5 \
    --write profiles trees -o out/
```

`coevolve` writes `species_tree.nwk`, `drivers.tsv` (per-node driver presence) and
`drivers_manifest.tsv` (the effect sizes and rates). `--driver-speciation`/`--driver-extinction` are
the per-driver log-rate effects; `--driver-transfer`/`--driver-loss`/`--driver-origination` the driver
dynamics. **What it recovers:** a successful key innovation both spreads across the tips
(high `tip_prevalence`) and leaves the clades that carry it more speciose — a genomic cause of a
diversification rate shift, rather than a trait one.

![Key-innovation diversification. **A**, the model: a lineage carrying the driver (D+) speciates faster than one without it (D−, the fork width is the speciation rate); the driver is gained by origination and by frequency-dependent transfer, and lost. **B**, one realization: the driver *originates* on a single branch (the +), and from there the carrier (heavy) clade radiates into a speciose group while the non-carrier (light) lineages stay sparse. Here it reaches 69% of the extant tips — a genomic cause of a diversification-rate shift.](figures/key_innovation.pdf){width=100%}

## Joint models: ClaSSE and co-diversification

Each into-S edge has a reverse that fires at speciation (the next chapter), and switching **both** on
gives a joint model in which the branching *itself* reshapes the driver.

**ClaSSE** (traits $\leftrightarrow$ species) turns on `traits:species` *and* `species:traits`: the
trait shapes the tree *and* is kicked by its own branching. Because one arrow still points into S, the
tree is again an output.

```python
from zombi2 import Cladogenesis
res = simulate_sse(BiSSE(1, 3, 0.05, 0.05, 0.05, 0.05),
                   n_tips=200, cladogenesis=Cladogenesis(shift=0.3), seed=3)
```

```bash
zombi2 coevolve --couple traits:species --couple species:traits \
    --lambda0 1 --lambda1 3 --q01 0.05 --q10 0.05 --clado-shift 0.3 \
    --tips 200 --seed 3 -o out/
```

**Co-diversification** (species $\leftrightarrow$ genomes) is the genomic twin: the *same* driver
families both set the diversification rates and are reshuffled by a cladogenetic burst at every
speciation. Because a burst can hand one daughter a key innovation and not its sister, speciation
seeds the rate heterogeneity that then plays out along the branches.

```python
from zombi2.coevolve import simulate_co_diversification

res = simulate_co_diversification(
    GeneDiversification(3, lambda0=1.0, mu0=0.15, driver_speciation=1.0,
                        loss=0.0, origination=0.0, transfer=0.0, root_drivers=1,
                        cladogenetic_loss=0.15, cladogenetic_gain=0.2),
    n_tips=200, seed=5)
res.tip_prevalence()   # drivers still spread: the genomes:species signal survives
```

```bash
zombi2 coevolve --couple genomes:species --couple species:genomes \
    --drivers 3 --lambda0 1 --mu0 0.15 --driver-speciation 1.0 \
    --driver-loss 0 --driver-origination 0 --driver-transfer 0 --root-drivers 1 \
    --driver-clado-loss 0.15 --driver-clado-gain 0.2 --tips 200 --seed 5 -o out/
```

`--driver-clado-loss`/`--driver-clado-gain` are the per-driver drop/gain probabilities of the burst;
with both `0` the `species:genomes` arrow is off and this is plain `genomes:species`. **What it
recovers:** the diversification signal (a key innovation over-represented among the tips) *and* the
punctuational signal (sisters differing where a burst split them) at once. The reverse edges these
joints build on — `species:traits` and `species:genomes` — are the subject of the next chapter.
