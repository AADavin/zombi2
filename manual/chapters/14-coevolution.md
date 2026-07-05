# Coevolution

## The idea

ZOMBI2 normally simulates in a **pipeline**: it builds a species tree, then evolves a trait along it,
then evolves gene families along it. This works because the couplings between the three axes are
*one-directional* and the joint distribution factorises,

$$P(\text{tree}) \cdot P(\text{trait} \mid \text{tree}) \cdot P(\text{genes} \mid \text{tree}, \text{trait}),$$

so each stage can be drawn on the frozen output of the previous one. The moment a trait or a genome
feeds **back** into the process that generated the tree, this factorisation breaks: the tree, the
trait and the gene content must be grown *together*. The **`coevolve`** mode is for these coupled
scenarios.

There are three processes:

- **S** — species diversification (the birth–death process that grows the tree),
- **T** — a phenotypic trait (BM/OU/EB/Mk/threshold),
- **G** — gene-family content (the DTL process).

A **coupling is a directed edge** `driver -> target`: the driver's state modulates the target's
rates. A coevolution scenario is a **set of directed edges** on these three nodes. That single
abstraction covers everything from state-dependent diversification to full three-way feedback, and it
makes the one thing that matters — *direction* — explicit.

On the command line the coupling is one repeatable flag, `--couple driver:target`, where the order
reads as the arrow (driver first):

```bash
zombi2 coevolve --couple traits:species ...   # T->S: trait sets speciation (SSE)
zombi2 coevolve --couple species:traits ...   # S->T: speciation drives the trait
# both arrows at once = ClaSSE:
zombi2 coevolve --couple traits:species --couple species:traits ...
```

So `--couple species:traits` and `--couple traits:species` are deliberately **different models**, and
a bidirectional coupling is simply *both* edges. The `:` (rather than `->`) keeps the flag shell-safe.

The three processes give six directed edges, each a distinct model:

| Edge (`--couple`) | Direction | Model |
|---|---|---|
| `traits:species` | T $\to$ S | state-dependent diversification (SSE / ClaSSE) |
| `species:traits` | S $\to$ T | cladogenetic trait jumps at speciation |
| `genes:species` | G $\to$ S | key-innovation diversification |
| `species:genes` | S $\to$ G | punctuational (cladogenetic) genome |
| `traits:genes` | T $\to$ G | trait-linked gene families |
| `genes:traits` | G $\to$ T | gene-conditioned trait |

![The six couplings of `coevolve` mode. Each directed arrow driver $\to$ target is one model, selected with `--couple driver:target`. The two edges that point *into* S are drawn heavy: an arrow into S makes the tree depend on the coupled state, so the tree becomes an output (forward-only). The other four are overlays on a tree you supply.](figures/coevolve_modes.pdf){width=100%}

The one rule that governs their difficulty is: **does any active edge point into S?** If no edge
points into S, the tree is fixed — it is read from `-t/--tree` and every coupling is an *overlay* on a
frozen tree. If an edge does point into S (`traits:species` or `genes:species`), the tree topology
depends on the coupled state and cannot be drawn first: the tree becomes an **output**, and those runs
are forward-only and take no `-t`.

The rest of the chapter takes the three node-pairs in turn — species and traits, species and genes,
traits and genes — documenting both directions of each, and the combined model where both arrows are
on. The overlay edges all run on a species tree you provide; the examples build one with

```python
from zombi2.species import BirthDeath, simulate_species_tree, prune
tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=60, age=6.0, seed=1)
```

## Species and traits

### `traits:species` — state-dependent diversification (SSE)

A discrete or continuous trait drives speciation and extinction, and the tree is grown *jointly* with
the trait. Because the trait shapes the topology, this edge **produces** the tree (it takes no `-t`)
and a stopping condition (`n_tips` or `age`) instead. The driver is `simulate_sse`:

```python
from zombi2.coevolve import simulate_sse, BiSSE

# BiSSE: state 1 speciates three times faster, so it comes to dominate the standing tips
res = simulate_sse(BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.2, q01=0.1, q10=0.1),
                   n_tips=200, seed=1)
res.tree               # the complete tree (extinct lineages kept; prune() for the reconstructed one)
res.labeled_values()   # the trait at the extant tips
```

`BiSSE` is the **binary** state-dependent birth–death process [@maddison2007bisse]; `MuSSE` the
**k-state** variant [@fitzjohn2012diversitree]; `QuaSSE` the **continuous-trait** variant, whose
speciation and extinction are functions of the trait value [@fitzjohn2010quasse]; and `HiSSE` the
**hidden-state** model [@beaulieu2016hisse], where diversification is driven by an unobserved class
rather than the focal trait — the honest null against which a real association must be judged.

```python
import numpy as np
from zombi2.coevolve import MuSSE, QuaSSE, HiSSE

MuSSE(birth=[1, 3], death=[0.2, 0.2], Q=np.array([[-0.1, 0.1], [0.1, -0.1]]))     # k-state
QuaSSE(speciation=lambda x: 1 + 2 / (1 + np.exp(-x)),                             # continuous
       extinction=lambda x: 0.2, sigma2=0.5, rate_bound=5.0, x0=0.0)
HiSSE(classes=[BiSSE(0.5, 0.7, 0.2, 0.2, 0.1, 0.1),                               # hidden classes
               BiSSE(2.0, 3.0, 0.2, 0.2, 0.1, 0.1)], hidden_transition=0.1)
```

From the command line, `--sse-model` picks the flavour:

```bash
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --mu0 0.2 --mu1 0.2 --q01 0.1 --q10 0.1 \
    --tips 200 --seed 1 -o out/
```

It writes `species_tree.nwk` (the tree the trait's rates shaped), `traits.tsv` (every node — tips
*and* ancestral states) and `trait_tree.nwk`. `--sse-model musse` takes `--birth`/`--death` vectors
and a `--q-matrix` file; `--sse-model quasse` takes a sigmoidal speciation (`--spec-low/high/center/slope`)
and Brownian `--diffusion`. **What it recovers:** the fast-speciating state accumulates lineages, so
it dominates the standing tips — the diversification signal is written into the tree shape itself.

### `species:traits` — cladogenetic trait jumps

The reverse arrow makes the trait jump *at* each speciation rather than (or in addition to) drifting
along the branches — speciational, or *cladogenetic*, evolution: sister species differ because
something happened at their split. On its own this edge has no arrow into S, so it runs on a **given**
tree, with a `Cladogenesis` kernel layered on an ordinary anagenetic trait model. Setting the
anagenetic rates to zero gives *purely* speciational change:

```python
import numpy as np
from zombi2 import Cladogenesis, simulate_traits, Mk

# a purely speciational binary trait: no within-branch change (zero-rate Mk), a jump at each split
res = simulate_traits(tree, Mk(np.zeros((2, 2))), cladogenesis=Cladogenesis(shift=0.4), seed=2)
```

`Cladogenesis(shift=…)` is the per-daughter state-hop probability for a discrete trait;
`Cladogenesis(jump_sigma2=…)` is the Gaussian jump variance for a continuous one. On the command line:

```bash
zombi2 coevolve --couple species:traits -t species_tree.nwk \
    --sse-model bisse --q01 0 --q10 0 --clado-shift 0.4 --seed 2 -o out/
```

**What it recovers:** change concentrated at the nodes — closely related tips can differ sharply while
long unbranched stretches stay constant, the signature a purely-gradual model cannot produce.

### Both arrows: ClaSSE

Turn on **both** `traits:species` and `species:traits` and you get the full ClaSSE feedback: the trait
shapes the tree *and* is kicked by its own branching. Because one arrow points into S, the tree is
again an output:

```python
res = simulate_sse(BiSSE(1, 3, 0.05, 0.05, 0.05, 0.05),
                   n_tips=200, cladogenesis=Cladogenesis(shift=0.3), seed=3)
```

```bash
zombi2 coevolve --couple traits:species --couple species:traits \
    --lambda0 1 --lambda1 3 --q01 0.05 --q10 0.05 --clado-shift 0.3 \
    --tips 200 --seed 3 -o out/
```

## Species and genes

### `genes:species` — key-innovation diversification

Gene content can drive the tree: a small panel of binary **driver** ("key innovation") families whose
presence sets each lineage's speciation and extinction rate. Drivers are gained de novo (origination)
and — the interesting part — by **transfer**, which is *frequency-dependent*: a driver carried by more
of the live population is donated more often, so it spreads as the tree grows. That feedback is why the
tree and the gene content must grow together; this edge produces the tree. The neutral bulk genome does
not affect diversification, so it is overlaid afterward on the finished tree with the ordinary
`genomes` machinery (exact under independent families, not an approximation).

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
zombi2 coevolve --couple genes:species --drivers 2 --root-drivers 1 \
    --lambda0 1 --mu0 0.2 --driver-speciation 1.2 --driver-transfer 0.8 --driver-loss 0.3 \
    --tips 200 --seed 1 -o out/

# overlay the neutral genome on the resulting tree (the factorisation, made explicit):
zombi2 genomes -t out/species_tree.nwk --trans 1 --loss 0.5 --write profiles trees -o out/
```

`coevolve` writes `species_tree.nwk`, `drivers.tsv` (per-node driver presence) and
`drivers_manifest.tsv` (the effect sizes and rates). `--driver-speciation`/`--driver-extinction` are
the per-driver log-rate effects; `--driver-transfer`/`--driver-loss`/`--driver-origination` the driver
dynamics. **What it recovers:** a successful key innovation both spreads across the tips
(high `tip_prevalence`) and leaves the clades that carry it more speciose — a genomic cause of a
diversification rate shift, rather than a trait one.

### `species:genes` — punctuational genome

The reverse of `genes:species`, and the genomic twin of `species:traits`: here gene content does *not*
affect diversification, so it is an overlay on a **given** tree. A genome is evolved down the tree with
a **cladogenetic burst** of gene loss and gain at each speciation (a founder-effect upheaval), on top
of the usual gradual along-branch change:

```python
from zombi2.coevolve import simulate_cladogenetic_genome, CladogeneticGenome

# purely punctuational: gene content changes ONLY at speciations (no anagenetic loss/origination)
res = simulate_cladogenetic_genome(
    tree, CladogeneticGenome(initial_families=30, loss=0.0, origination=0.0,
                             cladogenetic_loss=0.15, cladogenetic_gain=3.0), seed=2)
res.genome_sizes()     # {node: family count}
res.profile_matrix()   # families × extant tips
```

```bash
zombi2 coevolve --couple species:genes -t species_tree.nwk \
    --genome-size 30 --clado-gene-loss 0.15 --clado-gene-gain 3 --seed 2 -o out/
```

`--clado-gene-loss` is the per-family drop probability at each speciation and `--clado-gene-gain` the
mean number of families gained (Poisson); `--gene-loss`/`--gene-origination` add gradual along-branch
change (both `0` = pure punctuation). It writes `Profiles.tsv`/`Presence.tsv` and `genome_sizes.tsv`.
**What it recovers:** the model's signature is that **sister tips differ** — change is injected at
their split, not spread evenly along the branches.

## Traits and genes

### `traits:genes` — trait-linked gene families

A lineage that becomes aerobic retains and acquires oxygen-using gene families; one that reverts sheds
them. This edge simulates that link: a trait is evolved down the tree, then a panel of gene families is
evolved whose **loss depends on the local trait value**. Writing the trait on a branch as $s$, a
*present* family with coupling weight $w_i$ is lost at rate

$$\text{loss}_i = \text{base\_loss} \cdot \exp(-\text{effect\_loss} \cdot w_i \cdot s),$$

so where the trait favours a responsive family it is retained, and where it does not it is purged
faster than baseline; inert families ($w_i = 0$) always lose at `base_loss`. Gain is a trait-blind
horizontal influx, and the trait-modulated loss then *selectively retains* it — the same
retention mechanism as the gene-family coupling model, which is what writes the trait–gene association
into the profiles. The driver takes the tree, a **pre-simulated trait**, and the coupling:

```python
import numpy as np
from zombi2 import simulate_traits, BrownianMotion
from zombi2.coevolve import simulate_trait_linked_genomes, TraitGeneCoupling

trait = simulate_traits(tree, BrownianMotion(sigma2=1.0), seed=1)
weights = np.zeros(40); weights[::3] = 1.0        # every third family responds to the trait
res = simulate_trait_linked_genomes(
    tree, trait, TraitGeneCoupling(n_families=40, weights=weights,
                                   effect_loss=2.0, base_loss=1.0, transfer=0.5), seed=1)
res.profiles           # ProfileMatrix: families × extant species — the trait-linked data
res.trait              # the TraitResult the genes were conditioned on
res.genomes()          # promote to a full Genomes (gene trees, reconciliations)
```

The command builds the trait for you from `--trait-model` (any `zombi2 trait` model) and selects the
responsive families with `--responsive` (a count, a fraction, or an id list):

```bash
zombi2 coevolve --couple traits:genes -t species_tree.nwk \
    --trait-model bm --sigma2 1 --panel 40 --responsive 0.3 \
    --loss 0.5 --trans 1 --effect-loss 2 --seed 1 -o out/
```

It writes the usual gene-family output plus `traits.tsv`/`trait_tree.nwk` and `coupling.tsv` (the
per-family weights, so the exact trait–gene linkage that generated the profiles is on record for
downstream inference). **What it recovers:** responsive families are present where the trait favours
them and absent where it does not, while inert families do not distinguish the states — the signal is
entirely in the responsive panel, which is what an inference method should be able to pick out.

### `genes:traits` — gene-conditioned trait

The reverse: here gene content conditions a **trait**. A binary *modifier* gene is gained and lost
along the tree, and its presence sets a continuous trait's **OU optimum** — a lineage that acquires the
gene is pulled toward a new adaptive peak (`theta_present`), one that loses it drifts back toward
`theta_absent`. "Gene presence enables a phenotypic shift." Also an overlay on a given tree:

```python
from zombi2.coevolve import simulate_gene_conditioned_trait, GeneConditionedTrait

res = simulate_gene_conditioned_trait(
    tree, GeneConditionedTrait(gene_gain=0.6, gene_loss=0.6, theta_absent=0.0,
                               theta_present=5.0, alpha=2.5, sigma2=0.4), seed=2)
res.gene_presence()    # {leaf: 0/1} — modifier presence at the tips
res.trait_values()     # {leaf: value} — carriers near 5, non-carriers near 0
```

```bash
zombi2 coevolve --couple genes:traits -t species_tree.nwk \
    --modifier-gain 0.6 --modifier-loss 0.6 --theta-absent 0 --theta-present 5 \
    --trait-alpha 2.5 --trait-sigma2 0.4 --seed 2 -o out/
```

`--modifier-gain`/`--modifier-loss` set the binary modifier's dynamics (`--root-modifier` to start it
present); `--theta-absent`/`--theta-present` are the two OU optima; `--trait-alpha` is the
mean-reversion strength (`0` = Brownian) and `--trait-sigma2` the diffusion. **What it recovers:** tips
carrying the modifier sit near `theta_present`, those without near `theta_absent` — a discrete genomic
event reading out as a shift in a continuous phenotype.

## A note on inference

::: warning
State-dependent diversification models are notorious for high false-positive rates in *inference*: a
neutral trait can look strongly diversification-linked [@maddison2015unreplicated]. ZOMBI2 is a
**forward generator**, so it is not itself fooled — but the whole point of generating these scenarios
is to test downstream inference, so it is worth generating the null (a zero-effect coupling, or the
coupling simply absent) as a control. Every edge in this chapter offers one: the `HiSSE` hidden-state
model for `traits:species`, and the inert families that `traits:genes` leaves untouched. Generate the
null alongside the signal and check that a method can tell them apart.
:::
